#!/usr/bin/env python3
"""
wrap_layers.py - flat-domain STL -> per-layer wrapped STLs for CAL web printing.

Pipeline: 1 orient -> 2 position in print volume -> 3 slice into layers -> 4 wrap -> 5 export

Flat frame (after orientation):
    X = u   along the web = arc length around the mandrel
    Y = w   across the web width (becomes the vial / rotation axis)
    Z = v   stack axis = height above the web (each layer's bottom sits on the web)

Wrapped frame (vial axis = Z, origin on the rotation axis):
    theta = u / R,   r = R + v,   Z = w      R = D/2 = vial OD/2 (web wraps on the vial, thickness 0)
    The seam (where the web ends meet) is centered on theta = 0.

Requires: pip install trimesh manifold3d numpy matplotlib
Example:  python wrap_layers.py vessel.stl -D 10 -W 12 -H 1 --gap 1
"""
import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import trimesh

AXES = ['x', 'y', 'z', '-x', '-y', '-z']
UNIT = {'x': (1, 0, 0), 'y': (0, 1, 0), 'z': (0, 0, 1)}


# ---------- 1. orientation ----------
def axis_vec(s):
    return (-1.0 if s.startswith('-') else 1.0) * np.array(UNIT[s[-1]], float)


def orient_matrix(along, up):
    """Proper rotation taking model axis `along` -> X (u) and `up` -> Z (v)."""
    a, z = axis_vec(along), axis_vec(up)
    if abs(a @ z) > 1e-9:
        raise ValueError(f'--along {along} and --up {up} must be different axes')
    M = np.eye(4)
    M[:3, :3] = np.vstack([a, np.cross(z, a), z])  # (a, z x a, z) is right-handed
    return M


def evaluate(mesh, M, vol):
    """Fit check for one orientation: (fits, n_layers, slack, extents)."""
    ext = np.ptp(mesh.vertices @ M[:3, :3].T, axis=0)
    fits = ext[0] <= vol.arc and ext[1] <= vol.wid
    n = int(np.ceil(ext[2] / vol.H - 1e-6))
    return fits, n, min(vol.arc - ext[0], vol.wid - ext[1]), ext


def auto_orient(mesh, vol):
    """Among the 24 axis-aligned orientations: must fit, then fewest layers, then most slack."""
    best = None
    for along, up in itertools.product(AXES, AXES):
        if along[-1] == up[-1]:
            continue
        fits, n, slack, _ = evaluate(mesh, orient_matrix(along, up), vol)
        key = (not fits, n, -slack)
        if best is None or key < best[0]:
            best = (key, along, up)
    return best[1], best[2]


def mesh_issues(m):
    out = []
    if not m.is_watertight:
        out.append('not watertight')
    if not m.is_winding_consistent:
        out.append('inconsistent winding')
    if m.is_watertight and m.volume <= 0:
        out.append('inverted normals')
    return out


# ---------- 2. print volume + positioning ----------
class Volume:
    def __init__(s, a):
        s.D, s.W, s.H, s.gap, s.margin = a.diameter, a.width, a.height, a.gap, a.margin
        s.R = s.D / 2
        s.circ = np.pi * s.D           # full web length around the mandrel
        s.arc = s.circ - s.gap         # usable u
        s.wid = s.W - 2 * s.margin     # usable w


def position(m, vol):
    """Center in usable arc (seam opposite the part), center on web width (Z=0), bottom on web."""
    lo, hi = m.bounds
    ext = hi - lo
    u0 = vol.gap / 2 + (vol.arc - ext[0]) / 2
    m.apply_translation([u0 - lo[0], -ext[1] / 2 - lo[1], -lo[2]])


# ---------- 3. slicing ----------
def slice_layers(m, vol, n):
    lo, hi = m.bounds
    pad = 1.0
    layers = []
    for i in range(n):
        z0 = i * vol.H if i > 0 else lo[2] - pad          # avoid coplanar cut at part bottom/top
        z1 = (i + 1) * vol.H if i < n - 1 else hi[2] + pad
        box = trimesh.creation.box(bounds=[[lo[0] - pad, lo[1] - pad, z0], [hi[0] + pad, hi[1] + pad, z1]])
        s = trimesh.boolean.intersection([m, box], engine='manifold')
        if s.is_empty:
            layers.append(None)
            continue
        s.apply_translation([0, 0, -i * vol.H])             # local v: 0 = web surface
        layers.append(s)
    return layers


# ---------- 4. wrapping ----------
def wrap(m, vol, max_edge):
    v, f = trimesh.remesh.subdivide_to_size(m.vertices, m.faces, max_edge=max_edge, max_iter=40)
    th, r = v[:, 0] / vol.R, vol.R + v[:, 2]
    # (u, w, v) -> (theta, Z, r) is a cyclic permutation of cylindrical (r, theta, Z): normals preserved
    return trimesh.Trimesh(np.column_stack([r * np.cos(th), r * np.sin(th), v[:, 1]]), f)


# ---------- preview ----------
def preview(m, wrapped0, vol, n, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Rectangle

    kw = dict(facecolor='tab:blue', edgecolor='none', alpha=0.15)
    tri = m.vertices[m.faces]
    fig, ax = plt.subplots(1, 3, figsize=(16, 5), gridspec_kw={'width_ratios': [1.3, 1.3, 1]})

    ax[0].add_collection(PolyCollection(tri[:, :, [0, 1]], **kw))
    ax[0].add_patch(Rectangle((0, -vol.W / 2), vol.circ, vol.W, fill=False, ec='k', label='web'))
    ax[0].add_patch(Rectangle((vol.gap / 2, -vol.wid / 2), vol.arc, vol.wid, fill=False, ec='g', ls='--', label='usable'))
    ax[0].set(title='Top: web unrolled', xlabel='u along web [mm]', ylabel='w across web [mm]', aspect='equal')
    ax[0].legend(loc='upper right', fontsize=8)

    ax[1].add_collection(PolyCollection(tri[:, :, [0, 2]], **kw))
    for i in range(n + 1):
        ax[1].axhline(i * vol.H, color='r', lw=0.6)
    ax[1].axhline(0, color='k', lw=2)
    ax[1].set(title=f'Side: {n} layers x {vol.H} mm (black = web)', xlabel='u along web [mm]', ylabel='v stack [mm]')

    if wrapped0 is not None:
        t = wrapped0.vertices[wrapped0.faces]
        ax[2].add_collection(PolyCollection(t[:, :, [0, 1]], **kw))
    th = np.linspace(0, 2 * np.pi, 400)
    ax[2].plot(vol.R * np.cos(th), vol.R * np.sin(th), 'k-', lw=1, label='web')
    ax[2].plot((vol.R + vol.H) * np.cos(th), (vol.R + vol.H) * np.sin(th), 'k--', lw=0.8, label='resin top')
    ax[2].plot([vol.R, vol.R + vol.H], [0, 0], 'r-', lw=2, label='seam')
    ax[2].set(title='Layer 0 wrapped (down vial axis)', xlabel='X [mm]', ylabel='Y [mm]', aspect='equal')
    ax[2].legend(loc='upper right', fontsize=8)

    for x in ax:
        x.autoscale_view()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------- main ----------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('stl', help='flat-domain part STL (matrix = solid, vessels = voids)')
    p.add_argument('-D', '--diameter', type=float, required=True, help='vial outer diameter = web diameter (web thickness 0) [mm]')
    p.add_argument('-W', '--width', type=float, required=True, help='web width [mm]')
    p.add_argument('-H', '--height', type=float, required=True, help='cast resin height = layer thickness [mm]')
    p.add_argument('--gap', type=float, default=1.0, help='unprinted seam where the web ends meet [mm]')
    p.add_argument('--margin', type=float, default=0.5, help='clearance from each web edge [mm]')
    p.add_argument('--along', choices=AXES, help='model axis that runs along the web (-> u)')
    p.add_argument('--up', choices=AXES, help='model axis to stack layers along (-> v)')
    p.add_argument('--scale', type=float, default=1.0, help='multiply STL coords (e.g. 1000 if exported in m)')
    p.add_argument('--tol', type=float, default=0.001, help='max chord error of wrapped facets [mm]')
    p.add_argument('-o', '--out', default='layers', help='output folder')
    p.add_argument('--no-preview', action='store_true')
    a = p.parse_args()

    vol = Volume(a)
    out = Path(a.out)
    (out / 'flat').mkdir(parents=True, exist_ok=True)
    (out / 'wrapped').mkdir(parents=True, exist_ok=True)

    # 1. load + sanity
    mesh = trimesh.load(a.stl, force='mesh')
    mesh.apply_scale(a.scale)
    issues = mesh_issues(mesh)
    if 'inverted normals' in issues:
        mesh.invert()
        print('! normals were inverted; fixed')
        issues.remove('inverted normals')
    if issues:
        sys.exit(f'x input mesh: {", ".join(issues)}. Repair before slicing.')
    if not 0.5 < mesh.extents.max() < 500:
        print(f'! extents {np.round(mesh.extents, 3)} look wrong for mm; check --scale')

    # 1. orient
    if bool(a.along) != bool(a.up):
        sys.exit('x give both --along and --up, or neither (auto)')
    along, up = (a.along, a.up) if a.along else auto_orient(mesh, vol)
    M = orient_matrix(along, up)
    fits, n, slack, ext = evaluate(mesh, M, vol)
    print(f'1 orientation: model {along} -> u (along web), model {up} -> v (stack)'
          f'{"" if a.along else "  [auto]"}')
    print(f'  part u x w x v = {ext[0]:.2f} x {ext[1]:.2f} x {ext[2]:.2f} mm')
    if not fits:
        sys.exit(f'x does not fit: needs u <= {vol.arc:.2f}, w <= {vol.wid:.2f} mm')
    mesh.apply_transform(M)

    # 2. position
    position(mesh, vol)
    print(f'2 volume: circumference {vol.circ:.2f} mm, usable {vol.arc:.2f} x {vol.wid:.2f} mm, '
          f'layer {vol.H} mm; slack {slack:.2f} mm')
    mesh.export(out / 'part_oriented_flat.stl')

    # 3-5. slice, wrap, export
    max_edge = np.sqrt(8 * vol.R * a.tol)
    manifest = dict(params=vars(a), orientation=dict(along=along, up=up), n_layers=n, layers=[])
    wrapped0 = None
    for i, layer in enumerate(slice_layers(mesh, vol, n)):
        rec = dict(index=i, z_part=[round(i * vol.H, 4), round(min((i + 1) * vol.H, ext[2]), 4)])
        if layer is None:
            rec['status'] = 'empty, skipped'
            print(f'3 layer {i:02d}: empty, skipped')
            manifest['layers'].append(rec)
            continue
        w = wrap(layer, vol, max_edge)
        if wrapped0 is None:
            wrapped0 = w
        name = f'layer_{i:02d}.stl'
        layer.export(out / 'flat' / name)
        w.export(out / 'wrapped' / name)
        probs = mesh_issues(w)
        rec.update(file=name, volume_mm3=round(layer.volume, 4), faces=len(w.faces),
                   issues=probs, status='ok' if not probs else 'check')
        manifest['layers'].append(rec)
        print(f'3-5 layer {i:02d}: {layer.volume:7.3f} mm3, {len(w.faces):7d} faces'
              f'{"  ! " + ", ".join(probs) if probs else ""}')

    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    if not a.no_preview:
        preview(mesh, wrapped0, vol, n, out / 'preview.png')
    print(f'done -> {out}/  (flat/, wrapped/, manifest.json{"" if a.no_preview else ", preview.png"})')


if __name__ == '__main__':
    main()
