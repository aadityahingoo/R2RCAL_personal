# Project: layered CAL printing of vascularized hydrogel matrices

## Goal
Print hydrogel bioresin layers with computed axial lithography (CAL) on a web that is
rigid in tension, then stack and bond the layers into a 3D vascularized tissue construct.
Vessels exist for nutrient and metabolite diffusion. Motivation: extrusion bioprinting is
slow; SLA damages fragile hydrogels.

Based on R2R CAL (Toombs, Li, Taylor et al., UC Berkeley): https://arxiv.org/abs/2402.10955
The photoresist on a web is "unwrapped" from the tomographic volume.

## Current setup
- The R2R line is not working, so we use **CAL6** (standard rotating vial).
- The web is wrapped around the outside of the vial; resin is cast on the outside of the web.
  Radially outward: vial -> web -> resin. Web thickness is taken as 0.
- Design is done in Onshape in the **flat (unwrapped) domain**, exported as STL, then transformed.
- The user runs this on a Dell XPS 15 with Python and VS Code. Keep compute light.
- The user wants **brief, crisp answers**.

## Conventions (do not change without asking)
Flat frame:
- X = u, along the web (arc length)
- Y = w, across the web (becomes the vial axis)
- Z = v, the stack axis (height above web; each layer's bottom sits on the web)

Wrapped frame (vial axis = Z, origin on the rotation axis):
- theta = u / R, r = R + v, Z = w
- R = D/2 = vial OD/2 (the web sits on the vial)
- The seam is centered at theta = 0
- The map assumes the web is the inextensible reference surface, so the gel strains by
  about v/R. Keep H/R small.

## Files
- `wrap_layers.py`: orient (manual `--along/--up` or auto over 24 orientations) ->
  position (centered in usable arc and width, bottom on web) -> slice into H-thick layers
  (manifold boolean) -> subdivide + wrap -> export `flat/`, `wrapped/`, `manifest.json`,
  `preview.png`.
- Inputs: `-D` (web diameter), `-W` (web width), `-H` (resin height = layer thickness),
  `--gap`, `--margin`, `--scale`, `--tol` (chord error, use 0.005 on laptop),
  `--no-preview`.
- Example:
  `python wrap_layers.py vessel.stl -D 10 -W 12 -H 1 --gap 1 --tol 0.005`
- Validated on a test block: layer volumes sum to the part volume, wrapped volume matches
  (1 + v/R) scaling, and all layers are watertight within r in [R, R+H].
- Dependencies: `pip install trimesh manifold3d matplotlib numpy`

## Known issues / next steps
- Wrapped STLs are large (about 20 MB per layer at tol 1 um).
- The print target for VAMToolbox should be voxels. Add voxel export per layer (an inverse
  map via scipy `map_coordinates` was sketched earlier).
- Curling and residual stress from curing in the wrapped state (seen in the paper) are
  likely worse in soft hydrogels. Consider compensating, or delaying polymerization onset.
- Layer interfaces should cut vessels cleanly, not along thin walls. Stack axis choice is
  a design decision.
- Future: registration features for stacking and bonding the layers.
