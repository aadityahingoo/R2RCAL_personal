#%% import

import numpy as np, trimesh
from scipy.ndimage import map_coordinates
from skimage.measure import marching_cubes


#%% set params and file

R, H, W, dx = 5.0, 1.0, 10.0, 0.02  # mm
file = 'Part Studio 1 - Part 1 (3).stl'

flat = trimesh.load(file).voxelized(dx).fill().matrix.astype(float)  # axes u,v,w

#%% 

n = int(2*(R+H)/dx) + 2
x = (np.arange(n) - n/2) * dx
X, Y, Z = np.meshgrid(x, x, np.arange(flat.shape[2])*dx, indexing='ij')
r, th = np.hypot(X, Y), np.arctan2(Y, X)
u, v = R*(th + np.pi), r - R  # arc length measured at the web surface

wrapped = map_coordinates(flat, [u/dx, v/dx, Z/dx], order=1, cval=0) > 0.5

vt, fc, _, _ = marching_cubes(wrapped.astype(np.uint8), 0.5, spacing=(dx,)*3)
trimesh.Trimesh(vt, fc).export('wrapped_voxel.stl')

