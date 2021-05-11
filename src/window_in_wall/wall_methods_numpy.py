from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from math import fabs

import compas
import numpy as np

def apply_x_dir(i, j, val, Kx, x_size, z_size):
    rhsX = np.zeros(x_size*z_size)
    
    rhsX = rhsX - Kx[n_from_ij(i,j, z_size),:] * val
    Kx[n_from_ij(i,j, z_size),:] = 0
    Kx[:,n_from_ij(i,j, z_size)] = 0
    Kx[n_from_ij(i,j, z_size),n_from_ij(i,j, z_size)] = 1.0
    rhsX[n_from_ij(i,j, z_size)] = val
    
    return rhsX

def apply_z_dir(i, j, val, Kz, x_size, z_size):
    rhsZ = np.zeros(x_size*z_size)

    rhsZ = rhsZ - Kz[n_from_ij(i,j, z_size),:] * val
    Kz[n_from_ij(i,j, z_size),:] = 0
    Kz[:,n_from_ij(i,j, z_size)] = 0
    Kz[n_from_ij(i,j, z_size),n_from_ij(i,j, z_size)] = 1.0
    rhsZ[n_from_ij(i,j, z_size)] = val

    return rhsZ

def compute_displacement(mesh, gate_points, x_size, z_size):

    Kx = np.zeros((x_size * z_size, x_size * z_size))
    Kz = np.zeros((x_size * z_size, x_size * z_size))
    
    n = x_size * z_size
    kPres = 1.0

    for i in range(n):
        Kx[i,i] = kPres * 4.0

        if (i)>0 and (i%z_size) != 0:
            Kx[i,(i-1)] = Kx[i,(i-1)] - kPres         
        if (i<n-1) and (i % z_size) != (z_size-1):
            Kx[i,(i+1)] = Kx[i,(i+1)]-kPres
        if (n-i)>z_size:
            Kx[i,i+z_size] = Kx[i,i+z_size]-kPres
        if i>(z_size-1):
            Kx[i,i-z_size] = Kx[i,i-z_size]-kPres
        if (i%z_size) == 0:
            Kx[i,(i+1)]=Kx[i,(i+1)]-kPres
        if (i%z_size) == z_size-1:
            Kx[i,(i-1)] = Kx[i,(i-1)]-kPres

    for i in range(n):
        Kz[i,i] = kPres * 4.0

        if (i)>0 and (i%z_size)!=0:
            Kz[i,(i-1)]=Kz[i,(i-1)]-kPres         
        if (i<n-1) and (i%z_size)!=(z_size-1):
            Kz[i,(i+1)] = Kz[i,(i+1)]-kPres
        if (n-i) > z_size:
            Kz[i,i+z_size] = Kz[i,i+z_size]-kPres
        if i > (z_size-1):
            Kz[i,i-z_size] = Kz[i,i-z_size]-kPres
        if i < z_size:
            Kz[i,(i+z_size)] = Kz[i,(i+z_size)]-kPres
        if i/z_size > x_size-1:
            Kz[i,(i-z_size)] = Kz[i,(i-z_size)]-kPres
    
    for p in (gate_points):
        rhsX = apply_x_dir(p[0], p[1], p[2], Kx, x_size, z_size)
        rhsZ = apply_z_dir(p[0], p[1], p[3], Kz, x_size, z_size)
   
    solX = np.linalg.solve(Kx,rhsX)
    solZ = np.linalg.solve(Kz,rhsZ)

    return (solX, solZ)

def n_from_ij(i,j,size):
    return (i*size+j)