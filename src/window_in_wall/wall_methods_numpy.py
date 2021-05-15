from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from math import fabs

import compas
import numpy as np

def apply_x_dir(i, j, val, Kx, x_size, z_size, rhsX):
    rhsX = rhsX - Kx[n_from_ijs(i,j, z_size),:] * val
    Kx[n_from_ijs(i,j, z_size),:] = 0
    Kx[:,n_from_ijs(i,j, z_size)] = 0
    Kx[n_from_ijs(i,j, z_size),n_from_ijs(i,j, z_size)] = 1.0
    rhsX[n_from_ijs(i,j, z_size)] = val
    return rhsX

def apply_z_dir(i, j, val, Kz, x_size, z_size, rhsZ):
    rhsZ = rhsZ - Kz[n_from_ijs(i,j, z_size),:] * val
    Kz[n_from_ijs(i,j, z_size),:] = 0
    Kz[:,n_from_ijs(i,j, z_size)] = 0
    Kz[n_from_ijs(i,j, z_size),n_from_ijs(i,j, z_size)] = 1.0
    rhsZ[n_from_ijs(i,j, z_size)] = val
    return rhsZ

def compute_displacement(gate_points, x_size, z_size, coeff_diffusion=1):
    """

    """
    # stiffness matrices (diffusion equation for quad mesh topology)
    Kx = np.zeros((x_size * z_size, x_size * z_size))
    Kz = np.zeros((x_size * z_size, x_size * z_size))

    # righthand side (displacement boundary conditions)
    rhsX = np.zeros(x_size*z_size)
    rhsZ = np.zeros(x_size*z_size)
    
    n = x_size * z_size

    for i in range(n):
        Kx[i,i] = coeff_diffusion * 4.0

        if (i)>0 and (i%z_size) != 0:
            Kx[i,(i-1)] = Kx[i,(i-1)] - coeff_diffusion         
        if (i<n-1) and (i % z_size) != (z_size-1):
            Kx[i,(i+1)] = Kx[i,(i+1)]-coeff_diffusion
        if (n-i)>z_size:
            Kx[i,i+z_size] = Kx[i,i+z_size]-coeff_diffusion
        if i>(z_size-1):
            Kx[i,i-z_size] = Kx[i,i-z_size]-coeff_diffusion
        if (i%z_size) == 0:
            Kx[i,(i+1)]=Kx[i,(i+1)]-coeff_diffusion
        if (i%z_size) == z_size-1:
            Kx[i,(i-1)] = Kx[i,(i-1)]-coeff_diffusion

    for i in range(n):
        Kz[i,i] = coeff_diffusion * 4.0

        if (i)>0 and (i%z_size)!=0:
            Kz[i,(i-1)]=Kz[i,(i-1)]-coeff_diffusion         
        if (i<n-1) and (i%z_size)!=(z_size-1):
            Kz[i,(i+1)] = Kz[i,(i+1)]-coeff_diffusion
        if (n-i) > z_size:
            Kz[i,i+z_size] = Kz[i,i+z_size]-coeff_diffusion
        if i > (z_size-1):
            Kz[i,i-z_size] = Kz[i,i-z_size]-coeff_diffusion
        if i < z_size:
            Kz[i,(i+z_size)] = Kz[i,(i+z_size)]-coeff_diffusion
        if i/z_size > x_size-1:
            Kz[i,(i-z_size)] = Kz[i,(i-z_size)]-coeff_diffusion

    for p in (gate_points):
        rhsX = apply_x_dir(p[0], p[1], p[2], Kx, x_size, z_size, rhsX)
        rhsZ = apply_z_dir(p[0], p[1], p[3], Kz, x_size, z_size, rhsZ)
   
    solX = np.linalg.solve(Kx,rhsX)
    solZ = np.linalg.solve(Kz,rhsZ)

    return (solX, solZ)

def n_from_ijs(i,j,size):
    return (i*size+j)