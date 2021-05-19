from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from math import fabs

import compas
import numpy as np

def apply_x_dir(glob_id, x_disp, Kx, rhsX):
    rhsX = rhsX - Kx[glob_id,:] * x_disp
    Kx[glob_id,:] = 0
    Kx[:,glob_id] = 0
    Kx[glob_id,glob_id] = 1.0
    rhsX[glob_id] = x_disp
    return rhsX

def apply_z_dir(glob_id, z_disp, Kz, rhsZ):
    rhsZ = rhsZ - Kz[glob_id,:] * z_disp
    Kz[glob_id,:] = 0
    Kz[:,glob_id] = 0
    Kz[glob_id,glob_id] = 1.0
    rhsZ[glob_id] = z_disp
    return rhsZ

def compute_displacement(mesh, gate_points, x_size, z_size, coeff_diffusion=1):
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

    for vertex in (gate_points):
        glob_id = mesh.vertex_attribute(vertex, "glob_id")
        x_disp = mesh.vertex_attribute(vertex, "x_disp")
        z_disp = mesh.vertex_attribute(vertex, "z_disp")
        rhsX = apply_x_dir(glob_id, x_disp, Kx, rhsX)
        rhsZ = apply_z_dir(glob_id, z_disp, Kz, rhsX)
   
    solX = np.linalg.solve(Kx,rhsX)
    solZ = np.linalg.solve(Kz,rhsZ)

    return (solX, solZ)

def n_from_ijs(i,j,size):
    return (i*size+j)