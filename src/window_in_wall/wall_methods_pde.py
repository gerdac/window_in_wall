from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from math import fabs

import compas
import numpy as np
import scipy as sp
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

def compute_displacement_x(mesh, gate_points, x_size, z_size, coeff_diffusion=1):
    n = z_size
    m = x_size
    nm= n*m
    rhs = np.zeros(nm) #RHS
    dia1 = np.ones(nm)*-4
    dia2 = np.ones(nm-1)
    dia2[np.arange(0,nm-1,n)]+=1
    dia2[np.arange(n-1,nm-1,n)]-=1
    dia3 = np.ones(nm-1)
    dia3[np.arange(n-2,nm-1,n)]+=1
    dia3[np.arange(n-1,nm-1,n)]-=1
    dia4 = np.ones(nm-n)
    dia5 = np.ones(nm-n)
    K = sp.sparse.diags([dia1, dia2, dia3, dia4, dia5], [0, 1, -1, n, -n], format='csc')
    #leftright
    BCdofs=np.concatenate((np.arange(0,n),np.arange(nm-n,nm)))
    #up#
    #BCdofs=np.append(BCdofs,np.arange(n-1,nm,n))
    BCvals=np.zeros(BCdofs.size)
    rhs[BCdofs]=BCvals
    for vertex in (gate_points):
        glob_id = mesh.vertex_attribute(vertex, "glob_id")
        x_disp = mesh.vertex_attribute(vertex, "x_disp")
        rhs[glob_id]=x_disp
        BCdofs = np.append(BCdofs,glob_id)
    diaInter=np.ones(nm)
    diaInter[BCdofs]=0.0
    Iinter = sp.sparse.diags([diaInter], [0], format='csc')
    diaBC=np.zeros(nm)
    diaBC[BCdofs]=1.0
    Ibc = sp.sparse.diags([diaBC], [0], format='csc')
    A_BC= Iinter * K * Iinter + Ibc
    rhs = rhs - Iinter * (K-(Ibc * K ))* rhs
    sol = scipy.sparse.linalg.spsolve(A_BC,rhs)
    return sol