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

    """calculate an x displacement for every non-gate non-boundary point, 
    given the pre-assigned displacements of the gate points. The left and 
    right edges of the wall are treated as fixed boundaries (Dirichlet) and
    top and bottom as zero-gradient (Neumann)
    """    
    n = z_size
    m = x_size
    nm= n*m
    rhs = np.zeros(nm)              #equation right hand side
    dia1 = np.ones(nm)*-4           #main diagonal of K matrix
    dia2 = np.ones(nm-1)            #second diagonal (up) corresponding to right neighbors
    dia2[np.arange(0,nm-1,n)]+=1    #Neumann BC using ghost-points
    dia2[np.arange(n-1,nm-1,n)]-=1  #Neumann BC using ghost-points
    dia3 = np.ones(nm-1)            #diagonal corresponding to left neighbors
    dia3[np.arange(n-2,nm-1,n)]+=1  #Neumann
    dia3[np.arange(n-1,nm-1,n)]-=1  #Neumann
    dia4 = np.ones(nm-n)            #diagonal corresponding to up neghbors
    dia5 = np.ones(nm-n)            #diagonal corresponding to bottom neighbors
    #sparse coef. (stiffness) matrix built out of 5 diagonals
    K = sp.sparse.diags([dia1, dia2, dia3, dia4, dia5], [0, 1, -1, n, -n], format='csc')
    #BCdofs is the vector of degrees of freedom with dirichlet (displacement) boundary condition
    #at first left and right sides are added to BCdofs
    BCdofs=np.concatenate((np.arange(0,n),np.arange(nm-n,nm)))
    #assigning zero displacement to left and right
    BCvals=np.zeros(BCdofs.size)
    #for applying the dirichlet, for dofs with displacement BC, the rhs is set to the according dispalcement
    rhs[BCdofs]=BCvals
    #loop over gate points as additional displacement BCs (non-zero, in contrast to left & right edges)
    for vertex in (gate_points):
        glob_id = mesh.vertex_attribute(vertex, "glob_id")
        x_disp = mesh.vertex_attribute(vertex, "x_disp")
        rhs[glob_id]=x_disp
        BCdofs = np.append(BCdofs,glob_id)
    #constructing an identitiy matrix (named Iinter) of size (nm x nm) with zero value on diagonals
    # of displacement BC dofs. If applied on (multiplied by) K, only non-BC (internal) dofs will remain 
    diaInter=np.ones(nm)
    diaInter[BCdofs]=0.0
    Iinter = sp.sparse.diags([diaInter], [0], format='csc')
    #constructing an identitiy matrix (named Ibc) of size (nm x nm) with zero value on diagonals
    # of internal (non-BC) dofs. If applied on (multiplied by) K, only BC (internal) dofs will remain 
    diaBC=np.zeros(nm)
    diaBC[BCdofs]=1.0
    Ibc = sp.sparse.diags([diaBC], [0], format='csc')
    #applying dirichlet on K, by zeroing out rows and columns of BC dofs and setting BC-diagonals to 1
    K_BC= Iinter * K * Iinter + Ibc
    #modifying the rhs for non-BC dofs to account for the eliminated dofs 
    # the operation below assignes -K_internal*x_BC to rhs_internal (and doesn't change rhs_BC) 
    rhs = rhs - Iinter * (K-(Ibc * K ))* rhs
    #solving the system
    sol = scipy.sparse.linalg.spsolve(K_BC,rhs)
    return sol