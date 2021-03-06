# -*- coding: utf-8 -*-
# 1d Poisson Solver

import Constants as cst
import numpy as np
from scipy.signal import savgol_filter

def calc_invA(mesh):
    """
    
         2, -1,
        -1, 2, -1,             
    A = [   -1,  2, ...         ]
                 ...
                              -1
                          -1,  2
    compute invert A
    """
    # construct A
    A = np.zeros((mesh.ncellx+1,mesh.ncellx+1),dtype=np.float)
    A[0,0]       =  2.0; A[1,0]       = -1.0;
    A[-2,mesh.ncellx] = -1.0; A[-1,mesh.ncellx] =  2.0;
    for i in range(1,mesh.ncellx):
        A[i-1,i] = -1.0
        A[  i,i] =  2.0
        A[i+1,i] = -1.0
    
    return np.linalg.inv(A)
    

def Poisson_solver_1d(mesh,den_chrg,bc,invA):
    """
    Solve 1d Poisson's equation, return potential
    Using finite diferences, inverse the matrix
    Poisson's equation in vacuum
    d2/dx2*phi(x) = -rho(x)/eps0
    discretize it with finite differences
    A*P(x)/dx2 = -rho(x)/eps0

         2, -1,
        -1, 2, -1,             
    A = [   -1,  2, ...         ]
                 ...
                              -1
                          -1,  2

    phi(x) = -invA*rho(x)/eps0*dx2 + ax + b
    a = (phi[-1] - phi[0])/width
    b = phi[0] + (invA*rho)[0]*/eps0*dx2
    :param ncellx: number of cells
    :param dx: cell size, in meter
    :param bc=(bc_l,bc_r): a 2-tuple, boundary conditions
    """
    den_chrg = den_chrg*cst.UNIT_CHARGE/cst.EPS0
    invA_den_chrg = np.matmul(invA,den_chrg)*mesh.dx*mesh.dx
    a = (bc[-1] - bc[0]+invA_den_chrg[-1]-invA_den_chrg[0])/mesh.width
    b = bc[0] + invA_den_chrg[0] - a*mesh.gridx[0]
    
    pot = invA_den_chrg + a*mesh.gridx + b
    pot = savgol_filter(pot, 11, 3) # window size 11, polynomial order 3
    efld = -(pot[1:] - pot[0:-1])/mesh.dx
    return [pot, efld]
