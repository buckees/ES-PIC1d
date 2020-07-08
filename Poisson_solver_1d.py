# -*- coding: utf-8 -*-
# 1d Poisson Solver

import Constants as cst
import numpy as np
import matplotlib.pyplot as plt

def Poisson_solver_1d(ncellx,width,den_chrg,bc):
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
    # construct grids
    gridx, dx = np.linspace(0.0,width,ncellx+1,retstep=True)
    den_chrg = den_chrg*cst.UNIT_CHARGE
    # construct A
    A = np.zeros((ncellx+1,ncellx+1),dtype=np.float)
    A[0,0]       =  2.0; A[1,0]       = -1.0;
    A[-2,ncellx] = -1.0; A[-1,ncellx] =  2.0;
    for i in range(1,ncellx):
        A[i-1,i] = -1.0
        A[  i,i] =  2.0
        A[i+1,i] = -1.0
    invA = np.linalg.inv(A)
    a = (bc[-1] - bc[0])/width
    b = bc[0] + np.matmul(invA,den_chrg)[0]/cst.EPS0*dx*dx
    pot = -np.matmul(invA,den_chrg)/cst.EPS0*dx*dx + a*gridx + b
    return pot

width = 1.0
ncellx = 100
#den_chrg = np.zeros((ncellx+1,),dtype=np.float)
gridx, dx = np.linspace(0.0,width,ncellx+1,retstep=True)
den_chrg = np.power(gridx,1)*1e11
den_chrg.shape

pot = Poisson_solver_1d(ncellx,width,den_chrg,(100.0,0.0))
pot.shape
# diagnostic plot
fig, (ax0,ax1) = plt.subplots(2,1, figsize=(8,4),
      constrained_layout=True)
ax0.plot(gridx,den_chrg)
ax1.plot(gridx,pot)
plt.show()
print(pot[0],pot[-1])