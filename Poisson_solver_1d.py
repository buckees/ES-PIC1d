# -*- coding: utf-8 -*-
# 1d Poisson Solver

import Constants as cst
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

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
    den_chrg = den_chrg*cst.UNIT_CHARGE/cst.EPS0
    # construct A
    A = np.zeros((ncellx+1,ncellx+1),dtype=np.float)
    A[0,0]       =  2.0; A[1,0]       = -1.0;
    A[-2,ncellx] = -1.0; A[-1,ncellx] =  2.0;
    for i in range(1,ncellx):
        A[i-1,i] = -1.0
        A[  i,i] =  2.0
        A[i+1,i] = -1.0
    
    invA = np.linalg.inv(A)
    invA_den_chrg = np.matmul(invA,den_chrg)*dx*dx
    a = (bc[-1] - bc[0]+invA_den_chrg[-1]-invA_den_chrg[0])/width
    b = bc[0] + invA_den_chrg[0] - a*gridx[0]
    
    pot = invA_den_chrg + a*gridx + b
    pot = savgol_filter(pot, 7, 2) # window size 10, polynomial order 3
    efld = [-(pot[i+1]-pot[i])/dx for i in range(ncellx)]
    return pot, efld

# plot diagnostics
# validated with function below
#phi_f = @(t) t.*cos(t); % inline function for the exact solution
#rho_f = @(t) 2.*sin(t) + t.*cos(t); % inline function for the exact right-hand-side
#
#width = 1.0
#ncellx = 100
#gridx, dx = np.linspace(0.0,width,ncellx+1,retstep=True)
###den_chrg = (1-np.power(gridx,2))*1e11
##den_chrg = np.sin(gridx) + np.multiply(gridx,np.cos(gridx))*1e12
##den_chrg = np.zeros((ncellx+1,),dtype=np.float)
#den_chrg = [1.0 if i < 20 or i > 80 else 0 for i in range(ncellx+1) ]
#den_chrg = np.asarray(den_chrg)*1e10
#
#pot, efld = Poisson_solver_1d(ncellx,width,den_chrg,(0.0,0.0))
## diagnostic plot
#fig, (ax0,ax1,ax2) = plt.subplots(1,3, figsize=(9,3),
#      constrained_layout=True)
#ax0.plot(gridx,den_chrg)
#ax0.set_title('charge distribution')
#ax1.plot(gridx,pot)
#ax1.set_title('potential distribution')
#ax2.plot(gridx[1:],efld)
#ax2.set_title('E-field distribution')
#plt.show()
#print(pot[0],pot[-1],efld[0],efld[-1])