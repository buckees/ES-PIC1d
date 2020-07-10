# -*- coding: utf-8 -*-
# This ES-PIC_run.py runs the PIC simulation
# All in SI unit

from datetime import timedelta 
import time
t0 = time.time()

# modules from *.py in the same folder
import Constants as cst
from ESPIC1d_mesh import Mesh
import ESPIC1d_init as init
import ESPIC1d_move as move
import ESPIC1d_out as out
from Species import Eon, Arp, Ar
import Poisson_solver_1d as ps1d

# python built-in modules
import numpy as np
from scipy.signal import savgol_filter

"""
Geometry and Mesh are created in ESPIC1d_mesh.py
"""

# Initial Conditions
press_mT = 100.0 # pressure, in mTorr
pressure = press_mT/1.0e3*cst.TORR2PA # in Pa, 1 Torr = 133.322 Pa
Ar_den_init = pressure/(cst.KB*Ar.tmpt*cst.EV2K) # in m-3, N/V = P/(kb*T) ideal gas law
Eon_temp = [Eon.tmpt] # initial temperature for Eon, in eV
Arp_temp = [Arp.tmpt] # initial temperature for Arp, in eV
bc = (0.0, 0.0) # left and right boundary conditions, in Volt 
Eon_den_init = 1.0e15 # in m-3, initial electron density
den_limit = 1.0e11 # in m-3, lower limit of denisty, avoid 0 density 

# Operation Parameters
num_ptcl = 100000 # number of particles, should be >> ncellx to reduce noise
dt = 1.0e-13 # in sec

den_per_ptcl = Eon_den_init/num_ptcl # density contained in one particle

# initialize the position and velcotiy in a dataframe
Eon_pv = init.init_data(num_ptcl, Eon) # a list contain [posn, vels]
Arp_pv = init.init_data(num_ptcl, Arp)
Eon_pv[0] = Eon_pv[0]*Mesh.width
Arp_pv[0] = Eon_pv[0].copy()

# assign charge densities to grid nodes, unit in UNIT_CHARGE
Eon_den = move.den_asgmt(Eon_pv[0], Mesh)
Arp_den = move.den_asgmt(Arp_pv[0], Mesh)
# add up all charge densities
chrg_den = (Eon.charge*Eon_den + 
            Arp.charge*Arp_den)*den_per_ptcl # unit in UNIT_CHARGE
chrg_den = savgol_filter(chrg_den, 11, 3) # window size 11, polynomial order 3

# update potential according to assigned charges to nodes
invA = ps1d.calc_invA(Mesh)
pe = ps1d.Poisson_solver_1d(Mesh, chrg_den, bc, invA) # pe contains [pot, efld]
# move particles
Eon_pv = move.move_ptcl(Mesh, Eon, Eon_pv, pe[1], dt)
Arp_pv = move.move_ptcl(Mesh, Arp, Arp_pv, pe[1], dt)

num_iter = 15001 # number of iterations
nout_iter = 1000
for i in range(num_iter):
    # assign charge densities to grid nodes, unit in UNIT_CHARGE
    Eon_den = move.den_asgmt(Eon_pv[0], Mesh)
    Arp_den = move.den_asgmt(Arp_pv[0], Mesh)
    # add up all charge densities
    chrg_den = (Eon.charge*Eon_den + 
                Arp.charge*Arp_den)*den_per_ptcl # unit in UNIT_CHARGE
    chrg_den = savgol_filter(chrg_den, 11, 3) # window size 11, polynomial order 3
    
    # update potential according to assigned charges to nodes
    pe = ps1d.Poisson_solver_1d(Mesh, chrg_den, bc, invA) # pe contains [pot, efld]
    # move particles
    Eon_pv = move.move_ptcl(Mesh, Eon, Eon_pv, pe[1], dt)
    Arp_pv = move.move_ptcl(Mesh, Arp, Arp_pv, pe[1], dt)
    num_ptcl = len(Eon_pv[0])
    if i % nout_iter == 0:
        print("iter = %d" % i, 
              "- time %s -" % str(timedelta(seconds=(int(time.time() - t0)))))
        # plot animation
        out.plot_diag(Mesh, Eon_pv, Arp_pv, chrg_den, pe, i, num_ptcl)

print("-total time %s -" % str(timedelta(seconds=(int(time.time() - t0)))))
print("-- total plasma time %d ns --"    % (dt*num_iter/1e-9))