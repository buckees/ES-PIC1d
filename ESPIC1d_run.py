# -*- coding: utf-8 -*-
# This ES-PIC_run.py runs the PIC simulation
# All in SI unit

from os import system as sys
sys('rm Figures/*.png')

from datetime import timedelta 
import time
t0 = time.time()

# modules from *.py in the same folder
import Constants as cst
from ESPIC1d_mesh import Mesh
import ESPIC1d_init as init
import ESPIC1d_leapfrog as frog
import ESPIC1d_out as out
import ESPIC1d_rct as rct
from Species import Eon, Arp, Ar
import Poisson_solver_1d as ps1d

# python built-in modules
import numpy as np
import math
from scipy.signal import savgol_filter

"""
Geometry and Mesh are created in ESPIC1d_mesh.py
"""

# Initial Conditions
press_mT = 100.0 # pressure, in mTorr
pressure = press_mT/1.0e3*cst.TORR2PA # in Pa, 1 Torr = 133.322 Pa
Ar_den_init = pressure/(cst.KB*Ar.tmpt*cst.EV2K) # in m-3, N/V = P/(kb*T) ideal gas law
Eon_temp = [Eon.tmpt] # initial temperature for Eon, in eV
Arp_temp = [Arp.tmpt] # initial temperature for Ar+, in eV
bc = [0.0, 0.0] # left and right boundary conditions, in Volt 
Eon_clct = [[], []] # collect Eon particles bombarding left and right surface
Arp_clct = [[], []] # collect Ar+ particles bombarding left and right surface
Eon_den_init = 1.0e15 # in m-3, initial electron density
den_limit = 1.0e11 # in m-3, lower limit of denisty, avoid 0 density 

# Operation Parameters
num_ptcl = 10000 # number of particles, should be >> ncellx to reduce noise
dt = 1.0e-12 # in sec

den_per_ptcl = Eon_den_init/num_ptcl # density contained in one particle

# initialize the position and velcotiy in a dataframe
Eon_pv = init.init_data(num_ptcl, Eon) # a list contain [posn, vels]
Arp_pv = init.init_data(num_ptcl, Arp)
Eon_pv[0] = Eon_pv[0]*Mesh.width
Arp_pv[0] = Eon_pv[0].copy()

# assign charge densities to grid nodes, unit in UNIT_CHARGE
Eon_den = frog.den_asgmt(Eon_pv[0], Mesh)
Arp_den = frog.den_asgmt(Arp_pv[0], Mesh)
# add up all charge densities
chrg_den = (Eon.charge*Eon_den + 
            Arp.charge*Arp_den)*den_per_ptcl # unit in UNIT_CHARGE
chrg_den = savgol_filter(chrg_den, 11, 3) # window size 11, polynomial order 3

# update potential according to assigned charges to nodes
invA = ps1d.calc_invA(Mesh)
pe = ps1d.Poisson_solver_1d(Mesh, chrg_den, bc, invA) # pe contains [pot, efld]
# using leapfrog algrithm to update position
Eon_pv, v_clct = frog.move_leapfrog1(Mesh, Eon, Eon_pv, pe[1], dt)
Eon_clct[0] += v_clct[0]; Eon_clct[1] += v_clct[1];
Arp_pv, v_clct = frog.move_leapfrog1(Mesh, Arp, Arp_pv, pe[1], dt)
Arp_clct[0] += v_clct[0]; Arp_clct[1] += v_clct[1];

# update charge density at t1
Eon_den = frog.den_asgmt(Eon_pv[0], Mesh)
Arp_den = frog.den_asgmt(Arp_pv[0], Mesh)
chrg_den = (Eon.charge*Eon_den + 
            Arp.charge*Arp_den)*den_per_ptcl # unit in UNIT_CHARGE
chrg_den = savgol_filter(chrg_den, 11, 3) # window size 11, polynomial order 3
# update potential and E-field at t1
pe = ps1d.Poisson_solver_1d(Mesh, chrg_den, bc, invA) # pe contains [pot, efld]
# update only velocity at t1
Eon_pv = frog.move_leapfrog2(Mesh, Eon, Eon_pv, pe[1], dt)
Arp_pv = frog.move_leapfrog2(Mesh, Arp, Arp_pv, pe[1], dt)

ptcl_rec = [[], [], []] # [0] = num_ptcl; [1] = ptcl_rm; [2] = ptcl_add
ergs_mean, ergs_max = [[], []], [[], []] # [0] = Eon mean erg; [1] = Arp mean erg;
num_iter = 5000001 # number of iterations
nout_iter = 3000
for i in range(num_iter):
    # using leapfrog algrithm to update position
    Eon_pv, v_clct = frog.move_leapfrog1(Mesh, Eon, Eon_pv, pe[1], dt)
    Eon_clct[0] += v_clct[0]; Eon_clct[1] += v_clct[1];
    ptcl_rec[1].append(len(v_clct[0]) + len(v_clct[1]))
    Arp_pv, v_clct = frog.move_leapfrog1(Mesh, Arp, Arp_pv, pe[1], dt)
    Arp_clct[0] += v_clct[0]; Arp_clct[1] += v_clct[1];
    
    # update charge density at t1
    Eon_den = frog.den_asgmt(Eon_pv[0], Mesh)
    Arp_den = frog.den_asgmt(Arp_pv[0], Mesh)
    chrg_den = (Eon.charge*Eon_den + 
                Arp.charge*Arp_den)*den_per_ptcl # unit in UNIT_CHARGE
    chrg_den = savgol_filter(chrg_den, 11, 3) # window size 11, polynomial order 3
    # update potential and E-field at t1
    bc[0] = 10.0*math.sin(i/10000*2.0*math.pi)
    pe = ps1d.Poisson_solver_1d(Mesh, chrg_den, bc, invA) # pe contains [pot, efld]
    # update only velocity at t1
    Eon_pv = frog.move_leapfrog2(Mesh, Eon, Eon_pv, pe[1], dt)
    Arp_pv = frog.move_leapfrog2(Mesh, Arp, Arp_pv, pe[1], dt)
    # calc eon impact ionization
    # convert vels to ergs
    Eon_ergs = np.power(Eon_pv[1]*cst.VEL2EV,2) 
    Arp_ergs = np.power(Arp_pv[1]*cst.VEL2EV,2) 
    Eon_pv, Arp_pv, ptcl_rec[2] = rct.ioniz(Eon_pv, Arp_pv, Eon_ergs, dt, 
                                         ptcl_rec[2])
    
    ptcl_rec[0].append(len(Eon_pv[0]))
    ergs_mean[0].append(np.mean(Eon_ergs))
    ergs_mean[1].append(np.mean(Arp_ergs))
    ergs_max[0].append(np.amax(Eon_ergs))
    ergs_max[1].append(np.amax(Arp_ergs))
    if i % nout_iter == 0:
        print("iter = %d" % i, 
              "- time %s -" % str(timedelta(seconds=(int(time.time() - t0)))))
        # plot animation
        out.plot_diag(Mesh, Eon_pv, Arp_pv, Eon_clct, Arp_clct, 
                      chrg_den, pe, i, ergs_mean, ergs_max, ptcl_rec)
    if ptcl_rec[0][-1] < num_ptcl*0.01: break

print("-total time %s -" % str(timedelta(seconds=(int(time.time() - t0)))))
print("-- total plasma time %d ns --"    % (dt*num_iter/1e-9))