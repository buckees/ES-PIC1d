# -*- coding: utf-8 -*-
# This ES-PIC_run.py runs the PIC simulation
# All in SI unit

from os import system as sys
sys('mkdir Figures')
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
from Species import Eon, Hp, H
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
H_den_init = pressure/(cst.KB*H.tmpt*cst.EV2K) # in m-3, N/V = P/(kb*T) ideal gas law
Eon_temp = [Eon.tmpt] # initial temperature for Eon, in eV
Hp_temp = [Hp.tmpt] # initial temperature for H+, in eV
bc = [0.0, 0.0] # left and right boundary conditions, in Volt 
Eon_clct = [[], []] # collect Eon particles bombarding left and right surface
Hp_clct = [[], []] # collect H+ particles bombarding left and right surface
Eon_den_init = 1.0e15 # in m-3, initial electron density
den_limit = 1.0e11 # in m-3, lower limit of denisty, avoid 0 density 
# 1 mTorr, mfp ~ 10 cm, col_freq ~ 5e6 s-1 
col_freq = 5.0e8 # in s-1
powE = 100 # in W, power is deposited only to Eon

# Operation Parameters
num_ptcl = 20000 # number of particles, should be >> ncellx to reduce noise
freq = 10.0e6 # frequncy, in Hz, period = 100 ns
perd = 1.0/freq # period = 100 ns
dt = perd/2000.0 # in sec

den_per_ptcl = Eon_den_init/num_ptcl # density contained in one particle

# initialize the position and velcotiy in a dataframe
Eon_pv = init.init_data(num_ptcl, Eon) # a list contain [posn, vels]
Hp_pv = init.init_data(num_ptcl, Hp)
Eon_pv[0] = Eon_pv[0]*Mesh.width
Hp_pv[0] = Eon_pv[0].copy()

# assign charge densities to grid nodes, unit in UNIT_CHARGE
Eon_den = frog.den_asgmt(Eon_pv[0], Mesh)
Hp_den = frog.den_asgmt(Hp_pv[0], Mesh)
# add up all charge densities
chrg_den = (Eon.charge*Eon_den + 
            Hp.charge*Hp_den)*den_per_ptcl # unit in UNIT_CHARGE
chrg_den = savgol_filter(chrg_den, 11, 3) # window size 11, polynomial order 3

# update potential according to assigned charges to nodes
invA = ps1d.calc_invA(Mesh)
pe = ps1d.Poisson_solver_1d(Mesh, chrg_den, bc, invA) # pe contains [pot, efld]
# using leapfrog algrithm to update position
Eon_pv, v_clct = frog.move_leapfrog1(Mesh, Eon, Eon_pv, pe[1], dt)
Eon_clct[0] += v_clct[0]; Eon_clct[1] += v_clct[1];
Hp_pv, v_clct = frog.move_leapfrog1(Mesh, Hp, Hp_pv, pe[1], dt)
Hp_clct[0] += v_clct[0]; Hp_clct[1] += v_clct[1];

# update charge density at t1
Eon_den = frog.den_asgmt(Eon_pv[0], Mesh)
Hp_den = frog.den_asgmt(Hp_pv[0], Mesh)
chrg_den = (Eon.charge*Eon_den + 
            Hp.charge*Hp_den)*den_per_ptcl # unit in UNIT_CHARGE
chrg_den = savgol_filter(chrg_den, 11, 3) # window size 11, polynomial order 3
# update potential and E-field at t1
pe = ps1d.Poisson_solver_1d(Mesh, chrg_den, bc, invA) # pe contains [pot, efld]
# update only velocity at t1
Eon_pv = frog.move_leapfrog2(Mesh, Eon, Eon_pv, pe[1], dt)
Hp_pv = frog.move_leapfrog2(Mesh, Hp, Hp_pv, pe[1], dt)

# Eon: [0] = num_ptcl ; [1] = ptcl_rm; [2] = ptcl_add
# Arp: [3] = num_ptcl ; [4] = ptcl_rm; [5] = ptcl_add
ptcl_rec = [[] for j in range(6)]
ergs_mean, ergs_max = [[], []], [[], []] # [0] = Eon mean erg; [1] = Arp mean erg;
num_iter = 5000001 # number of iterations, total 50 us
nout_iter = int(perd/dt) # output every xxx ns
for i in range(num_iter):
    # using leapfrog algrithm to update position
    Eon_pv, v_clct = frog.move_leapfrog1(Mesh, Eon, Eon_pv, pe[1], dt)
    Eon_clct[0] += v_clct[0]; Eon_clct[1] += v_clct[1];
    ptcl_rec[1].append(len(v_clct[0]) + len(v_clct[1]))
    Hp_pv, v_clct = frog.move_leapfrog1(Mesh, Hp, Hp_pv, pe[1], dt)
    Hp_clct[0] += v_clct[0]; Hp_clct[1] += v_clct[1];
    ptcl_rec[4].append(len(v_clct[0]) + len(v_clct[1]))
    
    # 2nd Eon emission, 10% from Ion currrent
#    Eon_pv[0] = np.append(Eon_pv[0], 
#                          np.ones(len(v_clct[0]))*Mesh.dx*0.1)
#    Eon_pv[1] = np.append(Eon_pv[1], 
#                          np.abs(v_clct[0]))
#    ptcl_rec[3].append(len(v_clct[0]))
    
    # update charge density at t1
    Eon_den = frog.den_asgmt(Eon_pv[0], Mesh)
    Hp_den = frog.den_asgmt(Hp_pv[0], Mesh)
    chrg_den = (Eon.charge*Eon_den + 
                Hp.charge*Hp_den)*den_per_ptcl # unit in UNIT_CHARGE
    chrg_den = savgol_filter(chrg_den, 11, 3) # window size 11, polynomial order 3
    # update potential and E-field at t1
    pe = ps1d.Poisson_solver_1d(Mesh, chrg_den, bc, invA) # pe contains [pot, efld]
    # update only velocity at t1
    Eon_pv = frog.move_leapfrog2(Mesh, Eon, Eon_pv, pe[1], dt)
    # elastic collision only applies to Eon
    Eon_pv[1] = frog.move_coll(Eon_pv[1], col_freq*dt)
    Hp_pv = frog.move_leapfrog2(Mesh, Hp, Hp_pv, pe[1], dt)

    # calc eon impact ionization
    Eon_pv, Hp_pv, ptcl_rec[2] = rct.ioniz(Eon_pv, Hp_pv, dt, 
                                         ptcl_rec[2])
    ptcl_rec[5] = ptcl_rec[2].copy()

    # convert vels to ergs
    Eon_ergs = np.power(Eon_pv[1]*cst.VEL2EV,2) 
    Hp_ergs = np.power(Hp_pv[1]*cst.VEL2EV,2) 

    # update Eon_pv due to power deposition
    Eon_pv[1] = rct.pow_dep(Eon_pv[1], Eon_ergs, powE/den_per_ptcl, dt)
    
    ptcl_rec[0].append(len(Eon_pv[0]))
    ptcl_rec[3].append(len(Hp_pv[0]))
    ergs_mean[0].append(np.mean(Eon_ergs))
    ergs_mean[1].append(np.mean(Hp_ergs))
    ergs_max[0].append(np.amax(Eon_ergs))
    ergs_max[1].append(np.amax(Hp_ergs))
    
    if (i+1) % nout_iter == 0:
        t = divmod(int(i*dt/1.0e-9), 1000)
        print("iter = %d - " % i,
              "plasma time = %d us %d ns - " % t, 
              " time %s -" % str(timedelta(seconds=(int(time.time() - t0)))),
              "# of ptcl = %d" % ptcl_rec[0][-1])
        # plot animation
        out.plot_diag(Mesh, Eon_pv, Hp_pv, Eon_clct, Hp_clct, 
                      chrg_den, pe, i, ergs_mean, ergs_max, ptcl_rec,
                      nout_iter)
    if (ptcl_rec[0][-1] < num_ptcl*0.10 or 
        ptcl_rec[0][-1] > num_ptcl*10.0):
        break

