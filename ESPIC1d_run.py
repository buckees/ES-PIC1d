# -*- coding: utf-8 -*-
# This ES-PIC_run.py runs the PIC simulation
# All in SI unit

import Constants as cst
import numpy as np
import ESPIC1d_init as init
import ESPIC1d_move as move
import Particle as ptcl
import matplotlib.pyplot as plt
import pandas as pd

# Geometry and Mesh
width = 0.01 # lenght of domain, in meter
ncellx = 100 # number of cells in x direction

# Operation Conditions
press_mT = 100.0 # pressure, in mTorr
pressure = press_mT/1.0e3*cst.TORR2PA # in Pa, 1 Torr = 133.322 Pa
nAr_init = pressure/(cst.KB*pctl.Ar.temp*cst.EV2K) # in m-3, N/V = P/(kb*T) ideal gas law
pctl.Eon.temp = 2.0 # in eV
pctl.Arp.temp = 0.1 # in eV
pot_l = 100.0 # in V, left boundary potential
pot_r = 0.0 # in V, right boundary potential
nEon_init = 1.0e14 # in m-3, initial electron density


den_limit = 1.0e11 # in m-3, lower limit of denisty, avoid 0 density 

# Model Parameters
num_ptcl = 10000 # number of particles, should be >> ncellx to reduce noise
time_step = 1.0e-12 # in sec
num_iter = 200 # number of iterations

den_per_ptcl = nEon_init/num_ptcl # density contained in one particle

# initialize the position and velcotiy in a dataframe
data_Eon = init.init_data(num_ptcl,ptcl.Eon.temp,ptcl.Eon.mass,width)
data_Arp = init.init_data(num_ptcl,ptcl.Arp.temp,ptcl.Arp.mass,width)

# create mesh
gridx, cell_center, dx = init.init_mesh(ncellx,width)

# initilize potential on cell boudaries
pot, efld = init.init_pot(ncellx,dx,pot_l,pot_r)

# assign charge densities to grid nodes, unit in UNIT_CHARGE
den_Eon = move.den_asgmt(data_Eon['position'],gridx,dx)*den_per_ptcl
den_Arp = move.den_asgmt(data_Arp['position'],gridx,dx)*den_per_ptcl
den_chrg = ptcl.Eon.charge*den_Eon + ptcl.Arp.charge*den_Arp



# diagnostic plot
fig, (ax0,ax1) = plt.subplots(2,1, figsize=(8,4))
#data_Eon.plot(x='position',y='index',ax=ax0,kind='scatter',xlim=(0.0,width),
#              c='blue')
#data_Arp.plot(x='position',y='index',ax=ax0,kind='scatter',xlim=(0.0,width),
#              c='red')
data_Eon.plot.hist(y='position',bins=ncellx,ax=ax0)
data_Arp.plot.hist(y='position',bins=ncellx,ax=ax0)
ax1.plot(gridx,den_Eon,'b-')
ax1.plot(gridx,den_Arp,'r-')
ax1.plot(gridx,den_chrg,'k-')
plt.show()