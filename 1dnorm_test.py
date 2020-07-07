# -*- coding: utf-8 -*-
# This .py test the 1d normal velocity distribution in ESPIC1d_init.py
# All in SI unit

import Constants as cst
import numpy as np
import ESPIC1d_init as init
import Particle as pctl
import matplotlib.pyplot as plt

# Geometry and Mesh
length = 0.01 # lenght of domain, unit in m


# Operation Conditions
press_mT = 100.0 # pressure, unit in mTorr
pressure = press_mT/1e3*133.3 # unit in Pa, 1 Torr = 133.322 Pa
eon_temp = 2.0 # unit in eV
pctl.Arp.temp = 100.1 # unit in eV

# Model Parameters
num_cell = 10 # number of cells
num_ptcl = 1000 # number of particles
time_step = 1e-12 # unit in s
num_iter = 200 # number of iterations

vels = init.norm_distribution(pctl.N2.temp,pctl.N2.mass,num_ptcl)

fig, (ax0,ax1) = plt.subplots(2,1)
ax0.hist(vels, density=True, histtype='stepfilled', alpha=0.5)
ax0.legend(loc='best', frameon=False)

ax1.plot(vels,range(len(vels)),'o')
plt.show()