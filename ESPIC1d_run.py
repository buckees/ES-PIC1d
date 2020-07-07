# -*- coding: utf-8 -*-
# This ES-PIC_run.py runs the PIC simulation
# All in SI unit

import Constants as cst
import numpy as np
import ESPIC1d_init as init
import Particle as pctl
import matplotlib.pyplot as plt
import pandas as pd

# Geometry and Mesh
width = 0.01 # lenght of domain, in meter
ncellx = 10 # number of cells in x direction

# Operation Conditions
press_mT = 100.0 # pressure, in mTorr
pressure = press_mT/1.0e3*cst.TORR2PA # in Pa, 1 Torr = 133.322 Pa
nAr_init = pressure/(cst.KB*pctl.Ar.temp*cst.EV2K) # in m-3, N/V = P/(kb*T) ideal gas law
pctl.Eon.temp = 2.0 # in eV
pctl.Arp.temp = 0.1 # in eV
pot_l = 100.0 # in V, left boundary potential
pot_r = 0.0 # in V, right boundary potential
nEon_init = 1e14 # in m-3, initial electron density

den_limit = 1e11 # in m-3, lower limit of denisty, avoid 0 density 

# Model Parameters
num_cell = 10 # number of cells
num_ptcl = 10 # number of particles
time_step = 1e-12 # in sec
num_iter = 200 # number of iterations

# initialize the position and velcotiy in a dataframe
data_Eon = init.init_data(num_ptcl,pctl.Eon.temp,pctl.Eon.mass,width)
data_Arp = init.init_data(num_ptcl,pctl.Arp.temp,pctl.Arp.mass,width)

# create mesh
gridx, cell_center, dx = init.init_mesh(ncellx,width)

# initilize potential on cell boudaries
pot, efld = init.init_pot(ncellx,dx,pot_l,pot_r)
