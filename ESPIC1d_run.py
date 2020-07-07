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
width = 0.01 # lenght of domain, unit in m


# Operation Conditions
press_mT = 100.0 # pressure, unit in mTorr
pressure = press_mT/1.0e3*cst.TORR2PA # unit in Pa, 1 Torr = 133.322 Pa
pctl.Eon.temp = 2.0 # unit in eV
pctl.Arp.temp = 0.1 # unit in eV

# Model Parameters
num_cell = 10 # number of cells
num_ptcl = 10 # number of particles
time_step = 1e-12 # unit in s
num_iter = 200 # number of iterations

# initialize the position and velcotiy in a dataframe
data_Eon = init.init_data(num_ptcl,pctl.Eon.temp,pctl.Eon.mass,width)
data_Arp = init.init_data(num_ptcl,pctl.Arp.temp,pctl.Arp.mass,width)
