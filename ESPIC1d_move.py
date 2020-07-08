# -*- coding: utf-8 -*-
# This ES-PIC1d_move.py moves the particles
#

import Constants as cst
import Particle as ptcl
from scipy.stats import norm
import numpy as np
import pandas as pd
import math

def den_asgmt(posn, gridx, dx):
    """
    Return the density distribution on grids
    Assign the density to each grid/node based on 
    the proximity of the particle to that node. e.g.,
    node 1 -------------------- Eon -------- node 2
    1/3 Eon density is assigned to node 1, while 2/3 is assigned to node 2
    :param posn: particle positions
    :param gridx: grid/cell_boundary coordinate in x direction
    :param den_per_ptcl: density contained in a particle
    """
    density = np.full((len(gridx),),1.0e-5)
    for i, p in enumerate(posn):
        p = p/dx
        frac, whole = math.modf(p)
        whole = int(whole)
        density[whole] = 1-frac
        density[whole+1] = frac
    return density

def move_ptcl(efld,data,particle,dt):
    """
    Update position and velocity in dataframe at t1 = t0 + dt
    :param efld: E-field within each cell, in V/m
    :param data: particle position and velocity at t0
    :param particle: type of particle, class particle
    :param dt: time step, in sec
    """
    mass = particle.mass*cst.AMU # mass in kg
    chrg = particle.charge*cst.UNIT_CHARGE # charge in Coloumb
    for i in range(len(data)):
        frac, whole = math.modf(data.position.iloc[i])
        whole = int(whole)
        ef = efld[whole] # E-filed on i_th particle
        accel = ef*chrg/mass # acceleration in m/s2
        data.velocity.iloc[i] += accel*dt
        data.position.iloc[i] += data.velocity.iloc[i]*dt
    return data