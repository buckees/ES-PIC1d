# -*- coding: utf-8 -*-
# This ES-PIC1d_move.py moves the particles
#

import Constants as cst
import Particle as ptcl
from scipy.stats import norm
import numpy as np
import pandas as pd
import math

def den_asgmt(posn, mesh):
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
    density = np.full((mesh.ncellx+1,),1.0e-5)
    for i, p in enumerate(posn):
        frac, whole = math.modf(p/mesh.dx)
        whole = int(whole)
        density[whole] = 1-frac
        density[whole+1] = frac
    return density

def move_ptcl(particle,posn,vels,efld,dt,mesh):
    """
    Update position and velocity in dataframe at t1 = t0 + dt
    :param efld: E-field within each cell, in V/m
    :param posn: particle positions at t0
    :param vels: particle velocities at t0
    :param particle: type of particle, class particle
    :param dt: time step, in sec
    """
    mass = particle.mass*cst.AMU # mass in kg
    chrg = particle.charge*cst.UNIT_CHARGE # charge in Coloumb
    posn_new, vels_new = [], []
    for p, v in zip(posn,vels):
        frac, whole = math.modf(p/mesh.dx)
        whole = int(whole)
        ef = efld[whole] # E-filed on i_th particle
        accel = ef*chrg/mass # acceleration in m/s2
        p += v*dt; posn_new.append(p)
        v += accel*dt; vels_new.append(v)
    posn_new, vels_new = check_bdry(posn_new,vels_new,mesh)
    return posn_new, vels_new

def check_bdry(posn,vels,mesh):
    """
    Remove particles which go beyond the domain, (p < 0.0  or p > width)
    :param posn: particle positions at t0
    :param vels: particle velocities at t0
    :param width: domain width
    """
    index_remove = []
    for i, p in enumerate(posn):
        if not (0.0 < p < mesh.width):
            index_remove.append(i)
    return np.delete(posn,index_remove), np.delete(vels,index_remove)
    