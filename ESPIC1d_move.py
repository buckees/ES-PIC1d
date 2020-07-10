# -*- coding: utf-8 -*-
# This ES-PIC1d_move.py moves the particles
#

import Constants as cst
import numpy as np

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
    den = np.full((mesh.ncellx+1,),1.0e-5)
    frac, whole = np.modf(posn/mesh.dx)
    whole = whole.astype(int)
    for w, f in zip(whole, frac):
        den[w] += 1 - f
        den[w+1] += f
    return den

def move_ptcl(mesh, sp, pv, efld, dt):
    """
    Update position and velocity in dataframe at t1 = t0 + dt
    :param mesh: all mesh info
    :param sp: type of particle, class particle
    :param pv: sp positions and velocities at t0
    :param efld: E-field within each cell, in V/m
    :param dt: time step, in sec
    """
    mass = sp.mass*cst.AMU # mass in kg
    chrg = sp.charge*cst.UNIT_CHARGE # charge in Coloumb
    posn_new, vels_new = [], []
    frac, whole = np.modf(pv[0]/mesh.dx)
    whole = whole.astype(int)
    for p, v, w in zip(pv[0], pv[1], whole):
        ef = efld[w] # E-filed on i_th particle
        accel = ef*chrg/mass # acceleration in m/s2
        p += v*dt; posn_new.append(p)
        v += accel*dt; vels_new.append(v)
    posn_new, vels_new = check_bdry(posn_new,vels_new,mesh)
    return [posn_new, vels_new]

def check_bdry(posn, vels, mesh):
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
    