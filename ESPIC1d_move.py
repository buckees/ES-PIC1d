# -*- coding: utf-8 -*-
# This ES-PIC1d_move.py moves the particles
#

import Constants as cst
import numpy as np
from numba import njit

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
    frac, whole = np.modf(posn/mesh.dx)
    whole = whole.astype(int)
    return make_asgmt(whole, frac, mesh.ncellx)

@njit
def make_asgmt(whole, frac, ncellx):
    den = np.full((ncellx+1,),1.0e-5)
    for w, f in zip(whole, frac):
        den[w] += 1.0 - f
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
    accel = efld*sp.charge/sp.mass*cst.UNIT_CHARGE/cst.AMU
    frac, whole = np.modf(pv[0]/mesh.dx)
    whole = whole.astype(int)
    pv, v_clct = make_move(pv[0], pv[1], whole, accel, dt, mesh.width)
    return pv, v_clct

@njit
def make_move(posn, vels, whole, accel, dt, width):
    posn_keep, vels_keep = [], []
    v_left, v_right = [], []
    for p, v, w in zip(posn, vels, whole):
        p += v*dt; 
        if (0.001 < p < width*0.999):
            v += accel[w]*dt; 
            posn_keep.append(p); vels_keep.append(v)
        else:
            if p <= 0.001: 
                v_left.append(v)
            else:
                v_right.append(v)
    posn_keep = np.asarray(posn_keep); vels_keep = np.asarray(vels_keep)
#    v_left = np.asarray(v_left); v_right = np.asarray(v_right)
    return [posn_keep, vels_keep], [v_left, v_right]
    