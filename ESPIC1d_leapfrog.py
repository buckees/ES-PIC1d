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

def move_leapfrog1(mesh, sp, pv, efld, dt):
    """
    Update position and velocity using leapfrog
    p(t1) = p(t0) + v(t0)*dt + 0.5*a(t0)*dt*dt 
    v(t1) = v(t0) + 0.5*(a(t0) + a(t1))*dt
    this func calcs part I: 
    p(t1) = p(t0) + v(t0)*dt + 0.5*a(t0)*dt*dt 
    v(t1) = v(t0) + 0.5*(a(t0)
    :param mesh: all mesh info
    :param sp: type of particle, class particle
    :param pv: sp positions and velocities at t0
    :param efld: E-field within each cell, in V/m
    :param dt: time step, in sec
    """
    accel = efld*sp.charge/sp.mass*cst.UNIT_CHARGE/cst.AMU
    frac, whole = np.modf(pv[0]/mesh.dx)
    whole = whole.astype(int)
    pv[0] += pv[1]*dt + 0.5*accel[whole]*dt*dt
    pv[1] += 0.5*accel[whole]*dt
    # update posn, p(t1) = p(t0) + v(t0)*dt + 0.5*a(t0)*dt*dt 
    pv, v_clct = chk_bdry(pv[0], pv[1], mesh.width)
    return pv, v_clct

@njit
def chk_bdry(posn, vels, width):
    posn_keep, vels_keep = [], []
    v_left, v_right = [], []
    for p, v in zip(posn, vels):
        if (width*0.001 < p < width*0.999):
            posn_keep.append(p); vels_keep.append(v)
        else:
            if p <= width*0.001: 
                v_left.append(v)
            else:
                v_right.append(v)
    posn_keep = np.asarray(posn_keep); vels_keep = np.asarray(vels_keep)
    return [posn_keep, vels_keep], [v_left, v_right]

def move_leapfrog2(mesh, sp, pv, efld, dt):
    """
    Update position and velocity using leapfrog
    p(t1) = p(t0) + v(t0)*dt + 0.5*a(t0)*dt*dt 
    v(t1) = v(t0) + 0.5*(a(t0) + a(t1))*dt
    this func calcs part II: 
    v(t1) += 0.5*a(t1)
    :param mesh: all mesh info
    :param sp: type of particle, class particle
    :param pv: sp positions and velocities at t0
    :param efld: E-field within each cell, in V/m
    :param dt: time step, in sec
    """
    accel = efld*sp.charge/sp.mass*cst.UNIT_CHARGE/cst.AMU
    frac, whole = np.modf(pv[0]/mesh.dx)
    whole = whole.astype(int)
    # update posn, p(t1) = p(t0) + v(t0)*dt + 0.5*a(t0)*dt*dt
    pv[1] += 0.5*accel[whole]*dt
    return pv

def move_coll(vels, prob):
    rand = np.random.uniform(0.0, 1.0, size=(len(vels), ))
    vels = vels*np.sign(rand - prob)
    return vels