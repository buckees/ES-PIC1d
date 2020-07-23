# -*- coding: utf-8 -*-
# This ES-PIC1d_init.py initialize the PIC simulation

import Constants as cst
from scipy.stats import norm
import numpy as np

def thermal_velocity(charge, tmpt, mass):
    """
    Most probable speed in the velocity distribution = sqrt(2eT/m)
    :param charge: particle charge in C
    :param tmpt: particle temp in eV. 1 Joule = 1 Volt * 1 Coulomb 
    :param mass: particle mass in kg
    :return: thermal velocity in m/s
    """
    return np.sqrt(2*abs(charge)*tmpt/mass)

def norm_distribution(num_ptcl, tmpt, mass):
    """
    Returns an array of velocities sampled from the Normal distribution
    1d maxwell velocity distribution is a normal distribution
    f(v) = sqrt(m/2/pi/Kb/T)*exp(-v**2*m/2/kb/T)
    norm distribution f(x) = sqrt(1/2/pi)*exp(-x**2/2)
    v = x*m/kb/T, so scale = 1(m/kb/T)
    :param tmpt: particle tmpt, in eV
    :param mass: particle mass, in AMU
    :param num_velocities: number of velocities to generate
    :returns: 1d array array of velocities
    """
    scale = np.sqrt((tmpt/mass)*(cst.KB*cst.EV2K/cst.AMU))
    return norm.rvs(loc=0.0, scale=scale, size=num_ptcl)
    
def debye_length(eon_temp, ne):
    """
    Calculate plasma debye length
    :param eps0: relative permittivity, no unit
    :param eon_temp: electron temp, in eV
    :param ne: electron # density, in m-3
    :return: debye length, in meter
    """
    return np.sqrt(cst.EPS0*eon_temp/ne/cst.UNIT_CHARGE)

def init_posn(num_ptcl=10):
    """
    Return the random initial positions
    :param num_ptcl: number of particles
    :param ptcl_type: the type of particles, must be pre-defined in Particle.py
    """
    return np.random.uniform(low=0.001, high=0.999, size=(num_ptcl,))

def init_data(num_ptcl, sp):
    """
    return a DataFrame containing the particle positon and velocity
    :param num_ptcl: number of particles
    :param tmpt: particle tmpt, in eV
    :param mass: particle mass, in AMU
    :param width: the width of domain/geometry
    """
    return [init_posn(num_ptcl), norm_distribution(num_ptcl, sp.tmpt, sp.mass)]

def init_mesh(ncellx=10,width=0.01):
    """
    return grid/cell_boundry coordiates, cell center coordinates and cell size
    :param ncellx: number of cells in x direction
    :param width: the width of domain/geometry
    """
    gridx, dx = np.linspace(0.0,width,ncellx+1,retstep=True)
    cell_center = [np.mean(gridx[i:i+2]) for i in range(ncellx)]
    return gridx, cell_center, dx

def sinusoidal(amplitude, period, phase, t):
    return np.sin(2*np.pi*(t/period) + phase) * amplitude