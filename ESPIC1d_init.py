# -*- coding: utf-8 -*-
# This ES-PIC1d_init.py initialize the PIC simulation

import Constants as cst
from scipy.stats import norm
import numpy as np
#import pandas as pd

class Mesh(object):
    """Stores all particles"""
    def __init__(self, ncellx, width, dx, gridx, cell_cnt):
        self.ncellx = ncellx # number of cells in x direction
        self.width = width # domain size in 1d, in meter
        self.dx = dx # cell size, in meter
        self.gridx = gridx # an array of grids in x direction, in meter
        self.cell_cnt = cell_cnt # an array of cell center in x direction, in meter

# default Mesh
gridx = np.linspace(0.0,1.0,101)
cell_cnt = (gridx[0:-1] + gridx[1:])/2.0
Mesh = Mesh(100,1.0,0.01,gridx,cell_cnt)

def thermal_velocity(charge, temperature, mass):
    """
    Most probable speed in the velocity distribution = sqrt(2eT/m)
    :param charge: particle charge in C
    :param temperature: particle temp in eV. 1 Joule = 1 Volt * 1 Coulomb 
    :param mass: particle mass in kg
    :return: thermal velocity in m/s
    """
    return np.sqrt(2*abs(charge)*temperature/mass)

def norm_distribution(num_ptcl,temperature,mass):
    """
    Returns an array of velocities sampled from the Normal distribution
    1d maxwell velocity distribution is a normal distribution
    f(v) = sqrt(m/2/pi/Kb/T)*exp(-v**2*m/2/kb/T)
    norm distribution f(x) = sqrt(1/2/pi)*exp(-x**2/2)
    v = x*m/kb/T, so scale = 1(m/kb/T)
    :param temperature: particle temperature, in eV
    :param mass: particle mass, in AMU
    :param num_velocities: number of velocities to generate
    :returns: 1d array array of velocities
    """
    mass = mass*cst.AMU
    tempK = temperature*cst.EV2K
    scale = np.sqrt(cst.KB*tempK/mass)
    vels = norm.rvs(loc=0.0, scale=scale, size=num_ptcl)
    
    return vels
    
def debye_length(eps0, eon_temp, ne):
    """
    Calculate plasma debye length
    :param eps0: relative permittivity
    :param eon_temp: electron temp in eV
    :param ne: electron density in particles/volume
    :return: debye length
    """
    return np.sqrt(eps0*cst.EPS0*eon_temp/ne/cst.UNIT_CHARGE)

def init_posn(num_ptcl=10):
    """
    Return the random initial positions
    :param num_ptcl: number of particles
    :param ptcl_type: the type of particles, must be pre-defined in Particle.py
    """
    return np.random.uniform(low=0.001, high=0.999, size=(num_ptcl,))

def init_data(num_ptcl,particle):
    """
    return a DataFrame containing the particle positon and velocity
    :param num_ptcl: number of particles
    :param temperature: particle temperature, in eV
    :param mass: particle mass
    :param width: the width of domain/geometry
    """
    posn = init_posn(num_ptcl) # initial position
    vels = norm_distribution(num_ptcl,particle.temp,particle.mass)
#    data = pd.DataFrame(np.vstack((posn, vels)).T,
#                        columns=['position','velocity'])
    return posn, vels

def init_mesh(ncellx=10,width=0.01):
    """
    return grid/cell_boundry coordiates, cell center coordinates and cell size
    :param ncellx: number of cells in x direction
    :param width: the width of domain/geometry
    """
    gridx, dx = np.linspace(0.0,width,ncellx+1,retstep=True)
    cell_center = [np.mean(gridx[i:i+2]) for i in range(ncellx)]
    return gridx, cell_center, dx

def init_pot(mesh,bc):
    """
    return initial potential and E-field
    :param ncellx: number of cells in x direction
    :param dx: the cell size in x direction
    :param pot_l: the potential at left boundary
    :param pot_r: the potential at right boundary
    """
    pot = np.linspace(bc[0],bc[-1],mesh.ncellx+1)
    efld = -(pot[1:] - pot[0:-1])/mesh.dx
    return pot, efld


def sinusoidal(amplitude, period, phase, t):
    return np.sin(2*np.pi*(t/period) + phase) * amplitude