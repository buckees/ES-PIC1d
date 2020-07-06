# -*- coding: utf-8 -*-
# This ES-PIC1d_init.py initialize the PIC simulation

import Constants as cst
from scipy.stats import norm
import numpy as np

def thermal_velocity(charge, temperature, mass):
    """
    Most probable speed in the velocity distribution = sqrt(2eT/m)
    :param charge: particle charge in C
    :param temperature: particle temp in eV. 1 Joule = 1 Volt * 1 Coulomb 
    :param mass: particle mass in kg
    :return: thermal velocity in m/s
    """
    return np.sqrt(2*abs(charge)*temperature/mass)

def norm_distribution(temperature,mass,num_ptcl):
    """
    Returns an array of velocities sampled from the Normal distribution
    1d maxwell velocity distribution is a normal distribution
    f(v) = sqrt(m/2/pi/Kb/T)*exp(-v**2*m/2/kb/T)
    norm distribution f(x) = sqrt(1/2/pi)*exp(-x**2/2)
    v = x*m/kb/T, so scale = 1(m/kb/T)
    :param temperature: particle temperature, in eV
    :param mass: particle mass
    :param num_velocities: number of velocities to generate
    :returns: 1d array array of velocities
    """
    tempK = temperature*cst.EV2K
    scale = 1.0/(mass/cst.KB/tempK)
    vels = norm.rvs(loc=0.0, scale=scale, size=num_ptcl)
    
    return vels
    

def maxwell_velocity_distribution(v_thermal, num_velocities):
    """
    Returns an array of velocities sampled from the Maxwell distribution
    :param v_thermal: thermal velocity (most probable) of the distribution
    :param num_velocities: number of velocities to generate
    :returns: 2 by num_velocities array of velocities
    """
    a = v_thermal / np.sqrt(2)  # shape parameter of distribution

    maxwell = stats.maxwell

    speeds = maxwell.rvs(loc=0, scale=a, size=num_velocities)  # generate speeds
    theta = np.random.rand(num_velocities) * 2 * np.pi  # select random angle

    x_vels = speeds * np.sin(theta)
    y_vels = speeds * np.cos(theta)

    return np.stack((x_vels, y_vels))

def debye_length(eps0, eon_temp, ne):
    """
    Calculate plasma debye length
    :param eps0: relative permittivity
    :param eon_temp: electron temp in eV
    :param ne: electron density in particles/volume
    :return: debye length
    """
    return np.sqrt(eps0*cst.EPS0*eon_temp/ne/cst.UNIT_CHARGE)

def sinusoidal(amplitude, period, phase, t):
    return np.sin(2*np.pi*(t/period) + phase) * amplitude