# -*- coding: utf-8 -*-
"""
Calc Particle generations due to reactions
"""
import numpy as np
import Constants as cst

def ioniz(Eon_pv, Ion_pv, dt, ptcl_add):
    """
    assume Ar ionization threshold is 10 eV
    eon energy > 100 eV rate coeff 1e-9 cm3/s = 1e-15 m3/s
    eon energy > 10 eV and < 100 eV, linear interpolation
    prob = Ar_den*rate_coeff*dt
    """
    Eon_ergs = np.power(Eon_pv[1]*cst.VEL2EV, 2)
    # assuming 0.1% prob for ionization during dt for eon erg > 20 eV
    # assuming erg loss = threshold after ionization
    # created eon has initial vels = 0.0 m/s
    idx = np.where(Eon_ergs >= 10.0)[0]
    prob = ioniz_prob(Eon_ergs[idx])
    rand = np.random.uniform(0.0, 1.0, len(idx))
    idx = idx[rand <= prob*dt/1e-12]

    posn_add = Eon_pv[0][idx]
    Eon_pv[1][idx] *= np.sqrt((Eon_ergs[idx]-10.0)/Eon_ergs[idx])
    
    Eon_pv[0] = np.append(Eon_pv[0], posn_add)
    Ion_pv[0] = np.append(Ion_pv[0], posn_add)
    Eon_pv[1] = np.append(Eon_pv[1], np.zeros(posn_add.size))
    Ion_pv[1] = np.append(Ion_pv[1], np.zeros(posn_add.size))
    ptcl_add.append(posn_add.size)

    return Eon_pv, Ion_pv, ptcl_add

def ioniz_prob(ergs):
    """
    convert ergs to prob
    10 eV ---- 5.5e-4
    110 eV ---- 0.082
    exponential interpolation
    ergs:  0 eV ---- 10 eV ---- 110 eV ---- infinity
    prob:  0.0  ----  exp(0.05x - 8.0) ---- 0.082
    """
    prob = np.exp(0.05*ergs - 8.0)
    return prob

def pow_dep(vels, Eon_ergs, powE, dt):
    """
    update Eon vels due to power deposition
    power can be splitted to Eons evenly or due to some distribution
    power correlates to mean Te: more eons, less Te; less eons, more Te.
    :param vels: stores the Eon vels
    :param Eon_ergs: stores the Eon energies
    :param powE: power deposited to Eons only
    :param dt: timestep
    """
    # instead of evenly splitted the power
    # power is depositted more to "hot" Eons
    powE_add = powE*dt*cst.J2EV*(Eon_ergs/np.sum(Eon_ergs))
    ergs_change = (Eon_ergs + powE_add)/Eon_ergs # rate of energy change
    vels = vels*np.power(ergs_change,0.5)
    return vels