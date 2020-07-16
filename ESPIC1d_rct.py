# -*- coding: utf-8 -*-
"""
Calc Particle generations due to reactions
"""
import numpy as np
import Constants as cst

def ioniz(Eon_pv, Arp_pv, dt, ptcl_add):
    """
    assume Ar ionization threshold is 10 eV
    eon energy > 100 eV rate coeff 1e-9 cm3/s = 1e-15 m3/s
    eon energy > 10 eV and < 100 eV, linear interpolation
    prob = Ar_den*rate_coeff*dt
    """
    Eon_ergs = np.power(Eon_pv[1]*cst.VEL2EV,2)
    # assuming 0.1% prob for ionization during dt for eon erg > 20 eV
    # assuming erg loss = threshold after ionization
    # created eon has initial vels = 0.0 m/s
    idx = np.where(Eon_ergs >= 20.0)[0]
    rand = np.random.uniform(0.0, 1.0, len(idx))
    idx = idx[rand <= 0.01*dt/1e-12]

    posn_add = Eon_pv[0][idx]
    Eon_pv[1][idx] *= np.sqrt((Eon_ergs[idx]-20.0)/Eon_ergs[idx])
    
    Eon_pv[0] = np.append(Eon_pv[0], posn_add)
    Arp_pv[0] = np.append(Arp_pv[0], posn_add)
    Eon_pv[1] = np.append(Eon_pv[1], np.zeros(posn_add.size))
    Arp_pv[1] = np.append(Arp_pv[1], np.zeros(posn_add.size))
    ptcl_add.append(posn_add.size)

    return Eon_pv, Arp_pv, ptcl_add
