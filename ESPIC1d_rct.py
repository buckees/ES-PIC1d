# -*- coding: utf-8 -*-
"""
Calc Particle generations due to reactions
"""
import numpy as np

def ioniz(Eon_pv, Arp_pv, ergs, dt, ptcl_add):
    """
    assume Ar ionization threshold is 10 eV
    eon energy > 100 eV rate coeff 1e-9 cm3/s = 1e-15 m3/s
    eon energy > 10 eV and < 100 eV, linear interpolation
    prob = Ar_den*rate_coeff*dt
    """
    # assuming 0.1% prob for ionization during dt for eon erg > 100 eV
    # assuming no erg loss after ionization
    # created eon has initial vels = 0.0 m/s
    posn_add = Eon_pv[0][ergs >= 200.0]
    if not posn_add.size:
        ptcl_add.append(posn_add.size)
        return Eon_pv, Arp_pv, ptcl_add
    
    rand = np.random.uniform(0.0, 0.1, len(posn_add))
    posn_add = posn_add[rand <= 0.01*dt/1e-12]
    if not posn_add.size:
        ptcl_add.append(posn_add.size)
        return Eon_pv, Arp_pv, ptcl_add
    
    Eon_pv[0] = np.append(Eon_pv[0], posn_add)
    Arp_pv[0] = np.append(Arp_pv[0], posn_add)
    Eon_pv[1] = np.append(Eon_pv[1], np.zeros(posn_add.size))
    Arp_pv[1] = np.append(Arp_pv[1], np.zeros(posn_add.size))
    ptcl_add.append(posn_add.size)
#    print('%d particles are added' % posn_add.size)
    return Eon_pv, Arp_pv, ptcl_add