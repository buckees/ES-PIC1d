# -*- coding: utf-8 -*-
# Define the class of Particle

import Constants as cst

class Particle(object):
    """Stores all particles"""
    def __init__(self, type, mass, charge,temp):
        self.type = type # str, 'E','Ion' and 'Neut'
        self.mass = mass # unit in kg
        self.charge = charge # unit in Coulomb
        self.temp = temp # unit in eV

Eon = Particle('E',    cst.EON_MASS, -1.0*cst.UNIT_CHARGE, 2.0)
Arp = Particle('Ion',  32.0*cst.AMU,  1.0*cst.UNIT_CHARGE, 0.1)
Ar  = Particle('Neut', 32.0*cst.AMU,  0.0*cst.UNIT_CHARGE, 0.025)
H   = Particle('Neut',  1.0*cst.AMU,  0.0*cst.UNIT_CHARGE, 0.025)
N2  = Particle('Neut', 28.0*cst.AMU,  0.0*cst.UNIT_CHARGE, 0.025)