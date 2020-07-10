# -*- coding: utf-8 -*-
# Define the class of Particle

import Constants as cst

class Particle(object):
    """Stores all particles"""
    def __init__(self, name, type, mass, charge, tmpt):
        self.name = name # str
        self.type = type # str, 'E','Ion','Neut' or 'Bkg'
        self.mass = mass # unit in kg
        self.charge = charge # unit in Coulomb
        self.tmpt = tmpt # unit in eV

Eon = Particle('E',   'E',     1.0, -1.0, 2.0)
Arp = Particle('Ar+', 'Ion',  32.0,  1.0, 0.1)
Ar  = Particle('Ar',  'Bkg',  32.0,  0.0, 0.025)
H   = Particle('H',   'Neut',  1.0,  0.0, 0.025)
N2  = Particle('N2',  'Bkg',  28.0,  0.0, 0.025)