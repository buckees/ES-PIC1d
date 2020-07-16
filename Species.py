# -*- coding: utf-8 -*-
# Define the class of SP

class SP(object):
    """Stores all SPs"""
    def __init__(self, name, type, mass, charge, tmpt):
        self.name = name # str
        self.type = type # str, 'E','Ion','Neut' or 'Bkg'
        self.mass = mass # unit in kg
        self.charge = charge # unit in Coulomb
        self.tmpt = tmpt # unit in eV

# proton-electron mass ratio = 1836.15
Eon = SP('E',   'E',    5.45e-4, -1.0, 5.0) 
Arp = SP('Ar+', 'Ion',  32.0,     1.0, 0.2)
Ar  = SP('Ar',  'Bkg',  32.0,     0.0, 0.025)
Hp  = SP('H+',  'Ion',   1.0,     1.0, 0.2)
H   = SP('H',   'Bkg',   1.0,     0.0, 0.025)
N2  = SP('N2',  'Bkg',  28.0,     0.0, 0.025)