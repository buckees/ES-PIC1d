# -*- coding: utf-8 -*-
"""
create mesh class and mesh informaiton
"""
import numpy as np

class Mesh(object):
    """Stores all mesh information"""
    def __init__(self, ncellx, width, dx, gridx, cell_cnt):
        self.ncellx = ncellx # number of cells in x direction
        self.width = width # domain size in 1d, in meter
        self.dx = dx # cell size, in meter
        self.gridx = gridx # an array of grids in x direction, in meter
        self.cell_cnt = cell_cnt # an array of cell center in x direction, in meter

ncellx = 100
width = 0.01
gridx, dx = np.linspace(0.0, width, ncellx+1, retstep=True)
cell_cnt = (gridx[0:-1] + gridx[1:])/2.0
Mesh = Mesh(ncellx, width, dx, gridx, cell_cnt)
