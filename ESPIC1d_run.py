# -*- coding: utf-8 -*-
# This ES-PIC_run.py runs the PIC simulation
# All in SI unit

import time
t0 = time.time()


# modules from *.py in the same folder
import Constants as cst
import ESPIC1d_init as init
import ESPIC1d_move as move
import ESPIC1d_out as out
import Particle as ptcl
import Poisson_solver_1d as ps1d

# python built-in modules
import numpy as np
import matplotlib.pyplot as plt

# Geometry and Mesh
init.Mesh.ncellx = 100 # number of cells in x direction
init.Mesh.width = 0.01 # lenght of domain, in meter
# init.Mesh contains all mesh information
init.Mesh.gridx, init.Mesh.dx = np.linspace(0.0, init.Mesh.width,
                                            init.Mesh.ncellx+1, retstep=True)
init.Mesh.cell_cnt = (init.Mesh.gridx[0:-1] + init.Mesh.gridx[1:])/2.0

# Operation Conditions
press_mT = 100.0 # pressure, in mTorr
pressure = press_mT/1.0e3*cst.TORR2PA # in Pa, 1 Torr = 133.322 Pa
nAr_init = pressure/(cst.KB*ptcl.Ar.temp*cst.EV2K) # in m-3, N/V = P/(kb*T) ideal gas law
ptcl.Eon.temp = 1.0 # in eV
ptcl.Arp.temp = 0.1 # in eV
bc = (0.0, 0.0) # left and right boundary conditions, in Volt 
nEon_init = 1.0e15 # in m-3, initial electron density
den_limit = 1.0e11 # in m-3, lower limit of denisty, avoid 0 density 

# Model Parameters
num_ptcl = 10000 # number of particles, should be >> ncellx to reduce noise
dt = 1.0e-10 # in sec

den_per_ptcl = nEon_init/num_ptcl # density contained in one particle

# initialize the position and velcotiy in a dataframe
posn_Eon, vels_Eon = init.init_data(num_ptcl,ptcl.Eon) #posn in unit length
posn_Arp, vels_Arp = init.init_data(num_ptcl,ptcl.Arp)
posn_Eon = posn_Eon*init.Mesh.width
posn_Arp = posn_Eon.copy()

# create mesh
#gridx, cell_center, dx = init.init_mesh(ncellx,width)

# initilize potential on cell boudaries
pot, efld = init.init_pot(init.Mesh,bc)

# assign charge densities to grid nodes, unit in UNIT_CHARGE
den_Eon = move.den_asgmt(posn_Eon,init.Mesh)
den_Arp = move.den_asgmt(posn_Arp,init.Mesh)
den_chrg = (ptcl.Eon.charge*den_Eon + ptcl.Arp.charge*den_Arp)*den_per_ptcl

# update potential according to assigned charges to nodes
pot, efld = ps1d.Poisson_solver_1d(init.Mesh,den_chrg,bc)
# move particles
posn_Eon, vels_Eon = move.move_ptcl(ptcl.Eon,posn_Eon,vels_Eon,efld,dt,init.Mesh)
posn_Arp, vels_Arp = move.move_ptcl(ptcl.Arp,posn_Arp,vels_Arp,efld,dt,init.Mesh)

num_iter = 600001 # number of iterations
nout_iter = 300
for i in range(num_iter):
    # assign charge densities to grid nodes, unit in UNIT_CHARGE
    den_Eon = move.den_asgmt(posn_Eon,init.Mesh)
    den_Arp = move.den_asgmt(posn_Arp,init.Mesh)
    den_chrg = (ptcl.Eon.charge*den_Eon + 
                ptcl.Arp.charge*den_Arp)*den_per_ptcl
    
    # update potential according to assigned charges to nodes
    pot, efld = ps1d.Poisson_solver_1d(init.Mesh,den_chrg,bc)
    # move particles
    posn_Eon, vels_Eon = move.move_ptcl(ptcl.Eon,posn_Eon,vels_Eon,
                                        efld,dt,init.Mesh)
    posn_Arp, vels_Arp = move.move_ptcl(ptcl.Arp,posn_Arp,vels_Arp,
                                        efld,dt,init.Mesh)
    num_ptcl = len(posn_Eon)
    if i % nout_iter == 0:
        print('iter = %d' % i)
        # plot animation
        out.plot_diag(init.Mesh,posn_Eon, vels_Eon, posn_Arp, vels_Arp,
                      den_chrg, pot, efld, i, num_ptcl)


from datetime import timedelta
print("-total time %s -" % str(timedelta(seconds=(int(time.time() - t0)))))