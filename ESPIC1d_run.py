# -*- coding: utf-8 -*-
# This ES-PIC_run.py runs the PIC simulation
# All in SI unit

# modules from *.py in the same folder
import Constants as cst
import ESPIC1d_init as init
import ESPIC1d_move as move
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
ptcl.Eon.temp = 2.0 # in eV
ptcl.Arp.temp = 0.1 # in eV
bc = (0.0, 0.0) # left and right boundary conditions, in Volt 
nEon_init = 1.0e14 # in m-3, initial electron density
den_limit = 1.0e11 # in m-3, lower limit of denisty, avoid 0 density 

# Model Parameters
num_ptcl = 10000 # number of particles, should be >> ncellx to reduce noise
dt = 1.0e-12 # in sec

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

fig, ax = plt.subplots(2,3, figsize=(15,6),
      constrained_layout=True)
#ax[0,0].plot(gridx,den_Eon,'b-')
#ax[0,0].plot(gridx,den_Arp,'r-')
ax[0,0].hist(posn_Eon,bins=20,histtype='step',color='blue')
ax[0,0].hist(posn_Arp,bins=20,histtype='step',color='red')
ax[0,1].plot(init.Mesh.gridx,pot,'k-')
ax_temp0 = ax[0,1].twinx()
ax_temp0.plot(init.Mesh.cell_cnt,efld,'g-')
ax[0,2].plot(posn_Eon,vels_Eon,'bo')
ax[0,2].plot(posn_Arp,vels_Arp,'ro')
fig.canvas.draw()
#
#ax[1,0].plot(gridx,den_Eon,'b-')
#ax[1,0].plot(gridx,den_Arp,'r-')
ax[1,0].hist(posn_Eon,bins=20,histtype='step',color='blue')
ax[1,0].hist(posn_Arp,bins=20,histtype='step',color='red')
ax[1,1].plot(init.Mesh.gridx,pot,'k-')
ax[1,1].set_title('iter = 0')
ax_temp1 = ax[1,1].twinx()
ax_temp1.plot(init.Mesh.cell_cnt,efld,'g-')
ax[1,2].plot(posn_Eon,vels_Eon,'bo')
ax[1,2].plot(posn_Arp,vels_Arp,'ro')
ax[1,2].set_title('remaining particles = %d' % num_ptcl)

num_iter = 101 # number of iterations
nout_iter = 10
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
    if i % 50 == 0: print('iter = %d' % i)
    if i % nout_iter == 0:
        # plot animation    
        for item in ax[1,:]:
            item.cla()
        ax_temp1.cla()
#        ax[1,0].plot(gridx,den_Eon,'b-')
#        ax[1,0].plot(gridx,den_Arp,'r-')
        ax[1,0].hist(posn_Eon,bins=20,histtype='step',color='blue')
        ax[1,0].hist(posn_Arp,bins=20,histtype='step',color='red')
        ax[1,1].plot(init.Mesh.gridx,pot,'k-')
        ax[1,1].set_title('iter = %d' % i)
        ax_temp1 = ax[1,1].twinx()
        ax_temp1.plot(init.Mesh.cell_cnt,efld,'g-')
        ax[1,2].plot(posn_Eon,vels_Eon,'bo')
        ax[1,2].plot(posn_Arp,vels_Arp,'ro')
        ax[1,2].set_title('remaining particles = %d' % num_ptcl)
#        fig.canvas.draw()
        fig.savefig('.\Figures\ITER_{:05}.png'.format(i))
#        plt.pause(0.1)

fig.savefig('ESPIC1d.png')
#plt.show()
# diagnostic plot
#fig, (ax0,ax1,ax2) = plt.subplots(1,3, figsize=(9,3),
#      constrained_layout=True)
#ax0.plot(gridx,den_chrg)
#ax0.set_title('charge distribution')
#ax1.plot(gridx,pot)
#ax1.set_title('potential distribution')
#ax2.plot(gridx[1:],efld)
#ax2.set_title('E-field distribution')
#plt.show()

#fig, (ax0,ax1) = plt.subplots(2,1, figsize=(8,4))
#data_Eon.plot(x='position',y='index',ax=ax0,kind='scatter',xlim=(0.0,width),
#              c='blue')
#data_Arp.plot(x='position',y='index',ax=ax0,kind='scatter',xlim=(0.0,width),
#              c='red')
#data_Eon.plot.hist(y='position',bins=ncellx,ax=ax0)
#data_Arp.plot.hist(y='position',bins=ncellx,ax=ax0)
#ax1.plot(gridx,den_Eon,'b-')
#ax1.plot(gridx,den_Arp,'r-')
#ax1.plot(gridx,den_chrg,'k-')
#plt.show()