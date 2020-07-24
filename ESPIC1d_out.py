# -*- coding: utf-8 -*-
# This ES-PIC1d_out.py output all the results and diagnostics
#
# Output results and diagnostics
#
import Constants as cst
from Species import Eon, Hp, H
import matplotlib.pyplot as plt
import numpy as np

def vels2ergs(v_clct, sp):
    v_clct = np.asarray(v_clct)
    sign = np.sign(v_clct)
    e_clct = 0.5*sp.mass*cst.AMU*np.power(v_clct, 2)*cst.J2EV
    e_clct *= sign
    return e_clct

def correct_ylimit(ax, x, y):   
   # ax: axes object handle
   #  x: data for entire x-axes
   #  y: data for entire y-axes
   # assumption: you have already set the x-limit as desired
   lims = ax.get_xlim()
   i = np.where( (x > lims[0]) &  (x < lims[1]) )[0]
   ax.set_ylim( y[i].min()*0.9, y[i].max()*1.1 ) 

def plot_diag(mesh, Eon_pv, Ion_pv, Eon_clct, Ion_clct, 
              chrg_den, pe, iteration, ergs_mean, ergs_max, ptcl_rec,
              nout_iter):
    fig, ax = plt.subplots(3,3, figsize=(15,9),
      constrained_layout=True)
    ax[0,0].hist(Eon_pv[0], bins=20, histtype='step', color='blue')
    ax[0,0].hist(Ion_pv[0], bins=20, histtype='step', color='red')
    ax[0,0].set_title('density')
    ax_temp00 = ax[0,0].twinx()
    ax_temp00.plot(mesh.gridx, chrg_den, 'k-')
    ax[0,1].plot(mesh.gridx, pe[0], 'k-')
    ax_temp01 = ax[0,1].twinx()
    ax_temp01.plot(mesh.cell_cnt, pe[1], 'g-')
    ax[0,2].plot(Eon_pv[0], vels2ergs(Eon_pv[1], Eon), 'bo')
    ax[0,2].plot(Ion_pv[0], vels2ergs(Ion_pv[1], Hp),  'ro')
    ax[0,2].set_title('nPtcl=%d' % ptcl_rec[0][-1])
    
    ax[1,0].hist(vels2ergs(Eon_clct[0], Eon), bins=20, 
                 histtype='step', color='blue')
    ax[1,0].hist(vels2ergs(Eon_clct[1], Eon), bins=20, 
                 histtype='step', color='red')
    ax[1,0].set_title('Energy count')
    ax_temp10 = ax[1,0].twinx()
    ax_temp10.hist(vels2ergs(Ion_clct[0], Hp), bins=20, 
                   histtype='bar', color='blue')
    ax_temp10.hist(vels2ergs(Ion_clct[1], Hp), bins=20, 
                   histtype='bar', color='red')
    
    temp_xmin = iteration-10*nout_iter
    ax[1,1].set_xlim(temp_xmin, iteration)
    tempx = range(iteration+1)
    tempy = np.asarray(ergs_mean[0])
    ax[1,1].plot(tempx, tempy, 'b-')
    correct_ylimit(ax[1,1], tempx, tempy)
    ax_temp11 = ax[1,1].twinx()
    tempy = np.asarray(ergs_max[0])
    ax_temp11.plot(tempx, tempy, 'r-')
    correct_ylimit(ax_temp11, tempx, tempy)
    
    ax[2,0].plot(ptcl_rec[0], 'b-')
    ax[2,0].plot(ptcl_rec[3], 'r-')
    ax_temp20 = ax[2,0].twinx()
    ax_temp20.plot([tempi - tempe for tempe, tempi in 
                    zip(ptcl_rec[0], ptcl_rec[3])], 'k--')
    
    
    ax[2,1].set_xlim(temp_xmin, iteration)
    tempx = range(0, iteration, nout_iter)
    tempy = np.add.reduceat(ptcl_rec[1], np.arange(0, len(ptcl_rec[1]), nout_iter))
    ax[2,1].plot(tempx[0:-1], tempy[0:-1], 'bo-')
    tempy = np.add.reduceat(ptcl_rec[4], np.arange(0, len(ptcl_rec[4]), nout_iter))
    ax[2,1].plot(tempx[0:-1], tempy[0:-1], 'ro-')
    correct_ylimit(ax[2,1], tempx, tempy)
    ax_temp21 = ax[2,1].twinx()
    tempy = np.add.reduceat(ptcl_rec[2], np.arange(0, len(ptcl_rec[4]), nout_iter))
    ax_temp21.plot(tempx[0:-1], tempy[0:-1], 'bo--', mfc='none')
    tempy = np.add.reduceat(ptcl_rec[5], np.arange(0, len(ptcl_rec[5]), nout_iter))
    ax_temp21.plot(tempx[0:-1], tempy[0:-1], 'ro--', mfc='none')
    correct_ylimit(ax_temp21, tempx, tempy)
    fig.savefig('.\Figures\ITER_{:08}.png'.format(iteration))
    plt.close(fig)


    
