# -*- coding: utf-8 -*-
# This ES-PIC1d_out.py output all the results and diagnostics
#
# Output results and diagnostics
#
import matplotlib.pyplot as plt
import numpy as np

def plot_diag(mesh, Eon_pv, Arp_pv, Eon_clct, Arp_clct, 
              chrg_den, pe, iteration, ergs_mean, ergs_max, ptcl_rec,
              nout_iter):
    fig, ax = plt.subplots(2,3, figsize=(15,6),
      constrained_layout=True)
    ax[0,0].hist(Eon_pv[0], bins=20, histtype='step', color='blue')
    ax[0,0].hist(Arp_pv[0], bins=20, histtype='step', color='red')
    ax[0,0].set_title('density')
    ax_temp00 = ax[0,0].twinx()
    ax_temp00.plot(mesh.gridx, chrg_den, 'k-')
    ax[0,1].plot(mesh.gridx, pe[0], 'k-')
    ax_temp01 = ax[0,1].twinx()
    ax_temp01.plot(mesh.cell_cnt, pe[1], 'g-')
    ax[0,2].plot(Eon_pv[0], Eon_pv[1], 'bo')
    ax[0,2].plot(Arp_pv[0], Arp_pv[1], 'ro')
    ax[0,2].set_title('nPtcl=%d' % ptcl_rec[0][-1])
    
    ax[1,0].hist(Eon_clct[0], bins=20, histtype='step', color='blue')
    ax[1,0].hist(Eon_clct[1], bins=20, histtype='step', color='red')
    ax[1,0].set_title('v_Eon count')
    ax_temp10 = ax[1,0].twinx()
    ax_temp10.hist(Arp_clct[0], bins=20, histtype='step', color='blue')
    ax_temp10.hist(Arp_clct[1], bins=20, histtype='step', color='red')
    ax_temp10.set_title('v_Ar+ count')
    
    ax[1,1].plot(ergs_mean[0], 'b-')
    ax_temp11 = ax[1,1].twinx()
    ax_temp11.plot(ergs_max[0], 'r-')
    
    ax[1,2].plot(ptcl_rec[0], 'k-')
    ax_temp12 = ax[1,2].twinx()
    temp = np.add.reduceat(ptcl_rec[1], np.arange(0, len(ptcl_rec[1]), nout_iter))
    ax_temp12.plot(range(0, len(ptcl_rec[1]),nout_iter), temp, 'g-')
    temp = np.add.reduceat(ptcl_rec[2], np.arange(0, len(ptcl_rec[2]), nout_iter))
    ax_temp12.plot(range(0, len(ptcl_rec[2]),nout_iter), temp, 'y-')
    fig.savefig('.\Figures\ITER_{:08}.png'.format(iteration))
    plt.close(fig)

    