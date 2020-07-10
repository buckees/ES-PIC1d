# -*- coding: utf-8 -*-
# This ES-PIC1d_out.py output all the results and diagnostics
#
# Output results and diagnostics
#
import matplotlib.pyplot as plt

def plot_diag(mesh, Eon_pv, Arp_pv,
              chrg_den, pe, iteration, num_ptcl):
    fig, ax = plt.subplots(2,3, figsize=(15,6),
      constrained_layout=True)
    ax[0,0].hist(Eon_pv[0], bins=20, histtype='step', color='blue')
    ax[0,0].hist(Arp_pv[0], bins=20, histtype='step', color='red')
    ax[0,0].set_title('density')
    ax_temp0 = ax[0,0].twinx()
    ax_temp0.plot(mesh.gridx, chrg_den, 'k-')
    ax[0,1].plot(mesh.gridx, pe[0], 'k-')
    ax_temp1 = ax[0,1].twinx()
    ax_temp1.plot(mesh.cell_cnt, pe[1], 'g-')
    ax[0,2].plot(Eon_pv[0], Eon_pv[1], 'bo')
    ax[0,2].plot(Arp_pv[0], Arp_pv[1], 'ro')
    ax[0,2].set_title('nPtcl=%d' % num_ptcl)
    fig.savefig('.\Figures\ITER_{:08}.png'.format(iteration))
    plt.close(fig)

    