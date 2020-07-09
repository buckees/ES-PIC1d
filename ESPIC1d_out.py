# -*- coding: utf-8 -*-
# This ES-PIC1d_out.py output all the results and diagnostics
#
# Output results and diagnostics
#
import matplotlib.pyplot as plt

def plot_diag(mesh,posn_Eon, vels_Eon, posn_Arp, vels_Arp,
              den_chrg, pot, efld, iteration):
    fig, ax = plt.subplots(2,3, figsize=(15,6),
      constrained_layout=True)
    ax[0,0].hist(posn_Eon,bins=20,histtype='step',color='blue')
    ax[0,0].hist(posn_Arp,bins=20,histtype='step',color='red')
    ax_temp0 = ax[0,0].twinx()
    ax_temp0.plot(mesh.gridx,den_chrg,'k-')
    ax[0,1].plot(mesh.gridx,pot,'k-')
    ax_temp1 = ax[0,1].twinx()
    ax_temp1.plot(mesh.cell_cnt,efld,'g-')
    ax[0,2].plot(posn_Eon,vels_Eon,'bo')
    ax[0,2].plot(posn_Arp,vels_Arp,'ro')
    fig.savefig('.\Figures\ITER_{:08}.png'.format(iteration))
    plt.close(fig)

    