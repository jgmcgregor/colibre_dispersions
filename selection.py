# Script to produce selection plots for COLIBRE galaxies in a given snapshot. Galaxies are selected according to M_star, DeltaMS, and kappa_corot
# JG McGregor
# March 2025 

import argparse
import numpy as np
import matplotlib.pyplot as plt
import unyt
import swiftsimio as sw

parser = argparse.ArgumentParser()
parser.add_argument("device", help="local or cosma")
parser.add_argument("length", help="length, LXXX")
parser.add_argument("mres", help="mass res, mX")
parser.add_argument("snap", help="4-digit snapshot number")
args = parser.parse_args()

device = args.device
length = args.length
mres = args.mres
snap = args.snap
run = length+"_"+mres

if device == "local":
    filename = "cosma_files/"+run+"/THERMAL_AGN_"+mres+"/SOAP/halo_properties_"+snap+".hdf5"
elif device == "cosma":
    filename = "/cosma8/data/dp004/colibre/Runs/"+run+"/THERMAL_AGN_"+mres+"/SOAP/halo_properties_"+snap+".hdf5"
else:
    print("device not recognised")
    exit()


catalogue = sw.load(filename)
meta = catalogue.metadata

z=meta.redshift
z_short = round(z,3)

# universe age provided in inverse Hubble constant units
t = meta.cosmology_raw['Universe age [internal units]'][0] * unyt.Mpc * unyt.s / unyt.km
t.convert_to_units(unyt.Gyr)

# uses 30kpc spheres, not entire subhalo
Mstar = catalogue.exclusive_sphere_30kpc.stellar_mass
SFR = catalogue.exclusive_sphere_30kpc.star_formation_rate
kappa_corot = catalogue.exclusive_sphere_30kpc.kappa_corot_stars
disk_fraction = catalogue.exclusive_sphere_30kpc.disc_to_total_stellar_mass_fraction
Mgas = catalogue.exclusive_sphere_30kpc.gas_mass
gas_fraction = Mgas/(Mstar+Mgas)
mu_gas = Mgas/Mstar

# based on Popesso+23
SFMS_masses = np.logspace(8.5,11.5,100) * unyt.Msun
m = np.log10(SFMS_masses/unyt.Msun) + 0.025
tdim = t/unyt.Gyr
MSMZ = (-0.034*tdim+4.722)*m - 0.1925*m**2 + (0.2*tdim-26.16)
SFMS_SFRs = 10**MSMZ * (unyt.Msun/unyt.yr)

# can adjust galaxy selection here
masscut = 10**9.5 * unyt.Msun
softmasscut = 10**8.5 * unyt.Msun # galaxies below this don't have an accurate SFMS
dMScut = -0.3 # based on scatter from Popesso+23
kappacut = 0.4

m = np.log10(Mstar/unyt.Msun) + 0.025
SFR_MS = 10**((-0.034*tdim+4.722)*m - 0.1925*m**2 + (0.2*tdim-26.16)) * (unyt.Msun/unyt.yr) # bestfit from Popesso+23
deltaMS = np.log10(SFR/SFR_MS)

selection = (Mstar>masscut) & (deltaMS>dMScut) & (kappa_corot>kappacut)
softselection = Mstar>softmasscut

mass_lo=1e8
mass_hi=1e12
SFR_lo=1e5
SFR_hi=1e12
dMS_lo=-3
dMS_hi=1
kappa_lo=0
kappa_hi=1
gasfrac_lo=1e-2
gasfrac_hi=1e0
mugas_lo=1e-3
mugas_hi=1e2

unyt.matplotlib_support.label_style = '[]'
with unyt.matplotlib_support:

    fig,ax = plt.subplots(1,1,figsize=(6,6))

    ax.scatter(Mstar[softselection],SFR[softselection], xunits='Msun',yunits='Msun/Gyr',c='k')
    ax.scatter(Mstar[selection],SFR[selection], xunits='Msun',yunits='Msun/Gyr',c='r')
    #ax.scatter(Mstar[target],SFR[target], xunits='Msun',yunits='Msun/Gyr',c='b',marker='*',s=200)
    ax.plot(SFMS_masses,SFMS_SFRs,'k')
    ax.plot(SFMS_masses,SFMS_SFRs*10**dMScut,'r--')
    ax.vlines(x=masscut,ymin=SFR_lo,ymax=SFR_hi,colors='r',linestyles='--')
    ax.set_xscale('log')
    ax.set_xlim(mass_lo,mass_hi)
    ax.set_xlabel(r'$M_*$ [M$_\odot$]')
    ax.set_yscale('log')
    ax.set_ylim(SFR_lo,SFR_hi)
    ax.set_ylabel(r'SFR [M$_\odot$/Gyr]')
    ax.set_title('z='+str(z_short))

    imgname = 'plots/selection_plots/SFMS_' + str(run) + 'z' + str(z_short) + '.png'

    plt.savefig(imgname,bbox_inches='tight')



unyt.matplotlib_support.label_style = '[]'
with unyt.matplotlib_support:

    fig,((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(18,18))

    ax1.scatter(Mstar[softselection],deltaMS[softselection], xunits='Msun',c='k')
    ax1.scatter(Mstar[selection],deltaMS[selection], xunits='Msun',c='r')
    #ax1.scatter(Mstar[target],deltaMS[target], xunits='Msun',c='b',marker='*',s=200)
    ax1.vlines(x=masscut,ymin=dMS_lo,ymax=dMS_hi,colors='r',linestyles='--')
    ax1.hlines(y=dMScut,xmin=mass_lo,xmax=mass_hi,colors='r',linestyles='--')
    ax1.set_xscale('log')
    ax1.set_xlim(mass_lo,mass_hi)
    ax1.set_xlabel('')
    #ax1.set_xlabel(r'$M_*$ [M$_\odot$]')
    ax1.set_ylim(dMS_lo,dMS_hi)
    ax1.set_ylabel(r'$\Delta$MS [dex]')

    ax2.set_title('z='+str(z_short))
    ax2.set_axis_off()

    ax3.set_axis_off()

    ax4.scatter(Mstar[softselection],kappa_corot[softselection], xunits='Msun',c='k')
    ax4.scatter(Mstar[selection],kappa_corot[selection], xunits='Msun',c='r')
    #ax4.scatter(Mstar[target],kappa_corot[target], xunits='Msun',c='b',marker='*',s=200)
    ax4.vlines(x=masscut,ymin=kappa_lo,ymax=kappa_hi,colors='r',linestyles='--')
    ax4.hlines(y=kappacut,xmin=mass_lo,xmax=mass_hi,colors='r',linestyles='--')
    ax4.set_xscale('log')
    ax4.set_xlim(mass_lo,mass_hi)
    ax4.set_xlabel('')
    #ax4.set_xlabel(r'$M_*$ [M$_\odot$]')
    ax4.set_ylim(kappa_lo,kappa_hi)
    ax4.set_ylabel(r'$\kappa_{corot,*}$')

    ax5.scatter(deltaMS[softselection],kappa_corot[softselection],c='k')
    ax5.scatter(deltaMS[selection],kappa_corot[selection],c='r')
    #ax5.scatter(deltaMS[target],kappa_corot[target],c='b',marker='*',s=200)
    ax5.hlines(y=kappacut,xmin=dMS_lo,xmax=dMS_hi,colors='r',linestyles='--')
    ax5.vlines(x=dMScut,ymin=kappa_lo,ymax=kappa_hi,colors='r',linestyles='--')
    ax5.set_xlim(dMS_lo,dMS_hi)
    ax5.set_xlabel('')
    #ax5.set_xlabel(r'$\Delta$MS [dex]')
    ax5.set_ylim(kappa_lo,kappa_hi)
    ax5.set_ylabel('')
    #ax5.set_ylabel(r'$\kappa_{corot,*}$')

    ax6.set_axis_off()

    ax7.scatter(Mstar[softselection],mu_gas[softselection], xunits='Msun',c='k')
    ax7.scatter(Mstar[selection],mu_gas[selection], xunits='Msun',c='r')
    #ax7.scatter(Mstar[target],mu_gas[target], xunits='Msun',c='b',marker='*',s=200)
    ax7.vlines(x=masscut,ymin=mugas_lo,ymax=mugas_hi,colors='r',linestyles='--')
    ax7.set_xscale('log')
    ax7.set_xlim(mass_lo,mass_hi)
    ax7.set_xlabel(r'$M_*$ [M$_\odot$]')
    ax7.set_yscale('log')
    ax7.set_ylim(mugas_lo,mugas_hi)
    ax7.set_ylabel(r'$M_{gas} / M_*$')

    ax8.scatter(deltaMS[softselection],mu_gas[softselection],c='k')
    ax8.scatter(deltaMS[selection],mu_gas[selection],c='r')
    #ax8.scatter(deltaMS[target],mu_gas[target],c='b',marker='*',s=200)
    ax8.vlines(x=dMScut,ymin=mugas_lo,ymax=mugas_hi,colors='r',linestyles='--')
    ax8.set_xlim(dMS_lo,dMS_hi)
    ax8.set_xlabel(r'$\Delta$MS [dex]')
    ax8.set_yscale('log')
    ax8.set_ylim(mugas_lo,mugas_hi)
    ax8.set_ylabel('')
    #ax8.set_ylabel(r'$M_{gas} / M_*$')

    ax9.scatter(kappa_corot[softselection],mu_gas[softselection],c='k')
    ax9.scatter(kappa_corot[selection],mu_gas[selection],c='r')
    #ax9.scatter(kappa_corot[target],mu_gas[target],c='b',marker='*',s=200)
    ax9.vlines(x=kappacut,ymin=mugas_lo,ymax=mugas_hi,colors='r',linestyles='--')
    ax9.set_xlim(kappa_lo,kappa_hi)
    ax9.set_xlabel(r'$\kappa_{corot,*}$')
    ax9.set_yscale('log')
    ax9.set_ylim(mugas_lo,mugas_hi)
    ax9.set_ylabel('')
    #ax9.set_ylabel(r'$M_{gas} / M_*$')

    plt.subplots_adjust(wspace=0.1,hspace=0.1)
    
    imgname = 'plots/selection_plots/triangle_' + str(run) + 'z' + str(z_short) + '.png'

    plt.savefig(imgname,bbox_inches='tight')



