# Script to produce selection plots for COLIBRE galaxies in a given snapshot. Galaxies are selected according to M_star, DeltaMS, and kappa_corot
# JG McGregor
# March 2025 

import sys
import numpy as np
import matplotlib.pyplot as plt
import unyt
import swiftsimio as sw

run = sys.argv[1]
snap = sys.argv[2]

filename = "SOAP/" + run + "/halo_properties_" + snap + ".hdf5"

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

# based on Popesso+23
SFMS_masses = np.logspace(8.5,11.5,100) * unyt.Msun
m = np.log10(SFMS_masses/unyt.Msun) + 0.025
tdim = t/unyt.Gyr
MSMZ = (-0.034*tdim+4.722)*m - 0.1925*m**2 + (0.2*tdim-26.16)
SFMS_SFRs = 10**MSMZ * (unyt.Msun/unyt.yr)

# can adjust galaxy selection here
masscut = 10**9.5 * unyt.Msun
softmasscut = 10**8.5 * unyt.Msun # galaxies below this don't have an accurate SFMS
dMScut = -0.09 # based on scatter from Popesso+23
kappacut = 0.4

m = np.log10(Mstar/unyt.Msun) + 0.025
SFR_MS = 10**((-0.034*tdim+4.722)*m - 0.1925*m**2 + (0.2*tdim-26.16)) * (unyt.Msun/unyt.yr) # bestfit from Popesso+23
deltaMS = np.log10(SFR/SFR_MS)

selection = (Mstar>masscut) & (deltaMS>dMScut) & (kappa_corot>kappacut)
softselection = Mstar>softmasscut

unyt.matplotlib_support.label_style = '[]'
with unyt.matplotlib_support:

    fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(18,6))

    ax1.scatter(Mstar,SFR, xunits='Msun',yunits='Msun/Gyr',c='k')
    ax1.scatter(Mstar[selection],SFR[selection], xunits='Msun',yunits='Msun/Gyr',c='r')
    ax1.plot(SFMS_masses,SFMS_SFRs,'k')
    ax1.plot(SFMS_masses,SFMS_SFRs*10**dMScut,'r--')
    ax1.vlines(x=masscut,ymin=min(SFR),ymax=10*max(SFR),colors='r',linestyles='--')
    ax1.set_xscale('log')
    ax1.set_xlim(1e6,3e11)
    ax1.set_yscale('log')
    ax1.set_ylim(1e5,1e10)

    ax2.scatter(Mstar[softselection],kappa_corot[softselection], xunits='Msun',c='k')
    ax2.scatter(Mstar[selection],kappa_corot[selection], xunits='Msun',c='r')
    ax2.vlines(x=masscut,ymin=min(kappa_corot),ymax=max(kappa_corot),colors='r',linestyles='--')
    ax2.hlines(y=kappacut,xmin=min(Mstar),xmax=10*max(Mstar),colors='r',linestyles='--')
    ax2.set_xscale('log')
    ax2.set_xlim(1e6,3e11)
    ax2.set_ylabel(r'$\kappa_{corot,*}$')
    ax2.set_ylim(0,0.8)
    ax2.set_title('z='+str(z_short))

    ax3.scatter(deltaMS[softselection],kappa_corot[softselection],c='k')
    ax3.scatter(deltaMS[selection],kappa_corot[selection],c='r')
    ax3.hlines(y=kappacut,xmin=-2.5,xmax=1,colors='r',linestyles='--')
    ax3.vlines(x=dMScut,ymin=0,ymax=1,colors='r',linestyles='--')
    ax3.set_xlabel(r'$\Delta$MS [dex]')
    ax3.set_xlim(-2.5,1)
    ax3.set_ylabel(r'$\kappa_{corot,*}$')
    ax3.set_ylim(0,0.8)

    imgname = 'plots/selection_plots/z' + str(z_short) + '.png'

    plt.savefig(imgname,bbox_inches='tight')






