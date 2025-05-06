# Script to produce selection plots for COLIBRE galaxies in a given snapshot. Galaxies are selected according to M_star, DeltaMS, and kappa_corot.
# Images of random selected galaxies can be produced
# JG McGregor
# May 2025 

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import unyt
import swiftsimio as sw
from swiftgalaxy import SWIFTGalaxy, SOAP
from swiftsimio.visualisation.projection import project_gas, project_pixel_grid
from swiftsimio.visualisation import generate_smoothing_lengths
from scipy.spatial.transform import Rotation

parser = argparse.ArgumentParser()
parser.add_argument("device", help="local, cosma, or hyades")
parser.add_argument("length", help="length, LXXX")
parser.add_argument("mres", help="mass res, mX")
parser.add_argument("snap", help="4-digit snapshot number")
parser.add_argument("image", help="zero, one, or all")
args = parser.parse_args()

device = args.device
length = args.length
mres = args.mres
snap = args.snap
run = length+"_"+mres
run_hyades = length+mres
image = args.image

def find_nearest_idx(x,arr): # does what it says
    diff_arr = np.absolute(arr-x)
    index = diff_arr.argmin()
    return index


# hits the correct SOAP catalogues on different machines
if device == "local":
    filename = "cosma_files/"+run+"/THERMAL_AGN_"+mres+"/SOAP/halo_properties_"+snap+".hdf5"
elif device == "cosma":
    filename = "/cosma8/data/dp004/colibre/Runs/"+run+"/THERMAL_AGN_"+mres+"/SOAP/halo_properties_"+snap+".hdf5"
elif device == "hyades":
    filename = "/mnt/su3-pro/colibre/"+run_hyades+"/THERMAL_AGN/SOAP/halo_properties_"+snap+".hdf5"
else:
    print("device not recognised")
    exit()

catalogue = sw.load(filename)
meta = catalogue.metadata
z=meta.redshift
z_short = round(z,3) # used for file naming

# uses 30kpc spheres, not entire subhalo
Mstar = catalogue.exclusive_sphere_30kpc.stellar_mass
SFR = catalogue.exclusive_sphere_30kpc.star_formation_rate
kappa_corot = catalogue.exclusive_sphere_30kpc.kappa_corot_stars
disk_fraction = catalogue.exclusive_sphere_30kpc.disc_to_total_stellar_mass_fraction
Mgas = catalogue.exclusive_sphere_30kpc.gas_mass
Mmolgas = catalogue.exclusive_sphere_30kpc.molecular_hydrogen_mass
mu_molgas = Mmolgas/Mstar

# based on Popesso+23
SFMS_masses = np.logspace(8.5,11.5,100) * unyt.Msun
#m = np.log10(SFMS_masses/unyt.Msun) + 0.025
#tdim = t/unyt.Gyr
#MSMZ = (-0.034*tdim+4.722)*m - 0.1925*m**2 + (0.2*tdim-26.16)
#SFMS_SFRs = 10**MSMZ * (unyt.Msun/unyt.yr)

data = np.load("./COLIBRE_SFMS_m6res.npz") # contains Evgenii's SFMS, at integer redshifts
logSFR = data["log10_SFR_Msun_yr"]
logMstar = data["log10_Mstar_Msun"]
redshifts = data["redshifts"]

zidx = find_nearest_idx(z,redshifts)
zSFMS = logSFR[:,zidx] # grab SFMS at appropriate redshift
valid = np.isfinite(zSFMS) # only fit non-nan data points

a,b,c = np.polyfit(logMstar[valid],zSFMS[valid],deg=2)
def evg_SFMS(mstar):
    lgmstar = np.log10(mstar/unyt.Msun)
    lgSFR = a*lgmstar**2 + b*lgmstar + c
    evg_SFR = 10**lgSFR * (unyt.Msun/unyt.yr)
    return(evg_SFR)
SFMS_SFRs = evg_SFMS(SFMS_masses)

# can adjust galaxy selection here
masscut = 10**9.5 * unyt.Msun
softmasscut = 10**8.5 * unyt.Msun # galaxies below this don't have an accurate SFMS
dMScut = -0.3 # based on typical SFMS scatter
kappacut = 0.4 # based on Correa17

#m = np.log10(Mstar/unyt.Msun) + 0.025
#SFR_MS = 10**((-0.034*tdim+4.722)*m - 0.1925*m**2 + (0.2*tdim-26.16)) * (unyt.Msun/unyt.yr) # bestfit from Popesso+23

SFR_MS = evg_SFMS(Mstar) # expected SFMS for a MS galaxy of this size
deltaMS = np.log10(SFR/SFR_MS)

selection = (Mstar>masscut) & (deltaMS>dMScut) & (kappa_corot>kappacut) # my selected galaxies
softselection = Mstar>softmasscut # for plotting purposes - only galaxies above a certain mass are plotted

candidates = np.argwhere(selection)
tablename = 'tables/targets_' + str(run) + '_z' + str(z_short)
txtname = tablename + '.txt'
np.save(tablename, candidates) # saves target IDs to a .npy file, to be dealt with later...
np.savetxt(txtname, candidates, fmt="%d")

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
mumolgas_lo=1e-3
mumolgas_hi=1e1

unyt.matplotlib_support.label_style = '[]'
with unyt.matplotlib_support:

    fig,ax = plt.subplots(1,1,figsize=(6,6))

    ax.scatter(Mstar[softselection],SFR[softselection], xunits='Msun',yunits='Msun/Gyr',c='k')
    ax.scatter(Mstar[selection],SFR[selection], xunits='Msun',yunits='Msun/Gyr',c='r')
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

    if image=="one":
        candidates = np.argwhere(selection) # finding the ID of our galaxies
        if len(candidates) > 0:
            target = candidates[0][0]
            ax.scatter(Mstar[target],SFR[target], xunits='Msun',yunits='Msun/Gyr',c='b',marker='*',s=200)

    imgname = 'plots/selection_plots/SFMS_' + str(run) + 'z' + str(z_short) + '.png'

    plt.savefig(imgname,bbox_inches='tight')



unyt.matplotlib_support.label_style = '[]'
with unyt.matplotlib_support:

    fig,((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(3,3,figsize=(18,18))

    ax1.scatter(Mstar[softselection],deltaMS[softselection], xunits='Msun',c='k')
    ax1.scatter(Mstar[selection],deltaMS[selection], xunits='Msun',c='r')
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
    ax5.hlines(y=kappacut,xmin=dMS_lo,xmax=dMS_hi,colors='r',linestyles='--')
    ax5.vlines(x=dMScut,ymin=kappa_lo,ymax=kappa_hi,colors='r',linestyles='--')
    ax5.set_xlim(dMS_lo,dMS_hi)
    ax5.set_xlabel('')
    #ax5.set_xlabel(r'$\Delta$MS [dex]')
    ax5.set_ylim(kappa_lo,kappa_hi)
    ax5.set_ylabel('')
    #ax5.set_ylabel(r'$\kappa_{corot,*}$')

    ax6.set_axis_off()

    ax7.scatter(Mstar[softselection],mu_molgas[softselection], xunits='Msun',c='k')
    ax7.scatter(Mstar[selection],mu_molgas[selection], xunits='Msun',c='r')
    ax7.vlines(x=masscut,ymin=mumolgas_lo,ymax=mumolgas_hi,colors='r',linestyles='--')
    ax7.set_xscale('log')
    ax7.set_xlim(mass_lo,mass_hi)
    ax7.set_xlabel(r'$M_*$ [M$_\odot$]')
    ax7.set_yscale('log')
    ax7.set_ylim(mumolgas_lo,mumolgas_hi)
    ax7.set_ylabel(r'$M_{molgas} / M_*$')

    ax8.scatter(deltaMS[softselection],mu_molgas[softselection],c='k')
    ax8.scatter(deltaMS[selection],mu_molgas[selection],c='r')
    ax8.vlines(x=dMScut,ymin=mumolgas_lo,ymax=mumolgas_hi,colors='r',linestyles='--')
    ax8.set_xlim(dMS_lo,dMS_hi)
    ax8.set_xlabel(r'$\Delta$MS [dex]')
    ax8.set_yscale('log')
    ax8.set_ylim(mumolgas_lo,mumolgas_hi)
    ax8.set_ylabel('')
    #ax8.set_ylabel(r'$M_{gas} / M_*$')

    ax9.scatter(kappa_corot[softselection],mu_molgas[softselection],c='k')
    ax9.scatter(kappa_corot[selection],mu_molgas[selection],c='r')
    ax9.vlines(x=kappacut,ymin=mumolgas_lo,ymax=mumolgas_hi,colors='r',linestyles='--')
    ax9.set_xlim(kappa_lo,kappa_hi)
    ax9.set_xlabel(r'$\kappa_{corot,*}$')
    ax9.set_yscale('log')
    ax9.set_ylim(mumolgas_lo,mumolgas_hi)
    ax9.set_ylabel('')
    #ax9.set_ylabel(r'$M_{gas} / M_*$')

    if image=="one":
        candidates = np.argwhere(selection) # finding the ID of our galaxies
        if len(candidates) > 0:
            target = candidates[0][0]
            ax1.scatter(Mstar[target],deltaMS[target], xunits='Msun',c='b',marker='*',s=200)
            ax4.scatter(Mstar[target],kappa_corot[target], xunits='Msun',c='b',marker='*',s=200)
            ax5.scatter(deltaMS[target],kappa_corot[target],c='b',marker='*',s=200)
            ax7.scatter(Mstar[target],mu_molgas[target], xunits='Msun',c='b',marker='*',s=200)
            ax8.scatter(deltaMS[target],mu_molgas[target],c='b',marker='*',s=200)
            ax9.scatter(kappa_corot[target],mu_molgas[target],c='b',marker='*',s=200)

    plt.subplots_adjust(wspace=0.1,hspace=0.1)
    
    imgname = 'plots/selection_plots/triangle_' + str(run) + 'z' + str(z_short) + '.png'

    plt.savefig(imgname,bbox_inches='tight')



if (image == "one") or (image == "all"):
    candidates = np.argwhere(selection)   #finding the ID of our galaxies
    if len(candidates)==0:
        print("no galaxies selected")
        exit()
    if image == "one":
        target = candidates[0]   #list containing one galaxy ID
    elif image == "all":
        target = [candidates[i][0] for i in range(len(candidates))]   #list containing multiple IDs
        print("image=all not yet tested, exiting")
        exit()
    

    if device=="cosma":
        colibre_base_path = "/cosma8/data/dp004/colibre/Runs"
        simulation_dir = run+"/THERMAL_AGN_"+mres
    elif device=="hyades":
        colibre_base_path = "/mnt/su3-pro/colibre/"
        simulation_dir = run_hyades+"/THERMAL_AGN"

    soap_path = "SOAP/halo_properties_"+snap+".hdf5"
    soap_catalogue_file = os.path.join(colibre_base_path, simulation_dir, soap_path)
    virtual_snapshot_path = "SOAP/colibre_with_SOAP_membership_"+snap+".hdf5"
    virtual_snapshot_file = os.path.join(colibre_base_path, simulation_dir, virtual_snapshot_path)


    for ID in target:
        sg = SWIFTGalaxy(
            virtual_snapshot_file,
            SOAP(
                soap_catalogue_file,
                soap_index=ID
            ),
        )

        sg.gas.masses.convert_to_physical()
        sg.gas.coordinates.convert_to_physical()
        sg.gas.smoothing_lengths.convert_to_physical()

        sg.stars.masses.convert_to_physical()
        sg.stars.coordinates.convert_to_physical()
        sg.stars.smoothing_lengths.convert_to_physical()

        disc_radius = 15.0*unyt.kpc
        disc_region = sw.objects.cosmo_array(
                [-1*disc_radius, disc_radius, -1*disc_radius, disc_radius],
                comoving=False,
                scale_factor=sg.metadata.a,
                scale_exponent=1,
        )

        R50 = sg.halo_catalogue.exclusive_sphere_30kpc.half_mass_radius_stars.squeeze()
        R50.convert_to_physical()

        def mkpatch():
            circle = plt.Circle((0,0), R50/unyt.kpc, color='k',fill=False)
            return circle

        

        # min and max log gas and star densities (for plotting purposes)
        gmin = 4.5
        gmax = 9
        smin = 4.5
        smax = 9

        # unrotated
        gas_map = project_gas(
                sg,
                resolution=256,
                project="masses",
                parallel=True,
                periodic=False,
                region=disc_region,
        )
        star_map = project_pixel_grid(
                data=sg.stars,
                resolution=256,
                project="masses",
                parallel=True,
                periodic=False,
                region=disc_region,
        )

        fig,ax = plt.subplots(1,1,figsize=(8,6))
        R50circle = mkpatch()
        mp = ax.imshow(
            np.log10(gas_map.to_value(unyt.solMass / unyt.kpc**2).T),
            cmap="viridis",
            extent=disc_region,
            origin="lower",
            vmin=gmin,
            vmax=gmax,
        )
        ax.set_xlabel(f"x' [{disc_radius.units}]")
        ax.set_ylabel(f"y' [{disc_radius.units}]")
        ax.add_patch(R50circle)
        cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{gas}}{M_\odot kpc^{-2}}$')
        imgname = "plots/galplots/"+run+'_z'+str(z_short)+'_'+str(ID)+'gas.png'
        plt.savefig(imgname, bbox_inches='tight')


        fig,ax = plt.subplots(1,1,figsize=(8,6))
        R50circle = mkpatch()
        mp = ax.imshow(
            np.log10(star_map.to_value(unyt.solMass / unyt.kpc**2).T),
            cmap="magma",
            extent=disc_region,
            origin="lower",
            vmin=smin,
            vmax=smax,
        )
        ax.set_xlabel(f"x' [{disc_radius.units}]")
        ax.set_ylabel(f"y' [{disc_radius.units}]")
        ax.add_patch(R50circle)
        cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{*}}{M_\odot kpc^{-2}}$')
        imgname = "plots/galplots/"+run+'_z'+str(z_short)+'_'+str(ID)+'star.png'
        plt.savefig(imgname, bbox_inches='tight')



        # moving to face-on with stellar angular momentum (see SG Colibre quickstart)
        Lstars = sg.halo_catalogue.exclusive_sphere_30kpc.angular_momentum_stars.squeeze() #10kpc, 30kpc?
        Lstars.convert_to_physical()
        zhat = (Lstars / np.sqrt(np.sum(Lstars**2))).to_value(
            unyt.dimensionless
        )  # we'll align L with the z-axis
        arb = np.ones(3) / np.sqrt(
            3
        )  # we have one degree of freedom, we'll fix it by projecting onto an arbitrarily chosen vector
        xvec = arb - arb.dot(zhat) * zhat
        xhat = xvec / np.sum(xvec**2)
        yhat = np.cross(zhat, xhat)  # orthogonal, right-handed and normalized
        rotmat = np.vstack((xhat, yhat, zhat)).T # transpose!!!!
        sg.rotate(Rotation.from_matrix(rotmat)) # hopefully this puts the galaxy face-on

        


        # rotated
        gas_map = project_gas(
                sg,
                resolution=256,
                project="masses",
                parallel=True,
                periodic=False,
                region=disc_region,
        )
        star_map = project_pixel_grid(
                data=sg.stars,
                resolution=256,
                project="masses",
                parallel=True,
                periodic=False,
                region=disc_region,
        )

        fig,ax = plt.subplots(1,1,figsize=(8,6))
        R50circle = mkpatch()
        mp = ax.imshow(
            np.log10(gas_map.to_value(unyt.solMass / unyt.kpc**2).T),
            cmap="viridis",
            extent=disc_region,
            origin="lower",
            vmin=gmin,
            vmax=gmax,
        )
        ax.set_xlabel(f"x' [{disc_radius.units}]")
        ax.set_ylabel(f"y' [{disc_radius.units}]")
        ax.add_patch(R50circle)
        cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{gas}}{M_\odot kpc^{-2}}$')
        imgname = "plots/galplots/"+run+'_z'+str(z_short)+'_'+str(ID)+'gas_fo.png'
        plt.savefig(imgname, bbox_inches='tight')


        fig,ax = plt.subplots(1,1,figsize=(8,6))
        R50circle = mkpatch()
        mp = ax.imshow(
            np.log10(star_map.to_value(unyt.solMass / unyt.kpc**2).T),
            cmap="magma",
            extent=disc_region,
            origin="lower",
            vmin=smin,
            vmax=smax,
        )
        ax.set_xlabel(f"x' [{disc_radius.units}]")
        ax.set_ylabel(f"y' [{disc_radius.units}]")
        ax.add_patch(R50circle)
        cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{*}}{M_\odot kpc^{-2}}$')
        imgname = "plots/galplots/"+run+'_z'+str(z_short)+'_'+str(ID)+'star_fo.png'
        plt.savefig(imgname, bbox_inches='tight')
