# Fun little script to produce cool images for a selected COLIBRE galaxy
# JG McGregor
# May 2025 

import os
import numpy as np
import argparse
import unyt
import swiftsimio as sw
from swiftgalaxy import SWIFTGalaxy, SOAP
from swiftsimio.visualisation.projection import project_gas, project_pixel_grid
from swiftsimio.visualisation import generate_smoothing_lengths
import matplotlib.pyplot as plt
import matplotlib.colors
from scipy.spatial.transform import Rotation

parser = argparse.ArgumentParser()
parser.add_argument("device", help="local, cosma, or hyades")
parser.add_argument("length", help="length, LXXX")
parser.add_argument("mres", help="mass res, mX")
parser.add_argument("snap", help="snapshot number", type=int)
parser.add_argument("galaxyID", help="SOAP index of galaxy", type=int)
args = parser.parse_args()

device = args.device
length = args.length
mres = args.mres
snap = args.snap
run = length+"_"+mres
ID = args.galaxyID

snap4 = str(snap).zfill(4) # adds leading 0s to make a length 4 string

if device == "cosma":
    colibre_base_path = "/cosma8/data/dp004/colibre"
    simulation_dir = "Runs/"+run+"/THERMAL_AGN_"+mres
    soap_path = "SOAP/halo_properties_"+snap4+".hdf5"
    soap_catalogue_file = os.path.join(colibre_base_path, simulation_dir, soap_path)
    virtual_snapshot_path = "SOAP/colibre_with_SOAP_membership_"+snap4+".hdf5"
    virtual_snapshot_file = os.path.join(colibre_base_path, simulation_dir, virtual_snapshot_path)
else:
    print("device not supported")
    exit()

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

R50 = sg.halo_catalogue.exclusive_sphere_30kpc.half_mass_radius_stars.squeeze()
R50.convert_to_physical()
R50.convert_to_units('kpc')
R50_gas = sg.halo_catalogue.exclusive_sphere_30kpc.half_mass_radius_gas.squeeze()
R50_gas.convert_to_physical()
R50_gas.convert_to_units('kpc')
def mkpatch():
    circle1 = plt.Circle((0,0), R50/unyt.kpc, color='k',fill=False)
    circle2 = plt.Circle((0,0), R50_gas/unyt.kpc, color='w',fill=False)
    return circle1, circle2

img_region = sw.objects.cosmo_array(
    [-3*R50, 3*R50, -3*R50, 3*R50],
    comoving=False,
    scale_factor=sg.metadata.a,
    scale_exponent=1,
)

# min and max log gas and star densities (for plotting purposes)
gmin = 4.5
gmax = 9
smin = 4.5
smax = 9

#face-on
Lstars = sg.halo_catalogue.exclusive_sphere_30kpc.angular_momentum_stars.squeeze()
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

# unrotated
gas_map = project_gas(
    sg,
    resolution=256,
    project="masses",
    parallel=True,
    periodic=False,
    region=img_region,
)

star_map = project_pixel_grid(
    sg.stars,
    resolution=256,
    project="masses",
    parallel=True,
    periodic=False,
    region=img_region,
)

fig,ax = plt.subplots(1,1,figsize=(8,6))
circle1, circle2 = mkpatch()
mp = ax.imshow(
    np.log10(gas_map.to_value(unyt.solMass / unyt.kpc**2).T),
    cmap="viridis",
    extent=img_region,
    origin="lower",
    vmin=gmin,
    vmax=gmax,
)
ax.set_xlabel(f"x' [{R50.units}]")
ax.set_ylabel(f"y' [{R50.units}]")
ax.add_patch(circle1)
ax.add_patch(circle2)
cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{gas}}{M_\odot kpc^{-2}}$')
imgname = "plots/extra_plots/"+run+'_z'+snap4+'_'+str(ID)+'gasfo.png'
plt.savefig(imgname, bbox_inches='tight')

fig,ax = plt.subplots(1,1,figsize=(8,6))
circle1, circle2 = mkpatch()
mp = ax.imshow(
    np.log10(star_map.to_value(unyt.solMass / unyt.kpc**2).T),
    cmap="magma",
    extent=img_region,
    origin="lower",
    vmin=smin,
    vmax=smax,
)
ax.set_xlabel(f"x' [{R50.units}]")
ax.set_ylabel(f"y' [{R50.units}]")
ax.add_patch(circle1)
ax.add_patch(circle2)
cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{*}}{M_\odot kpc^{-2}}$')
imgname = "plots/extra_plots/"+run+'_z'+snap4+'_'+str(ID)+'starfo.png'
plt.savefig(imgname, bbox_inches='tight')


#edge=on
rotmat = np.vstack(([1,0,0], [0,0,-1], [0,1,0]))
sg.rotate(Rotation.from_matrix(rotmat))

gas_map = project_gas(
    sg,
    resolution=256,
    project="masses",
    parallel=True,
    periodic=False,
    region=img_region,
)

star_map = project_pixel_grid(
    sg.stars,
    resolution=256,
    project="masses",
    parallel=True,
    periodic=False,
    region=img_region,
)

fig,ax = plt.subplots(1,1,figsize=(8,6))
circle1, circle2 = mkpatch()
mp = ax.imshow(
    np.log10(gas_map.to_value(unyt.solMass / unyt.kpc**2).T),
    cmap="viridis",
    extent=img_region,
    origin="lower",
    vmin=gmin,
    vmax=gmax,
)
ax.set_xlabel(f"x' [{R50.units}]")
ax.set_ylabel(f"y' [{R50.units}]")
ax.add_patch(circle1)
ax.add_patch(circle2)
cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{gas}}{M_\odot kpc^{-2}}$')
imgname = "plots/extra_plots/"+run+'_z'+snap4+'_'+str(ID)+'gaseo.png'
plt.savefig(imgname, bbox_inches='tight')

fig,ax = plt.subplots(1,1,figsize=(8,6))
circle1, circle2 = mkpatch()
mp = ax.imshow(
    np.log10(star_map.to_value(unyt.solMass / unyt.kpc**2).T),
    cmap="magma",
    extent=img_region,
    origin="lower",
    vmin=smin,
    vmax=smax,
)
ax.set_xlabel(f"x' [{R50.units}]")
ax.set_ylabel(f"y' [{R50.units}]")
ax.add_patch(circle1)
ax.add_patch(circle2)
cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{*}}{M_\odot kpc^{-2}}$')
imgname = "plots/extra_plots/"+run+'_z'+snap4+'_'+str(ID)+'stareo.png'
plt.savefig(imgname, bbox_inches='tight')



