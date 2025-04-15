# Modified from Kyle Oman's SWIFTGalaxy_Colibre_QuickStart.ipynb

import os
import numpy as np
import unyt
import swiftsimio as sw
from swiftgalaxy import SWIFTGalaxy, SOAP
from swiftsimio.visualisation.projection import project_gas, project_pixel_grid
from swiftsimio.visualisation import generate_smoothing_lengths
import matplotlib.pyplot as plt
import matplotlib.colors
from scipy.spatial.transform import Rotation

colibre_base_path = "/cosma8/data/dp004/colibre/Runs"
simulation_dir = "L012_m6/THERMAL_AGN_m6"
soap_catalogue_file = os.path.join(colibre_base_path, simulation_dir, "SOAP/halo_properties_0127.hdf5")
virtual_snapshot_file = os.path.join(colibre_base_path, simulation_dir, "SOAP/colibre_with_SOAP_membership_0127.hdf5")

chosen_halo_index = 13797 #11956 # chosen with my selection fns - can automate this better later

sg = SWIFTGalaxy(
    virtual_snapshot_file,
    SOAP(
        soap_catalogue_file,
        soap_index=chosen_halo_index
    ),
)

#gas_coord = sg.gas.coordinates

#fig,ax = plt.subplots(figsize=(10,10))
#ax.scatter(gas_coord[:,0],gas_coord[:,1],s=0.1,c='k')
#plt.savefig('plots/galtest_xy.png',bbox_inches='tight')

#fig,ax = plt.subplots(figsize=(10,10))
#ax.scatter(gas_coord[:,0],gas_coord[:,2],s=0.1,c='k')
#plt.savefig('plots/galtest_xz.png',bbox_inches='tight')

#fig,ax = plt.subplots(figsize=(10,10))
#ax.scatter(gas_coord[:,1],gas_coord[:,2],s=0.1,c='k')
#plt.savefig('plots/galtest_yz.png',bbox_inches='tight')

#disc_radius = 15*unyt.kpc
#halo_radius = 200*unyt.kpc

if not hasattr(sg.dark_matter, "smoothing_lengths"):
    sg.dark_matter.smoothing_lengths = generate_smoothing_lengths(
            (sg.dark_matter.coordinates + sg.centre) % sg.metadata.boxsize,
            sg.metadata.boxsize,
            kernel_gamma=1.8,
            neighbours=57,
            speedup_fac=2,
            dimension=3,
     )

#disc_region = [-disc_radius, disc_radius, -disc_radius, disc_radius]
#halo_region = [-halo_radius, halo_radius, -halo_radius, halo_radius]

disc_radius = 15.0 #kpc
halo_radius = 200.0 #kpc


disc_region = sw.objects.cosmo_array(
        [-1*disc_radius, disc_radius, -1*disc_radius, disc_radius],
        unyt.kpc,
        comoving=False,
        scale_factor=sg.metadata.a,
        scale_exponent=1,
)

halo_region = sw.objects.cosmo_array(
        [-1*halo_radius, halo_radius, -1*halo_radius, halo_radius],
        unyt.kpc,
        comoving=False,
        scale_factor=sg.metadata.a,
        scale_exponent=1,
)

gas_map = project_gas(
        sg,
        resolution=256,
        project="masses",
        parallel=True,
        region=disc_region,
)

gas_map.convert_to_units(unyt.msun / unyt.kpc**2)

star_map = project_pixel_grid(
        data=sg.stars,
        resolution=256,
        project="masses",
        parallel=True,
        region=disc_region,
)

star_map.convert_to_units(unyt.msun / unyt.kpc**2)

#fig,ax = plt.subplots(1,1,figsize=(6,6))
#ax.imshow(matplotlib.colors.LogNorm()(gas_map.value), cmap="viridis", extent=disc_region)
#imgname = "plots/galplots/L012_m6_example2_gas.png"
#plt.savefig(imgname, bbox_inches='tight')

fig,ax = plt.subplots(1,1,figsize=(8,6))
#mp = ax.imshow(matplotlib.colors.LogNorm()(gas_map.value), cmap="viridis", extent=disc_region)
mp = ax.imshow(np.log10(gas_map.value), cmap="viridis", extent=disc_region)
cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{gas}}{M_\odot kpc^{-2}}$')
imgname = "plots/galplots/L012_m6_example2_gas_cb.png"
plt.savefig(imgname, bbox_inches='tight')


#fig,ax = plt.subplots(1,1,figsize=(6,6))
#ax.imshow(matplotlib.colors.LogNorm()(star_map), cmap="magma", extent=disc_region)
#imgname = "plots/galplots/L012_m6_example2_star.png"
#plt.savefig(imgname, bbox_inches='tight')

fig,ax = plt.subplots(1,1,figsize=(8,6))
#mp = ax.imshow(matplotlib.colors.LogNorm()(gas_map.value), cmap="viridis", extent=disc_region)
mp = ax.imshow(np.log10(star_map.value), cmap="magma", extent=disc_region)
cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{*}}{M_\odot kpc^{-2}}$')
imgname = "plots/galplots/L012_m6_example2_star_cb.png"
plt.savefig(imgname, bbox_inches='tight')

Lstars = sg.halo_catalogue.exclusive_sphere_30kpc.angular_momentum_stars.squeeze()
zhat = (Lstars / np.sqrt(np.sum(Lstars**2))).to_value(
    unyt.dimensionless
)  # we'll align L with the z-axis
arb = np.ones(3) / np.sqrt(
    3
)  # we have one degree of freedom, we'll fix it by projecting onto an arbitrarily chosen vector
xvec = arb - arb.dot(zhat) * zhat
xhat = xvec / np.sum(xvec**2)
yhat = np.cross(zhat, xhat)  # orthogonal, right-handed and normalized
rotmat = np.vstack((xhat, yhat, zhat))
sg.rotate(Rotation.from_matrix(rotmat)) # hopefully this puts the galaxy face-on

gas_map = project_gas(
        sg,
        resolution=256,
        project="masses",
        parallel=True,
        region=disc_region,
)
gas_map.convert_to_units(unyt.msun / unyt.kpc**2)

star_map = project_pixel_grid(
        data=sg.stars,
        resolution=256,
        project="masses",
        parallel=True,
        region=disc_region,
)
star_map.convert_to_units(unyt.msun / unyt.kpc**2)

#fig,ax = plt.subplots(1,1,figsize=(6,6))
#ax.imshow(matplotlib.colors.LogNorm()(gas_map.value), cmap="viridis", extent=disc_region)
#imgname = "plots/galplots/L012_m6_example2_gas_fo.png"
#plt.savefig(imgname, bbox_inches='tight')

fig,ax = plt.subplots(1,1,figsize=(8,6))
#mp = ax.imshow(matplotlib.colors.LogNorm()(gas_map.value), cmap="viridis", extent=disc_region)
mp = ax.imshow(np.log10(gas_map.value), cmap="viridis", extent=disc_region)
cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{gas}}{M_\odot kpc^{-2}}$')
imgname = "plots/galplots/L012_m6_example2_gas_focb.png"
plt.savefig(imgname, bbox_inches='tight')


#fig,ax = plt.subplots(1,1,figsize=(6,6))
#ax.imshow(matplotlib.colors.LogNorm()(star_map), cmap="magma", extent=disc_region)
#imgname = "plots/galplots/L012_m6_example2_star_fo.png"
#plt.savefig(imgname, bbox_inches='tight')

fig,ax = plt.subplots(1,1,figsize=(8,6))
#mp = ax.imshow(matplotlib.colors.LogNorm()(gas_map.value), cmap="viridis", extent=disc_region)
mp = ax.imshow(np.log10(star_map.value), cmap="magma", extent=disc_region)
cb = fig.colorbar(mp, ax=ax,label=r'$\log \frac{\Sigma_{*}}{M_\odot kpc^{-2}}$')
imgname = "plots/galplots/L012_m6_example2_star_focb.png"
plt.savefig(imgname, bbox_inches='tight')


