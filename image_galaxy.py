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

colibre_base_path = "/cosma8/data/dp004/colibre/Runs"
simulation_dir = "L012_m6/THERMAL_AGN_m6"
soap_catalogue_file = os.path.join(colibre_base_path, simulation_dir, "SOAP/halo_properties_0127.hdf5")
virtual_snapshot_file = os.path.join(colibre_base_path, simulation_dir, "SOAP/colibre_with_SOAP_membership_0127.hdf5")

chosen_halo_index = 11956 # chosen with my selection fns - can automate this better later

sg = SWIFTGalaxy(
    virtual_snapshot_file,
    SOAP(
        soap_catalogue_file,
        soap_index=chosen_halo_index
    ),
)




def myvis(sg, disc_radius=15 * unyt.kpc, halo_radius=200 * unyt.kpc, fignum=1):
    # Before visualising the galaxy, we need to initialise some smoothing
    # lengths for the dark matter particles (this bit of code is taken
    # from the `swiftsimio` visualisation documentation).
    if not hasattr(sg.dark_matter, "smoothing_lengths"):
        sg.dark_matter.smoothing_lengths = generate_smoothing_lengths(
            (sg.dark_matter.coordinates + sg.centre) % sg.metadata.boxsize,
            sg.metadata.boxsize,
            kernel_gamma=1.8,
            neighbours=57,
            speedup_fac=2,
            dimension=3,
        )
    disc_region = [-disc_radius, disc_radius, -disc_radius, disc_radius]
    halo_region = [-halo_radius, halo_radius, -halo_radius, halo_radius]
    gas_map = project_gas(
        sg,
        resolution=256,
        project="masses",
        parallel=True,
        region=disc_region,
    )
    dm_map = project_pixel_grid(
        data=sg.dark_matter,
        boxsize=sg.metadata.boxsize,
        resolution=256,
        project="masses",
        parallel=True,
        region=halo_region,
    )
    star_map = project_pixel_grid(
        data=sg.stars,
        boxsize=sg.metadata.boxsize,
        resolution=256,
        project="masses",
        parallel=True,
        region=disc_region,
    )

    fig = plt.figure(fignum, figsize=(10, 3))
    sp1, sp2, sp3 = [fig.add_subplot(1, 3, i) for i in range(1, 4)]
    sp1.imshow(
        matplotlib.colors.LogNorm()(gas_map.value), cmap="viridis", extent=disc_region
    )
    sp1.set_xlabel(f"x' [{disc_radius.units}]")
    sp1.set_ylabel(f"y' [{disc_radius.units}]")
    sp1.text(
        0.9, 0.9, "gas", color="white", ha="right", va="top", transform=sp1.transAxes
    )
    sp2.imshow(
        matplotlib.colors.LogNorm()(dm_map),
        cmap="inferno",
        extent=halo_region,
    )
    sp2.plot(
        [-disc_radius, -disc_radius, disc_radius, disc_radius, -disc_radius],
        [-disc_radius, disc_radius, disc_radius, -disc_radius, -disc_radius],
        "-k",
    )
    sp2.set_xlabel(f"x' [{halo_radius.units}]")
    sp2.set_ylabel(f"y' [{halo_radius.units}]")
    sp2.text(
        0.9, 0.9, "DM", ha="right", va="top", color="white", transform=sp2.transAxes
    )
    sp3.imshow(
        matplotlib.colors.LogNorm()(star_map),
        cmap="magma",
        extent=disc_region,
    )
    sp3.set_xlabel(f"x' [{disc_radius.units}]")
    sp3.set_ylabel(f"y' [{disc_radius.units}]")
    sp3.text(0.9, 0.9, "stars", ha="right", va="top", transform=sp3.transAxes)
    sp2.set_title(f"soap_index={sg.halo_catalogue.soap_index}")
    fig.subplots_adjust(wspace=0.4)
    return fig

