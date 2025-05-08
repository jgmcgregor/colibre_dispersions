# Measuring turbulence in COLIBRE disk galaxies

This is a hodge-podge collection of notebooks and scripts I've written as part of my first PhD project. I'll neaten it up once I've built it out a bit more.

## Important scripts
* *selection.py*: takes a subhalo catalogue and selects target galaxies based on their stellar mass, distance from the main sequence, and 'diskiness', produces tables of target galaxy IDs, some simple plots, and can also produce simple images of the galaxies.
  * *multiselect_cosma.sh*, *multiselect_hyades.sh*: sbatch scripts which run the selection script across multiple redshifts. Different redshifts run in parallel, but it's pretty quick regardless.
* *mergercheck.py*: takes a list of galaxies and records how many mergers (above a given mass ratio) each galaxy has had recently.
  * *multimerger_cosma.sh*: sbatch script for the merger checking. Each redshift looks at different HBT files, so I run them in parallel. No Hyades script, as that computer doesn't have the HBT data.
* *turbulence.py* [UNDER CONSTRUCTION]: calculates velocity dispersion, disk thickness, and any other quantities I feel like calculating for the targeted galaxies.
  * *multiturb.sh* [UNDER CONSTRUCTION]: sbatch script for *turbulence.py*. It may be worth parallelising by galaxy, since the particle data for each will be opened separately.

## Other junk
* *colibre_requirements.txt* contains a non-exhaustive list of the different Python package dependencies, which can be used to set up an environment. (I may have missed a few...)
* *______.ipynb*: notebooks I've been messing around in. I'll move them to *testing* once I don't need them anymore



:)
