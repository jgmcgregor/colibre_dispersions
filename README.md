# Measuring turbulence in COLIBRE disk galaxies

This is a relatively hodge-podge collection of notebooks and scripts I've written as part of my first PhD project. I'll neaten it up once I've built it out a bit more.
*colibre_requirements.txt* contains a non-exhaustive list of the different Python package dependencies, which can be used to set up an environment. (I may have missed a few...)

## Important scripts
* *selection.py*: takes a subhalo catalogue and selects target galaxies based on their stellar mass, distance from the main sequence, and 'diskiness', produces tables of target galaxy IDs, some simple plots, and can also produce simple images of the galaxies.
* *multiselect_cosma.sh*, *multiselect_hyades.sh*: sbatch scripts which run the selection script across multiple redshifts. Different redshifts run in parallel, but it's pretty quick regardless.
* *mergercheck.py* [UNDER CONSTRUCTION]: takes targeted galaxies and checks if they've had any major or minor mergers recently. Haven't decided whether it runs redshifts or galaxies in parallel.
* *multimerger.sh* [UNDER CONSTRUCTION]
* *turbulence.py* [UNDER CONSTRUCTION]: calculates velocity dispersion, disk thickness, and any other quantities I feel like calculating for the targeted galaxies. Haven't decided whether it runs redshifts or galaxies in parallel.
* *multiturb.sh* [UNDER CONSTRUCTION]

## Other junk
* *______.ipynb*: notebooks I've been messing around in. I'll move them to *testing* once I don't need them anymore



:)
