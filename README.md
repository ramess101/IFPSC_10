# IFPSC_10

Directory description:

Scripts- Simulation and post-simulation analysis files

IOCTANE_REFPROP- REFPROP data for isooctane

Ethane- Ethane simulation input files and results


Currently the scripts in Ethane are for simulating 200 different epsilon/sigma parameter sets using a Mie lambda-6 potential. These parameter sets are determined using Markov Chain Monte Carlo.

The Green-Kubo analysis uses bootstrapping to assign uncertainties. These uncertainties account for both the fitting of the Green-Kubo integral as well as the uncertainty in the force field parameters.
