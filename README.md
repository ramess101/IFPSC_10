# IFPSC_10

Directory description:

Scripts- Simulation and post-simulation analysis files

IOCTANE_REFPROP- REFPROP data for isooctane

Ethane- Ethane simulation input files and results

C3H8- Propane simulation input files and results


Currently the scripts in Ethane are for simulating 200 different epsilon/sigma parameter sets using a Mie lambda-6 potential. These parameter sets are determined using Markov Chain Monte Carlo.

The Green-Kubo analysis uses bootstrapping to assign uncertainties. These uncertainties account for both the fitting of the Green-Kubo integral as well as the uncertainty in the force field parameters.


nAlkanes_BayesianSaturatedViscosity is the script that should be used to simulate all n-alkanes larger than ethane. 

eps_sig_lam_MCMC are the epsilon, sigma, lambda sets for each interaction site (CH3, CH2, CH, C) as determined with a Bayesian Markov Chain Monte Carlo approach.


