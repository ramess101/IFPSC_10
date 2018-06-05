#!/bin/bash
# This small script automates rerunning the Green Kubo analysis for
# however many rhos are needed. 
# Specify the output path by hand.
Compound=C4H10
rho_low=0  # First rho to be rerun
rho_high=0  # Second rho to be rerun
low_lim=0  # Low limit of number of replicates to pass to GK
high_lim=3  # High limit of number of replicates to pass to GK
scripts_path=/Scripts  # Where to find GK
output_path="$PWD"  # Where to find the data

cd "$output_path" || exit 1  # Enter analysis directory or quit

for ((rho="$rho_low"; rho <= "$rho_high"; rho++))
do
python "$scripts_path"/GreenKubo_analyze.py --ilow "$low_lim" --ihigh "$high_lim" --nReps 1 --irho "$rho" --sat
# Rerun for all MCMC parameter sets
done
