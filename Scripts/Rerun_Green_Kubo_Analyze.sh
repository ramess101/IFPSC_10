#!/bin/bash
# This small script automates rerunning the Green Kubo analysis for
# however many rhos are needed. 
# Specify the output path by hand.
Compound=IC8H18
Model=TAMie
Conditions_type=Saturation
Nmol=400
rho_low=0  # First rho to be rerun
rho_high=4  # Last rho to be rerun
low_lim=0  # Low limit of number of replicates to pass to GK
high_lim=59  # High limit of number of replicates to pass to GK
scripts_path=~/Scripts  # Where to find GreenKubo_analyze
# Where to find the data
output_path=~/"$Compound"/Gromacs/"$Conditions_type"_Viscosity/"$Model"_N"$Nmol" #
#output_path="$PWD" 
#output_path=~/"$Compound"/Gromacs/"$Model"_N"$Nmol"
#output_path=~/"$Compound"/Gromacs/"$Model"

cd "$output_path" || exit 1  # Enter analysis directory or quit

for ((rho="$rho_low"; rho <= "$rho_high"; rho++))
do
python "$scripts_path"/GreenKubo_analyze.py --ilow "$low_lim" --ihigh "$high_lim" --nReps 1 --irho "$rho" --sat
# Rerun for all MCMC parameter sets
done
