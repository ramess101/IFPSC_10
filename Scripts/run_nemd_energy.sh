#!/bin/bash
# I don't know why this happened.
# Don't ask supercomputers why things happen.
Compound=C4H10
Model=TraPPE
Conditions_type=T293highP # ie T293highP
BondType=LINCS  #harmonic (flexible) or LINCS (fixed)
jlim=1  # number of condition sets to run
batches=1  # Number of batches to run
NREPS=10 #Run NREPS replicates in parallel/# in a batch
#Set the number of molecules
Nmol=400

#Specify the path location for files
scripts_path=~/Scripts
conditions_path=~/"$Conditions_type"_Conditions  
output_path=~/"$Compound"/Gromacs/"$Conditions_type"_Viscosity/"$Model"_N"$Nmol"_"$BondType"_NEMD

jlim=$((jlim - 1)) # Zero indexed
total=$((batches * NREPS - 1))  # Total number of repeats (0 indexed)
nRep=0

for j in $(seq 0 $jlim)
do
# For every MCMC that exists in this j
for iMCMC in $(seq 0 $total)
do
cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep"/NVT_eq/NVT_prod/NVT_vis || error_report "Unable to change to directory NVT_vis" "$j" "$iMCMC" "NVT viscosity"  #start fresh for do cycle instead of multiple "cd .."'s
echo 34 | gmx energy -f nvt_vis.edr -s nvt_vis.tpr > vis_out 2>> vis_out   # Get 1/NEMD_VISCOSITY
awk -f "$scripts_path"/vis_arrange.awk < energy.xvg > visco_temp.xvg  # Convert to NEMD_VISCOSITY
awk -f "$scripts_path"/runavg.awk < visco_temp.xvg > visco.xvg  # Perform running average
vis_avg=$(<visco_avg.txt)  # Read in average
echo "${cos_acc[iMCMC]}	$vis_avg" >> "$output_path"/cosacc_vs_vis.txt
rm visco_temp.xvg  # Remove bulky files
#rm energy.xvg
done
done  # For each j
