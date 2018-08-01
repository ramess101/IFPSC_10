#!/bin/bash
####
#
# This code exists to correct viscosity determined by previous programs
# which were using buggy code (ie wrong volume was passed to the gmx 
# energy program)
Compound=3MPentane
Model=Potoff
Conditions_type=Saturation # ie T293highP
BondType=LINCS  #harmonic (flexible) or LINCS (fixed)
jlim=3  # number of condition sets to run
batches=3  # Number of batches to run
NREPS=20 #Run NREPS replicates in parallel/# in a batch
#Set the number of molecules
Nmol=400

#Specify the path location for files
scripts_path=~/Scripts
conditions_path=~/"$Conditions_type"_Conditions  
output_path="$PWD"
#output_path=~/"$Compound"/Gromacs/"$Conditions_type"_Viscosity/"$Model"_N"$Nmol"_"$BondType"

jlim=$((jlim - 1)) # Zero indexed
total=$((batches * NREPS - 1))  # Total number of repeats (0 indexed)

# Begin actual work
for j in $(seq 0 0)
do
# For every MCMC that exists in this j
for iMCMC in $(seq 0 $total)
do
cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep0/NPT_eq/NPT_prod || exit 1

if [ -f Lbox_NPT_ave ]
then
Lbox=$(<Lbox_NPT_ave)
echo "$Lbox" >> "$output_path"/Lbox_all
else
echo "Warning: LBOX not found for j $j iMCMC $iMCMC"
fi

done  # For each j

done

exit 0

#######

