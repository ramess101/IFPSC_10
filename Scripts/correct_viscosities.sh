#!/bin/bash
####
#
# This code exists to correct viscosity determined by previous programs
# which were using buggy code (ie wrong volume was passed to the gmx 
# energy program)
Compound=IC4H10
Model=TAMie
Conditions_type=T293highP # ie T293highP
BondType=LINCS  #harmonic (flexible) or LINCS (fixed)
jlim=4  # number of condition sets to run
batches=3  # Number of batches to run
NREPS=20 #Run NREPS replicates in parallel/# in a batch
#Set the number of molecules
Nmol=400

#Specify the path location for files
scripts_path=~/Scripts
conditions_path=~/"$Conditions_type"_Conditions  
output_path=~/"$Compound"/Gromacs/"$Conditions_type"_Viscosity/"$Model"_N"$Nmol"_"$BondType"

jlim=$((jlim - 1)) # Zero indexed
total=$((batches * NREPS - 1))  # Total number of repeats (0 indexed)

### Read the box size, that depends on number of molecules
while read line
do 
liquid_box+=("$line")
done < "$conditions_path"/"$Compound"_liquid_box_N"$Nmol"

# Begin actual correction work
for j in $(seq 4 4)
do
# For every MCMC that exists in this j
for iMCMC in $(seq 0 $total)
do
cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep0/NVT_eq/NVT_prod/NVT_vis || exit 1

# Test for the presence of the NPT file and warn if it isn't available
if [ -f ../../../NPT_eq/NPT_prod/Lbox_NPT_ave ]
then
Lbox=$(<../../../NPT_eq/NPT_prod/Lbox_NPT_ave)
wrong_box=${liquid_box[j]}
else
echo "Warning: Lbox_NPT_ave file not found for j $j iMCMC $iMCMC. No correction applied"
echo "Warning: Lbox_NPT_ave file not found for j $j iMCMC $iMCMC. No correction applied" >> "$output_path"/correction_warnings.txt
Lbox=${liquid_box[j]}
wrong_box=${liquid_box[j]}
fi

Vol=$(echo $Lbox|awk '{print $1*$1*$1}')
wrong_vol=$(echo $wrong_box|awk '{print $1*$1*$1}')
echo "Right box $Lbox wrong box $wrong_box right vol $Vol wrong vol $wrong_vol for iMCMC $iMCMC j $j"

# Now that we have those in hand... do the corrections
# first off, we want to make sure we don't run this script over again
if [ -f wrong_visco.xvg ]
then
echo "WARNING: wrong_visco.xvg exists in iMCMC $iMCMC j $j : has this script been run before?"
echo "Type N to continue, anything else to exit."
read -t 300 cont  # Times out after five minutes
if [ "$cont" = "N" ]
then
echo "Continuing"
# Do nothing. Continue
else
exit 0  # Abort due to duplicate file
fi
fi

mv visco.xvg wrong_visco.xvg
# Back up the file and run the correction
awk -f "$scripts_path"/correct.awk -v first="$Vol" -v second="$wrong_vol" < wrong_visco.xvg > visco.xvg

# Append this informaiton in case it becomes needed
echo "# Right box $Lbox wrong box $wrong_box right vol $Vol wrong vol $wrong_vol for iMCMC $iMCMC j $j" >> wrong_visco.xvg
done

cd "$output_path" || exit 1
###GreenKubo_analyze for all MCMC parameter sets
python "$scripts_path"/GreenKubo_analyze.py --ilow 0 --ihigh $total --nReps 1 --irho "$j" --sat

done  # For each j

cd "$output_path" || exit
###Create plots to compare with experimental data and correlations
if [ "$Conditions_type" = "Saturation" ]
then
python "$scripts_path"/compare_TDE_REFPROP.py --comp "$Compound" --nrhomax $((j+1)) --sat
else
python "$scripts_path"/compare_TDE_REFPROP.py --comp "$Compound" --nrhomax $((j+1)) --"$Conditions_type"
fi

exit 0

#######

