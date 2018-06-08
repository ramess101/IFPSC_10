#!/bin/bash
####
#
#This code submits a set of CH3 and CH2 parameter sets (determined by Bayesian MCMC with a different code)
#Simulations are performed at 5 REFPROP saturation conditions
#Temperatures are performed sequentially, batches of NREPS parameter sets are submitted at a time
#Green-Kubo analysis is used after each temperature loop is completed
#
#
####
clean() {   # Make it so everything is killed on an interrupt
local pids=$(jobs -pr)
echo "On exit sending kill signal to: $pids"
[ -n "$pids" ] && kill $pids
exit 1
}
trap "clean" SIGINT SIGTERM EXIT SIGQUIT  # Call cleanup when asked to

job_date=$(date "+%Y_%m_%d_%H_%M_%S")

Compound=C4H10
Model=TraPPE
Conditions_type=T293highP
BondType=harmonic  #Harmonic (glexible) or LINCS (fixed)
Temp=293
jlim=2  # number of condition sets to run
batches=1  # Number of batches to run
NREPS=2 #Run NREPS replicates in parallel/# in a batch
pin0=28  # Default pinoffset
nt_eq=1  # Thread number during equilibration
nt_vis=1  # This thread number will serve in production and viscosity runs
NPT=YES  # YES indicates NPT runs should be carried out prior to NVT runs
#Set the number of molecules
Nmol=400

# Runtiming control: override default runtimes=YES, otherwise use NO
OVERRIDE_STEPS=YES
# If OVERRIDE_STEPS=YES, arrays must be of length j 
# indicating how many equilibration and production steps,
# respectively, to perform for each j
equil_steps=(100 500 500 500 500)
prod_steps=(500 100 500 500 500)


#Specify the path location

scripts_path=~/Scripts
conditions_path=~/"$Conditions_type"_Conditions  
exp_data_path=~/TDE_REFPROP_values
input_path=~/"$Compound"/Gromacs/Gromacs_input
output_path=~/"$Compound"/Gromacs/"$Conditions_type"_Viscosity/"$Model"_N"$Nmol"_"$BondType" #Code assumes that a folder already exists with this directory and the eps_sig_lam_MCMC file in it

jobfile="$output_path"/"$Compound"_job_"$job_date" 
cp "$scripts_path"/AlkanesViscosity.sh "$jobfile" #Keep track of what jobs have been submitted
cat "$scripts_path"/nAlkanesNPTsteps >> "$jobfile"  # Inellegant append
cat "$scripts_path"/run_single.sh >> "$jobfile"  # Another inellegant append
touch "$output_path"/warnings_job_"$job_date"

cd "$output_path"

touch press_all_log
touch Lbox_all

### Read the box size, that depends on number of molecules

while read line
do 
liquid_box+=("$line")
done < "$conditions_path"/"$Compound"_liquid_box_N"$Nmol"

### Read the pressure file into the array 'press'
while read line
do
press+=("$line")
done < "$conditions_path"/"$Compound"_press
echo Pressures ${press[@]} read from "$conditions_path"/"$Compound"_press

### Determine Temperature. In some cases temperature should be read from a file
# IF a temperature file exists, it will use that 
echo Looking for temperatures at "$conditions_path"/"$Compound"_Temp
if [ -e "$conditions_path"/"$Compound"_Temp ]
then
while read line
do 
temps+=("$line")
done < "$conditions_path"/"$Compound"_Temp

# Otherwise, fill the array with copies of the given temperature
# Make the array the same length as the pressure
echo Temperatures ${temps[@]} read from "$conditions_path"/"$Compound"_Temp
else
for ((x=0; x < ${#press[@]}; x++))  # C style loop 
do
temps+=("$Temp")
done
echo Temperatures ${temps[@]} from default temperature used
fi

###Equilibration time is different for certain compounds
# 0.002 ps time step

if [ "$Compound" = 'C3H8' ]
then

nsteps_eq=50000 # 0.1 ns total

elif [ "$Compound" = 'C4H10' ]
then

nsteps_eq=100000 # 0.2 ns total

else

nsteps_eq=500000 # 1 ns total

fi

# We're going to ignore normal step sizes
if [ "$OVERRIDE_STEPS" = YES ]
then
echo Overriding step numbers using equil=${equil_steps[*]}, prod=${prod_steps[*]}
# Use predefined step sizes
else
for ((x=0; x < jlim; x++))
do
equil_steps[$x]=$nsteps_eq
prod_steps[$x]=500000
echo Step numbers : equil=${equil_steps[*]}, prod=${prod_steps[*]}
done
fi



###Cut-off distance is different for some force fields
if [ "$Model" = 'Potoff' ]
then

rvdw=1.0 # [nm]

else

rvdw=1.4 # [nm]

fi

### Assign force field parameters
if [ "$Model" = 'TraPPE' ]
then

#TraPPE

bondlength_CH3=0.1540

epsCH3_sim=98. # (K)
sigCH3_sim=0.375 # (nm)
lamCH3_sim=12.0

epsCH2_sim=46. # (K)
sigCH2_sim=0.395 # (nm)
lamCH2_sim=12.0

epsCH_sim=10. # (K)
sigCH_sim=0.468 # (nm)
lamCH_sim=12.0

epsC_sim=0.5 # (K)
sigC_sim=0.640 # (nm)
lamC_sim=12.0

lam_sim=12.0

elif [ "$Model" = 'TraPPEAUA' ]
then

#TraPPE_AUA

bondlength_CH3=0.192

epsCH3_sim=134.5 # (K)
sigCH3_sim=0.352 # (nm)
lamCH3_sim=12.0

#Borrow the rest from TraPPE-UA

epsCH2_sim=46. # (K)
sigCH2_sim=0.395 # (nm)
lamCH2_sim=12.0

epsCH_sim=10. # (K)
sigCH_sim=0.468 # (nm)
lamCH_sim=12.0

epsC_sim=0.5 # (K)
sigC_sim=0.640 # (nm)
lamC_sim=12.0

lam_sim=12.0

elif [ "$Model" = 'Potoff' ]
then

#Potoff

bondlength_CH3=0.1540

epsCH3_sim=121.25 # (K)
sigCH3_sim=0.3783 # (nm)
lamCH3_sim=16.0

epsCH2_sim=61. # (K)
sigCH2_sim=0.399 # (nm)
lamCH2_sim=16.0

# Transferable sites
#epsCH_sim=15.00 # (K)
#sigCH_sim=0.460 # (nm)
#lamCH_sim=16.0

#epsC_sim=0.5 # (K)
#sigC_sim=0.610 # (nm)
#lamC_sim=16.0

#S/L
#Short (4 or fewer C backbone)
epsCH_sim=15.00 # (K)
sigCH_sim=0.470 # (nm)
lamCH_sim=16.0

epsC_sim=1.45 # (K)
sigC_sim=0.610 # (nm)
lamC_sim=16.0

lam_sim=16.0

elif [ "$Model" = 'TAMie' ]
then

#TAMie

bondlength_CH3=0.174

epsCH3_sim=136.318 # (K)
sigCH3_sim=0.36034 # (nm)
lamCH3_sim=14.0

epsCH2_sim=52.9133 # (K)
sigCH2_sim=0.404 # (nm)
lamCH2_sim=14.0

epsCH_sim=14.5392 # (K)
sigCH_sim=0.43656 # (nm)
lamCH_sim=14.0

epsC_sim=0.0 # (K)
sigC_sim=0.0 # (nm)
lamC_sim=14.0

lam_sim=14.0

elif [ "$Model" = 'AUA4' ]
then

#AUA4

bondlength_CH3=0.1751

epsCH3_sim=120.15 # (K)
sigCH3_sim=0.3607 # (nm)
lamCH3_sim=12.0

epsCH2_sim=86.29 # (K)
sigCH2_sim=0.3431 # (nm)
lamCH2_sim=12.0

epsCH_sim=50.98 # (K)
sigCH_sim=0.3363 # (nm)
lamCH_sim=12.0

epsC_sim=15.04 # (K)
sigC_sim=0.244 # (nm)
lamC_sim=12.0

lam_sim=12.0

#Standard CC bondlength for AUA4 = 0.1535 nm
#delta_CH3=0.0216 nm
#delta_CH2=0.0384 nm
#delta_CH=0.0646 nm

fi

if [ "$lam_sim" = 12.0 ]
then

echo "Using native Lennard-Jones potential"
mdp_path=~/LennardJonesVar
#Originally I tried to have a tab_flag variable, but it was not working for tabulated because
#using "$tab_flag" in mdrun was different than echo. Fortunately, I discovered that Gromacs
#will just ignore the provided table if the vdwtype is cut-off instead of user.
#tab_flag=""

else

echo "Using tabulated potential"
mdp_path=~/Tabulated
#This also required copying tab_it into the output_path whereas it is typically output_path/MCMC_iMCMC
#tab_flag="-table $output_path/tab_it.xvg "

fi

### Copy experimental data files to output directory

cp "$exp_data_path"/"$Compound"_REFPROP_eta_"$Conditions_type".txt "$output_path"/REFPROP_eta_"$Conditions_type".txt

### Create file and header for state point conditions if not already created
if [ ! -e "$output_path"/"$Conditions_type"Settings.txt ]
then
echo "NMol" "Length (nm)" "Temp (K)" > "$output_path"/"$Conditions_type"Settings.txt	
fi

### Loop through temperatures and densities		
jlim=$((jlim - 1))
batches=$((batches - 1))  # Limits for 0 indexed loops later

for j in $(seq 0 $jlim) # Number of conditions to run

do

### Record state points
echo "$Nmol" "${liquid_box[j]}" "${temps[j]}" >> "$output_path"/"$Conditions_type"Settings.txt

nRep=0

NREP_low=0
NREP_high=$((NREPS-1))

####

for iRep in $(seq 0 $batches) # Number of batches to run (0 9), 10 batches with 20 replicates is 200 parameter sets

do

for iMCMC in $(seq $NREP_low $NREP_high)
do

cd "$output_path" || exit

echo Run Saturated Viscosity for epsCH3 = "$epsCH3_sim" sigCH3 = "$sigCH3_sim" lamCH3 = "$lam_sim" epsCH2 = "$epsCH2_sim" sigCH2 = "$sigCH2_sim" lamCH2 = "$lam_sim" epsCH = "$epsCH_sim" sigCH = "$sigCH_sim" lamCH = "$lam_sim" epsC = "$epsC_sim" sigC = "$sigC_sim" lamC = "$lam_sim"

####

mkdir -p MCMC_"$iMCMC"
cd MCMC_"$iMCMC" || exit

###Create files with force field parameters
"$scripts_path"/force_field_params "$Compound" "$input_path" "$scripts_path" "$lam_sim" "$epsCH3_sim" "$sigCH3_sim" "$epsCH2_sim" "$sigCH2_sim" "$epsCH_sim" "$sigCH_sim" "$epsC_sim" "$sigC_sim" "$bondlength_CH3" "$Nmol"

### Copy tab_it.xvg to previous directory, assumes that all MCMC parameter sets use the same value of lambda
#Not needed anymore
#cp tab_it.xvg ../tab_it.xvg

### Record state points, include header if creating file
### Moved this outside of loop
#if [ -e "$output_path"/MCMC_"$iMCMC"/"$Conditions_type"Settings.txt ]
#then
#echo "NMol" "Length (nm)" "Temp (K)" > "$output_path"/MCMC_"$iMCMC"/"$Conditions_type"Settings.txt
#fi

# Initialize the folders, copy files, and insert variables

cd "$output_path"/MCMC_"$iMCMC" || exit #presumes this dir was made previously

mkdir -p Saturated
cd Saturated || exit

mkdir -p rho"$j"
cd    rho"$j" || exit

#echo "$Nmol" "${liquid_box[j]}" "${Temp}" >> "$output_path"/MCMC_"$iMCMC"/"$Conditions_type"Settings.txt

mkdir -p Rep"$nRep"  
cd    Rep"$nRep" || exit

gmx insert-molecules -ci ../../../"$Compound".gro -nmol "$Nmol" -try 500 -box "${liquid_box[j]}" "${liquid_box[j]}" "${liquid_box[j]}" -o "$Compound"_box.gro > insertout 2>> insertout

#Copy the minimization files
echo mdp_path "$mdp_path"
cp "$mdp_path"/em_steep.mdp em_steep.mdp
cp "$mdp_path"/em_l-bfgs.mdp em_l-bfgs.mdp

sed -i -e s/some_rvdw/"$rvdw"/ em_steep.mdp
sed -i -e s/some_rvdw/"$rvdw"/ em_l-bfgs.mdp

mkdir -p NVT_eq
cd NVT_eq || exit

# Copy the equilibration files and edit the temperature

cp "$mdp_path"/nvt_eq_no_output_"$BondType".mdp nvt_eq.mdp
sed -i -e s/some_temperature/"${temps[j]}"/ nvt_eq.mdp
sed -i -e s/some_nsteps/"${equil_steps[j]}"/ nvt_eq.mdp
sed -i -e s/some_rvdw/"$rvdw"/ nvt_eq.mdp

# Still creating NVT_prod directory so that our codes are backwards compatible with data analysis (i.e. same directory hierarchy)

mkdir -p NVT_prod
cd NVT_prod || exit

# Create new directory for viscosity run

mkdir -p NVT_vis
cd NVT_vis || exit

# Copy the viscosity files and edit the temperature

cp "$mdp_path"/nvt_vis_no_xv_"$BondType".mdp nvt_vis.mdp
sed -i -e s/some_temperature/"${temps[j]}"/ nvt_vis.mdp
sed -i -e s/some_nsteps/"${prod_steps[j]}"/ nvt_vis.mdp
sed -i -e s/some_rvdw/"$rvdw"/ nvt_vis.mdp

done # for loop over iMCMC

### Run the NPT steps to determine the box sizes
# Skip this if NPT != YES
if [ "$NPT" = "YES" ]
then
bash "$scripts_path"/nAlkanesNPTsteps "$Compound" "$Nmol" "${liquid_box[j]}" "$mdp_path" "$BondType" "${temps[j]}" "${press[j]}" "${equil_steps[j]}" "$rvdw" "$NREP_low" "$NREP_high" "$output_path" "$j" "$nRep" "$scripts_path" "$pin0" "$nt_eq" "$nt_vis" "${prod_steps[j]}"
fi

###First energy minimization

pinoffset=$pin0

for iMCMC in $(seq $NREP_low $NREP_high)
do

echo pinoffset = "$pinoffset"

cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep" || exit  #start fresh for do cycle instead of multiple "cd .."'s

gmx grompp -f em_steep.mdp -c "$Compound"_box.gro -p ../../../"$Compound".top -o em_steep.tpr > gromppout 2>> gromppout
#We now use a different approach for assigning nodes
#gmx mdrun -table "$output_path"/MCMC_"$iMCMC"/tab_it.xvg -pin on -pinoffset "$pinoffset" -pinstride 1 -ntomp 1 -nt 1 -nb cpu -deffnm em_steep > runout 2>> runout &
gmx mdrun -table "$output_path"/MCMC_"$iMCMC"/tab_it.xvg -nt 1 -nb cpu -pme cpu -deffnm em_steep > runout 2>> runout &
cur_pid=$!
taskset -cp "$pinoffset" $cur_pid > /dev/null 2>&1
min_3_pids[${iMCMC}]=$cur_pid  # This is the third such step, hence 3

pinoffset=$((pinoffset+1))

done #for iMCMC

echo "Waiting for em_steep.tpr: Energy Minimization Part1"

for pid in ${min_3_pids[*]}  # Loop over array
do 
wait $pid
done

###Second energy minimization

pinoffset=$pin0

for iMCMC in $(seq $NREP_low $NREP_high)
do

echo pinoffset = "$pinoffset"

cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep" || exit  #start fresh for do cycle instead of multiple "cd .."'s

gmx grompp -f em_l-bfgs.mdp -c em_steep.gro -p ../../../"$Compound".top -o em_l_bfgs.tpr -maxwarn 1 >> gromppout 2>> gromppout
#gmx mdrun -table "$output_path"/MCMC_"$iMCMC"/tab_it.xvg -pin on -pinoffset "$pinoffset" -pinstride 1 -ntomp 1 -nt 1 -nb cpu -deffnm em_l_bfgs >> runout2 2>> runout2 &
gmx mdrun -table "$output_path"/MCMC_"$iMCMC"/tab_it.xvg -nt 1 -nb cpu -pme cpu -deffnm em_l_bfgs > runout2 2>> runout2 &
cur_pid=$!
taskset -cp "$pinoffset" $cur_pid > /dev/null 2>&1
min_4_pids[${iMCMC}]=$cur_pid

pinoffset=$((pinoffset+1))

done #for iMCMC

echo "Waiting for second energy minimization"

for pid in ${min_4_pids[*]}  # Iterate over 4th minimization array
do
wait $pid
done

###Equilibration period

pinoffset=$pin0

for iMCMC in $(seq $NREP_low $NREP_high)
do

echo pinoffset = "$pinoffset"

cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep"/NVT_eq || exit  #start fresh for do cycle instead of multiple "cd .."'s

if ls ../step*.pdb 1> /dev/null 2>&1 #Remove these warning files
then
echo some energy minimizations might have failed for "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep" >> "$output_path"/warnings_job_"$job_date"
rm ../step*.pdb
fi

gmx grompp -f nvt_eq.mdp -c ../em_l_bfgs.gro -p ../../../../"$Compound".top -o nvt_eq.tpr > gromppout 2>> gromppout
#gmx mdrun -table "$output_path"/MCMC_"$iMCMC"/tab_it.xvg -pin on -pinoffset "$pinoffset" -pinstride 1 -nt "$nt_eq" -nb cpu -deffnm nvt_eq > runout 2>> runout &
Tempbox=${liquid_box[j]}  # This is the default for no NPT
if [ $NPT = YES ]
then
Tempbox=$(<../NPT_eq/NPT_prod/Lbox_NPT_ave)  # This is the default for NPT
fi
# Pass the proper box size to the subscript in case it needs to restart it
"$scripts_path"/run_single.sh "$output_path"/MCMC_"$iMCMC"/tab_it.xvg "$nt_eq" cpu cpu nvt_eq "$pinoffset" "$j" "$nRep" "$output_path" "$NREP_low" "$NREP_high" "$Compound" "$Nmol" "$Tempbox" nvt &

pinoffset=$((pinoffset+nt_eq))

done #for iMCMC

echo Waiting for "$Conditions_type" equilibration

nit=0
maxit=3000000
ndone=$(cat "$output_path"/MCMC_*/Saturated/rho"$j"/Rep"$nRep"/NVT_eq/runout | grep -c "GROMACS reminds you")
while [ $ndone -lt $((NREP_high+1)) ] && [ $nit -lt $maxit ]
do
nit=$((nit+1))
sleep 100s 
echo Waiting for "$Conditions_type" equilibration
ndone=$(cat "$output_path"/MCMC_*/Saturated/rho"$j"/Rep"$nRep"/NVT_eq/runout | grep -c "GROMACS reminds you")
done

###Removed production period

###Viscosity period

pinoffset=$pin0
for iMCMC in $(seq $NREP_low $NREP_high)
do

echo pinoffset = "$pinoffset"

cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep"/NVT_eq/NVT_prod/NVT_vis || exit  #start fresh for do cycle instead of multiple "cd .."'s

gmx grompp -f nvt_vis.mdp -c ../../nvt_eq.gro -p ../../../../../../"$Compound".top -o nvt_vis.tpr > gromppout 2>> gromppout
#gmx mdrun -table "$output_path"/MCMC_"$iMCMC"/tab_it.xvg -pin on -pinoffset "$pinoffset" -pinstride 1 -nt "$nt_vis" -nb cpu -deffnm nvt_vis > runout 2>> runout & #Can use more cores in liquid phase since vapor phase will have already finished
gmx mdrun -table "$output_path"/MCMC_"$iMCMC"/tab_it.xvg -nt "$nt_vis" -nb cpu -pme cpu -deffnm nvt_vis > runout 2>> runout &
taskset -cp "$pinoffset"-"$((pinoffset+nt_vis-1))" $! > /dev/null 2>&1

pinoffset=$((pinoffset+nt_vis))

done #for iMCMC


echo Waiting for "$Conditions_type" viscosities

nit=0
maxit=3000000
ndone=$(cat "$output_path"/MCMC_*/Saturated/rho"$j"/Rep"$nRep"/NVT_eq/NVT_prod/NVT_vis/runout | grep -c "GROMACS reminds you")
while [ $ndone -lt $((NREP_high+1)) ] && [ $nit -lt $maxit ]
do
nit=$((nit+1))
sleep 100s #00s 
echo Waiting for "$Conditions_type" viscosities
ndone=$(cat "$output_path"/MCMC_*/Saturated/rho"$j"/Rep"$nRep"/NVT_eq/NVT_prod/NVT_vis/runout | grep -c "GROMACS reminds you")
done

###Data analysis

echo "Waiting for post processing viscosity data"

for iMCMC in $(seq $NREP_low $NREP_high)
do

cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep"/NVT_eq/NVT_prod/NVT_vis || exit  #start fresh for do cycle instead of multiple "cd .."'s

# Don't need trajectories of eq and prod with viscosity
# Modified .mdp so that trajectories are not produced, i.e. this is no longer needed
#rm ../nvt_prod.trr
#rm ../../nvt_eq.trr

### Analyze trajectories for TCAF method
# No longer using TCAF to obtain viscosities
#echo 0 | gmx tcaf -f nvt_vis.trr -s nvt_vis.tpr > tcaf_out 2>> tcaf_out &

Lbox="${liquid_box[j]}"
Vbox=$(echo $Lbox|awk '{print $1*$1*$1}')

### Analyze Green-Kubo and Einstein

echo "$Vbox" | gmx energy -f nvt_vis.edr -s nvt_vis.tpr -vis > vis_out 2>> vis_out &

done #for iMCMC

nit=0
maxit=60 #Don't want this too high, so it will finish eventually
ndone=$(cat "$output_path"/MCMC_*/Saturated/rho"$j"/Rep"$nRep"/NVT_eq/NVT_prod/NVT_vis/vis_out | grep -c "GROMACS reminds you")
while [ $ndone -lt $((NREP_high+1)) ] && [ $nit -lt $maxit ]
do
nit=$((nit+1))
sleep 10s
echo "Still post processing viscosity data"
ndone=$(cat "$output_path"/MCMC_*/Saturated/rho"$j"/Rep"$nRep"/NVT_eq/NVT_prod/NVT_vis/vis_out | grep -c "GROMACS reminds you")
done

echo "Removing large viscosity output files"

for iMCMC in $(seq $NREP_low $NREP_high)
do

cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep"/NVT_eq/NVT_prod/NVT_vis || exit  #start fresh for do cycle instead of multiple "cd .."'s

if [ -e vis_out ]
then
#rm nvt_vis.trr
rm energy.xvg
rm nvt_vis.edr
rm enecorr.xvg
rm -f \#*
else
echo WARNING: No vis_out file for "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep" >> "$output_path"/warnings_job_"$job_date"
fi


### Compile the pressure values; Due to different storage for the bond types, requires nested ifs
if [ "$BondType" = harmonic ]
then

if [ "$Compound" = Ethane ]
then
sed -e '1,/A V E R A G E S/d' nvt_vis.log | grep -m 1 -A1 'Pressure' | tail -n 1 | awk '{print $1}' > press_log
elif [ "$Compound" = C3H8 ] || [ "$Compound" = IC4H10 ] || [ "$Compound" = NEOC5H12 ]
then
sed -e '1,/A V E R A G E S/d' nvt_vis.log | grep -m 1 -A1 'Pressure' | tail -n 1 | awk '{print $2}' > press_log
else
sed -e '1,/A V E R A G E S/d' nvt_vis.log | grep -m 1 -A1 'Pressure' | tail -n 1 | awk '{print $3}' > press_log
fi

else

if [ "$Compound" = Ethane ]
then
sed -e '1,/A V E R A G E S/d' nvt_vis.log | grep -m 1 -A1 'Pressure' | tail -n 1 | awk '{print $5}' > press_log
elif [ "$Compound" = C3H8 ] || [ "$Compound" = IC4H10 ] || [ "$Compound" = NEOC5H12 ]
then
sed -e '1,/A V E R A G E S/d' nvt_vis.log | grep -m 1 -A1 'Pressure' | tail -n 1 | awk '{print $1}' > press_log
else
sed -e '1,/A V E R A G E S/d' nvt_vis.log | grep -m 1 -A1 'Pressure' | tail -n 1 | awk '{print $2}' > press_log
fi
fi

cat "$output_path"/press_all_log press_log > "$output_path"/press_all_temp
cp "$output_path"/press_all_temp "$output_path"/press_all_log
rm "$output_path"/press_all_temp
###


### Compile the box sizes
cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep"/NPT_eq/NPT_prod || exit 1
cat "$output_path"/Lbox_all Lbox_NPT_ave > "$output_path"/Lbox_all_temp
cp "$output_path"/Lbox_all_temp "$output_path"/Lbox_all
rm "$output_path"/Lbox_all_temp

done #for iMCMC

cd "$output_path"/MCMC_"$iMCMC"/Saturated || exit 1

#rm -f rho"$j"/Rep"$nRep"/NVT_eq/NVT_prod/NVT_vis/tcaf_all.xvg

NREP_low=$((NREP_low+NREPS))
NREP_high=$((NREP_high+NREPS))

echo iRep = "$iRep"
echo NREP_low = "$NREP_low"
echo NREP_high = "$NREP_high"

done # for iRep

cd "$output_path" || exit
echo "In $PWD"

###GreenKubo_analyze for all MCMC parameter sets
python "$scripts_path"/GreenKubo_analyze.py --ilow 0 --ihigh $((NREP_low-1)) --nReps 1 --irho "$j" --sat

done

cd "$output_path" || exit
###Create plots to compare with experimental data and correlations
python "$scripts_path"/compare_TDE_REFPROP.py --comp "$Compound" --nrhomax $((j+1)) --"$Conditions_type"

exit 0

#######
