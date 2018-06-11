#!/bin/bash
# Will run a computation with the exponential 6 forcefield for
# the enumerated mollecule.
# LennardJones mdp files will be used, Buckingham topology files
# and .gro files will be employed.
# Condition sets are performed bin sequential batches as in all similar scripts
clean() {   # Make it so everything is killed on an interrupt
local pids=$(jobs -pr)
echo "Kill on exit: $pids"
[ -n "$pids" ] && kill $pids
exit 1
}
trap "clean" SIGINT SIGTERM EXIT SIGQUIT  # Call cleanup when asked to

job_date=$(date "+%Y_%m_%d_%H_%M_%S")

Compound=C4H10
Model=Buckingham
Conditions_type=T293highP
BondType=LINCS  #Harmonic (glexible) or LINCS (fixed)
Temp=273

#Set the number of molecules
Nmol=400

#Specify the path location

scripts_path=~/Scripts
conditions_path=~/"$Conditions_type"_Conditions
topology_path=~/Buckingham  # Where to find topology/gro files
exp_data_path=~/TDE_REFPROP_values
output_path=~/"$Compound"/Gromacs/"$Conditions_type"_Viscosity/"$Model"_N"$Nmol"_"$BondType" #Code assumes that a folder already exists with this directory and the eps_sig_lam_MCMC file in it


cp "$scripts_path"/nAlkanes_Buckingham.sh "$output_path"/"$Compound"_job_"$job_date" #Keep track of what jobs have been submitted
touch "$output_path"/warnings_job_"$job_date"

cd "$output_path" || exit

touch press_all_log

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

nsteps_eq=500 #00 # 0.1 ns total 

elif [ "$Compound" = 'C4H10' ]
then

nsteps_eq=100 #000 # 0.2 ns total

else

nsteps_eq=500 #000 # 1 ns total

fi

mdp_path=~/Buckingham

### Loop through temperatures and densities		

for j in $(seq 0 0) #4) # Number of temperatures to run (0 4)

do

### Record state points
echo "$Nmol" "${liquid_box[j]}" "${temps[j]}" >> "$output_path"/"$Conditions_type"Settings.txt

nRep=0

NREPS=2 #20 #Run 20 replicates in parallel
NREP_low=0
NREP_high=$((NREPS-1))

####

for iRep in $(seq 0 0) #2) # Number of batches to run (0 9), 10 batches with 20 replicates is 200 parameter sets

do

for iMCMC in $(seq $NREP_low $NREP_high)
do

cd "$output_path" || exit

echo Run "$Conditions_type" Viscosity for epsCH3 = "$epsCH3_sim" sigCH3 = "$sigCH3_sim" lamCH3 = "$lam_sim" epsCH2 = "$epsCH2_sim" sigCH2 = "$sigCH2_sim" lamCH2 = "$lam_sim"

####

mkdir -p MCMC_"$iMCMC"
cd MCMC_"$iMCMC" || exit

# Put the number of moles into the topology file
cp "$topology_path"/"$Compound".top "$Compound".top 
cp "$topology_path"/"$Compound".gro "$Compound".gro 
sed -i -e s/some_Nmol/"$Nmol"/ "$Compound".top

cd "$output_path"/MCMC_"$iMCMC" || exit #presumes this dir was made previously

mkdir -p Saturated
cd Saturated || exit

mkdir -p rho"$j"
cd    rho"$j" || exit

mkdir -p Rep"$nRep"  
cd Rep"$nRep" || exit

gmx insert-molecules -ci ../../../"$Compound".gro -nmol "$Nmol" -try 500 -box "${liquid_box[j]}" "${liquid_box[j]}" "${liquid_box[j]}" -o "$Compound"_box.gro > insertout 2>> insertout

#Copy the minimization files
echo mdp_path "$mdp_path"
cp "$mdp_path"/em_steep.mdp em_steep.mdp
cp "$mdp_path"/em_l-bfgs.mdp em_l-bfgs.mdp

rvdw=1.4
sed -i -e s/some_rvdw/"$rvdw"/ em_steep.mdp
sed -i -e s/some_rvdw/"$rvdw"/ em_l-bfgs.mdp

mkdir -p NVT_eq
cd NVT_eq || exit

# Copy the equilibration files and edit the temperature

cp "$mdp_path"/nvt_eq_no_output_"$BondType".mdp nvt_eq.mdp
sed -i -e s/some_temperature/"${temps[j]}"/ nvt_eq.mdp
sed -i -e s/some_nsteps/"$nsteps_eq"/ nvt_eq.mdp
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
sed -i -e s/some_rvdw/"$rvdw"/ nvt_vis.mdp

done # for loop over iMCMC

### Run the NPT steps to determine the box sizes
bash "$scripts_path"/nAlkanesNPTbsteps "$Compound" "$Nmol" "${liquid_box[j]}" "$mdp_path" "$BondType" "${temps[j]}" "${press[j]}" "$nsteps_eq" "$rvdw" "$NREP_low" "$NREP_high" "$output_path" "$j" "$nRep" "$scripts_path"


pinoffset=28

for iMCMC in $(seq $NREP_low $NREP_high)
do

echo pinoffset = "$pinoffset"

cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep" || exit  #start fresh for do cycle instead of multiple "cd .."'s

gmx grompp -f em_steep.mdp -c "$Compound"_box.gro -p ../../../"$Compound".top -o em_steep.tpr > gromppout 2>> gromppout
#We now use a different approach for assigning nodes
gmx mdrun -nt 1 -nb cpu -pme cpu -deffnm em_steep > runout 2>> runout &
min_1_pids[${iMCMC}]=$!  # Record the PID from that process
taskset -cp "$pinoffset" $! > /dev/null 2>&1

pinoffset=$((pinoffset+1))

done #for iMCMC

echo "Waiting for em_steep.tpr: Energy Minimization Part1"
for pid in ${min_1_pids[*]}   # For all PID's, wait for completion
do
wait $pid
done

###Second energy minimization

pinoffset=28

for iMCMC in $(seq $NREP_low $NREP_high)
do

echo pinoffset = "$pinoffset"

cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep" || exit  #start fresh for do cycle instead of multiple "cd .."'s

gmx grompp -f em_l-bfgs.mdp -c em_steep.gro -p ../../../"$Compound".top -o em_l_bfgs.tpr -maxwarn 1 >> gromppout 2>> gromppout
gmx mdrun -nt 1 -nb cpu -pme cpu -deffnm em_l_bfgs > runout2 2>> runout2 &
min_2_pids[${iMCMC}]=$!
taskset -cp "$pinoffset" $! > /dev/null 2>&1

pinoffset=$((pinoffset+1))

done #for iMCMC

echo "Waiting for second energy minimization"

for pid in ${min_2_pids[*]}  # Wait until all of these minimizations are completed
do
wait $pid
done

###Equilibration period

pinoffset=28
ntomp_eq=1
nt_eq=2 

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
# Run an NVT run with restart.
"$scripts_path"/run_single_buckingham.sh "$output_path"/MCMC_"$iMCMC"/tab_it.xvg "$nt_eq" cpu cpu nvt_eq "$pinoffset" "$j" "$nRep" "$output_path" "$NREP_low" "$NREP_high" "$Compound" "$Nmol" "${liquid_box[j]}" nvt &

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

###Viscosity period

pinoffset=28
ntomp_vis=1
nt_vis=2

for iMCMC in $(seq $NREP_low $NREP_high)
do

echo pinoffset = "$pinoffset"

cd "$output_path"/MCMC_"$iMCMC"/Saturated/rho"$j"/Rep"$nRep"/NVT_eq/NVT_prod/NVT_vis || exit  #start fresh for do cycle instead of multiple "cd .."'s

gmx grompp -f nvt_vis.mdp -c ../../nvt_eq.gro -p ../../../../../../"$Compound".top -o nvt_vis.tpr > gromppout 2>> gromppout
gmx mdrun -nt "$nt_vis" -nb cpu -pme cpu -deffnm nvt_vis > runout 2>> runout &
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


done #for iMCMC

cd "$output_path"/MCMC_"$iMCMC"/Saturated || exit

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
# python "$scripts_path"/compare_TDE_REFPROP_test.py --comp "$Compound" --nrhomax $((j+1)) --"$Conditions_type"

exit 0

#######

