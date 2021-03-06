#!/bin/bash
source /usr/local/gromacs/bin/GMXRC
# This is designed to restart an equilibration, either NVT or NPT, run.
# Script to run a single gmx mdrun job and restart it if it fails
# This script will move right into a restart without any user input
# if a failure is detected.
# This will attempt a maximum number of restarts as specified in
# the while loop control.

# Gives a more informative error when something goes wrong
# with the script.
error_report() {
echo "Error $1 on j $2, iMCMC $3, stage $4"
exit 1
}

# Assemble the parameters required for the run and potential restart
clean() {   # Make it so that everything is killed on an interrupt
local pids=$(jobs -pr)
echo "On exit sending kill signal to: $pids"
[ -n "$pids" ] && kill $pids
exit 1
}
trap "clean" SIGINT SIGTERM EXIT SIGQUIT  # Call clean for any such signal

table="$1"
nt="$2"
nb="$3"
pme="$4"
deffnm="$5"
pinoffset="$6"
j="$7"
nRep="$8"
output_path="${9}"
NREP_low="${10}"  # Present for historiacal reasons
NREP_high="${11}"
Compound="${12}"
Nmol="${13}"
liquid_box="${14}"
ensemble="${15}"

# Run given gmx command- if the command fails, the program will not exit but will 
# rather move on to the restart phase
gmx mdrun -table "$table" -nt "$nt" -nb "$nb" -pme $pme -deffnm "$deffnm" > runout 2>> runout &  # UNDO pme
pid=$!  # Record PID of launched process
taskset -cp "$pinoffset"-"$((pinoffset+nt-1))" $pid > /dev/null 2>&1

wait $pid && exit 0 # Wait until the process is finished; exit if job is successful

# Otherwise, the job has failed.
cd .. || error_report "Unable to change directory to .. for restart" "$j" "?" "$ensemble"   # Backup to the proper location to restart
# Necessary for restart
if [ "$ensemble" = NVT ] || [ "$ensemble" = nvt ]
then
ENSEMBLE=NVT
ensemble=nvt
elif [ "$ensemble" = NPT ] || [ "$ensemble" = npt ]
then
ENSEMBLE=NPT
ensemble=npt
else
echo Provide the correct ensemble, either NVT or NPT
exit 1
fi

tries=0

while [ "$tries" -lt 4 ]  # While we have not tried more than three times.
do

tries=$((tries+1))  # Increment the loop counter; we don't want to do this too many times

echo "Failure to run $ensemble reading from table $table. Restart # $tries."

echo "$PWD is working directory"
cp "$ENSEMBLE"_eq/runout "$ENSEMBLE"_eq/runout_fatal_error
cp "$ENSEMBLE"_eq/gromppout "$ENSEMBLE"_eq/gromppout_fatal_error
cp "$ENSEMBLE"_eq/"$ensemble"_eq.log "$ENSEMBLE"_eq/"$ensemble"_eq_fatal_error.log

gmx insert-molecules -ci ../../../"$Compound".gro -nmol "$Nmol" -try 100 -box "$liquid_box" "$liquid_box" "$liquid_box" -o "$Compound"_box.gro > insertout 2>> insertout || echo Insert molecules failed.

echo Inserted molecules for restart

# Now, run the first energy minimization
gmx grompp -f em_l-bfgs.mdp -c "$Compound"_box.gro -p ../../../"$Compound".top -o em_l-bfgs.tpr -maxwarn 1 > gromppout 2>> gromppout
gmx mdrun -table "$table" -nt 1 -nb "$nb" -pme $pme -deffnm em_l-bfgs > runout 2>> runout &
pid=$!
taskset -cp "$pinoffset" $pid > /dev/null 2>&1

wait $pid || echo First energy minimization failed.  # Report on status
echo Finished first energy minimization attempt for restart.

# Run the second energy minimization
gmx grompp -f em_steep.mdp -c em_l-bfgs.gro -p ../../../"$Compound".top -o em_steep.tpr >> gromppout 2>> gromppout
gmx mdrun -table "$table" -nt 1 -nb "$nb" -pme $pme -deffnm em_steep > runout2 2>> runout2
pid=$!
taskset -cp "$pinoffset" $pid > /dev/null 2>&1

wait $pid || echo Second energy minimization failed.
echo Finished second energy minimization attempt for restart.

cd "$ENSEMBLE"_eq || error_report "Unable to change directory to $ENSEMBLE_eq for restart" "$j" "?" "$ensemble"   # Transfer into subdirectory to complete run

# Do the actual restart run. Be more tolerant of warnings from grompp this time around.
gmx grompp -f "$ensemble"_eq.mdp -c ../em_steep.gro -p ../../../../"$Compound".top -o "$ensemble"_eq.tpr -maxwarn 5 > gromppout 2>> gromppout
gmx mdrun -table "$table" -nt "$nt" -nb "$nb" -pme $pme -deffnm "$ensemble"_eq > runout 2>> runout &
pid=$!
taskset -cp "$pinoffset"-"$((pinoffset+nt-1))" $pid > /dev/null 2>&1

echo Restarted "$ensemble" run from table "$table".

wait $pid && exit 0  # GMX ran if the && statement is executed. Exit.

cd .. || error_report "Unable to change directory to .. for restart" "$j" "?" "$ensemble"  # Back off again to repeat calculations.

done

echo Warning: exceeded maximum restart attempts for "$ensemble" job in "$PWD" 
echo User intervention required. Exiting.

exit 1  # GMX never ran


