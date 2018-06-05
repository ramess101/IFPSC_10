#!/bin/bash
nmol=400  # Number of molecules
Script_Dir=~/IFPSC_10/Scripts

for condition in T293highP Near_Crit_HighP sat
do
  for compound in Propane Butane Octane Isobutane Isopentane Isohexane Isooctane Neopentane
  do
    echo "BEGINNING $condition for $compound" >> gen_out
    python "$Script_Dir"/simulation_conditions.py --comp "$compound" --Nmol "$nmol" --"$condition" >> gen_out 2>>gen_out || echo "Error on generation of $condition for $compound"
  done
done
