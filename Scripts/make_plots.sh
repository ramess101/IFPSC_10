#!/bin/bash
# These are the nams of the scripts as posted on Git.
cff=~/Scripts/compare_force_fields_mca_mod.py
ctr=~/Scripts/compare_TDE_REFPROP_force_fields_mca_mod.py

# Short normal
python "$cff" --mod Potoff TAMie TraPPE --comp C3H8 C4H10 C8H18 --settings 4 --devPlot --uncer --validation Potoff TAMie TraPPE
mv ~/Viscosity_compare_force_fields/compare_force_fields.pdf ~/Viscosity_compare_force_fields/short_norm.pdf
# Long normal
python "$cff" --mod Potoff TAMie TraPPE --comp C12H26 C16H34 C22H46 --settings 5 --devPlot --uncer --validation Potoff TAMie TraPPE
mv ~/Viscosity_compare_force_fields/compare_force_fields.pdf ~/Viscosity_compare_force_fields/long_norm.pdf
# Short branched
python "$cff" --mod Potoff TAMie TraPPE AUA4 --comp IC4H10 IC5H12 23DMButane NEOC5H12 --settings 1 --devPlot --uncer --validation Potoff TAMie TraPPE
mv ~/Viscosity_compare_force_fields/compare_force_fields.pdf ~/Viscosity_compare_force_fields/short_branch.pdf
# Long branched
python "$cff" --mod Potoff TAMie TraPPE --comp 3MPentane IC8H18 IC6H14 --settings 2 --devPlot --uncer --validation Potoff TAMie TraPPE
mv ~/Viscosity_compare_force_fields/compare_force_fields.pdf ~/Viscosity_compare_force_fields/long_branch.pdf
# C3H8
python "$ctr" --comp C3H8 --T293highP --devPlot --uncer --mod Potoff TAMie TraPPE
# C4H10
python "$ctr" --comp C4H10 --T293highP --devPlot --uncer --mod Potoff TAMie TraPPE
# IC4H10
python "$ctr" --comp IC4H10 --T293highP --uncer --mod Potoff TAMie TraPPE
# IC5H12
python "$ctr" --comp IC5H12 --T293highP --uncer --mod Potoff TAMie TraPPE
# IC8H18
python "$ctr" --comp IC8H18 --T293highP --uncer --mod Potoff TraPPE --nrhomax 4
# C8H18
python "$ctr" --comp C8H18 --T293highP --uncer --mod Potoff TAMie TraPPE --nrhomax 4
# 3MPentane
python "$ctr" --comp 3MPentane --T293highP --uncer --mod Potoff TAMie TraPPE

