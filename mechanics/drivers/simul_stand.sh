#!/bin/bash

# STANDARD SIMULATION PROTOCOL
# Example use: >> drivers/simul_stand.sh rep_full_shear_aster_actin 30
# Used for following simulations/parameter sweeps:
#   rep_full_shear_aster_actin      (numRep=30)
#   rep_full_tens_aster_actin       (numRep=30)
#   rep_full_shear_aster_otherRhos  (numRep=30)
#   modulus_distr_shear             (numRep=100)
#   modulus_distr_tens              (numRep=100)
#   rf_shear_aster_eq_actin         (numRep=30)
#   rf_shear_aster_eq_actin_otherRhos   (numRep=30)
#   crosslinker_dens_shear          (numRep=15)
#   crosslinker_stiff_shear         (numRep=15)
#   kbend_shear                     (numRep=10)
#   krot_shear                      (numRep=10)

cd ~/Documents/AstralMikadoCYM
export myroot=$PWD
export rundir=$PWD/runs
# Assumes cytosim root is in an "adjacent" folder to AstralMikadoCYM
cd ../cytosim-bcb
export cymroot=$PWD

export runName=$1

# Specify number of sweep replicates (each with fixed seed)
export numRep=$2

cd $rundir
# prepare directory for each individual run (last expanded parameter = force value)
mkdir -p $runName
cd $runName
# clean any files/folders created from previous runs
rm -rf *
# create config files in $rundir/$runName for every force value specified in template
$cymroot/python/run/preconfig.py $numRep $myroot/templates/${runName}.cym.tpl
# place them into their own subdirectories of $runName
$cymroot/python/look/collect.py run%04i/config.cym ${runName}????.cym
# run the simulations in parallel
$cymroot/python/look/scan.py $cymroot/bin/sim nproc=24 run????

# # create an .html file with final-frame snapshots from each simulation
# $cymroot/python/look/scan.py - "$cymroot/bin/play image frame=30" run????
# $cymroot/python/look/make_page.py tile=6 width=100 final_snapshots.html run????

# # for visual inspection, also generate snapshots at initial times
# $cymroot/python/look/scan.py - "$cymroot/bin/play image frame=0,5" run????

# `report` calls relocated to within .cym.tpl file (end of file)
# allows discarding of .cmo files to save storage space
$cymroot/python/look/scan.py - "rm *.cmo" run????