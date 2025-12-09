#!/bin/bash

cd ~/Documents/AstralMikadoCYM
export myroot=$PWD
export rundir=$PWD/runs
# Assumes cytosim root is in an "adjacent" folder to AstralMikadoCYM
cd ../cytosim-bcb
export cymroot=$PWD

export runName="rf_shear_aster_actin_physkT"

# Specify number of sweep replicates (each with fixed seed)
export numRep=20

cd $rundir
# prepare directory for each individual run (force value)
mkdir -p $runName
cd $runName
# clean any files/folders created from previous runs
rm -rf *
# create config files in $rundir/rep_full_shear_aster_actin for every force value specified in template
$cymroot/python/run/preconfig.py $numRep $myroot/templates/${runName}.cym.tpl
# place them into their own subdirectories of rep_full_shear_aster_actin
$cymroot/python/look/collect.py run%04i/config.cym ${runName}????.cym
# run the simulations in parallel
$cymroot/python/look/scan.py $cymroot/bin/sim nproc=24 run????

# create an .html file with final-frame snapshots from each simulation
$cymroot/python/look/scan.py - "$cymroot/bin/play image frame=48" nproc=24 run????
$cymroot/python/look/make_page.py tile=6 width=100 final_snapshots.html run????

# for visual inspection, also generate snapshots at initial times
$cymroot/python/look/scan.py - "$cymroot/bin/play image frame=0,14" nproc=24 run????

# report position of fibers
# initial position frames: 5-14 (41..50 sec)
# final position frames: 39-48 (291..300 sec)
$cymroot/python/look/scan.py "$cymroot/bin/reportF fiber:point root=pos" nproc=24 run????
# delete unneeded position data (may cause rm error (argument list too long))
rm -rf run*/pos{0000..0004}.txt 
rm -rf run*/pos{0015..0038}.txt
