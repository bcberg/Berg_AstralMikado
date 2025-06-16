#!/bin/bash

# Adjust below to indicate location of Berg_AstralMikado/ on your machine
cd ~/Documents/Berg_AstralMikado/mechanics
export myroot=$PWD
export rundir=$PWD/runs
# Adjust below to indicate location of cytosim-bcb/ on your machine
cd ~/Documents/cytosim-bcb
export cymroot=$PWD

export runName="visuals_crslnkdens"

cd $rundir
# prepare directory for each individual run (force value)
mkdir -p $runName
cd $runName
# clean any files/folders created from previous runs
rm -rf *
# create config files in $rundir/$runName for every force value specified in template
$cymroot/python/run/preconfig.py $myroot/templates/${runName}.cym.tpl
# place them into their own subdirectories of $runName
$cymroot/python/look/collect.py run%04i/config.cym ${runName}????.cym
# run the simulations in parallel
$cymroot/python/look/scan.py $cymroot/bin/sim nproc=12 run????

# import display parameters
$cymroot/python/look/scan.py "cp $myroot/displays/${runName}.cyp display.cyp" run????
# create an .html file with final-frame snapshots from each simulation
$cymroot/python/look/scan.py - "$cymroot/bin/play display.cyp image zoom=0.9 frame=5,45" run????
$cymroot/python/look/make_page.py tile=3 width=200 final_snapshots.html run????