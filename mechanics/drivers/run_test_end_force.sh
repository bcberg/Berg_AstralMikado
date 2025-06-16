#!/bin/bash

# Adjust below to indicate location of Berg_AstralMikado/ on your machine
cd ~/Documents/Berg_AstralMikado/mechanics
export myroot=$PWD
export rundir=$PWD/runs
# Adjust below to indicate location of cytosim-bcb/ on your machine
cd ~/Documents/cytosim-bcb
export cymroot=$PWD

cd $rundir
# prepare directory for each individual run (force value)
mkdir -p test_end_force
cd test_end_force
# clean any files/folders created from previous runs
rm -rf *
# create config files in $rundir/test_end_force for every force value specified in template
$cymroot/python/run/preconfig.py $myroot/templates/test_end_force.cym.tpl
# place them into their own subdirectories of test_end_force
$cymroot/python/look/collect.py run%04i/config.cym test_end_force????.cym
# run the simulations in parallel
$cymroot/python/look/scan.py $cymroot/bin/sim nproc=4 run????

# create an .html file with final-frame snapshots from each simulation
$cymroot/python/look/scan.py - "$cymroot/bin/play image frame=100" run????
$cymroot/python/look/make_page.py tile=5 width=200 final_snapshots.html run????

# report positions of fiber ends
$cymroot/python/look/scan.py - "$cymroot/bin/report fiber:end frame=100 > fiber_ends.txt" run????