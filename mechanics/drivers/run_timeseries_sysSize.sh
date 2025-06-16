#!/bin/bash

# Adjust below to indicate location of Berg_AstralMikado/ on your machine
cd ~/Documents/Berg_AstralMikado/mechanics
export myroot=$PWD
export rundir=$PWD/runs
# Adjust below to indicate location of cytosim-bcb/ on your machine
cd ~/Documents/cytosim-bcb
export cymroot=$PWD

export simName="timeseries_sysSize"

cd $rundir
# prepare directory for each individual run
mkdir -p $simName
cd $simName
# clean any files/folders created from previous runs
rm -rf *
# create config files in $rundir/timeseries_sysSize
$cymroot/python/run/preconfig.py $myroot/templates/${simName}.cym.tpl
# place them into their own subdirectories of timeseries_sysSize
$cymroot/python/look/collect.py run%04i/config.cym ${simName}????.cym
# run the simulations in parallel
$cymroot/python/look/scan.py $cymroot/bin/sim nproc=4 run????

# report, framewise, network position over entire simulation
$cymroot/python/look/scan.py "$cymroot/bin/reportF fiber:point root=pos" nproc=4 run????
