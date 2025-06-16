#!/bin/bash

# Adjust below to indicate location of Berg_AstralMikado/ on your machine
cd ~/Documents/Berg_AstralMikado/mechanics
export myroot=$PWD
export rundir=$PWD/runs
# Adjust below to indicate location of cytosim-bcb/ on your machine
cd ~/Documents/cytosim-bcb
export cymroot=$PWD

export simName="crosslinker_dens_shear"

# Specify number of sweep replicates (each with fixed seed)
export numRep=10

cd $rundir
# prepare directory for each individual run (force value)
mkdir -p $simName
cd $simName
# clean any files/folders created from previous runs
rm -rf *
# create config files in $rundir/$simName for every parameter combo specified in template
$cymroot/python/run/preconfig.py $numRep $myroot/templates/${simName}.cym.tpl
# place them into their own subdirectories of $simName
$cymroot/python/look/collect.py run%04i/config.cym ${simName}????.cym
# run the simulations in parallel
$cymroot/python/look/scan.py $cymroot/bin/sim nproc=24 run????

# create an .html file with final-frame snapshots from each simulation
$cymroot/python/look/scan.py - "$cymroot/bin/play image frame=45" run????
$cymroot/python/look/make_page.py tile=8 width=100 final_snapshots.html run????

# `report` calls relocated to within .cym.tpl file (end of file)
# allows discarding of .cmo files to save storage space
$cymroot/python/look/scan.py - "rm *.cmo" run????
