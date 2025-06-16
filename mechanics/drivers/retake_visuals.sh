#!/bin/bash

# Script to regenerate snapshots of specified runs
# Arguments: $1 = name of run (e.g. visuals_crslnkdens)
#            $2 = frames to generate pictures of (e.g. 5,45)


# Adjust below to indicate location of Berg_AstralMikado/ on your machine
cd ~/Documents/Berg_AstralMikado/mechanics
export myroot=$PWD
export rundir=$PWD/runs
# Adjust below to indicate location of cytosim-bcb/ on your machine
cd ~/Documents/cytosim-bcb
export cymroot=$PWD

export runName=$1

cd $rundir/$runName

# import display parameters
$cymroot/python/look/scan.py "cp $myroot/displays/${runName}.cyp display.cyp" run????
# create an .html file with final-frame snapshots from each simulation
$cymroot/python/look/scan.py - "$cymroot/bin/play display.cyp image zoom=0.9 frame=${2} size=1024,1024" run????
$cymroot/python/look/make_page.py tile=3 width=200 final_snapshots.html run????