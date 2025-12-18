#!/bin/bash

export repoLocation="~/Documents/Berg_AstralMikado"
export myroot="${repoLocation}/mechanics"
export rundir="${myroot}/runs"
# Assumes Cytosim root is in an "adjacent" folder to Berg_AstralMikado
export cymroot="~/Documents/cytosim-bcb"

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

# nonastral display parameters for runs0000..0002
for runIdx in {0..2}
do
    runCode=$(printf "%04d" "$runIdx")
    cp $myroot/displays/visuals_nonastral.cyp run${runCode}/display.cyp
done
# regular displays for other runs
for runIdx in {3..11}
do
    runCode=$(printf "%04d" "$runIdx")
    cp $myroot/displays/${runName}.cyp run${runCode}/display.cyp
done
# create an .html file with final-frame snapshots from each simulation
$cymroot/python/look/scan.py - "$cymroot/bin/play display.cyp image zoom=0.9 frame=5,30" run????
$cymroot/python/look/make_page.py tile=3 width=200 final_snapshots.html run????