#!/bin/bash

export repoLocation="~/Documents/Berg_AstralMikado"
export myroot="${repoLocation}/mechanics"
export rundir="${myroot}/runs"
# Assumes Cytosim root is in an "adjacent" folder to Berg_AstralMikado
export cymroot="~/Documents/cytosim-bcb"

export runName="perc_in_cytosim"

# Specify number of sweep replicates (each with fixed seed)
export numRep=100

cd $rundir
# prepare directory for each individual run
mkdir -p $runName
cd $runName
# clean any files/folders created from previous runs
rm -rf *
# create config files in $rundir/$runName for every parameter value specified in template
# (last expanded parameter = crosslinker bind rate)
# 5-wide digits necessary for 35*4*100 = 14000 sims.
$cymroot/python/run/preconfig.py -5 $numRep $myroot/templates/${runName}.cym.tpl
# place them into their own subdirectories of $runName (note 5-wide digits)
$cymroot/python/look/collect.py run%05d/config.cym ${runName}?????.cym
# run the simulations in parallel
# `report` calls are within .cym.tpl file (end of file), so can immediately delete .cmo files
$cymroot/python/look/scan.py "$cymroot/bin/sim; rm *.cmo" nproc=28 run?????
