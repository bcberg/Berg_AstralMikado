#!/bin/bash

# Script intended to run from AstralMikadoCYM
# Assumes cytosim root is in an "adjacent" folder to AstralMikadoCYM
cd ~/Documents/AstralMikadoCYM
export myroot=$PWD
export rundir=$PWD/runs
cd ../cytosim-bcb
export cymroot=$PWD
export simName="timeseries_astralNum_physkT"

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
$cymroot/python/look/scan.py $cymroot/bin/sim nproc=8 run????

# report, framewise, network position over entire simulation
$cymroot/python/look/scan.py "$cymroot/bin/reportF fiber:point root=pos" nproc=8 run????

# import display parameters
# nonastral display parameters for run0000 and run0004
for runIdx in {0,4}
do
    runCode=$(printf "%04d" "$runIdx")
    cp $myroot/displays/visuals_nonastral.cyp run${runCode}/display.cyp
done
# regular displays for other runs
for runIdx in {1,2,3,5,6,7}
do
    runCode=$(printf "%04d" "$runIdx")
    cp $myroot/displays/timeseries_astralNum.cyp run${runCode}/display.cyp
done
# generate .png of each frame for movies
$cymroot/python/look/scan.py "$cymroot/bin/play display.cyp movie zoom=0.8 size=1024,1024" nproc=8 run????

# example command using ffmpeg to make a movie out of image%04d.png files
# -r specifies frames per second
# -crf specifies quality (between 0 (better) and 51 (worse))
# >> ffmpeg -r 10 -i image%04d.png -pix_fmt yuv420p -c:v libx264 -crf 23 movie.mp4