#!/bin/bash

# user input
SPECFEMDIR=$1

echo
echo " start: `date`"
echo " setting up enviroment"
echo

# remove trash
rm -rf bin
rm -rf DATABASES_MPI*
rm -rf OUTPUT_FILES*
rm -rf SYNTHETICS

rm -f NOISE_TOMOGRAPHY.py
rm -f PetersonNoiseModel.m

rm -rf DATA/topo_bathy
rm -rf DATA/s20rts
rm -rf DATA/crust1.0
rm -f NOISE_TOMOGRAPHY/S_squared

# make directories
mkdir bin
mkdir DATABASES_MPI
mkdir OUTPUT_FILES

# link executables
cd ./bin
ln -s $SPECFEMDIR/bin/* .
cd ..

# link earth model and topography
cd ./DATA
ln -s $SPECFEMDIR/DATA/topo_bathy .
ln -s $SPECFEMDIR/DATA/crust1.0 .
ln -s $SPECFEMDIR/DATA/s20rts .
cd ..

# link utilities
ln -s $SPECFEMDIR/EXAMPLES/noise_examples/NOISE_TOMOGRAPHY.py .
ln -s $SPECFEMDIR/EXAMPLES/noise_examples/PetersonNoiseModel.m .

# generate noise source time function
./run_generate_S_squared.sh
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo " finished setting up enviroment"
echo " end: `date`"
echo
