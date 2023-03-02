#!/bin/bash

# user input
NPROCTOT_VAL=$1

echo
echo " start: `date`"
echo " running mesher on $NPROCTOT_VAL cores"
echo

mpirun -np $NPROCTOT_VAL ./bin/xmeshfem3D
if [[ $? -ne 0 ]]; then exit 1; fi

mkdir -p OUTPUT_FILES_MESHER
mv OUTPUT_FILES/* OUTPUT_FILES_MESHER

mkdir -p DATABASES_MPI_MESHER
mv DATABASES_MPI/* DATABASES_MPI_MESHER

echo
echo " finished mesher"
echo " end: `date`"
echo
