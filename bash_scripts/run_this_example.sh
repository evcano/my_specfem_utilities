#!/bin/bash

# ENVIROMENT VARIABLES
# ====================
WORKDIR=`pwd`
SPECFEMDIR="/home/valeroe/packages/specfem3d_globe"
BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `

# read Par_file and compute total number of cores
NPROC_XI=`grep ^NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep ^NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep ^NCHUNKS DATA/Par_file | cut -d = -f 2 `
NPROCTOT_VAL=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))


# SETUP ENVIROMENT
# ==================
./scripts/setup_enviroment.sh $SPECFEMDIR
if [[ $? -ne 0 ]]; then exit 1; fi

# copy noise distribution
cp NOISE_TOMOGRAPHY/noise_distribution/* DATABASES_MPI/

# RUN MESHER
# =================
./scripts/run_mesher.sh $NPROCTOT_VAL
if [[ $? -ne 0 ]]; then exit 1; fi

# FIRST CONTRIBUTION STEPS 1 AND 2
# ================================
./scripts/change_master_receiver.sh 1

cp DATABASES_MPI_MESHER/* DATABASES_MPI
cp OUTPUT_FILES_MESHER/* OUTPUT_FILES

./scripts/run_solver_noise_step.sh $NPROCTOT_VAL 1
if [[ $? -ne 0 ]]; then exit 1; fi

./scripts/run_solver_noise_step.sh $NPROCTOT_VAL 2
if [[ $? -ne 0 ]]; then exit 1; fi

rm DATABASES_MPI/proc*_noise_surface_movie.bin
mv DATABASES_MPI DATABASES_MPI_REC1
mv OUTPUT_FILES OUTPUT_FILES_REC1

# SECOND CONTRIBUTION STEPS 1 AND 2
# ================================
./scripts/change_master_receiver.sh 2

mv DATABASES_MPI_MESHER DATABASES_MPI
mv OUTPUT_FILES_MESHER OUTPUT_FILES

./scripts/run_solver_noise_step.sh $NPROCTOT_VAL 1
if [[ $? -ne 0 ]]; then exit 1; fi

./scripts/run_solver_noise_step.sh $NPROCTOT_VAL 2
if [[ $? -ne 0 ]]; then exit 1; fi

rm DATABASES_MPI/proc*_noise_surface_movie.bin
mv DATABASES_MPI DATABASES_MPI_REC2
mv OUTPUT_FILES OUTPUT_FILES_REC2
