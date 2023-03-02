#!/bin/bash

# get user input
NPROCTOT_VAL=$1
STEP=$2

if [ "$STEP" == "" ]; then echo "usage: ./run_solver_noise_step.sh NPROCTOT_VAL STEP[1/2/3]"; exit 1; fi

# modify Par_file
case $STEP in
1) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 1:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" DATA/Par_file
   ;;
2) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 1:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 2:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .true.:" DATA/Par_file
   ;;
3) sed -i "s:^SIMULATION_TYPE .*:SIMULATION_TYPE = 3:" DATA/Par_file
   sed -i "s:^NOISE_TOMOGRAPHY .*:NOISE_TOMOGRAPHY = 3:" DATA/Par_file
   sed -i "s:^SAVE_FORWARD .*:SAVE_FORWARD = .false.:" DATA/Par_file
   ;;
*) echo "step not recognized: $step"; echo "please use as step number 1, 2 or 3"; exit 1 ;;
esac

# run solver
echo
echo " start: `date`"
echo " running noise simulation step $STEP on $NPROCTOT_VAL cores"
echo

mpirun -np $NPROCTOT_VAL ./bin/xspecfem3D
if [[ $? -ne 0 ]]; then exit 1; fi

cp DATA/Par_file OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

echo
echo " finished noise simulation step $STEP"
echo " end: `date`"
echo
