#!/bin/bash

# ==============
# SLURM FLAGS
# ==============
#SBATCH --job-name=specfem_test
#SBATCH --output=my_test.%J.out
#SBATCH --error=my_test.%J.err
#SBATCH --gres=gpu:v100:1
#SBATCH --time=00:35:00

# ================
# ENV VARIABLES
# ================
MESHER_BIN=./bin/xmeshfem3D
DATABASE_BIN=./bin/xgenerate_databases
SOLVER_BIN=./bin/xspecfem3D

# ========
# STEPS
# ========
./run_this_example.sh
