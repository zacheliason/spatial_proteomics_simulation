#!/bin/bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10240M   # memory per CPU core
#SBATCH -J "sc_sim_enhance"   # job name
#SBATCH --mail-user=jhgirald@byu.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load miniconda3
conda activate #TODO: Add env name

Rscript BayesSpaceEnhanceScript.R #TODO: Add dir name to line on python script
