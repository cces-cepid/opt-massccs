#!/bin/bash
#SBATCH -N 4
#SBATCH --time=900:00:00
#SBATCH --exclusive
#SBATCH --partition cpu

module purge
module load python/3.11.2 

export APPTAINER_UNSHARE_IPC=1

python diff_evol_ompc.py CCS_exp_CO2.dat bounds.dat

