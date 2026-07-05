#!/bin/bash

#SBATCH --job-name=null_model
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=9:59:59
#SBATCH --mem=20000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.username@institution.edu

module purge
module load r/4.5.0
module load intel/2025.0.1
export R_LIBS_USER=/path-to-folder/R-4.5.0-MKL 
/path-to-folder/R-4.5.0/bin/R --vanilla --args ${SLURM_ARRAY_TASK_ID} < $1 > "${1}.${SLURM_ARRAY_TASK_ID}.out" 2> "${1}.${SLURM_ARRAY_TASK_ID}.err"