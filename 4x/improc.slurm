#!/bin/bash

#SBATCH -J ImProc_4x

#SBATCH -t 02:00:00
#SBATCH -N 1 # node
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4096M
#SBATCH --array=0-999

#SBATCH --mail-user=qcaudron@princeton.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL



python batch.py $SLURM_ARRAY_TASK_ID



# submit using sbatch
# scontrol show job 1373263
# squeue -u qcaudron
# scancel 1423435