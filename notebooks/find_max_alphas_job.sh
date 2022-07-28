#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --job-name=find_max_alphas
#SBATCH --account=rrg-dcurtin
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --output=/project/d/dcurtin/jpbarron/ADM/find_max_alphas.out
cd $SLURM_SUBMIT_DIR
module load gcc
module load mkl
module load openmpi
source $HOME/ADM/class_adm_env/bin/activate
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
python find_max_alphas.py
