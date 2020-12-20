#!/bin/bash
#SBATCH --job-name=stability
#SBATCH --time=00:03:00 # hh:mm:ss
#SBATCH --ntasks=1 #number of MPI processes
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3970 # megabytes
#SBATCH --partition=batch
#SBATCH --output=output_stability

module load OpenMPI/3.1.4-GCC-8.3.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd $SLURM_SUBMIT_DIR

echo "calling the python script"
python3 stability.py
