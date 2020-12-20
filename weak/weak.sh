#!/bin/bash
#SBATCH --job-name=weak
#SBATCH --time=00:03:00 # hh:mm:ss
#SBATCH --ntasks=2 #number of MPI processes
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3970 # megabytes
#SBATCH --partition=batch
#SBATCH --output=output_weak

module load OpenMPI/3.1.4-GCC-8.3.0
cd $SLURM_SUBMIT_DIR

export reso_scheme=0
export param_file=weak.txt
export nbProcess=$SLURM_NTASKS
echo "calling the python script"
python3 weak.py $SLURM_NTASKS
