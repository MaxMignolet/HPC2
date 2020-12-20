#!/bin/bash
#SBATCH --job-name=w2_expl
#SBATCH --time=00:10:00 # hh:mm:ss
#SBATCH --ntasks=2 #number of MPI processes
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3970 # megabytes
#SBATCH --partition=batch
#SBATCH --output=w2_expl

module load OpenMPI/3.1.4-GCC-8.3.0
cd $SLURM_SUBMIT_DIR

export reso_scheme=0
export param_file=w2_expl.txt
export nbProcess=$SLURM_NTASKS
echo "calling the python script"
python3 weak.py $SLURM_NTASKS $param_file $reso_scheme
