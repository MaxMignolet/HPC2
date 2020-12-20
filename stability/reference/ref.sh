#!/bin/bash
#SBATCH --job-name=stab_ref
#SBATCH --time=00:30:00 # hh:mm:ss
#SBATCH --ntasks=4 #number of MPI processes
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3970 # megabytes
#SBATCH --partition=batch
#SBATCH --output=stab_ref

module load OpenMPI/3.1.4-GCC-8.3.0
cd $SLURM_SUBMIT_DIR

declare -a reso_scheme=1
declare -a param_file=ref.txt

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

mpirun -np $SLURM_NTASKS --bind-to none ./main $param_file $reso_scheme;
