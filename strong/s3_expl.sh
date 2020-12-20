#!/bin/bash
#SBATCH --job-name=s3_expl
#SBATCH --time=00:15:00 # hh:mm:ss
#SBATCH --ntasks=3 #number of MPI processes
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3970 # megabytes
#SBATCH --partition=batch
#SBATCH --output=strong_expl_3node

module load OpenMPI/3.1.4-GCC-8.3.0
cd $SLURM_SUBMIT_DIR

declare -a n_threads=(24 20 16 12 8 4 2 1)
declare -a reso_scheme=0
declare -a param_file=strong_expl.txt

for ((i=0; i<${#n_threads[@]}; i++ ));
do
	export OMP_NUM_THREADS=${n_threads[$i]};
	echo "-The number of threads is ${n_threads[$i]}-";
	mpirun -np $SLURM_NTASKS --bind-to none ./main $param_file $reso_scheme;
done
