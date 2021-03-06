#!/bin/bash
#SBATCH --job-name=test_1
#SBATCH --time=00:06:00 # hh:mm:ss
#SBATCH --ntasks=1 #number of MPI processes
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=3970 # megabytes
#SBATCH --partition=batch
#SBATCH --output=output_expl_1node

module load OpenMPI/3.1.4-GCC-8.3.0
cd $SLURM_SUBMIT_DIR

declare -a n_threads=(24 16 12 4 1)
declare -a reso_scheme=0
declare -a param_file=test_expl.txt

for ((i=0; i<${#n_threads[@]}; i++ ));
do
	export OMP_NUM_THREADS=${n_threads[$i]};
	echo "-The number of threads is ${n_threads[$i]}-";
	mpirun -np $SLURM_NTASKS ./main $param_file $reso_scheme;
done
