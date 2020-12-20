#!/bin/bash
echo $nbProcess
echo $param_file
echo $reso_scheme
echo $OMP_NUM_THREADS
mpirun -np $nbProcess ./main $param_file $reso_scheme;
