#!/bin/bash
declare -a reso_scheme=0
declare -a param_file=stability.txt

mpirun -np 1 ./main $param_file $reso_scheme;
