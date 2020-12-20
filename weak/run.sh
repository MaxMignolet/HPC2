#!/bin/bash
mpirun -np $nbProcess ./main $param_file $reso_scheme;
