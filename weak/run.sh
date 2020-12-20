#!/bin/bash
mpirun --bind-to none -np $nbProcess ./main $param_file $reso_scheme;
