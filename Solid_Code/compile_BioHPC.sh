#!/bin/bash
clear
module purge
module add shared slurm gcc/8.3.0

gcc -O2 -o Solid_Solver.o Solid_Solver.c -lm -Wall -fopenmp

echo '----------------------------------------------------------------------'
