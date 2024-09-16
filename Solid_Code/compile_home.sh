#!/bin/bash
clear

echo '----------------------------------------------------------------------'
echo ' Compiling on Home PC'

gcc -O2 -o Solid_Solver.o Solid_Solver.c -fopenmp -lm -Wall
gcc -O2 -o Solid_Solver_Debug.o Solid_Solver.c -fopenmp -lm -Wall -g

echo ' Compiled Solid_Solver code'

echo ' Done.'
echo '-----------------------------------------------------------------------'

