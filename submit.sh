#!/bin/bash -l
#SBATCH -A uppmax2024-2-9 -n 8 -t 01:10:00
#SBATCH -M snowy 

module load gcc/12.2.0 openmpi/4.1.4
make sum
time -p mpirun -n 8 ./sum 28
time -p mpirun -n 4 ./sum 28
time -p mpirun -n 2 ./sum 28
time -p mpirun -n 1 ./sum 28