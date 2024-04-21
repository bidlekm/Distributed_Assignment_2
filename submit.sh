#!/bin/bash -l
#SBATCH -A uppmax2024-2-9 -n 8 -t 01:10:00
#SBATCH -M snowy 
#SBATCH -J Rafitos_job

module load gcc/12.2.0 openmpi/4.1.4
make
mpirun -n 8 ./stencil input96.txt output96_1.txt 1
mpirun -n 8 ./stencil input96.txt output96_2.txt 2