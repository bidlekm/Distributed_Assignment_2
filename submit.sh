#!/bin/bash -l
#SBATCH -A uppmax2024-2-9 -n 16 -t 01:10:00
#SBATCH -M snowy 
#SBATCH -J Rafitos_job

module load gcc/12.2.0 openmpi/4.1.4
make

for np in $(seq 1 16); do
    mpirun --bind-to none -n $np ./stencil /proj/uppmax2024-2-9/A2/input4000000.txt output120_2.txt 4 >> timing_4M_4.txt
done