#!/bin/bash

#SBATCH -J sim
#SBATCH --account=carney-dbadre-condo
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -o data/sweep-00_out/sim-%j.out
#SBATCH -e data/sweep-00_out/sim-%j.out

module load r
Rscript src/sim.R