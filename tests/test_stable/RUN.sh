#!/bin/bash

#SBATCH --job-name="STABLE"
#SBATCH --workdir=.
#SBATCH --output=macroc_%j.out
#SBATCH --error=macroc_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=NP
#SBATCH --time=00:40:00

EXEC="../../macroc"

srun -n NP ${EXEC} -da_grid_x NX -da_grid_y NY -da_grid_z NZ -ts 10
