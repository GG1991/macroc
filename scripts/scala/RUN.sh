#!/bin/bash

#SBATCH --job-name="TRACE"
#SBATCH --workdir=.
#SBATCH --output=macroc_%j.out
#SBATCH --error=macroc_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=NP
#SBATCH --qos=xlarge
#SBATCH --time=00:10:00

mpirun -np NP ./macroc -da_grid_x NX -da_grid_y NY -da_grid_z NZ -ts 10
#mpirun -np NP hpcrun -t ./macroc -da_grid_x 100 -da_grid_y 100 -da_grid_z NZ -ts 10
#mpirun -np 1440 hpcprof -S macroc.hpcstruct -I ../src/'*' hpctoolkit-macroc-measurements/
