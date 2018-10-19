#!/bin/bash

#SBATCH --job-name="TRACE"
#SBATCH --workdir=.
#SBATCH --output=macro_%j.out
#SBATCH --error=macro_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1440
#SBATCH --qos=xlarge
#SBATCH --time=00:10:00

NX=100
NY=100
NZ=100

mpirun -np 1440 hpcrun -t ./macroc -da_grid_x $NX -da_grid_y $NY -da_grid_z $NZ

#mpirun -np 1440 hpcprof -S macroc.hpcstruct -I ../src/'*' hpctoolkit-macroc-measurements/
