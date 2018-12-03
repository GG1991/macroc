#!/bin/bash

# Generates lot of scripts in a folder saved with the date

path="/gpfs/projects/bsc21/bsc21774"
today=$(date +%d-%m-%Y)
folder="${path}/macroc-${today}"

EXEC="~/GIT/macroc/build/macroc"

mkdir -p ${folder}

NPROC=4
NX=50
NY=3
NZ=50
MICRON=10
VTU_FREQ=-1
HS=02
MIN=00

#QUEUE="bsc_case"
QUEUE="interactive"
#QUEUE="debug"
#QUEUE="xlarge"

subfolder="NN_$((NX * NY * NZ))_MICRON_$((MICRON * MICRON * MICRON))_NP_${NPROC}_${QUEUE}"
execfolder="${folder}/${subfolder}"

mkdir -p ${execfolder}

echo -e \
"#!/bin/bash

#SBATCH --job-name=\"MAC-${NPROC}\"
#SBATCH --workdir=${execfolder}
#SBATCH --output=macro_%j.out
#SBATCH --error=macro_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=${NPROC}
#SBATCH --qos=${QUEUE}
#SBATCH --time=${HS}:${MIN}:00

srun -n ${NPROC} ${EXEC} \\
	-da_grid_x $NX -da_grid_y $NY -da_grid_z $NZ \\
	-vtu_freq ${VTU_FREQ} \\
	-new_its 4 \\
	-ts 10000 \\
        -dt 0.001 \\
	-micro_n ${MICRON} \\
	-micro_mat_1 1.0e7,0.25,1.0e4,1.0e4 \\
	-micro_mat_2 1.1e7,0.25,1.0e4,1.0e7
" > "${execfolder}/job.sh"

chmod 766 "${execfolder}/job.sh"

sbatch "${execfolder}/job.sh"

