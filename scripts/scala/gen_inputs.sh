#!/bin/bash

NX=100
NY=100
NZ=100

nodes=(30 35 40 50 55 60 65 70)

function generate {
for i in ${nodes[@]}; do 
        rm -rf $i
	mkdir $i;
	NP=$(($i * 48))
	sed  "\
		s/\(.*\)NX\(.*\)/\1${NX}\2/ ; \
		s/\(.*\)NY\(.*\)/\1${NY}\2/ ; \
		s/\(.*\)NZ\(.*\)/\1${NY}\2/ ; \
		s/\(.*\)NP\(.*\)/\1${NP}\2/ ; \
		" RUN.sh > $i/RUN.sh 
done
}

function launch {
for i in ${nodes[@]}; do 
	cd $i
	sbatch RUN.sh
	cd ..
done
}

#generate
launch
