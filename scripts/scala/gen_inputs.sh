#!/bin/bash

NX=100
NY=100
NZ=100

nodes=(2 4 6 8 10 12 14 16)

function generate {
for i in ${nodes[@]}; do 
        rm -rf $i
	mkdir $i;
	NP=$(($i * 48))
	sed  "\
		s/\(.*\)NX\(.*\)/\1${NX}\2/ ; \
		s/\(.*\)NY\(.*\)/\1${NY}\2/ ; \
		s/\(.*\)NZ\(.*\)/\1${NZ}\2/ ; \
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

function take_times {
rm -rf times.dat time_aux.dat
for i in ${nodes[@]}; do 
	awk -v p=$i '/Elapsed/{print p"\t"$4}' $i/macroc_*.out >> times.dat
done
awk 'NR==1 {t1=($1 * $2);}{print $1"\t"$2"\t"t1/$2}' times.dat >> times_aux.dat
mv times_aux.dat times.dat
}

#generate
#launch
#take_times
