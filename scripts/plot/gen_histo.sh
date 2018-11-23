#!/bin/bash

file=$1
time=$2

awk -v line=$time 'BEGIN{for(i=0;i<line;i++) getline;}{for(i=1;i<=NF;i++) print i"\t"$i;exit;}' $file > data.dat
./plot.gpl
