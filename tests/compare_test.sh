#!/bin/bash

nodes=(10 20 30 40 50)

for np in ${nodes[@]}; do

diff test_stable/${np}/macroc_*.out test_develop/${np}/macroc_*.out

done
