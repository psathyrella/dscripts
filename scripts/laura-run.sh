#!/bin/bash

bin=./datascripts/run.py
extra_str=v20

batchargs="--n-procs 20 --n-max-jobs 8"
for action in partition; do #cache-parameters seed-partition annotate-seed-partitions; do
    for study in mg505-synth; do  #kate-qrs laura-mb laura-mb-2 qa255-synth; do
	common="--study $study --extra-str $extra_str"  #  --write-yaml"  # --no-merged
	common="$common --n-random-queries 50000"  # --samples BF520-g-merged" #
	# common="$common --get-naive-prob"
	echo $bin $action $common $batchargs --logstr=50k --n-random-subsets 3  # --count-parameters
    done
done
