#!/bin/bash

bin=./datascripts/run.py
extra_str=nabl-v1

batchargs="--n-procs 20 --n-max-jobs 5" # --no-slurm"
for action in cache-parameters; do #cache-parameters seed-partition annotate-seed-partitions; do
    for study in bf520-synth bg505-synth mg505-synth qa255-synth qb850-synth qa013-synth; do
	common="--study $study --extra-str $extra_str"  #  --write-yaml"  # --no-merged
	# common="$common --samples BF520-g-merged" # --n-random-queries 1000"
	# common="$common --get-naive-prob"
	echo $bin $action $common $batchargs --dry
	break
	echo ""
	
    done
done
