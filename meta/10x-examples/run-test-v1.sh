#!/bin/bash

label=v0
n_procs=10
echo ./datascripts/run.py cache-parameters --study 10x-examples --version $label --paired-loci --no-slurm --n-procs $n_procs
echo ./datascripts/run.py partition --study 10x-examples --version $label --paired-loci --no-slurm --n-procs $n_procs #--extra-args="--count-parameters"  # --view-ascii --write-to-log-file
