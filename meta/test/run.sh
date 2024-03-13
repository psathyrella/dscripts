#!/bin/bash

label=v0
n_procs=10
study=test
for action in cache-parameters partition seed-partition; do  # NOTE cache-parameters has to finish before you can run anything else
    echo ./datascripts/run.py $action --study $study --version $label --paired-loci --n-procs $n_procs --dry-run
done
