#!/bin/bash

study=vanwinkle-170

pvsn=v3

cmd=./datascripts/meta/$study/process.py
bd=/fh/fast/matsen_e/data/$study
# echo "$cmd --version $pvsn >$bd/processed/$pvsn.log"
# exit 0

bin=./datascripts/run.py
ovsn=v2

batchargs="--n-procs 20"
for action in cache-parameters partition; do
# :all-bulk-merged-100k-igh
    common="--study $study --version $ovsn"
    # common=" $common --plot-partitions --extra-args=\"--meta-info-key-to-color tissues --meta-info-bkg-val pbmc\" --ex-out-version tissues"
    # common=" $common --plot-partitions --extra-args=\"--meta-info-key-to-color timepoint-tissues --meta-info-bkg-val None-pbmc\" --ex-out-version timepoint-tissues"
    echo $bin $action $common $batchargs --dry # --loci igh:igk:igl  # bulk-only samples
done
