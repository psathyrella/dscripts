#!/bin/bash

pre_vsn=v2
prod_vsn=v0 #v1
max_subjects=${1:-10}
vwdir=/fh/fast/matsen_e/processed-data/partis/vanwinkle-170

# find subjects that have all three loci
subjects=$(for locus in h k l; do
    ls /fh/fast/matsen_e/data/vanwinkle-170/processed/$pre_vsn/ig$locus/*.tsv 2>/dev/null | xargs -I{} basename {} .tsv
done | sort | uniq -c | awk '$1==3{print $2}')
# subjects=110042935 #-igk/plots/multi-hmm/

isubj=0
for subj in $subjects; do
    # if [[ $isubj -ge $max_subjects ]]; then
    #     break
    # fi
    isubj=$((isubj + 1))
    for locus in h k l; do

        # compaare germlines+usage for v0 (germline inference) vs v1 (per-subject genomic germlines)
        cf-germlines.py --names v0:v1 --locus ig$locus $vwdir/v0/$subj-ig$locus/hmm $vwdir/v1/$subj-ig$locus/hmm >$vwdir/v1/$subj-ig$locus/cf-germlines.log &

        # # compare mutations from changeo vs partis (it turns out they're essentially identical)
        # pbase=$vwdir/$prod_vsn
        # ./datascripts/meta/vanwinkle-170/cf-identity-vs-mfreqs.py --tsv /fh/fast/matsen_e/data/vanwinkle-170/processed/$pre_vsn/ig$locus/$subj.tsv \
        #                                --mfreq-plotdirs $pbase/$subj-ig$locus/plots/hmm/mute-freqs/overall:$pbase/partitions/$subj-ig$locus/plots/multi-hmm/mute-freqs/overall \
        #                                --names changeo:hmm:multi-hmm \
        #                                --outdir $fs/tmp/mfreq-test-v2/$subj-ig$locus
    done
    sleep 0.5
done
