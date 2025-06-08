#!/bin/bash

#!/bin/bash

bd=/fh/fast/matsen_e/data/parul-vrc01gh

# NOTE species is wrong[ish], since these seqs are a mix of human and mouse germlines, but seems to work fine

# echo ./datascripts/meta/parul-vrc01gh/process.py  # NOTE edit in/out version by hand

# ivsn=v1
# ovsn=v2
# # NOTE first run process.py (above) to get timepoint meta file
# echo ./bin/split-loci.py $bd/processed/$ivsn/all-seqs.fa --outdir $bd/processed/$ovsn/split-loci --guess-pairing-info --droplet-id-separators - --droplet-id-indices 0:1:2 --input-metafname $bd/processed/$ivsn/meta.yaml # --germline-dir maybe?
# NOTE run make-pseudo-samples.py here (in which you have to update the various versions by hand)
# exit 0

dvsn=v5

bin=./datascripts/run.py
for action in cache-parameters partition; do
    echo $bin $action --study parul-vrc01gh --version $dvsn --paired-loci --no-slurm --dry \
         --plot-partitions --ex-out-version no-mut-labels --extra-args=\"--add-unpaired-seqs-to-fake-paired-annotations\"
         # --plot-partitions --ex-out-version vrc01-muts --extra-args=\"--add-unpaired-seqs-to-fake-paired-annotations --branch-color-key vrc01-muts\"
         # --plot-partitions --ex-out-version all-mut-labels --extra-args=\"--add-unpaired-seqs-to-fake-paired-annotations --mutation-label-cfg all:mut-strs\"
         # --plot-partitions --ex-out-version label-loci --extra-args=\"--add-unpaired-seqs-to-fake-paired-annotations --label-leaf-nodes --node-label-regex=[hk]\"
         # --plot-partitions --ex-out-version no-unpaired-seqs --extra-args=\"--tree-inference-subdir no-unpaired-seqs --mutation-label-cfg all\" #:mut-strs
done
exit 0

# not sure if this linearham stuff is still needed, maybe it's old?
# # ./test/linearham-run.py --outdir /fh/fast/matsen_e/processed-data/partis/parul-vrc01gh/no-dups/d1/linearham --partis-outdir /fh/fast/matsen_e/processed-data/partis/parul-vrc01gh/no-dups/d1 --n-sim-events 10 --docker --n-procs 15 --dry
# ./test/linearham-run.py --outdir /fh/fast/matsen_e/processed-data/partis/parul-vrc01gh/v3/d1/linearham --partis-outdir /fh/fast/matsen_e/processed-data/partis/parul-vrc01gh/v3/d1 --docker --n-procs 15 --local-docker-image --ignore-unmutated-seqs
# for i in 0; do ./bin/partis plot-partitions --outfname /fh/fast/matsen_e/processed-data/partis/parul-vrc01gh/v3/d1/linearham/with-inferred-ancestors/itree-$i/partition-igk.yaml --plotdir /fh/fast/matsen_e/dralph/partis/tmp-plots-v3-igk/$i --meta-info-key-to-color timepoints --label-mutations --label-tree-nodes; done

svsn=v5 #def-v0 #v6-$ird #

# # make boosted and non-boosted simulation
# for ird in 0 1 2; do # 1 2 3 4; do
# mkey=timepoints # gc-rounds
# bd=/fh/fast/matsen_e/processed-data/partis/parul-vrc01gh/simulation/$svsn-$ird # /fh/fast/matsen_e/dralph/partis/tmp/gcr-sim-$svsn
# mkdir -p $bd
# bin="./bin/bcr-phylo-run.py --paired-loci --all-inference-plots --carry-cap 1500 --n-procs 10 --dont-observe-common-ancestors --meta-info-key-to-color $mkey --seed $ird --n-sim-events 5 --restrict-to-single-naive-seq --tree-inference-method iqtree --actions cache-parameters" #  --actions simu
# echo $bin --base-outdir $bd/default --obs-times 1:10:25:50:100:150:200:300:400 --n-sim-seqs-per-generation 150 --sequence-sample-time-fname datascripts/meta/parul-vrc01gh/default-sample.yaml #>$bd/default.log &
# echo $bin --n-gc-rounds 2 --n-reentry-seqs 30 --base-outdir $bd/boost --obs-times 1,5,10,20,50,75,100,150,200,250,300:1,30,50,75,100,200,300 --n-sim-seqs-per-generation 150 --sequence-sample-time-fname datascripts/meta/parul-vrc01gh/boost-sample.yaml #>$bd/boost.log &
# break
# done
# exit 0

# make a bunch of different types of comparison plots
pvsn=v0
isim=0
tp=w25
bd=/fh/fast/matsen_e/processed-data/partis/parul-vrc01gh
od=$bd/cf-plots/$pvsn
bin=./bin/compare-plotdirs.py
# for ltmp in igh igk; do
# $bin --outdir $od/$ltmp --file-glob-str="*-$tp.csv" --file-replace-str=-$tp.csv:all_timepoints_ \
#      --plotdirs $bd/simulation/$svsn-$isim/default/selection/partis/single-chain/plots/$ltmp/parameters/hmm/mute-freqs/overall:$bd/simulation/$svsn-$isim/boost/selection/partis/single-chain/plots/$ltmp/parameters/hmm/mute-freqs/overall:$bd/$dvsn/d1-mutd/single-chain/plots/$ltmp/parameters/hmm/mute-freqs/overall \
#      --names "no@boost:boost:data" --colors '#006600:#2b65ec:black' --linewidths 5:3:2 --plottitle="$tp $ltmp" --alphas 0.45:0.45:0.6
# done
# $bin --outdir $od/subtree-purity/tp-shuffled --file-glob-str="*-$tp.csv" --file-replace-strs=-$tp.csv:-iclust-0 --swarm-meta-key timepoints \
#      --plotdirs $bd/$dvsn/d1-mutd/partitions/inferred/subtree-purity:$bd/$dvsn/d1-mutd-shuf/partitions/inferred/subtree-purity \
#      --names "data:shuffled data" --colors 'black:#cc0000' --linewidths 2:3 --plottitle="$tp" --alphas 0.6:0.45 # --xtitle-list subtree@size:mean@dist.@to@ancestor:asdf
for itmpsim in 0 1 2; do
$bin --outdir $od/subtree-purity/simu-boost-vs-no-$itmpsim --file-glob-str="*-$tp.csv" --file-replace-strs=-$tp.csv:-iclust-0 --swarm-meta-key timepoints \
     --plotdirs $bd/simulation/$svsn-$itmpsim/boost/selection/partis/partitions/inferred/subtree-purity:$bd/simulation/$svsn-$itmpsim/default/selection/partis/partitions/inferred/subtree-purity \
     --names "simulation (boost):simulation (no boost)" --colors 'black:#cc0000' --linewidths 2:3 --plottitle="$tp" --alphas 0.6:0.45 # --xtitle-list subtree@size:mean@dist.@to@ancestor:asdf
done

# $bin --outdir $od/subtree-purity/simu --file-glob-str="*-$tp.csv" --file-replace-strs=-$tp.csv:-iclust-0 \
#      --plotdirs $bd/simulation/$svsn-$isim/default/selection/partis/partitions/inferred/subtree-purity:$bd/simulation/$svsn-$isim/boost/selection/partis/partitions/inferred/subtree-purity:$bd/$dvsn/d1-mutd/partitions/inferred/subtree-purity \
#      --names "no@boost:boost:data" --colors '#006600:#2b65ec:black' --linewidths 5:3:2 --plottitle="$tp" --alphas 0.45:0.45:0.6 # --xtitle-list subtree@size:mean@dist.@to@ancestor:asdf

# for ltmp in igh igk igl; do
#     for tp in w2 w23 w25; do
#         subd=single-chain/plots/$ltmp/parameters/hmm/mute-freqs/overall
#         $bin --outdir $od/shm/$ltmp/$tp --file-glob-str="*-$tp.csv" --file-replace-strs=-$tp.csv:all_timepoints_ \
#              --plotdirs $bd/v5/d1/$subd:$bd/v5/d1-mutd/$subd \
#              --names "all@seqs:h+l@mutd@seqs" --colors '#006600:#cc0000' --linewidths 5:3 --plottitle="$tp $ltmp" --alphas 0.45:0.45 --square-bins --no-errors --translegend 0:0.15 --normalize
# done
# done
