#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
import sys
import csv
import os
import numpy
import json
import colored_traceback.always
import itertools
from io import open
from six.moves import zip

# NOTE run from inside main partis dir
partis_dir = os.getcwd()  # os.path.dirname(os.path.realpath(__file__)).replace('/datascripts/meta/goo-dengue-10x', '')
sys.path.insert(1, partis_dir)
import python.utils as utils
import paircluster

raise Exception('if you rerun this, you need to write .fa files to paired subdir (igh+igk/) at bottom below, for now i\'m just making links by hand')
ivsn = 'v2'  # input version
pvsn = 'v5'  # partis output version
ovsn = 'v0'  # output version from this script
ibd = '/fh/fast/matsen_e/data/parul-vrc01gh'
indir = '%s/processed/%s/split-loci' % (ibd, ivsn)
pdir = '/fh/fast/matsen_e/processed-data/partis/parul-vrc01gh/%s/d1' % pvsn
odir = '%s/pseudo/%s' % (ibd, ovsn)
timepoints = ['w2', 'w23', 'w25']
lpairs = [['igh', 'igk']]

def fnfcn(l, lpair=None): return paircluster.paired_fn(pdir, l, lpair=lpair,  actstr='partition', suffix='.yaml')
lp_infos = paircluster.read_lpair_output_files(lpairs, fnfcn, debug=True)
antn_pairs = []
for lpair in lpairs: #[lpk for lpk in utils.locus_pairs[ig_or_tr] if tuple(lpk) in lp_infos]:
    antn_pairs += paircluster.find_cluster_pairs(lp_infos, lpair) #, min_cluster_size=min_cluster_size)  # , required_keys=['tree-info']

with open('%s/meta.yaml'%indir) as mfile:
    metafos = json.load(mfile)
all_timepoints = []
for mfo in metafos.values():
    all_timepoints.append(mfo['timepoint'])
print('    read timepoint info from meta file (timepoint count): %s' % '  '.join('%s %d' % (tp, len(list(group))) for tp, group in itertools.groupby(sorted(all_timepoints))))

mutated_uids, mutated_seqfos = set(), []
shuffled_seqfos = []
n_mut_info = {s : 0 for s in [l for lp in lpairs for l in lp] + ['h+l']}
debug = False
if debug:
    print('         N muts')
    print('        h   l   tot')
for iclust, (h_atn, l_atn) in enumerate(antn_pairs):
    for tln, oln in zip((h_atn, l_atn), reversed((h_atn, l_atn))):
        if debug:
            print('  --->')
        for iseq, (uid, tseq) in enumerate(zip(tln['unique_ids'], tln['input_seqs'])):
            if uid in mutated_uids:
                continue
            sfos = [{'name' : uid, 'seq' : tseq}]
            n_muts = tln['n_mutations'][iseq]
            pid = None  # just for debug
            if len(tln['paired-uids'][iseq]) > 0:
                pid = utils.get_single_entry(tln['paired-uids'][iseq])
                assert pid in oln['unique_ids']
                assert pid not in mutated_uids  # if it is in there, <uid> should also be, so we shouldn't be able to get here
                sfos.append({'name' : pid, 'seq' : utils.per_seq_val(oln, 'input_seqs', pid)})
                n_muts += utils.per_seq_val(oln, 'n_mutations', pid)
            if n_muts > 0:
                for sfo in sfos:
                    assert sfo['name'] not in mutated_uids
                    mutated_seqfos.append(sfo)
                    mutated_uids.add(sfo['name'])
                # print tln['loci'][iseq] if pid is None else 'h+l', 1 if pid is None else 2
                n_mut_info[tln['loci'][iseq] if pid is None else 'h+l'] += 1 # if pid is None else 2
            if debug:
                print('   %s %3d  %3s  %3s  %25s  %25s' % ('  ' if n_muts > 0 else utils.color('blue', '->'), tln['n_mutations'][iseq], ' ' if pid is None else n_muts - tln['n_mutations'][iseq], ' ' if pid is None else n_muts, uid, ' ' if pid is None else pid))

            new_timepoint = numpy.random.choice(all_timepoints)
            metafos[uid]['timepoint'] = new_timepoint
            if pid is not None:
                metafos[pid]['timepoint'] = new_timepoint
    # utils.print_reco_event(tln, extra_print_keys=['paired-uids'])
    # utils.print_reco_event(oln)
    # break
assert len(mutated_uids) == len(mutated_seqfos)
print('  found %d seqs with at least one h+l mutation: %s %d (pairs)   %s %d   %s %d'  % (len(mutated_seqfos), 'h+l', n_mut_info['h+l'], 'igh', n_mut_info['igh'], 'igk', n_mut_info['igk']))

mdir = '%s/mutated' % odir
print('  writing mutated seqs to %s' % mdir)
utils.mkdir(mdir)
for ltmp in sorted(set([l for lp in lpairs for l in lp])):
    outfos = utils.read_fastx('%s/%s.fa' % (indir, ltmp))
    n_before = len(outfos)
    outfos = [o for o in outfos if o['name'] in mutated_uids]
    print('  %s %d --> %d' % (ltmp, n_before, len(outfos)))
    utils.write_fasta('%s/%s.fa' % (mdir, ltmp), outfos)
utils.simplerun('cp %s/meta.yaml %s/' % (indir, mdir))

shdir = '%s/shuffled' % odir
print('  writing timepoint-shuffled seqs to %s' % shdir)
utils.mkdir(shdir)
for ltmp in sorted(set([l for lp in lpairs for l in lp])):
    utils.simplerun('cp %s/%s.fa %s/' % (indir, ltmp, shdir))
with open('%s/meta.yaml'%shdir, 'w') as mfile:
    json.dump(metafos, mfile)

msdir = '%s/mutd-shufd' % odir
print('  writing mutated + timepoint-shuffled seqs to %s' % msdir)
utils.mkdir(msdir)
for ltmp in sorted(set([l for lp in lpairs for l in lp])):
    utils.simplerun('cp %s/%s.fa %s/' % (mdir, ltmp, msdir))
utils.simplerun('cp %s/meta.yaml %s/' % (shdir, msdir))
