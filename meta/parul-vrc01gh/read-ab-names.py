#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
import sys
import csv
import os
import numpy
import json
import yaml
from yaml import CLoader as Loader, CDumper as Dumper
import colored_traceback.always
from io import open

# NOTE run from inside main partis dir
partis_dir = os.getcwd()  # os.path.dirname(os.path.realpath(__file__)).replace('/datascripts/meta/goo-dengue-10x', '')
sys.path.insert(1, partis_dir)
import python.utils as utils

# reads fasta input files and ab name file (abs that were made+tested), outputs yaml meta file
#  - need to read fasta input files because they didn't name the seqs in a consistent way (especially the chain stuff)
#  - example: 033021-P1-C11-HC-igh

tloci = ['igh', 'igk']
subjects = ['d1']
metadir = 'datascripts/meta/parul-vrc01gh'
cfn = '%s/ab-names.csv' % metadir
ivsn = 'v2'
fadir = '/fh/fast/matsen_e/data/parul-vrc01gh/processed/%s/split-loci' % ivsn

all_uids = set()
for ltmp in tloci:
    seqfos = utils.read_fastx('%s/%s.fa' % (fadir, ltmp))
    for sfo in seqfos:
        all_uids.add(sfo['name'])

missing_ids, n_tot = [], 0
seedfos = {l : [] for l in tloci}  # maybe not actually seeds any more
with open(cfn) as cfile:
    reader = csv.DictReader(cfile)
    for line in reader:  # loop over abs that were chosen/tested
        n_tot += 1
        basename = '%s-%s-%s' % (line['date'], line['plate'].upper(), line['well'].upper())
        bids = [u for u in all_uids if basename+'-' in u]  # find uids from input fasta with the same basename
        if len(bids) != 2:
            missing_ids.append(basename)
            continue
        for bid, other_id in zip(bids, reversed(bids)):
            ltmp = bid.split('-')[-1]
            if ltmp not in tloci:
                raise Exception('unexpected locus \'%s\'' % ltmp)
            seedfos[ltmp].append({'uid' : bid, 'seq' : None, 'paired-uid' : other_id, 'alternate-uid' : line['name']})

if len(missing_ids) > 0:
    print('    %s couldn\'t find uids for %d / %d seed base names: %s' % (utils.wrnstr(), len(missing_ids), n_tot, ' '.join(sorted(missing_ids))))

# sfn = '%s/seeds.yaml' % metadir
# print('  writing seeds to %s' % sfn)
# jfo = {'subjects' : subjects, 'common' : seedfos}
# utils.jsdump(sfn, jfo)

mfn = '%s/meta-chosen.yaml' % metadir
print('  writing chosen meta for %d abs to %s' % (len(seedfos[tloci[0]]), mfn))
jfo = {s['uid'] : {'chosen' : True, 'alternate-uid' : s['alternate-uid']} for lsfos in seedfos.values() for s in lsfos}
utils.jsdump(mfn, jfo)
