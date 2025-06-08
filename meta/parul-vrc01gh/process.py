#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
import sys
import csv
import os
import numpy
import json
import colored_traceback.always
from io import open

# NOTE run from inside main partis dir
partis_dir = os.getcwd()  # os.path.dirname(os.path.realpath(__file__)).replace('/datascripts/meta/goo-dengue-10x', '')
sys.path.insert(1, partis_dir)
import python.utils as utils

dvsn = '2022-06-16' #'2022-08-05'
ovsn = 'v3'
bdir = '/fh/fast/matsen_e/data/parul-vrc01gh'
odir = '%s/processed/%s' % (bdir, ovsn)
timepoints = ['w2', 'w23', 'w25']

# ----------------------------------------------------------------------------------------
def getname(iname):
    namestrs = [n.upper() for nstr in iname.rstrip(':').split('_') for x in nstr.split('-') for n in x.split('.')]
    for ins, nstr in enumerate(namestrs):
        if nstr in ['GLVRC01', 'GLVRCO1', 'GLVCR01', 'KAPPA', 'MVKAPPA', '5MVKAPPA', 'AB1']:  # yes two are just typos
            namestrs[ins] = ''  # will remove after loop
        # if nstr in ['hc', 'lc']:
    revstr = None
    if 'R' in namestrs:
        revstr = 'R'
        namestrs[namestrs.index('R')] = ''
    namestrs = [n for n in namestrs if n!='']
    if len(namestrs) != 4:
        raise Exception('len not 4: %s' % namestrs)
    datestr = namestrs[0]
    if datestr[0] != '0':
        datestr = '0' + datestr
    if namestrs[1] in ['HC', 'LC']:  # there's two different orderings
        chainstr = namestrs[1]
        platenum = namestrs[2]
        wellnum = namestrs[3]
    elif namestrs[3] in ['HC', 'LC']:
        platenum = namestrs[1]
        wellnum = namestrs[2]
        chainstr = namestrs[3]
    elif namestrs[3] == 'K':
        platenum = namestrs[1]
        wellnum = namestrs[2]
        chainstr = namestrs[3]
    else:
        raise Exception('couldn\'t parse %s (%s not in {HL}C' % (iname, namestrs))
    rname = '-'.join([datestr, platenum, wellnum, chainstr])
    if revstr is not None:
        rname += '-%s' % revstr
    return rname

# ----------------------------------------------------------------------------------------
translations = {}
all_seqfos, all_metafos = [], {}
for tp in timepoints:
    ifn = '%s/data/%s/%s.fa' % (bdir, dvsn, tp)
    seqfos = utils.read_fastx(ifn)
    for sfo in seqfos:
        new_name = getname(sfo['name'])
        if new_name in all_metafos:
            print('already there: %s' % translations[new_name])
            print('now: %s' % sfo['name'])
            raise Exception()
        translations[new_name] = sfo['name']
        sfo['name'] = new_name
        if '-LC' in new_name or '-K' in new_name:
            sfo['seq'] += 'TTCGGTGGAGGCACCAAGCTGGAAATCAAAC'  # in light chain, they only sequence to end of v, so add j1*01 to each sequence: TTCGGTGGAGGCACCAAGCTGGAAATCAAAC (i.e. j 5' deletion is GTGGACG)
        all_metafos[sfo['name']] = {'timepoint' : tp}
    all_seqfos += seqfos

print('  writing %d timepoints to %s' % (len(timepoints), odir))
utils.write_fasta('%s/all-seqs.fa'%odir, all_seqfos)
with open('%s/meta.yaml'%odir, 'w') as mfile:
    json.dump(all_metafos, mfile)
