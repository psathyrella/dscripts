#!/usr/bin/env python
import sys
import operator
import os
import csv
csv.field_size_limit(sys.maxsize)
import subprocess
import collections
from Bio.Seq import Seq

sys.path.insert(1, './python')
import utils
import glutils
from clusterpath import ClusterPath

min_size = 50

# ----------------------------------------------------------------------------------------
def get_gldir(sample):
    return '/fh/fast/matsen_e/processed-data/partis/%s/%s/%s/hmm/germline-sets' % (study, version, sample)

# ----------------------------------------------------------------------------------------
def get_outfname(cdr3_seq=None):
    cmd = './datascripts/run.py %s --study %s --sample %s --extra-str=%s --logstr=%s --only-merged --logfnames' % (action, study, sample, version, logstr)
    if cdr3_seq is not None:
        cmd += ' --seed-uids %s' % cdr3_seq
    logfname = subprocess.check_output(cmd.split()).strip()
    outfname = logfname.replace('.log', '.csv')
    if not os.path.exists(outfname):
        print '  output file %s don\'t exist' % outfname
        assert False
    return outfname, cmd

# ----------------------------------------------------------------------------------------
def get_naive_seq(requested_aa_cdr3_seq, clusterfo, fname, glfo, cmd, debug=False):
    requested_uids = clusterfo['unique_ids'].split(':')
    cpath = ClusterPath()
    cpath.readfile(fname)
    overlapping_clusters = []
    # print '        size     overlap'
    for cluster in cpath.partitions[cpath.i_best]:
        if len(cluster) < min_size:
            continue
        overlap_ids = set(cluster) & set(requested_uids)
        if len(overlap_ids) == 0:
            continue
        # print '    %5d    %3d' % (len(cluster), len(overlap_ids))
        overlapping_clusters.append((len(overlap_ids), cluster))

    if len(overlapping_clusters) == 0:
        print '  didn\'t find an overlappy cluster'
        return 'XXX'
    overlapping_clusters = sorted(overlapping_clusters, reverse=True)
    if debug:
        print '      %4d  %s                  %s' % (len(requested_uids), glfo['locus'], cmd)
        for n_overlap, cluster in overlapping_clusters:
            print '                %4d  / %d' % (n_overlap, len(cluster))
    most_overlappy_cluster = overlapping_clusters[0][1]

    annotation = None
    with open(outfname.replace('.csv', '-cluster-annotations.csv')) as afile:
        reader = csv.DictReader(afile)
        for line in reader:
            if line['unique_ids'].split(':') == most_overlappy_cluster:
                annotation = line
                break
    if annotation is None:
        print '  couldn\'t find annotation'
        return 'XXX'

    utils.process_input_line(annotation)
    utils.add_implicit_info(glfo, annotation)

    tmpx, _ = utils.subset_sequences(annotation, iseq=0, restrict_to_region='cdr3')
    naive_aa_cdr3_seq = Seq(tmpx).translate()
    if naive_aa_cdr3_seq != requested_aa_cdr3_seq:
        colored_result = utils.color_mutants(requested_aa_cdr3_seq, naive_aa_cdr3_seq, amino_acid=True)
        print '                              cdr3s don\'t match: %s  %s' % (requested_aa_cdr3_seq, colored_result)
    if annotation['naive_seq'] != clusterfo['naive_seq']:
        colored_result = utils.color_mutants(annotation['naive_seq'], clusterfo['naive_seq'], align=True)
        print '                              naive seqs don\'t match: %s' % annotation['naive_seq']
        print '                                                      %s' % colored_result

    return annotation['naive_seq']

# logstr = 'isub-0'; action = 'partition'
logstr = 'computer-seeds'; action = 'seed-partition'
version = 'v18'
study = 'bf520-synth'
infname = 'datascripts/seeds/%s/families.csv' % study

familyfo = collections.OrderedDict()
with open(infname) as infile:
    reader = csv.DictReader(infile)
    for line in reader:
        assert line['study'] == study
        if line['sample'] not in familyfo:
            familyfo[line['sample']] = collections.OrderedDict()
        familyfo[line['sample']][line['cdr3_seq']] = line  # assumes only one instance of this cdr3 seq for each sample... which is probably ok

print '   requested'
print '      uids     overlap / (found total)'
for sample, samplefo in familyfo.items():
    if action == 'partition':  # same for all of 'em
        outfname = get_outfname()

    for cdr3_seq, clusterfo in samplefo.items():
        # if clusterfo['locus'] != 'igh':
        #     continue
        outfname, cmd = get_outfname(cdr3_seq)
        glfo = glutils.read_glfo(get_gldir(sample), locus=clusterfo['locus'])
        clusterfo['naive_seq'] = get_naive_seq(cdr3_seq, clusterfo, outfname, glfo, cmd, debug=True)
        # sys.exit()
