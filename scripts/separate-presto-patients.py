#!/usr/bin/env python
import argparse
import subprocess
import glob
import csv
import numpy
import operator
import contextlib
import pysam
import os
import sys
import itertools

indir = '/fh/fast/matsen_e/data' #/jason-mg-2017-02-01/processed'
hstrs = 'adgme'

parser = argparse.ArgumentParser()
parser.add_argument('study')
args = parser.parse_args()

if args.study != 'jason-influenza':
    raise Exception('undhandled study (actually maybe jason-mg is ok as well?)')

# ----------------------------------------------------------------------------------------
def get_locus(prcons):
    if 'jason-mg' in args.study:  # in the mg samples this is Ig[ADGMEKL]
        pr = prcons.lower().replace('ig', '')
        if pr in hstrs:
            return 'igh'
        else:
            assert pr == 'k' or pr == 'l'
            return 'ig' + pr
    elif 'jason-influenza' in args.study:
        # but in the influenza it's like
        #    heavy: p5-hIGHC-{A1,A2,D,E,G,M}_bs-0
        #    light: p5-hIG{K,L}C_bs-0
        if prcons[:6] == 'p5-hIG':
            return 'ig' + prcons[6].lower()
        else:
            assert False
    else:
        assert False

# ----------------------------------------------------------------------------------------
def write_sample(sample, info):
        strlist = sample.split('-')
        if len(strlist) == 2:
            subject, locusstr, tpstr = strlist + [None]
        else:
            subject, locusstr, tpstr = strlist
        print '  %5s %5s  %5s   %d lines' % (subject, locusstr, tpstr if tpstr is not None else '', len(info))
        outfname = indir + '/' + args.study + '/processed/patients/' + sample + '.csv'
        if os.path.exists(outfname):
            raise Exception('wtf %s' % outfname)
        subprocess.check_call(['mkdir', '-p', os.path.dirname(outfname)])
        with open(outfname, 'w') as outfile:
            writer = csv.DictWriter(outfile, header)
            writer.writeheader()
            for line in info:
                writer.writerow(line)

# ----------------------------------------------------------------------------------------
for fname in glob.glob(indir + '/' + args.study + '/processed/*.tsv'):
    samples = {}
    subj_samples = {}  # same, but with timepoints merged together
    header = None
    print fname
    with open(fname) as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        for line in reader:
            if line['PRCONS'] == 'NA':
                continue
            if header is None:
                header = line.keys()
            lstr = get_locus(line['PRCONS'])  # I don't know what PRCONS stands for, but it's the code for locus, isotype, plus other stuff
            if 'HeavyChain' in fname and lstr != 'igh':  # there's like one igl in the heavy chain file, no I don't fucking no why
                print 'fuck it'
                continue
            sample = line['SUBJECT'] + '-' + lstr
            if 'TIME_POINT' in line:
                tpstr = line['TIME_POINT'].lower().replace('-', 'm').replace('+', 'p')
                sample += '-' + tpstr
            if sample not in samples:
                samples[sample] = []
            samples[sample].append(line)

            subj_sample = line['SUBJECT'] + '-' + lstr
            if subj_sample not in subj_samples:
                subj_samples[subj_sample] = []
            subj_samples[subj_sample].append(line)

#             if len(samples[sample]) > 100:
#                 break

    # write each individiual time point
    for sample, info in samples.items():
        write_sample(sample, info)

    # write each subject with all timepoints together
    for subj_sample, info in subj_samples.items():
        write_sample(subj_sample, info)
