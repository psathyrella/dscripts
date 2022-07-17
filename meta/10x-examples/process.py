#!/usr/bin/env python
import sys
import csv
import os
import numpy
import glob

# NOTE run from inside main partis dir
partis_dir = os.getcwd()
sys.path.insert(1, partis_dir + '/python')
import utils

# ----------------------------------------------------------------------------------------
version = 'v1'
bdir = '/fh/fast/matsen_e/data/10x-examples'
odir = '%s/processed-data/%s' % (bdir, version)

for ifname in glob.glob('%s/data/*.fasta'%bdir):
    dset = utils.getprefix(os.path.basename(ifname))
    mfname = '%s/%s/meta.yaml' % (odir, dset)
    spstr = ' --species mouse' if 'balbc' in ifname else ''
    utils.simplerun('./bin/extract-pairing-info.py %s %s' % (ifname, mfname), logfname='%s/%s/extract-pairing.log'%(odir, dset)) #, dryrun=True)
    utils.simplerun('./bin/split-loci.py %s --overwrite --outdir %s/%s --input-metafname %s%s' % (ifname, odir, dset, mfname, spstr), logfname='%s/%s/split-loci.log'%(odir, dset)) #, dryrun=True)
