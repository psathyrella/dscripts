#!/usr/bin/env python3
import sys
import argparse
import csv
import os
import numpy
import json
import yaml
from yaml import Loader, Dumper
import colored_traceback.always
import collections
from io import open
import copy
import itertools
import re
import subprocess
import hashlib

import partis.utils as utils

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clean_molecules import clean_donor_tsv, load_reference_genes, read_exclusions, split_master_tsv


# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--version', default='test', help='processed/ subdir we write to')
parser.add_argument('--source-version', default='v2', help='processed/ subdir that says which (locus, subject) pairs exist. Has to be a version that already exists: testing --version for this skips everything on a fresh tree.')
parser.add_argument('--sharing-table', help='csv from data-central scripts/vanwinkle_sharing.py, one row per (locus, donor) with a keep/drop decision. Requires --shared-molecules.')
parser.add_argument('--shared-molecules', help='csv from the same script, one row per copy a kept donor holds of a shared molecule. Requires --sharing-table.')
parser.add_argument('--write-meta', action='store_true', help='write samples.yaml pointing at --version, rather than building the tree. Run after the tree exists.')
parser.add_argument('--only-samples', help='colon-separated subject ids to restrict to, for a smoke test before committing to the whole cohort')
parser.add_argument('--overwrite', action='store_true')
args = parser.parse_args()
if (args.sharing_table is None) != (args.shared_molecules is None):
    raise Exception('--sharing-table and --shared-molecules have to be set together: the first says which donors to keep, the second which of a kept donor\'s rows to delete, and neither is any use alone')

base_dir = '/fh/fast/matsen_e/data/vanwinkle-170'
outdir = '%s/processed/%s' % (base_dir, args.version)
glfo_dir = '%s/germlines/processed' % base_dir  # NOTE deleted by datascripts 6ecbde1 while keeping the use below, so this crashed with a NameError

sample_ids = ['110040277', '110040504', '110040811', '110040851', '110040887', '110040958', '110040975', '110042592', '110042606', '110042713', '110042813', '110042935', '110044118', '110044155', '110044355', '110044376', '110044782', '110044786', '110044845', '110044978', '110046119', '110046337', '110046617', '110046912', '222693519', '22301641N', '888174491', '888179544', '888214522', '888262983', '888286853', '888287118', '888313858', '888333863', '888402542', '888407672', '888408987', '888436895', '888508776', '888527032', '888559379', '888581247', '888624327', '888634432', '888640612', '888646440', '888647920', '888669369', '888696188', '888703260', '888760853', '888777004', '888816126', '888921048', '888927334', '888943370', '888945988', '888956204', '888976935', '888997753', '889005836', '889018176', 'CE0005500', 'CE0005604', 'CE0005662', 'CE0005779', 'CE0005786', 'CE0005797', 'CE0005802', 'CE0005876', 'CE0005924', 'CE0005988', 'CE0006103', 'CE0006205', 'CE0006233', 'CE0006371', 'CE0006394', 'CE0006435', 'CE0006481', 'CE0006524', 'CE0006547', 'CE0006623', 'CE0006705', 'CE0006781', 'CE0006829', 'CE0006943', 'CE0007007', 'CE0007070', 'CE0007129', 'CE0007408', 'CE0007424', 'CE0007624', 'CE0007765', 'CE0007769', 'CE0007814', 'CE0007821', 'CE0007855', 'CE0007865', 'CE0007908', 'CE0007925', 'CE0008114', 'CE0008136', 'CE0008202', 'CE0008264', 'CE0008375', 'CE0008411', 'CE0008484', 'CE0008542', 'CE0008593', 'CE0008661', 'CE0008682', 'CE0008817', 'CE0008924', 'CE0009275', 'CE0009991', 'RG1083', 'RG1193', 'RG1200', 'RG1234', 'RG1265', 'RG1272', 'RG1288', 'RG1289', 'RG1367', 'RG1443', 'RG1565', 'RG1578', 'RG1791', 'RG1802', 'RG1822', 'RG1823', 'RG1849', 'RG1889', 'RG1895', 'RG1936', 'RG2037', 'RG2105', 'RG2115', 'RG2118', 'RG2155', 'RG2195', 'RG2207', 'RG2226', 'RG2251', 'RG2284', 'RG2350', 'RG2379', 'RG2459', 'RG2512', 'RG2617', 'RG2674', 'RG2685', 'RG2838', 'RG2903', 'RG3105', 'RG3309', 'RG3355', 'RG3444', 'RG3549', 'RG3556', 'RG3578', 'RG3666', 'RG3758', 'RG3783', 'RG3844', 'RG3945', 'RG3948', 'RG3956', 'RG3987', 'RG3995', 'RG4094', 'RG4543', 'SCT2601343', 'SCT2794524', 'SCT4168329', 'SCT6222960', 'SCT6338383', 'SCT8789722', 'SCT8995047', 'SCT9160998']

if args.sharing_table is None:
    print('  %s no --sharing-table, so writing every donor and every row unfiltered' % utils.color('yellow', 'note'))
    dropped, delete_keys = set(), None
else:
    dropped, delete_keys = read_exclusions(args.sharing_table, args.shared_molecules)

metafos = collections.OrderedDict()
partition_args = ['--fast']
for ltmp in ['igh', 'igk', 'igl']:
    print(utils.color('blue_bkg', ltmp))
    tfn = '%s/raw/master_%s_changeo.tsv' % (base_dir, ltmp.upper())

    # # command to get all sample ids (this might need updating, but in any case I just ran it once and copied the results to <sample_ids> above so I don't have to wait on every file):
    # sample_ids, errstr = utils.simplerun('cut -f1 %s | sort | uniq' % tfn, shell=True, return_out_err=True, debug=False, dryrun=True)
    # sample_ids = [s.strip() for s in sample_ids.strip().split('\n')]

    cmdfos = []
    to_split = collections.OrderedDict()  # subject id -> per-subject tsv to write
    handled = set()  # every kept (locus, donor) has to end up in here: three separate continues below can drop one, and doing that silently would quietly shrink the cohort
    for sid in sample_ids:
        if args.only_samples is not None and sid not in args.only_samples.split(':'):
            continue
        sfn = '%s/%s/%s.tsv' % (outdir, ltmp, sid)
        sampledir = '%s/%s/%s' % (outdir, ltmp, sid)
        ffn = '%s/seqs.fa' % sampledir
        srcfn = '%s/processed/%s/%s/%s.tsv' % (base_dir, args.source_version, ltmp, sid)
        if not os.path.exists(srcfn):  # skip missing input (some subjects are missing loci). NOTE has to test --source-version: testing <sfn> (under --version) skips every sample on a fresh tree, and the script then exits 0 having written nothing
            continue
        if (ltmp, sid) in dropped:  # over the sharing threshold for this locus
            continue
        if args.write_meta:
            if os.path.exists(ffn):
                metafos['%s-%s'%(sid, ltmp)] = {'infname' : ffn, 'locus' : ltmp, 'species' : 'human', 'subject' : sid,
                                                'extra-args' : {'partition' : partition_args, 'cache-parameters' : ['--leave-default-germline', '--initial-germline-dir', '%s/germline-sets' % sampledir, '--crash-on-duplicate-uids']}}  # NOTE clean_molecules makes every uid unique, so a rename here would mean that silently failed (its --help says 'not intended for direct use', which is about the paired-loci parent, not about using it as a guard)
            else:
                print('    skipping missing file %s' % sfn)
            continue
        if utils.output_exists(args, ffn):
            handled.add(sid)  # written by a previous run
            continue
        to_split[sid] = sfn

    if len(to_split) > 0:  # one pass over the master tsv for all of them, rather than one pass per subject
        for sfn in to_split.values():
            utils.mkdir(sfn, isfile=True)
        split_master_tsv(tfn, to_split, load_reference_genes(glfo_dir, ltmp))
    for sid, sfn in to_split.items():
        if delete_keys is not None:
            clean_donor_tsv(sfn, delete_keys[(ltmp, sid)])
        sampledir = '%s/%s/%s' % (outdir, ltmp, sid)
        utils.mkdir(sampledir)
        handled.add(sid)
        cmdfos.append({
            'cmd_str' : 'parse-output.py %s %s --airr-input --glfo-dir %s --locus %s' % (sfn, '%s/seqs.fa' % sampledir, glfo_dir, ltmp),
            'outfname' : '%s/seqs.fa' % sampledir,
            'workdir' : sampledir,
        })
    if delete_keys is not None and not args.write_meta and args.only_samples is None:
        expected = set(d for l, d in delete_keys if l == ltmp)
        if handled != expected:
            raise Exception('%s: %d kept donors in the sharing table but %d written (missing %s)' % (ltmp, len(expected), len(handled), ' '.join(sorted(expected - handled)) or 'none'))
    if len(cmdfos) > 0:
        utils.run_cmds(cmdfos, n_max_procs=utils.auto_n_procs(), debug='write')

if args.write_meta:
    with open('datascripts/meta/vanwinkle-170/samples.yaml', 'w') as mfile:
        yaml.dump(metafos, mfile, Dumper=Dumper, default_flow_style=False)
