#!/usr/bin/env python
import math
import numpy
import yaml
import random
import tempfile
import colored_traceback.always
import collections
import argparse
import sys
import csv
csv.field_size_limit(sys.maxsize)
import glob
import time
import sys
from subprocess import check_call, check_output, Popen
import os

script_dir = os.path.dirname(os.path.realpath(__file__))

# ----------------------------------------------------------------------------------------
example_str = '\n    '.join(['merge-samples: by default merges all samples (timepoints and replicates) of each isotype for each subject, but see --dont-merge-timepoints',
                             'example usage:',
                             './run.py cache-parameters --study kate-qrs --samples 1g:2l:4k',
                             './run.py seed-partition --study kate-qrs --samples 1g',
                             './run.py partition --study katie --samples vmo:11303'])

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, epilog=example_str)
partis_actions = {  # correspondence between run.py actions and partis actions (default: same)
    'seed-partition' : 'partition',
}

def pact(act):
    return partis_actions.get(act, act)
all_actions = ['cache-parameters', 'partition', 'seed-partition', 'check-seeds-in-input-files', 'process-vlad-data', 'merge-samples', 'convert-csv-meta-to-yaml', 'simulate']  # NOTE args like --write-seed-summary-table and --get-naive-probabilities below are kind of like actions
all_actions += ['alternate-seed-naive-seqs']  # DEPRECATED
parser.add_argument('action', choices=all_actions)

# ----------------------------------------------------------------------------------------
def iostr(io):
    return ('--paired-%sdir'%io) if args.paired_loci else '--%sfname'%io

# ----------------------------------------------------------------------------------------
# action-like arguments:
def are_we_reading_existing_output(args):  # are we actually running new stuff, or just running on existing output?
    return args.view_ascii or args.plot_partitions or args.merge_paired_partitions or args.get_naive_probabilities or args.write_cluster_summary_tables or args.get_subst_probabilities or args.get_selection_metrics or args.view_alternative_annotations or args.update_meta_info
parser.add_argument('--write-cluster-summary-tables', action='store_true', help='write two csv files to the output dir, each summarizing clusters (family size and rank, gene calls, shm, etc.). One details the seed cluster for each seed sequence, while the other describes the top clusters by size. Use with either \'partition\' or \'seed-partition\', and note that the results mean *very* different things for each.')  # for search: print_seed_summary_table print seed summary
parser.add_argument('--get-naive-probabilities', action='store_true', help='runs partis/bin/get-naive-probabiliites.py (see that script for details)')
parser.add_argument('--get-subst-probabilities', action='store_true', help='write substitution probabilites for each posiiton in the seed sequence plus any other sequences specified with --inferred-intermediate-fname')
parser.add_argument('--plot-partitions', action='store_true')
parser.add_argument('--merge-paired-partitions', action='store_true')
parser.add_argument('--get-selection-metrics', action='store_true', help='get tree metrics on existing partition output (if you instead want to get them *while* partitioning, add --get-selection-metrics to \'extra-args\' in the the sample\'s yaml metafo)')
parser.add_argument('--view-alternative-annotations', action='store_true')
parser.add_argument('--view-alternative-naive-seqs', action='store_true', help='DEPRECATED')
parser.add_argument('--update-meta-info', action='store_true', help='')

# ----------------------------------------------------------------------------------------
parser.add_argument('--study', required=True)
parser.add_argument('--n-procs', type=int, default=1)
parser.add_argument('--samples')
parser.add_argument('--subjects')
parser.add_argument('--loci')
parser.add_argument('--isotypes')
parser.add_argument('--timepoints')
parser.add_argument('--only-umid', action='store_true')  # if you want only non-umid, you'll have to specify the samples by hand
parser.add_argument('--version', default='latest', help='string to distinguish entire runs, i.e. stuff is all put in a new subdir with name <--version>')
parser.add_argument('--logstr', default='', help='extra string to append (with prepended \'-\') to log and output names to distinguish different runs with e.g. different command line args. Generally assumes that you want to use the original (no --logstr) parameter dir, although if set *during* cache-parameters, will also be added to parameter dir (although i think not subsequently used for inference, which yes is kind of silly, but it\'s not atm intended to be used that way). When running merge-samples, this gets appended to the merged sample name')
parser.add_argument('--ex-out-version', default='', help='extra string for different runs of existing output actions, e.g. selection metrics and partition plotting. Feeds through to e.g. log file and chosen ab file name.')
parser.add_argument('--dry-run', action='store_true', help='print commands and exit without running (note that you can also set --extra-args="--dry-run" to dry-run partis paired loci commands')
parser.add_argument('--rm', action='store_true')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--force-csv-output', action='store_true')
parser.add_argument('--rm-zero-length', action='store_true')
parser.add_argument('--logfnames', action='store_true')
parser.add_argument('--no-slurm', action='store_true')
parser.add_argument('--force-default-initial-germline-dir', action='store_true')
parser.add_argument('--start-n-max-and-exit', action='store_true')
parser.add_argument('--check', action='store_true')
parser.add_argument('--view-ascii', action='store_true')
parser.add_argument('--only-print-seed-cluster', action='store_true')
parser.add_argument('--only-print-queries-to-include-clusters', action='store_true')
parser.add_argument('--naive-prob-config-fname', help='yaml configuration file for naive probability calculations')
parser.add_argument('--inferred-intermediate-fname', help='fasta file for --get-subst-probabilities with inferred sequences from the phylo program for which we\'ll need to get annotations (file should also contain the sed sequence and naive [?])')
parser.add_argument('--write-yaml', action='store_true')
parser.add_argument('--n-max-jobs', type=int)
parser.add_argument('--sleep', type=int)
parser.add_argument('--n-random-queries', type=int, help='see partis help for --n-random-queries')
parser.add_argument('--n-random-subsets', type=int, help='run this many independent runs (you must also set --n-random-queries), incrementing --random-seed by one each time, and adding the corresponding string to each --logstr')
parser.add_argument('--n-max-queries', type=int)
parser.add_argument('--n-jobs', type=int, default=0, help='for internal use only, don\'t set this')  # not really an argument
parser.add_argument('--random-seed-seqs', type=int, help='choose this many seed sequences at random from the input file')
parser.add_argument('--random-seed', type=int, help='passed to partis')
parser.add_argument('--no-indels', action='store_true', help='set partis --no-indels')
parser.add_argument('--extra-args', help='add extra partis arguments by hand here')
parser.add_argument('--print-width', type=int, default=300, help='set to 0 for infinite')
parser.add_argument('--partis-dir', default=os.getcwd(), help='path to main partis install dir')  # os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
parser.add_argument('--no-merged', action='store_true')
parser.add_argument('--only-merged', action='store_true')
parser.add_argument('--no-simu', action='store_true')
parser.add_argument('--dont-merge-timepoints', action='store_true', help='when running \'merge-samples\', only merge samples with the same timepoint, i.e. don\'t merge different timepoints')
parser.add_argument('--dont-merge-single-samples', action='store_true', help='when running \'merge-samples\', by default we (now) make a merged sample even when there\'s only one component sample (we used to skip it, and setting this goes back to doing that)')
parser.add_argument('--seed-origin')
parser.add_argument('--other-method', choices=['tigger-default', 'igdiscover'])
parser.add_argument('--small-clusters-to-ignore')
parser.add_argument('--test', action='store_true', help='only used for timepoint merging NOTE: use this if you want to write merged meta info without rewriting the merged files')
parser.add_argument('--write-meta', action='store_true', help='only used for preprocessing NOTE merged meta info is written automatically when running merge-samples')
parser.add_argument('--synth-component-studies', help='only used when action is process-vlad-data for {}-synth studies (all laura studies: kate-qrs:laura-mb:laura-mb-2:dana-qrs:qa255-jul-18:bf520-jul-18:qa013-mar-8, all laura subjects: BF520:BG505:MG505:QA255:QB850:QA013)')
parser.add_argument('--no-plots', action='store_true', help='turn off plotting, i.e. don\'t pass in a --plotdir')
parser.add_argument('--seed-uids')
parser.add_argument('--paired-loci', action='store_true', help='run with paired locus input + output (see partis arg of same name)')
parser.add_argument('--simulate-from-scratch', action='store_true')
parser.add_argument('--extra-str', help='DEPRECATED')
parser.add_argument('--write-seed-summary-table', action='store_true', help='DEPRECATED use --write-cluster-summary-tables')

args = parser.parse_args()
if args.write_seed_summary_table:
    raise Exception('use new arg --write-cluster-summary-tables')

if args.extra_str is not None:
    print '  note: transferring value of deprecated arg --extra-str to --version'
    args.version = args.extra_str
    args.extra_str = None
if args.paired_loci and args.action != 'simulate' and args.logstr != '':
    raise Exception('can\'t set --logstr if --paired-loci is set (since the whole idea of --logstr is to run different versions on the same parameters, but with --paired-loci they all go in the same dir). Best way is probably to mkdir and ln the parameters by hand, then run a new --version')
if args.action == 'alternate-seed-naive-seqs':
    raise Exception('action \'%s\' is deprecated, use --view-alternative-annotations instead' % args.action)
if args.view_alternative_naive_seqs:
    print '  note: transferring deprecated option --view-alternative-naive-seqs to --view-alternative-annotations'
    args.view_alternative_naive_seqs = args.view_alternative_annotations
    delattr(args, 'view_alternative_naive_seqs')

if args.start_n_max_and_exit and args.n_max_jobs is None:
    raise Exception('have to set --n-max-jobs if you\'re setting --start-n-max-and-exit')

if not os.path.exists(args.partis_dir):
    print 'WARNING current script dir %s doesn\'t exist, so python path may not be correctly set' % args.partis_dir
sys.path.insert(1, args.partis_dir + '/python')
import utils
import seqfileopener
from clusterpath import ClusterPath
import processargs
import paircluster

import heads
from yamlwriter import YamlWriter
import preprocess

args.metafo = None
if not args.write_meta:
    try:
        args.metafo = heads.read_metadata(args.study)
        args.seedfos = heads.read_seed_info(args.study)
    except IOError as e:
        print e
        raise Exception('couldn\'t find meta and/or seed info (see previous line) -- maybe need to run process-vlad-data with --write-meta set?')

if args.n_jobs != 0:
    raise Exception('don\'t set this, it\'s for internal use only')
if args.action == 'merge-samples':
    args.no_merged = True
args.samples = utils.get_arg_list(args.samples)
if args.samples is None:
    if args.no_merged:
        args.samples = [ds for ds in args.metafo if args.metafo[ds]['timepoint'] != 'merged']
    elif args.only_merged:
        args.samples = [ds for ds in args.metafo if '-merged' in args.metafo[ds]['sample']]
    elif args.write_meta:  # uh... maybe
        args.samples = None
    else:
        args.samples = args.metafo.keys()
if args.samples is not None:  # will only be None if args.write_meta... but oh well
    for isample in range(len(args.samples)):
        if args.samples[isample] in args.metafo:  # actual full name was specified on the command line
            continue
        args.samples[isample] = heads.full_sample(args.metafo, args.samples[isample])
args.subjects = utils.get_arg_list(args.subjects)
if args.subjects is not None and not args.write_meta:
    args.samples = [s for s in args.samples if args.metafo[s]['subject'] in args.subjects]

args.loci = utils.get_arg_list(args.loci)
if args.loci is not None:
    if len(set(args.loci) - set(utils.loci)) > 0:
        raise Exception('unexpected --loci \'%s\'' % ' '.join(set(args.loci) - set(utils.loci)))
    args.samples = [s for s in args.samples if args.metafo[s]['locus'] in args.loci]

args.isotypes = utils.get_arg_list(args.isotypes)
if args.isotypes is not None:
    if len(set(args.isotypes) - set(utils.isotypes)) > 0:
        raise Exception('unexpected --isotypes \'%s\'' % ' '.join(set(args.isotypes) - set(utils.isotypes)))
    args.samples = [s for s in args.samples if args.metafo[s].get('isotype') in args.isotypes]

args.timepoints = utils.get_arg_list(args.timepoints)
if args.timepoints is not None:
    args.samples = [s for s in args.samples if args.metafo[s]['timepoint'] in args.timepoints]

if args.only_umid:
    args.samples = [s for s in args.samples if args.metafo[s].get('umi')==True]

if not args.write_meta and len(args.samples) == 0:
    raise Exception('args.samples is empty')

args.seed_uids = utils.get_arg_list(args.seed_uids)

if args.print_width == 0:
    args.print_width = 99999

baseoutdir = heads.get_datadir(args.study, 'processed', extra_str=args.version)
if args.write_yaml:
    if args.version == 'latest':
        raise Exception('--version must be explicitly specified in order to write yamls (look in %s)' % baseoutdir.replace('/latest', ''))
    args.yamlwriter = YamlWriter(args, baseoutdir)

args.synth_component_studies = utils.get_arg_list(args.synth_component_studies)

# if args.ex_out_version != '':
#     assert args.get_selection_metrics or args.plot_partitions  # would have to implement, but really you always want to run get-selection-metrics as a separate step (not when partitioning) anyway

if args.paired_loci:
    assert args.action in ['cache-parameters', 'partition', 'seed-partition', 'simulate']  # others still need to be checked/implemented
    assert not args.get_naive_probabilities
    assert not args.view_alternative_annotations

# ----------------------------------------------------------------------------------------
if args.action == 'process-vlad-data':
    if '-synth' in args.study and args.synth_component_studies is None:
        raise Exception('specify which studies to be merged into the synth study with --synth-component-studies')
    preprocess.process_vlad_data(args)  # , add_seed_seqs=True)  # turn on for old (pre-'qa013-mar-8') samples
    sys.exit(0)
elif args.action == 'merge-samples':
    if args.subjects is None:
        subjects = set([mfo['subject'] for mfo in args.metafo.values()])
    else:
        subjects = args.subjects
    for subject in subjects:
        if args.dry_run:
            print '  %s \n      --dry-run not implemented for timepoint merging, use --test instead' % subject  # not implemented
            continue
        preprocess.merge_samples(args, args.study, subject)  # , add_seed_seqs=True)  # turn on for old (pre-'qa013-mar-8') samples
    sys.exit(0)
elif args.action == 'convert-csv-meta-to-yaml':
    # heads.parse_laura_neut_csv(args.study, args.seedfos)  # adds info to <args.seedfos>
    heads.write_yaml_metafo(args.study, args.metafo, args.seedfos)
    sys.exit(0)

# ----------------------------------------------------------------------------------------
def get_output_suffix(baseprefix):  # figure out whether we want to use old-style .csv output, or new-style .yaml
    if args.force_csv_output or os.path.exists(baseprefix + '.csv'):
        assert not os.path.exists(baseprefix + '.yaml')  # whole point of this function is that we're not mixing .csv and .yaml output files in the same directories
        return '.csv'
    else:
        assert not os.path.exists(baseprefix + '.csv')  # whole point of this function is that we're not mixing .csv and .yaml output files in the same directories (ok, so this is also check in the if above...)
        return '.yaml'

# ----------------------------------------------------------------------------------------
def pad_lines(linestr, padwidth=12):
    lines = [padwidth * ' ' + l for l in linestr.split('\n')]
    return '\n'.join(lines).rstrip()

# ----------------------------------------------------------------------------------------
def check_seeds_in_files(mfo, sample):
    subfos = heads.subset_seed_info(args.seedfos, mfo['subject'], mfo['locus'], seed_uids=args.seed_uids, seed_origin=args.seed_origin)
    if len(subfos) == 0:
        print '  no seeds'
        return
    ref_seeds = set(subfos)
    infname = heads.get_infname(args.study, sample, mfo)  # NOTE not checked after changing from: studies[args.study]['infname-fcn'](heads.get_datadir(args.study, 'raw'), sample)
    fastafo = utils.read_fastx(infname, n_max_queries=len(subfos))  # aw, screw it, just require that they're at the top of the file, and don't worry about there being extras further down
    fasta_subfos = collections.OrderedDict([(l['name'], {'uid' : l['name'], 'seq' : l['seq']}) for l in fastafo])
    seeds_from_file = set(fasta_subfos)
    if len(ref_seeds - seeds_from_file) > 0:
        raise Exception('missing %d seeds from %s:\n    %s' % (len(ref_seeds - seeds_from_file), infname, ' '.join(ref_seeds - seeds_from_file)))
    if ref_seeds != seeds_from_file:
        print '  ref - file: %s' % ' '.join(ref_seeds - seeds_from_file)
        print '  file - ref: %s' % ' '.join(seeds_from_file - ref_seeds)
        raise Exception('wtf %s' % infname)
    for seed in subfos:
        if subfos[seed]['seq'] != fasta_subfos[seed]['seq']:
            fasta_seq, meta_seq = utils.color_mutants(subfos[seed]['seq'], fasta_subfos[seed]['seq'], return_ref=True, align=True)
            raise Exception('different sequences for %s in %s\n   meta: %s\n  fasta: %s' % (seed, infname, meta_seq, fasta_seq))
    print '    %s   %2d seeds in %s' % (utils.color('green', 'ok'), len(subfos), infname)

# ----------------------------------------------------------------------------------------
def read_logs(logfname):  # NOTE this duplicates code in utils.process_out_err()/utils.finish_process()
    print '    log/err tails:'
    if os.path.exists(logfname):
        print '        %s tail (%s)' % (utils.color('green', 'log'), check_output(['ls', '-l', logfname]).strip())
        print pad_lines(check_output(['tail', '-n5', logfname]))
    errfname = logfname.replace('.log', '.err')
    if os.path.exists(errfname) and os.stat(errfname).st_size > 0:
        print '        %s tail (%s)' % (utils.color('red', 'err'), errfname)
        print pad_lines(check_output(['tail', '-n5', errfname]))
    if not os.path.exists(logfname) and not os.path.exists(errfname):
        print ''

# ----------------------------------------------------------------------------------------
def other_method_outfname(parameter_dir, method, metafo):
    sys.path.insert(1, './python')
    import glutils
    return glutils.get_fname(parameter_dir + '/' + method, metafo['locus'], 'v')

# ----------------------------------------------------------------------------------------
def run_other_method(args, parameter_dir, mfo, infname, outfname):
    if args.other_method not in ['tigger-default', 'igdiscover']:  # really just to make it easier to search for this fcn
        assert False
    if utils.output_exists(args, outfname, are_we_reading_existing_output(args)):
        return
    if utils.getsuffix(infname) not in  ['.fa', '.fasta']:
        if utils.getsuffix(infname) != '.csv':
            raise Exception('need to write an .fa file for %s (only happens automatically for csv, but input file is %s)' % (args.other_method, infname))
        fasta_infname = infname.replace('.csv', '.fa')
        seq_column = None
        if 'extra-args' in mfo and '--seq-column' in mfo['extra-args']['all']:
            tmpvar = mfo['extra-args']['all']
            seq_column = tmpvar[tmpvar.index('--seq-column') + 1]
        if not os.path.exists(fasta_infname):
            print '\nmissing input fasta'
            utils.csv_to_fasta(infname, outfname=fasta_infname, overwrite=False, remove_duplicates=True, seq_column=seq_column, name_column=None, debug=True)  # if <name_column> is None, it uses a hash of the sequence as the name
        infname = fasta_infname

    cmd = './test/%s-run.py' % args.other_method.replace('-default', '')  # kinda ugly, only effects tigger... but whatever
    cmd += ' --infname ' + infname
    cmd += ' --outfname ' + outfname
    cmd += ' --n-procs ' + str(args.n_procs)
    # cmd += ' --glfo-dir <uh>'
    if args.other_method == 'igdiscover':
        cmd += ' --glfo-dir ' + args.partis_dir + '/data/germlines/' + mfo['species']  # the partis methods have this as the default internally, but we want/have to set it explicitly here
        if mfo['species'] != 'human':
            cmd += ' --species ' + mfo['species']
    else:  # for now we're saving all the igdiscover output/intermediate files, so we write them to an output dir
        cmd += ' --workdir ' + parameter_dir + '/' + args.other_method + '/work'
    if args.n_random_queries is not None:
        assert args.n_random_subsets is None  # would need to be implemented
        cmd += ' --n-random-queries %d' % args.n_random_queries
    cmd += ' --gls-gen'  # just so it doesn't think we're running the non-gls-gen simu plots
    cmd += ' --locus ' + mfo['locus']
    if args.overwrite:
        cmd += ' --overwrite'

    if args.dry_run:
        utils.simplerun(cmd, dryrun=True)
    else:  # should really add a dry run option to utils.run_cmds()
        utils.run_cmds([{'cmd_str' : cmd,
                         'workdir' : parameter_dir + '/' + args.other_method,
                         'outfname' : outfname}],
                       n_max_tries=1,
                       debug='write')

# ----------------------------------------------------------------------------------------
def get_naive_prob_config_fname(args, mfo):
    if args.naive_prob_config_fname is None:
        return '%s/%s/naive-prob-config.yaml' % (heads.metadir, args.study)
    else:
        return args.naive_prob_config_fname

# ----------------------------------------------------------------------------------------
def do_stuff_with_existing_output(cmdstr, logfname, args, mfo, padwidth=16):
    # ----------------------------------------------------------------------------------------
    def get_single_cmd(actstr, outfname):
        singlecmd = '%s %s' % (cmdlist[0], actstr)
        if args.paired_loci:
            singlecmd += ' --paired-loci --paired-outdir %s' % utils.get_val_from_arglist(cmdlist, '--paired-outdir')
        else:
            singlecmd += ' --outfname %s --locus %s' % (outfname, utils.get_val_from_arglist(cmdlist, '--locus'))
        if args.extra_args is not None:
            singlecmd += ' %s' % args.extra_args
        if utils.getsuffix(outfname) == '.csv':  # old-style output files need germline info from somewhere
            singlecmd += ' --parameter-dir ' + utils.get_val_from_arglist(cmdlist, '--parameter-dir')
        if actstr in ['plot-partitions', 'merge-paired-partitions', 'get-selection-metrics'] and not args.no_plots:
            pdir = utils.get_val_from_arglist(cmdlist, '--plotdir')
            if args.ex_out_version != '':
                if pdir == 'paired-outdir':
                    print '    %s --ex-out-version not supported for paired plotting' % utils.color('yellow', 'warning')
                else:
                    pdir += '-' + args.ex_out_version
            singlecmd += ' --plotdir %s' % pdir
        if actstr == 'get-selection-metrics':  # NOTE you can also add --get-selection-metrics to 'extra-args' in the sample's metafo
            singlecmd += ' --only-print-best-partition'  # causes it to also only get tree metrics on best
            chastr = 'chosen-abs'
            if args.ex_out_version != '':
                chastr += '-' + args.ex_out_version
            singlecmd += ' --chosen-ab-fname %s' % (('%s/%s-%s.csv'%(outfname, mfo['subject'], chastr)) if utils.getsuffix(outfname) == '' else utils.insert_before_suffix('-'+chastr, outfname))  # the suffix is empty if it's a dir (i.e. paired non-seed clustering)
            abcfname = heads.getmetafname(args.study, ab_choice=True, sample=mfo['sample'])
            if os.path.exists(abcfname):
                singlecmd += ' --ab-choice-cfg %s' % abcfname
            # singlecmd += ' --ete-path None'
        if actstr == 'update-meta-info':
            singlecmd += ' %s %s' % (iostr('in'), utils.get_val_from_arglist(cmdlist, iostr('in')))
        if actstr == 'view-output':
            singlecmd += ' --abbreviate --only-print-best-partition'  # --abbreviate only affects the partition ascii art
            if args.only_print_seed_cluster:
                singlecmd += ' --only-print-seed-cluster'
            if args.only_print_queries_to_include_clusters:
                singlecmd += ' --only-print-queries-to-include-clusters'
        if args.action == 'seed-partition':
            singlecmd += ' --seed-unique-id %s' % utils.get_val_from_arglist(cmdlist, '--seed-unique-id')
            if args.paired_loci:
                singlecmd += ' --seed-loci %s' % utils.get_val_from_arglist(cmdlist, '--seed-loci')
        if '--queries-to-include-fname' in cmdlist:
            write_queries_to_include_fname(cmdstr, args, mfo)
            singlecmd += ' --queries-to-include-fname ' + utils.get_val_from_arglist(cmdlist, '--queries-to-include-fname')
        if 'extra-args' in mfo:
            for eact in set(['all', actstr]) & set(mfo['extra-args']):
                singlecmd += ' ' + ' '.join(mfo['extra-args'][eact])
        return singlecmd
    # ----------------------------------------------------------------------------------------
    def run_print_single(actstr, outfname, write_log_files=True):  # 'single' because we used to need to run more than one step
        logfname = '%s-%s.log' % (utils.getprefix(outfname), actstr)
        if args.ex_out_version != '':
            logfname = utils.insert_before_suffix('-'+args.ex_out_version, logfname)
        if 'seed-' in args.action and args.paired_loci:  # ok this is super ugly, but i only want to check for the existence of the igh output file (since don't want to have to figure out which other files should be there), but then we don't want it to look like the log file is specific only to igh
            assert logfname.count('partition-igh') == 1
            logfname = logfname.replace('partition-igh', 'partition')
        if write_log_files:
            print '                                              writing log to %s' % logfname
        if args.view_ascii and args.view_alternative_annotations:
            if not os.path.exists(logfname):
                raise Exception('log file %s doesn\'t exist: have to run without --view-ascii first in order to write log file' % logfname)
            with open(logfname) as logfile:
                print ''.join(logfile.readlines())
            return

        singlecmd = get_single_cmd(actstr, outfname)  # also (potentially) writes queries_to_include file
        outstr, errstr = utils.simplerun(singlecmd, dryrun=args.dry_run, return_out_err=True, extra_str=' '*padwidth)  # if you use check_call(), it gets out of sync with the stuff you're printing directly from this script
        if args.dry_run:  # ok, yes, this is weird to have after the actual fcn call, but we do actually want both of them
            return
        if write_log_files:
            with open(logfname, 'w') as logfile:
                logfile.write(singlecmd + '\n')
                logfile.write(outstr)
            if errstr.strip() != '':
                errfname = logfname.replace('.log', '.err')
                with open(errfname, 'w') as errfile:
                    errfile.write(errstr)
            print '    wrote out%s to %s%s' % ('' if errstr.strip() == '' else '/err', logfname, '' if errstr.strip() == '' else ('   ' + errfname))
        else:
            print utils.pad_lines(outstr, padwidth=padwidth)
    # ----------------------------------------------------------------------------------------
    if args.get_naive_probabilities:  # NOTE has to match parameter out dir in main loop
        assert not args.paired_loci  # would need to be updated
        with tempfile.NamedTemporaryFile() as outfile:
            tmpcmd = '%s/bin/get-naive-probabilities.py %s/parameters/all-probs.csv --debug --config-fname %s --outfname %s' % (args.partis_dir, os.path.dirname(logfname), get_naive_prob_config_fname(args, mfo), outfile.name)
            utils.simplerun(tmpcmd, extra_str=' '*padwidth, dryrun=args.dry_run)
            naive_probs.append(yaml.load(outfile, Loader=yaml.Loader))  # add it to the global variable (temporary hack to get the info out of this function)
        return

    cmdlist = cmdstr.split()
    if args.paired_loci and ('seed-' in args.action or args.write_cluster_summary_tables):
        outfname = paircluster.paired_fn(os.path.dirname(logfname), 'igh', suffix='.yaml', actstr='partition')
    else:
        outfname = utils.get_val_from_arglist(cmdlist, iostr('out'))
    if not os.path.exists(outfname):
        print '  output file doesn\'t exist, skipping:  %s' % outfname
        return

    if args.write_cluster_summary_tables:
        heads.write_cluster_summary_tables(outfname, heads.subset_seed_info(args.seedfos, mfo['subject'], mfo['locus']), args)
        return

    if args.get_subst_probabilities:
        assert args.action == 'seed-partition'  # at least for now? not sure what it'd mean to do it otherwise
        tmp_seedfos = heads.subset_seed_info(args.seedfos, mfo['subject'], mfo['locus'], seed_uids=[utils.get_val_from_arglist(cmdlist, '--seed-unique-id')])
        heads.write_subst_prob_table(outfname, utils.get_val_from_arglist(cmdlist, '--parameter-dir'), mfo['locus'], tmp_seedfos, args, debug=True)
        return

    if args.view_ascii and not args.view_alternative_annotations:
        run_print_single('view-output', outfname)
    elif args.plot_partitions:
        run_print_single('plot-partitions', outfname)
    elif args.merge_paired_partitions:
        run_print_single('merge-paired-partitions', outfname)
    elif args.get_selection_metrics:
        run_print_single('get-selection-metrics', outfname)
    elif args.view_alternative_annotations:
        run_print_single('view-alternative-annotations', outfname)
    elif args.update_meta_info:
        run_print_single('update-meta-info', outfname)
    else:
        raise Exception('not sure what to do!')

# ----------------------------------------------------------------------------------------
def write_queries_to_include_fname(cmdstr, args, mfo):  # this is super hackey, but I can't figure out a better way
    subfos = heads.subset_seed_info(args.seedfos, mfo['subject'], mfo['locus'])  # dont apply --seed-uids or --seed-origin here (since we want *all* available seeds to be included here)
    if len(subfos) == 0:
        return
    outfos = []
    for sfo in subfos.values():
        tmpfo = {k : v for k, v in sfo.items() if k != 'uid'}
        tmpfo['name'] = sfo['uid']  # this sucks, but the datascripts seed file uses 'uid', whereas all the partis seqfo stuff uses 'name'
        outfos.append(tmpfo)
    if args.paired_loci:
        heads.add_seed_pairing_info(subfos, outfos)
    qtifname = heads.get_qti_fname(args.study, mfo)
    if os.path.exists(qtifname):  # make sure it's the same as the one we would write (if it isn't, we just rewrite it, presumably/hopefully it just means the file is old)
        existing_vals = utils.read_seqfos(qtifname)
        new_vals = outfos
        if existing_vals == new_vals:
            return

    utils.write_seqfos(qtifname, outfos)  # each job will rewrite this, which is messy, but it's an easy way to make sure it's up to date

# ----------------------------------------------------------------------------------------
def run(cmdstr, logfname, args, mfo):
    cmdstr += ' --n-procs ' + str(args.n_procs)

    if are_we_reading_existing_output(args):
        do_stuff_with_existing_output(cmdstr, logfname, args, mfo)
        return

    print '%s %s%s' % (utils.color('red', 'run'), cmdstr[:args.print_width], '' if len(cmdstr) < args.print_width else '[...]')
    if args.dry_run:
        return
    write_queries_to_include_fname(cmdstr, args, mfo)
    # if heads.has_input_metafo_in_seedfos(mfo, args.seedfos) and '--input-metafname' in cmdstr.split():  # write once when caching parameters, then after that you should be getting the affinities from the partis output files
    #     heads.write_neut_info(utils.get_val_from_arglist(cmdstr.split(), '--input-metafname'), args.study, args.seedfos, locus=mfo['locus'], debug=True)

    # this duplicates some code in utils.run_cmds(), but whatever, it's not very much, and it's enough different it's not worth merging
    if not os.path.exists(os.path.dirname(logfname)):
        os.makedirs(os.path.dirname(logfname))
    errfname = logfname.replace('.log', '.err')
    logfile, errfile = open(logfname, 'w'), open(errfname, 'w')
    logfile.write(cmdstr + '\n')
    Popen(cmdstr.split(), stdout=logfile, stderr=errfile)  # should this be using the fcns in utils? not sure
    if args.sleep is not None:
        time.sleep(args.sleep)

# ----------------------------------------------------------------------------------------
def get_seed_cluster(partition, seed_uid):
    seed_clusters = [cl for cl in partition if seed_uid in cl]
    assert len(seed_clusters) == 1
    return seed_clusters[0]

# ----------------------------------------------------------------------------------------
def run_single_seed(logfname, cmd, sample, seedstr, queries=None, extra_logstr=None):
    if args.logfnames:
        print logfname
        return

    if args.paired_loci:  # it's just one of many files, but it's about the last one to be written
        outfname = paircluster.paired_fn(os.path.dirname(logfname), 'igh', suffix='.yaml', actstr='partition')
    else:
        outfname = utils.getprefix(logfname) + get_output_suffix(utils.getprefix(logfname))

    if args.write_yaml:
        assert not args.paired_loci  # needs to be implemented
        args.yamlwriter.edit(args, outfname, sample, seedstr=seedstr, extra_logstr=extra_logstr)
        return

    if heads.output_exists(args, outfname, are_we_reading_existing_output(args)):
        return

    if args.check:
        read_logs(logfname)
        return

    if not args.paired_loci:
        cmd += ' --outfname %s' % outfname

    if args.random_seed_seqs is None:
        cmd += ' --seed-unique-id ' + seedstr
        if args.paired_loci:
            assert seedstr.count(':') == 1
            loci = [s.split('-')[-1] for s in seedstr.split(':')]
            assert all(l in utils.loci for l in loci)  # already checked in heads.collect_paired_seeds()
            cmd += ' --seed-loci %s' % ':'.join(loci)
    else:
        assert '--seed' not in cmd
        cmd += ' --random-seed-seq --seed ' + seedstr  # if --random-seed-seqs is set, <seedstr> will be str(i) for i in range(--random-seed-seqs)

    # cmd += ' --calculate-alternative-annotations'

    run(cmd, logfname, args, args.metafo[sample])

# ----------------------------------------------------------------------------------------
def check_number_of_jobs():
    if args.start_n_max_and_exit:
        args.n_jobs += 1
        if args.n_jobs > args.n_max_jobs:
            if not args.logfnames:
                print 'finished starting %d jobs' % args.n_max_jobs
            sys.exit()
    else:
        if args.n_max_jobs is None or args.dry_run or args.check or args.logfnames or args.view_ascii or args.write_yaml:  # not really an exhaustive list
            return
        utils.limit_procs('bin/partis', n_max_procs=args.n_max_jobs, sleep_time=30)

# ----------------------------------------------------------------------------------------
def run_seed_things(mfo, sample, cmd, seedstr, extra_logstr):  # seedstr is str(iseed) if random seed seqs, otherwise just the seed uid str
    def seed_logstr():
        rstr = seedstr if args.random_seed_seqs is None else 'iseed-' + seedstr
        if args.paired_loci:
            rstr = rstr.replace(':', '+')
        return rstr

    check_number_of_jobs()
    if args.paired_loci:
        outdir = baseoutdir + '/' + sample + heads.logstr_str(args, extra_logstr) + '/seeds/' + seed_logstr()  # NOTE heads.logstr_str() is where args.logstr gets used
    else:
        outdir = baseoutdir + '/seeds/' + seed_logstr() + '/' + sample + heads.logstr_str(args, extra_logstr)
    logfname = outdir + '/' + pact(args.action) + '.log'

    if not args.logfnames:
        print '   %s    ' % utils.color('green', seed_logstr(), width=13),

    cmd += ' --plotdir %s' % ('paired-outdir' if args.paired_loci else (os.path.dirname(logfname)+'/plots'))

    # run seed-partition or alternate naive seqs
    if args.action == 'seed-partition':
        run_single_seed(logfname, cmd, sample, seedstr, extra_logstr=extra_logstr)
        return

# ----------------------------------------------------------------------------------------
def ldummy(ig_or_tr='ig'):  # ick (it would be a lot better to check all loci, but at least now we're checking the last one to be run)
    return utils.sub_loci(ig_or_tr)[-1]

# ----------------------------------------------------------------------------------------
def run_sample(sample, extra_logstr=None, seed_increment=0):
    check_number_of_jobs()
    mfo = args.metafo[sample]
    infname = heads.get_infname(args.study, sample, mfo)

    if mfo.get('simu'):  # i think this should really go somewhere else (maybe during meta info reading?)
        if args.action == 'simulate' or args.no_simu:
            return

    if not args.logfnames:
        print '%s %s %s   %s' % (utils.color('blue', mfo['subject'], width=5, padside='right'),
                                 utils.color('purple', mfo['timepoint'], width=7, padside='right'),
                                 utils.color('yellow', mfo['locus']),
                                 ('%-' + str(shorthand_width) + 's') % mfo.get('shorthand', mfo['sample'])),

    cmd = args.partis_dir + '/bin/partis %s %s %s' % (pact(args.action), iostr('in'), infname)
    cmd += ' --print-git-commit'
    if not args.paired_loci:
        cmd += ' --locus ' + mfo['locus']
    if mfo['species'] != 'human':
        cmd += ' --species ' + mfo['species']
    if 'extra-args' in mfo:
        # if len(set(mfo['extra-args']) - set(['all'] + all_actions)) > 0:
        #     raise Exception('unexpected action/key in extra-args in meta info: %s\n   meta info: \'%s\'\n   should be like: \'%s\'' % (' '.join(set(mfo['extra-args']) - set(['all'] + all_actions)), mfo['extra-args'], 'extra-args: {all: [--foop, --woop, doop]}'))
        for eact in set(['all', args.action]) & set(mfo['extra-args']):
            cmd += ' ' + ' '.join(str(a) for a in mfo['extra-args'][eact])  # have to convert to str since the yaml reader auto-converts to e.g. int

    parameter_dir = heads.get_parameter_dir(args.study, sample, extra_str=args.version)
    if args.action == 'cache-parameters':  # assume that parameters were cached once with the full sample, and we always want to use those (i.e. the general use case is --logstr is only set for non-parameter-caching options)
        parameter_dir += heads.logstr_str(args, extra_logstr)
    else:
        cmd += ' --refuse-to-cache-parameters'

    assert args.n_max_queries is None or args.n_random_queries is None  # doesn't make sense to specify both
    if args.n_random_queries is not None:
        cmd += ' --n-random-queries ' + str(args.n_random_queries)
    elif args.n_max_queries is not None:
        cmd += ' --n-max-queries ' + str(args.n_max_queries)

    cmd += ' --random-seed %d' % (utils.non_none([args.random_seed, 0]) + seed_increment)

    if not args.no_slurm:
        cmd += ' --batch-system slurm'

    if args.extra_args is not None:
        cmd += ' %s' % args.extra_args

    if args.paired_loci:
        cmd += ' --paired-loci --paired-outdir ' + parameter_dir
    else:
        cmd += ' --parameter-dir ' + parameter_dir
        cmd += ' --sw-cachefname %s/sw-cache%s' % (parameter_dir, get_output_suffix(parameter_dir + '/sw-cache'))

    if mfo.get('simu'):
        cmd += ' --is-simu'

    # cmd += ' --workdir %s' % processargs.get_workdir(batch_system=None if args.no_slurm else 'slurm')
    if args.no_indels:
        cmd += ' --no-indels'

    subfos = heads.subset_seed_info(args.seedfos, mfo['subject'], mfo['locus'])
    if len(subfos) > 0:  # dont apply --seed-uids or --seed-origin here (since we want *all* available seeds to be included here)
        cmd += ' --queries-to-include-fname %s' % heads.get_qti_fname(args.study, mfo)  # don't write this file until we get to run()

    if args.action == 'check-seeds-in-input-files':
        raise Exception('shouldn\'t have seeds anymore in input files')
        check_seeds_in_files(mfo, sample)
        return

    if args.random_seed_seqs is None:  # these are (only) the seed uids we're actually going to run
        subseedfos = heads.subset_seed_info(args.seedfos, mfo['subject'], mfo['locus'], seed_uids=args.seed_uids, seed_origin=args.seed_origin)
        seed_uids = heads.collect_paired_seeds(subseedfos) if args.paired_loci else subseedfos.keys()
    else:
        assert not args.paired_loci
        seed_uids = [str(i) for i in range(args.random_seed_seqs)]
    if 'seed-' in args.action:  # have to loop through all the seeds, so for now there's some code duplication in a separate fcn
        print ''
        for seedstr in seed_uids:
            run_seed_things(mfo, sample, cmd, seedstr, extra_logstr)
        return

    if args.action == 'cache-parameters':
        assert args.logstr == ''  # it just doesn't make that much sense, so i think it's better to forbid it for now
        logfname = baseoutdir + '/' + sample + heads.logstr_str(args, extra_logstr) + '.log'
        if args.other_method is None:
            outpath = '%s%s/hmm/hmms' % (parameter_dir, ('/parameters/%s'%ldummy()) if args.paired_loci else '')
        else:
            outpath = other_method_outfname(parameter_dir, args.other_method, mfo)
        if not args.no_plots:
            cmd += ' --plotdir %s' % ('paired-outdir' if args.paired_loci else (parameter_dir+'/plots'))

        cmd += ' --debug-allele-finding'
        if 'find-new-alleles' in mfo:  # was only ever there so it could be set to false
            raise Exception('I think I don\'t need this any more, but think about how well the germline inference will behave on short reads before taking it out altogether')

        if not args.force_default_initial_germline_dir and 'initial-gl-dirs' in mfo and mfo['locus'] in mfo['initial-gl-dirs']:
            cmd += ' --initial-germline-dir ' + heads.metadir + '/' + args.study + '/' + mfo['initial-gl-dirs'][mfo['locus']]
            cmd += ' --default-initial-germline-dir ' + heads.metadir + '/' + args.study + '/' + mfo['initial-gl-dirs'][mfo['locus']]  # arg, complicatedness

        # NOTE commenting this because everything except the neut processing should happen/get written to the queries-to-include file, and the neut will only come up in old old samples that I should just rejigger to work in the new system
        # if heads.has_input_metafo_in_seedfos(mfo, args.seedfos):  # NOTE this is checking for old-style input metafo in <args.seedfos>, whereas now we only use input metafo from a separate file (see exception text)
        #     if '--input-metafname' in cmd:
        #         raise Exception('can\'t merge auto-generated input meta file (from info in seed yaml file) with another input meta file (you should really use the latter -- the input meta info from seed files was really just a hack to avoid doing ic50 calculations)')
        #     cmd += ' --input-metafname %s/input-meta.yaml' % parameter_dir  # write once when caching parameters, then after that you should be getting the affinities from the partis output files

    elif pact(args.action) == 'partition':
        logfname = '%s%s/%s%s/partition.log' % (baseoutdir, '' if args.paired_loci else '/partitions', sample, heads.logstr_str(args, extra_logstr))
        if args.paired_loci:
            outpath = paircluster.paired_fn(os.path.dirname(logfname), 'igh', suffix='.yaml', actstr='partition')  # just checks for one file, but it is one of the last to be written
        else:
            outpath = utils.getprefix(logfname) + get_output_suffix(utils.getprefix(logfname))

        if args.small_clusters_to_ignore is not None:
            cmd += ' --small-clusters-to-ignore ' + args.small_clusters_to_ignore  # NOTE this is just a string, but in bin/partis it gets converted to/with a list/range
        if not args.paired_loci:
            cmd += ' --outfname %s' % outpath
        if not args.no_plots:
            cmd += ' --plotdir %s' % ('paired-outdir' if args.paired_loci else (os.path.dirname(logfname)+'/plots'))
        if args.extra_args is not None and '--count-parameters' in args.extra_args and not args.paired_loci:  # for paired loci they automatically get written to the paired outdir
            cmd += ' --parameter-out-dir %s/parameters' % os.path.dirname(logfname)  # NOTE has to match parameter out dir in do_stuff_with_existing_output()
    elif pact(args.action) == 'simulate':
        logfname = '%s/%s/simulate%s/simulate.log' % (baseoutdir, sample, heads.logstr_str(args, extra_logstr))
        clist = cmd.split()
        if args.simulate_from_scratch:
            clist += ['--simulate-from-scratch']
        else:
            if args.paired_loci:
                clist += ['--parameter-dir', '%s/parameters' % utils.get_val_from_arglist(clist, iostr('out'))]
        if args.paired_loci:
            utils.replace_in_arglist(clist, iostr('out'), os.path.dirname(logfname))
            outpath = paircluster.paired_fn(os.path.dirname(logfname), 'igh', suffix='.yaml')  # just checks for one file, but it is one of the last to be written
        else:
            outpath = utils.replace_suffix(logfname, '.yaml')
        utils.remove_from_arglist(clist, iostr('in'), has_arg=True)
        utils.remove_from_arglist(clist, '--input-metafname', has_arg=True)
        utils.remove_from_arglist(clist, '--input-metafnames', has_arg=True)
        utils.remove_from_arglist(clist, '--refuse-to-cache-parameters')
        cmd = ' '.join(clist)
        if not args.paired_loci:
            cmd += ' --outfname %s' % outpath
        if not args.paired_loci:
            raise Exception('i think non-paired simulation is ok, but the cmd line manipulation above needs checking')
    else:
        raise Exception('unexpected action: \'%s\'' % args.action)

    if args.logfnames:
        print logfname
        return
    if args.write_yaml:
        args.yamlwriter.edit(args, outpath, sample, extra_logstr=extra_logstr)
        return
    if heads.output_exists(args, outpath, are_we_reading_existing_output(args)):
        return
    if args.check:
        read_logs(logfname)
        return

    if args.other_method is not None:
        run_other_method(args, parameter_dir, mfo, infname, outpath)  # NOTE doesn't use <cmd>
    else:
        run(cmd, logfname, args, mfo)

# ----------------------------------------------------------------------------------------
naive_probs = []  # ok global variables are terrible, but this avoids having to rewrite the entire script just now
shorthand_width = max([len(args.metafo[sam].get('shorthand', sam)) for sam in args.samples])
for sample in args.samples:
    if args.n_random_subsets is None:  # default, just run once for each sample
        run_sample(sample)
    else:  # run a bunch of random subsets, chosen with different random seeds
        if args.random_seed is None:  # this is kind of weird to have separately, but it makes it clear what happens
            args.random_seed = 0
        for isub in range(args.n_random_subsets):
            run_sample(sample, extra_logstr='isub-%d' % isub, seed_increment=isub)

# ----------------------------------------------------------------------------------------
if args.get_naive_probabilities:
    print '  calc\'d probs over %d subsets:' % len(naive_probs)
    for probfo in naive_probs:
        print '     %6d / %-6d = %.4f +/- %.4f'  % (probfo['counts'], probfo['total'], probfo['fraction'], probfo['frac_err'])
    err_over_subsets = numpy.std([pfo['fraction'] for pfo in naive_probs], ddof=1) / math.sqrt(len(naive_probs))
    print '  mean over subsets:   %.4f +/- %.4f' % (numpy.mean([pfo['fraction'] for pfo in naive_probs]), err_over_subsets)
