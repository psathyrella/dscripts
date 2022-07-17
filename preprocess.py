import glob
import random
import copy
import string
import subprocess
import csv
import os
import sys
import collections
import json

partis_dir = os.getcwd()  # os.path.dirname(os.path.realpath(__file__)).replace('/bin', '')
sys.path.insert(1, partis_dir + '/python')
import heads
import utils

# ----------------------------------------------------------------------------------------
per_seq_metafo_keys = ['unique_id', 'timepoint', 'multiplicity']

letter_translations = {
    # NOTE can also be dmXXX for negative time points
    'laura-mb-2': {
        'A' : ['BG505', 'W14'],
        'B' : ['MG505', 'P31'],
        'C' : ['BF520', 'W1'],
        'D' : ['BF520', 'M9'],
        'E' : ['QA255', 'D462'],
    },
    'laura-mb': {
        'A' : ['BG505', 'w14'],
        'B' : ['MG505', 'p31'],
        'C' : ['BF520', 'w1'],
        'D' : ['BF520', 'm9'],
    },
    'kate-qrs': {
        '1': ['QB850', '240dpi'],
        '2': ['QA255', '462dpi'],
        '3': ['QA255', '1512dpi'],
        '4': ['QB850', '1586dpi'],
    },
}

loci = {
    'g' : 'igh',
    'm' : 'igh',
    'k' : 'igk',
    'l' : 'igl',
}

# ----------------------------------------------------------------------------------------
def get_per_seq_metafo_fname(seqfname):
    return seqfname.replace(utils.getsuffix(seqfname), '-per-seq-meta.yaml')

# ----------------------------------------------------------------------------------------
def get_vlad_paths(study):
    glob_str = heads.vlad_path_glob_strs[study]
    vpaths = sorted(glob.glob(heads.get_datadir(study, 'raw') + '/' + glob_str))
    vnames = [os.path.basename(vpath).split('.')[0] for vpath in vpaths]
    print '  found %d vlad paths for %s' % (len(vpaths), study)
    return zip(vnames, vpaths)

# ----------------------------------------------------------------------------------------
def get_multiplicity(study, seqfo):
    if study in ['crotty-fna', 'cap256-vrc26']:  # ick
        return 1
    uid = seqfo['name']
    if '_contig_' in uid:  # 10x data (no multiplicity info)
        multstr = '1'
    elif len(uid.split('-')) == 2:
        rank, multstr = uid.split('-')
    elif 'infostrs' in seqfo and 'size' in seqfo['infostrs']:
        multstr = seqfo['infostrs']['size']
    elif 'size=' in uid:
        tmpfo = collections.OrderedDict(tstr.split('=') for tstr in uid.split(';') if '=' in tstr)
        if 'size' in tmpfo:
            multstr = tmpfo['size']
        else:
            raise Exception('couldn\'t parse \'size\' from %s' % uid)
    else:
        raise Exception('couldn\'t get multiplicity from uid \'%s\'' % uid)
    return int(multstr)

# ----------------------------------------------------------------------------------------
def write_per_seq_metafo(per_seq_metafo, seqfname):
    if set(per_seq_metafo[0].keys()) != set(per_seq_metafo_keys):
        raise Exception('unexpected per seq meta info keys: %s' % per_seq_metafo[0].keys())
    mfname = get_per_seq_metafo_fname(seqfname)
    print '    writing per-seq metafo to %s' % mfname
    with open(mfname, 'w') as outfile:
        outfo = {mfo['unique_id'] : {k : v for k, v in mfo.items() if k!='unique_id'} for mfo in per_seq_metafo}  # this will break if i want more metafo keys, but oh well
        json.dump(outfo, outfile)

# ----------------------------------------------------------------------------------------
def get_standardized_timepoint(tpstr):  # NOTE see <letter_translations> above, which is I think where most of these come from (well, originally from vlad's paths, and in retrospect I should have changed those, sigh) NOTE dana-qrs isn't there, oh right, I think because it's only for the older samples where I had to translate the darn HS names
    tpstr = tpstr.lower()
    if 'dpi' in tpstr:
        tpstr = 'd' + tpstr.replace('dpi', '')
    return tpstr

# ----------------------------------------------------------------------------------------
def get_sample_merge_info(metafo, study, subject, args, dont_merge_timepoints=False):
    creg_key = None
    for tkey in ['isotype', 'locus']:  # merge by isotype if they all have it, otherwise locus (which i think is required)
        if all(tkey in mfo for mfo in metafo.values()):
            creg_key = tkey
            print '    using \'%s\' to merge samples (chose from isotype/locus)' % creg_key
            break
    merge_info = {}  # map from <keystr> (string uniquely identifying the group of samples we'll merge together) : those samples
    for sample, mfo in metafo.items():
        if '-merged' in mfo['sample']:
            continue
        if subject is not None and mfo['subject'] != subject:
            continue
        if sample not in args.samples:
            continue
        if study ==  'kate-qrs' or study == 'laura-mb':
            keystr = mfo['subject'] + '-' + mfo['locus'][2] + '-' + 'Ig' + mfo['isotype'].upper()
        elif study ==  'crotty-fna':
            keystr = mfo['subject']
        elif study ==  'scherer-hpv':
            keystr = mfo['subject'] + '-merged'
        else:
            if creg_key is None:
                raise Exception('couldn\'t choose a constant region key (isotype vs locus)')
            keystr = mfo['subject'] + '-' + mfo[creg_key] + '-merged'
        if dont_merge_timepoints:
            assert mfo['timepoint'] != ''
            tpstr = get_standardized_timepoint(mfo['timepoint'])
            if '-merged' in keystr:
                keystr = keystr.replace('-merged', '-%s-merged' % tpstr)
            else:
                keystr += '%s-merged' % tpstr  # not really sure that this is what I want, but it's probably ok
        if args.logstr != '':  # NOTE see heads.logstr_str() (but i think it's ok to use args.logstr by itself here)
            keystr = '%s-%s' % (keystr, args.logstr)

        if keystr not in merge_info:
            merge_info[keystr] = []
        merge_info[keystr].append(sample)

    if len(merge_info) == 0:
        raise Exception('no time point info')
    return merge_info

# ----------------------------------------------------------------------------------------
def merge_samples(args, study, subject, add_seed_seqs=False):
    # ----------------------------------------------------------------------------------------
    def remove_input_metafname(tmfo):
        if 'extra-args' not in tmfo:
            return
        for astr in tmfo['extra-args'].keys():
            if '--input-metafname' in tmfo['extra-args'][astr]:
                utils.remove_from_arglist(tmfo['extra-args'][astr], '--input-metafname', has_arg=True)
            if len(tmfo['extra-args'][astr]) == 0:
                del tmfo['extra-args'][astr]
        return tmfo
    # ----------------------------------------------------------------------------------------
    def update_merged_mfo(merged_sample, merged_mfo, sub_mfo, merged_outfname):
        if merged_mfo is None:
            merged_mfo = copy.deepcopy(mfo)
            remove_input_metafname(merged_mfo)  # need to remove these from the merged metafo, since it'd be redundant, and it'll have the wrong uids
            merged_mfo['sample'] = merged_sample
            merged_mfo['infname'] = merged_outfname
            if 'extra-args' not in merged_mfo:
                merged_mfo['extra-args'] = {}
            if 'all' not in merged_mfo['extra-args']:
                merged_mfo['extra-args']['all'] = []
            merged_mfo['extra-args']['all'] += ['--input-metafname', get_per_seq_metafo_fname(merged_outfname)]
            if args.dont_merge_timepoints:
                merged_mfo['timepoint'] = get_standardized_timepoint(merged_mfo['timepoint'])
            for k in ['shorthand', 'vlad-fname']:
                if k in merged_mfo:
                    del merged_mfo[k]
        else:
            sub_mfo = remove_input_metafname(copy.deepcopy(sub_mfo))  # ok this is hackey, but whatever
            for key in (set(merged_mfo) & set(sub_mfo)) - set(['sample', 'infname']):
                subval = sub_mfo[key] if key != 'timepoint' else get_standardized_timepoint(sub_mfo[key])
                if merged_mfo[key] != subval:
                    if key in heads.mandatory_metafo:
                        print '  %s mandatory meta value for \'%s\' differs among subsamples: %s %s' % (utils.color('red', 'warning'), key, merged_mfo[key], subval)
                    if key == 'extra-args':
                        for subkey in set(merged_mfo[key]) | set(sub_mfo[key]):
                            merged_mfo[key][subkey] = utils.merge_arg_lists(merged_mfo[key].get(subkey, []), sub_mfo[key].get(subkey, []))
                        continue
                    merged_mfo[key] = 'merged'
        return merged_mfo

    merge_info = get_sample_merge_info(args.metafo, study, subject, args, dont_merge_timepoints=args.dont_merge_timepoints)
    if args.test:
        tmpdir = '/tmp/merge-timepoint-test/subd'
        print 'writing test to %s' % tmpdir
        if os.path.exists(tmpdir):
            for fn in glob.glob(tmpdir + '/*.fasta'):
                os.remove(fn)
        else:
            os.makedirs(tmpdir)

    for merged_name, samples in merge_info.items():
        if len(samples) < 2 and args.dont_merge_single_samples:
            print '    too few samples (%d) for %s' % (len(samples), merged_name)
            continue
        outfname = None
        infnames, outfo = [], []
        if add_seed_seqs:
            added_seeds = set()
        per_seq_metafo = []  # if we don't collect this here, then we have to reconstruct the time point and multiplicity info from the merged unique_id, which would kind of suck
        merged_sample_mfo = None

        # read through the input files
        input_lengths = {}
        for sample in samples:
            mfo = args.metafo[sample]
            if args.loci is not None and mfo.get('locus') not in args.loci:
                continue
            if args.isotypes is not None and mfo.get('isotype') not in args.isotypes:
                continue
            fname = mfo['infname']
            if args.test:
                testfname = tmpdir + '/' + os.path.basename(fname)
                utils.write_fasta(testfname, utils.read_fastx(fname, n_max_queries=100))
                fname = testfname
            infnames.append(fname)
            if outfname is None:
                mergedir = os.path.abspath(os.path.dirname(fname) + '/' + os.path.pardir) + '/merged'
                outfname = mergedir + '/' + merged_name + os.path.splitext(fname)[1]
                print '  %s: merging %d sample%s' % (merged_name, len(samples), utils.plural(len(samples)))
                if os.path.exists(outfname):  # ok, this is a weird place to put this, but it works ok
                    if not os.path.exists(get_per_seq_metafo_fname(outfname)):
                        raise Exception('processed file %s exists but per-seq metafo file %s does not (probably should delete entire directory and start over)' % (outfname, get_per_seq_metafo_fname(outfname)))
                    print '  %s merged output file already exists, exiting: %s' % (utils.color('yellow', 'warning'), outfname)
                    return
            merged_sample_mfo = update_merged_mfo(merged_name, merged_sample_mfo, mfo, outfname)
            if add_seed_seqs:
                expected_seeds = heads.subset_seed_info(args.seedfos, mfo['subject'], mfo['locus'])  # it's kind of weird to have this inside the loop, but the seeds are the same for all the samples, and we want to get the subject and locus while we're in the loop, so...
                found_seeds = set()
            input_info = utils.read_fastx(fname, n_max_queries=args.n_max_queries, n_random_queries=args.n_random_queries)
            input_lengths[fname] = len(input_info)
            input_metafo = {}
            if 'extra-args' in mfo:
                mfnames = [utils.get_val_from_arglist(eargs, '--input-metafname') for eargs in mfo['extra-args'].values() if '--input-metafname' in eargs]
                if len(mfnames) > 0:  # i guess it's ok if there's more than one (would have to use seqfileopener fcn), although i'm not sure why there would be
                    input_metafo = utils.read_json_yaml(utils.get_single_entry(mfnames))  # maybe should use seqfile.read_input_metafo()? but i think this is ok
            this_psmfo = []
            def namefcn(u): return u if add_seed_seqs and u in expected_seeds else '%s-%s' % (mfo.get('shorthand', mfo['sample']), u)
            for seqfo in input_info:
                new_name = namefcn(seqfo['name'])
                if add_seed_seqs and seqfo['name'] in expected_seeds:
                    found_seeds.add(seqfo['name'])
                    if seqfo['name'] in added_seeds:  # make sure to only add the seed once
                        continue
                    added_seeds.add(seqfo['name'])  # we need both <found_seeds> and <added_seeds> because the former keeps track of just this sample (to make sure we find every seed in each sample), while the latter is over all samples to be merged (to make sure we only add each seed once)
                    multiplicity = 1
                else:
                    outfo.append({'name' : new_name, 'seq' : seqfo['seq']})
                    multiplicity = get_multiplicity(study, seqfo)
                mdict = {'unique_id' : new_name, 'timepoint' : mfo['timepoint'], 'multiplicity' : multiplicity}
                if 'paired' in mfo:
                    mdict['paired'] = True
                if seqfo['name'] in input_metafo:
                    imfo = input_metafo[seqfo['name']]
                    for tk in set(imfo) & set(mdict):
                        if imfo[tk] != mdict[tk]:
                            raise Exception('existing input metafo from %s doesn\'t agree with new values for %s: \'%s\' vs \'%s\'' % (', '.join(mfnames), tk, imfo[tk], mdict[tk]))
                    # gdmit, can't do this since this gives the prefix as *this* sample, whereas it needs to be the paired id's sample, and we have no real way of getting that here
                    # if 'paired-uids' in imfo:
                    #     imfo['paired-uids'] = [namefcn(u) for u in imfo['paired-uids']]
                    mdict.update(imfo)
                this_psmfo.append(mdict)
            if add_seed_seqs and len(set(expected_seeds) - found_seeds) > 0:
                raise Exception('missing seeds %s in %s' % (' '.join(set(expected_seeds) - found_seeds), infnames[-1]))
            per_seq_metafo += this_psmfo
        if add_seed_seqs and len(set(expected_seeds) - added_seeds)> 0:
            raise Exception('missing seeds %s from files %s' % (' '.join(set(expected_seeds) - set(seedfo.keys())), ' '.join(infnames)))
        if outfname is None:
            continue

        print '  merged extra-args: %s' % merged_sample_mfo['extra-args']
        args.metafo[merged_name] = merged_sample_mfo

        utils.mkdir(outfname, isfile=True)
        with open(outfname, 'w') as outfile:
            if add_seed_seqs:
                for seqfo in expected_seeds.values():  # first write the seeds
                    outfile.write('>%s\n%s\n' % (seqfo['uid'], seqfo['seq']))
            random.shuffle(outfo)  # has to be in random order so we can use --n-max-queries UPDATE actually I think I'm always using --n-random-queries now, but it's probably still better to have the different time points shuffled together
            for seqfo in outfo:  # then write everybody else
                outfile.write('>%s\n%s\n' % (seqfo['name'], seqfo['seq']))
        write_per_seq_metafo(per_seq_metafo, outfname)
        print '    input file seqs%s%s:' % ('' if args.n_max_queries is None else ', --n-max-queries %d'%args.n_max_queries, '' if args.n_random_queries is None else ', --n-random-queries %d'%args.n_random_queries)
        for ifn, count in sorted(input_lengths.items()):  #, key=opertator.itemgetter(1), reverse=True
            print '     %4d %s' % (count, ifn)
        print '     %4d total' % sum(input_lengths.values())

        print '    output:'
        out, _ = utils.simplerun('grep \'>\' %s | wc -l' % outfname, shell=True, return_out_err=True, debug=False)
        print '     %4d %s' % (int(out), outfname)

    heads.write_yaml_metafo(args.study, args.metafo, mfname=('%s/meta.yaml'%tmpdir) if args.test else None)

# ----------------------------------------------------------------------------------------
def translate_sample_name(study, vname, vpath, add_seed_seqs=False):
    """ translate Hs-<crap> name to more readable format, e.g. MG505-g-P31 """
    chunks = vname.split('-')
    mfo = {'species' : 'human'}
    # NOTE conflation of constant_region and isotype. IDGAF, it's convenient
    if study == 'qa255-synth':
        assert False  # shouldn't be used any more
        if '-LN2-' in vname:
            return translate_sample_name('kate-qrs', vname)
        elif '-LNE-' in vname:
            return translate_sample_name('laura-mb-2', vname)
        elif 'Hs-QA255-' in vname:
            return translate_sample_name('dana-qrs', vname)
        else:
            raise Exception('unexpected vname for qa255-synth: %s' % vname)
    elif study == 'bf520-synth':
        assert False  # shouldn't be used any more
        if '-LN-D-' in vname:
            return translate_sample_name('laura-mb', vname)
        elif '-LND-' in vname:
            return translate_sample_name('laura-mb-2', vname)
        else:
            raise Exception('unexpected vname for bf520-synth: %s' % vname)
    elif study == 'mg505-synth':
        assert False  # shouldn't be used any more
        if '-LN-B-' in vname:
            return translate_sample_name('laura-mb', vname)
        elif '-LNB-' in vname:
            return translate_sample_name('laura-mb-2', vname)
        else:
            raise Exception('unexpected vname for bf520-synth: %s' % vname)
    elif study == 'qa255-jul-18':
        assert chunks[0] == 'Hs'
        mfo['subject'] = chunks[1]
        mfo['timepoint'] = chunks[2]
        assert len(chunks[4]) == 3 and chunks[4][:2] == 'Ig'
        mfo['isotype'] = chunks[4][2].lower()
        mfo['locus'] = loci[mfo['isotype']]
        assert '5RACE' in chunks[3]
        mfo['umi'] = 'UMI' in chunks[3]
        mfo['sample'] = '-'.join([mfo['subject'], mfo['isotype'], mfo['timepoint'], 'j18'])
        if mfo['umi']:
            mfo['sample'] += '-umi'
    elif study == 'qa013-mar-8':  # Hs-QA013-d765-UMI5RACE-IgM
        assert chunks[0] == 'Hs'
        assert chunks[1] == 'QA013'
        mfo['subject'] = chunks[1]
        mfo['timepoint'] = chunks[2]
        assert chunks[3] == 'UMI5RACE'
        assert len(chunks[4]) == 3 and chunks[4][:2] == 'Ig'
        mfo['isotype'] = chunks[4][2].lower()
        mfo['locus'] = loci[mfo['isotype']]
        mfo['umi'] = True
        mfo['sample'] = '-'.join([mfo['subject'], mfo['isotype'], mfo['timepoint'], 'm8'])
    elif study == 'bf520-jul-18':  # Hs-BF250-UMI5RACE-IgM.igblast.prod.scrub.clon.fasta
        assert chunks[0] == 'Hs'
        mfo['subject'] = chunks[1]
        assert '5RACE' in chunks[2]
        mfo['umi'] = 'UMI' in chunks[2]
        mfo['timepoint'] = 'M9'
        assert len(chunks[3]) == 3 and chunks[3][:2] == 'Ig'
        mfo['isotype'] = chunks[3][2].lower()
        mfo['locus'] = loci[mfo['isotype']]
        mfo['sample'] = '-'.join([mfo['subject'], mfo['isotype'], mfo['timepoint'], 'j18'])
        if mfo['umi']:
            mfo['sample'] += '-umi'
    elif study == 'dana-qrs':
        assert chunks[0] == 'Hs'
        mfo['subject'] = chunks[1]
        mfo['timepoint'] = chunks[2]
        mfo['isotype'] = chunks[4][2].lower()
        mfo['locus'] = loci[mfo['isotype']]
        assert chunks[3] == '5RACE'
        assert chunks[4] == 'Ig' + mfo['isotype'].upper()
        mfo['sample'] = '-'.join([mfo['subject'], mfo['isotype'], mfo['timepoint']])
    elif study == 'laura-mb-2':
        assert chunks[0] == 'Hs'
        assert chunks[1][:2] == 'LN'
        letter = chunks[1][2]
        mfo['subject'], mfo['timepoint'] = letter_translations[study][letter]
        mfo['isotype'] = chunks[2].lower()
        mfo['locus'] = loci[mfo['isotype']]
        assert chunks[3] == '5RACE'
        assert chunks[4] == 'Ig' + mfo['isotype'].upper()
        mfo['sample'] = '-'.join([mfo['subject'], mfo['isotype'], mfo['timepoint']])
    elif study == 'laura-mb':  # wasn't actually using the nice new readable names at this point
        assert chunks[0] == 'Hs'
        assert chunks[1] == 'LN'
        letter = chunks[2]
        mfo['subject'], mfo['timepoint'] = letter_translations[study][letter]
        assert chunks[3] == '5RACE'
        assert chunks[4][:2] == 'Ig'
        mfo['isotype'] = chunks[4][2].lower()
        mfo['locus'] = loci[mfo['isotype']]
        mfo['sample'] = vname  # oh, right... i guess i can erase the above, except for getting the shorthand below'-'.join([subject, isotype])
        mfo['shorthand'] = letter.lower() + mfo['isotype']
    elif study == 'kate-qrs':  # wasn't actually using the nice new readable names at this point
        assert chunks[0] == 'Hs'
        assert chunks[1][:2] == 'LN'
        number = chunks[1][2]
        mfo['subject'], mfo['timepoint'] = letter_translations[study][number]
        assert chunks[2] == '5RACE'
        assert chunks[3][:2] == 'Ig'
        mfo['isotype'] = chunks[3][2].lower()
        mfo['locus'] = loci[mfo['isotype']]
        mfo['shorthand'] = number + mfo['isotype']
        mfo['sample'] = vname
    else:
        raise Exception('study %s not yet handled' % study)

    # NOTE for qa255-synth we don't make it down here in the first/outer call, but do in the nested call

    mfo['vlad-fname'] = vpath
    mfo['infname'] = '%s/vlad-processed-data%s/%s.fasta' % (heads.get_datadir(study, 'raw'), '-with-seed-seqs' if add_seed_seqs else '', mfo['sample'])
    mfo['extra-args'] = {'all' : heads.laura_base_args, 'seed-partition' : heads.laura_seed_args}

    return mfo

# ----------------------------------------------------------------------------------------
def process_vlad_data(args, add_seed_seqs=False):  # if <add_seed_seqs> is set, we do it the old way where we actually add the seed sequences to the fasta files. The new way is using the new partis option --queries-to-include-fname, so we don't have to put the seed sequences into the main input fastas.
    if args.test:
        raise Exception('that\'s only for timepoint merging')

    if args.write_meta:
        if args.synth_component_studies is not None:  # combining existing samples from existing studies
            metafos = collections.OrderedDict()
            for sub_study in args.synth_component_studies:
                sub_metafos = heads.read_metadata(sub_study, subjects=args.subjects)
                for sample in sub_metafos:
                    if sample in metafos:
                        raise Exception('sample %s already in metafos' % sample)
                    sub_metafos[sample]['infname'] = sub_metafos[sample]['infname'].replace(sub_study, args.study)
                    sub_metafos[sample]['sub-study'] = sub_study
                metafos.update(sub_metafos)  # meta info corresponding to that sub study, but modified to be a part of the synth study
        else:
            metafos = make_vlad_meta(args, args.study, add_seed_seqs=add_seed_seqs)
        heads.write_yaml_metafo(args.study, metafos)
    else:
        write_processed_vlad_files(args, args.study, add_seed_seqs=add_seed_seqs)

# ----------------------------------------------------------------------------------------
def make_vlad_meta(args, study, add_seed_seqs=False):  # NOTE do _not_ use <args.study> in here
    metafos = collections.OrderedDict()
    for vname, vpath in get_vlad_paths(study):  # uses glob.glob to get file names
        mfo = translate_sample_name(study, vname, vpath, add_seed_seqs=add_seed_seqs)
        if mfo['sample'] in metafos:
            raise Exception('sample %s already in metafos' % mfo['sample'])
        metafos[mfo['sample']] = mfo
    return metafos

# ----------------------------------------------------------------------------------------
# writes per-seq meta info, adds seeds sequences (if <add_seed_seqs>), removes chimeras
def write_processed_vlad_files(args, study, add_seed_seqs=False):  # NOTE do _not_ use <args.study> in here (well, at this point i think it's always the same, but oh well)
    # NOTE if making a -synth dataset, some files will be in vlad-processed-data-with-seed-seqs, some will be in vlad-processed-data, which I think is what we want, since we want the paths to be the same in the -synth and non-synth metafo
    metafos = heads.read_metadata(study)
    seedfos = heads.read_seed_info(study)

    first_time_through = True
    for sample, mfo in metafos.items():
        if first_time_through:
            print '  using vlad paths like: %s' % mfo['vlad-fname']
            first_time_through = False
        if args.subjects is not None and mfo['subject'] not in args.subjects:
            continue
        if mfo['timepoint'] == 'merged':
            continue
        outfname = mfo['infname']
        print '  %-20s  %20s       %s' % (mfo.get('sub-study', ''), mfo['sample'], outfname)
        if args.dry_run:
            continue
        if not os.path.exists(os.path.dirname(outfname)):
            os.makedirs(os.path.dirname(outfname))
        if os.path.exists(outfname):
            if not os.path.exists(get_per_seq_metafo_fname(outfname)):
                raise Exception('processed file exists but per-seq metafo file does not (probably should delete entire directory and start over)')
            print '  processed vlad files already exist %s, returning' % outfname
            return

        with open(outfname, 'w') as outfile:
            per_seq_metafo = []

            seedfo = heads.subset_seed_info(seedfos, mfo['subject'], mfo['locus'])
            print '     read %d seeds' % len(seedfo)
            sys.stdout.flush()
            for sfo in seedfo.values():
                if add_seed_seqs:
                    outfile.write('>%s\n%s\n' % (sfo['uid'], sfo['seq']))
                per_seq_metafo.append({'unique_id' : sfo['uid'], 'timepoint' : None, 'multiplicity' : 1})

            print '     writing fasta lines'
            sys.stdout.flush()
            remove_chimeras = '.igblast.' in mfo['vlad-fname']  # if we're using vlad's igblast'd output, we want to remove sequences that he's marked as chimeric NOTE chimera removal also happens in get_per_seq_metafo()  UPDATE: wtf? I don't see this fcn, I guess it doesn't exist any more
            n_chimeras = 0
            for seqfo in utils.read_fastx(mfo['vlad-fname']):
                if remove_chimeras and 'chimera' in seqfo['name']:
                    n_chimeras += 1
                    continue
                if ';' in seqfo['name']:  # stuff to the right of here just has to do with vlad's clonal family clustering
                    seqfo['name'] = seqfo['name'].split(';')[0]
                per_seq_metafo.append({'unique_id' : seqfo['name'], 'timepoint' : mfo['timepoint'], 'multiplicity' : get_vlad_multiplicity(seqfo)})
                outfile.write('>%s\n%s\n' % (seqfo['name'], seqfo['seq']))

        if remove_chimeras:
            print '    removed %d chimeras' % n_chimeras

        write_per_seq_metafo(per_seq_metafo, outfname)

    # don't screw up this indenting again if you add something here
