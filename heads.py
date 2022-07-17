import glob
import copy
import os
import collections
import sys
import csv
import yaml
import operator
import random
import traceback
import numpy
import itertools

metadir = os.path.dirname(os.path.realpath(__file__)) + '/meta'
try:
    sys.path.insert(1, './python')
    import utils
    import glutils
    import indelutils
except:
    pass

isotypes = ['m', 'g', 'k', 'l']  # non-exhaustive
loci = {'igh' : 'vdj',  # damnit, this is already in utils, but I can't figure out a good way to import that
        'igk' : 'vj',
        'igl' : 'vj',
        'tra' : 'vj',
        'trb' : 'vdj',
        'trg' : 'vj',
        'trd' : 'vdj',
}

# NOTE if you modify either of these, you need to actually write code in read_csv_metafo() (EDIT hm, maybe not any more?)
mandatory_metafo = ['sample', 'infname', 'species', 'subject', 'locus']  # would be nice to add 'vlad-name' at some point, but it will be a bit involed
optional_metafo = ['timepoint', 'tissue', 'shorthand', 'isotype', 'umi', 'extra-args', 'vlad-fname']  # not using this at the moment I think

# ----------------------------------------------------------------------------------------
def logstr_str(args, extra_logstr):
    return_str = ''
    if args.logstr != '':
        return_str += '-' + args.logstr
    if extra_logstr is not None:
        return_str += '-' + extra_logstr
    return return_str

# ----------------------------------------------------------------------------------------
def read_yaml_metafo(metafile, metafo, study):
    cfg = yaml.load(metafile, Loader=yaml.Loader)
    for sample in cfg:
        samplefo = cfg[sample]
        samplefo['sample'] = sample
        if 'timepoint' not in samplefo:
            samplefo['timepoint'] = ''
        if 'shorthand' in samplefo and samplefo['shorthand'] == samplefo['sample']:  # clean up old files
            del samplefo['shorthand']
        for mname in mandatory_metafo:
            if mname not in samplefo:
                raise Exception('required meta info \'%s\' not present for %s' % (mname, study))
        metafo[sample] = samplefo

# ----------------------------------------------------------------------------------------
def read_csv_metafo(metafile, metafo, study):
    reader = csv.DictReader(metafile)
    for line in reader:
        if line[reader.fieldnames[0]][0] == '#':
            continue
        if line['species'] == '':
            raise Exception('species not specified for study %s' % study)
        if 'timepoint' not in line:
            line['timepoint'] = ''
        if 'shorthand' in line and line['shorthand'] == line['sample']:  # clean up old files
            del line['shorthand']
        assert 'locus' in line and line['locus'] != ''  # isotype is optional, locus is not
        for key in line:
            if line[key] in ['True', 'False']:
                line[key] = utils.useful_bool(line[key])
        if 'isotype' not in line:
            samplelist = line['sample'].split('-')
            if study == 'kate-qrs' or study == 'laura-mb':
                isostr = samplelist[-1]
                assert isostr[:2] == 'Ig'  # and isostr[2].lower() in ['m', 'g']:
                line['isotype'] = isostr[2].lower()
            elif study == 'laura-mb-2':
                assert len(samplelist) == 3
                line['isotype'] = samplelist[1]
            else:
                line['isotype'] = ''
            assert line['isotype'] == '' or line['isotype'] in 'mgkl'

        metafo[line['sample']] = line

# ----------------------------------------------------------------------------------------
def get_all_csv_seedfnames(study):
    return glob.glob('%s/%s/seeds/*-*.csv' % (metadir, study))

# ----------------------------------------------------------------------------------------
def getseedfname(study, yaml=False, subject=None, locus=None):
    if not yaml:
        fname = '%s/%s/seeds/%s-%s.csv' % (metadir, study, subject, locus)
    else:  # NOTE doesn't check if it exists, cause we want to do that in get_seeds()
        fname = '%s/%s/seeds.yaml' % (metadir, study)
    return fname

# ----------------------------------------------------------------------------------------
def getmetafname(study, force_yaml=False, ab_choice=False, sample=None):
    path = '%s/%s' % (metadir, study)
    fname = path + '/meta.csv'
    if force_yaml or ab_choice or not os.path.exists(fname):
        fname = '%s/%s.yaml' % (path, 'ab-choice' if ab_choice else 'samples')
    if ab_choice and os.path.isdir('%s/ab-choice'%path):
        if sample is not None and os.path.exists('%s/ab-choice/%s.yaml' % (path, sample)):
            fname = '%s/ab-choice/%s.yaml' % (path, sample)
        else:
            fname = '%s/ab-choice/default.yaml' % path
    return fname

# ----------------------------------------------------------------------------------------
def get_qti_fname(study, mfo):
    return '%s/%s/_tmp/%s-%s-queries-to-include.yaml' % (metadir, study, mfo['subject'], mfo['locus'])

# ----------------------------------------------------------------------------------------
def write_yaml_metafo(study, metafo, seedfos=None, mfname=None):  # e.g. if you have meta/seed csvs, and you want to switch to yamls for this study UPDATE i think this comment is out of date?
    if study in studies:  # convert the lambda functions to their return value
        for sample, mfo in metafo.items():
            for key in mfo.keys():
                if key.split('-')[-1] == 'fcn':
                    mfo['-'.join(key.split('-')[:-1])] = mfo[key](get_datadir(study, 'raw'), sample)
                    del mfo[key]
    if mfname is None:
        mfname = getmetafname(study, force_yaml=True)
    print '  writing yaml sample info to %s' % mfname
    if not os.path.exists(os.path.dirname(mfname)):
        os.makedirs(os.path.dirname(mfname))
    with open(mfname, 'w') as yfile:
        yaml.dump(metafo, yfile, default_flow_style=False)  # if you cast <meta> to dict(), it is easier to read the file by hand, but you lose sample order

    if seedfos is not None:
        yamlfo = seedfos
        for subject in yamlfo:
            for locus in yamlfo[subject]:
                yamlfo[subject][locus] = yamlfo[subject][locus].values()  # avoid writing an ordered dict, so the file's more human readable
        print '  writing yaml seed info to %s' % getseedfname(study, yaml=True)
        with open(getseedfname(study, yaml=True), 'w') as yfile:
            yaml.dump(yamlfo, yfile, default_flow_style=False)

# ----------------------------------------------------------------------------------------
def read_metadata(study, subjects=None):
    metafo = collections.OrderedDict()
    fname = getmetafname(study)
    with open(fname) as metafile:
        if '.yaml' in fname:
            read_yaml_metafo(metafile, metafo, study)
        else:
            read_csv_metafo(metafile, metafo, study)
    if study in studies:  # and add in the old-style info from <studies> (below)
        for key in studies[study]:
            for sample in metafo:
                metafo[sample][key] = studies[study][key]

    if subjects is not None:
        metafo = collections.OrderedDict([(sample, mfo) for sample, mfo in metafo.items() if mfo['subject'] in subjects])

    return metafo

# # ----------------------------------------------------------------------------------------
# def has_input_metafo_in_seedfos(mfo, seedfos):  # NOTE this is "input metafo" (i.e. affinity, multiplicity, etc, per-sequence stuff), not "metafo" as in the fcns above
#     if len(seedfos) == 0:
#         return False
#     imkeys = ['ic50_vals'] + utils.input_metafile_keys.keys()
#     return any(len(set(sfo) & set(imkeys)) > 0 for sfo in seedfos[mfo['subject']][mfo['locus']].values())  # if any of the seeds have input meta info

# ----------------------------------------------------------------------------------------
def subset_seed_info(seedfos, subject, locus, seed_uids=None, seed_origin=None):
    if subject not in seedfos:
        return {}
    if locus is None or locus == 'None':  # paired h/l
        subfos = collections.OrderedDict((u, s) for l in seedfos[subject] for u, s in seedfos[subject][l].items())
    else:
        if locus not in seedfos[subject]:
            return {}
        subfos = seedfos[subject][locus]

    if seed_uids is not None:
        uids_before = subfos.keys()
        subfos = collections.OrderedDict([(uid, sfo) for uid, sfo in subfos.items() if uid in seed_uids])
        if len(subfos) == 0:
            print '\n  %s no seedfos passing restriction %s (started with %s)' % (utils.color('yellow', 'warning'), seed_uids, uids_before)
    if seed_origin is not None:
        subfos = collections.OrderedDict([(uid, sfo) for uid, sfo in subfos.items() if sfo['origin'] == seed_origin])
    return subfos

# ----------------------------------------------------------------------------------------
def read_seed_info(study, debug=False):
    seedfos = {}

    seedfnames = get_all_csv_seedfnames(study)
    if len(seedfnames) == 0:  # look for yaml if no csvs exist
        seedfnames = [getseedfname(study, yaml=True)]
        if not os.path.exists(seedfnames[0]):
            return seedfos

    if debug:
        print '\n  reading seeds from %s' % ' '.join(seedfnames)
    for fname in seedfnames:   # shitty loop name, but I don't want to break backwards compatibility by changing the arg name
        if os.path.splitext(fname)[1] == '.csv':
            subject, locus = os.path.basename(fname).split('.')[0].split('-')
            if subject not in seedfos:
                seedfos[subject] = {}
            if locus not in seedfos[subject]:
                seedfos[subject][locus] = collections.OrderedDict()
            with open(fname) as seedfile:
                reader = csv.DictReader(seedfile)
                for line in reader:
                    suid = line['seed_uid'].strip()  # everywhere else it's just 'uid' now, but I don't want to change all the csv files
                    seedfos[subject][locus][suid] = {'uid' : suid, 'seq' : line['seq'].strip().upper()}
                    if 'origin' in line:
                        seedfos[subject][locus][suid]['origin'] = line['origin']
        elif os.path.splitext(fname)[1] == '.yaml':
            with open(fname) as seedfile:
                seedfos = yaml.load(seedfile, Loader=yaml.Loader)
            for subject in seedfos:
                for locus in seedfos[subject]:
                    seedfos[subject][locus] = collections.OrderedDict([(sfo['uid'], sfo) for sfo in seedfos[subject][locus]])
        else:
            assert False

    if debug:
        for subject in seedfos:
            print subject
            for locus in seedfos[subject]:
                print '  %s' % locus
                max_len = max([len(s) for s in seedfos[subject][locus]])
                for uid in seedfos[subject][locus]:
                    print ('   %' + str(max_len) + 's  %s') % (uid, seedfos[subject][locus][uid]['seq'])

    return seedfos

# ----------------------------------------------------------------------------------------
def collect_paired_seeds(seedfos):
    spairs = []

    n_before = len(seedfos)
    seedfos = [s for s in seedfos.values() if s['uid'].count('-')>0 and s['uid'].split('-')[-1] in utils.loci]
    if len(seedfos) < n_before:
        print '\n    %s removed %d/%d seedfos with uids that we couldn\'t parse for locus (should be e.g. name-stuff-igh)' % (utils.color('yellow', 'warning'), n_before - len(seedfos), n_before)
    def keyfunc(x):
        xl = x['uid'].split('-')
        if xl[-1] in utils.loci:  # if it has the locus at the end, remove it
            xl = xl[:-1]
        return '-'.join(xl)
    for ustr, sfos in itertools.groupby(sorted(seedfos, key=keyfunc), key=keyfunc):
        sfos = list(sfos)
        if len(sfos) > 2:
            print '    too many seeds for %s: %s' % (ustr, ' '.join(s['uid'] for s in sfos))
        elif len(sfos) <= 1:
            print '    not enough seeds for %s: %s' % (ustr, ' '.join(s['uid'] for s in sfos))
        if '-igh' in sfos[1]['uid']:
            sfos.reverse()
        spairs.append(':'.join(s['uid'] for s in sfos))

    return spairs

# ----------------------------------------------------------------------------------------
def add_seed_pairing_info(seedfos, outfos):
    spairs = collect_paired_seeds(seedfos)
    for sfo in outfos:
        spair = utils.get_single_entry([j for j in spairs if sfo['name'] in j.split(':')])
        pid = utils.get_single_entry([u for u in spair.split(':') if u != sfo['name']])
        sfo['paired-uids'] = [pid]
        sfo['locus'] = sfo['name'].split('-')[-1]
        if sfo['locus'] not in utils.loci:
            raise Exception('couldn\'t parse locus from %s (got %s)' % (sfo['name'], sfo['locus']))

# ----------------------------------------------------------------------------------------
laura_base_args = ['--skip-unproductive']
# NOTE these are *only* used for the samples below! newer samples use the meta info files
laura_seed_args = [
    '--min-largest-cluster-size', '50',  # NOTE used to be --n-final-clusters 1
    ' --write-additional-cluster-annotations', '0:99',
    '--n-partitions-to-write', '30',
    # '--max-cluster-size', 'XXX',  # would (probably) reduce memory usage
]
laura_extra_args = {'all' : laura_base_args, 'seed-partition' : laura_seed_args}

vlad_path_glob_strs = {
    'kate-qrs' : 'processed_data/*/04_igblast_out/*.igblast.prod.scrub.clon.fasta',  # changed from 03_fastx_out jul 16 (when making new qa255-synth using qa255-jul-18)
    'laura-mb' : 'processed_data/*/04_igblast_out/*.igblast.prod.scrub.clon.fasta',  # changed from 03_fastx_out jul 16 (when making new qa255-synth using qa255-jul-18)
    'laura-mb-2' : 'processed_data/*/04_igblast_out/*.igblast.prod.scrub.clon.fasta',  # changed from 03_fastx_out jul 16 (when making new qa255-synth using qa255-jul-18)
    'dana-qrs' : 'processed_data/*/04_igblast_out/*.igblast.prod.scrub.clon.fasta',
    'qa255-jul-18' : 'processed_data/*/04_igblast_out/*.igblast.prod.scrub.clon.fasta',
    'bf520-jul-18' : '*.igblast.prod.scrub.clon.fasta',
    'qa013-mar-8' : 'Hs*/04_igblast_out/*.igblast.prod.scrub.clon.fasta',
}

# NOTE 'extra-args' keys correspond to datascripts/run.py action, not partis/pact() action
studies = {  # DEPRECATED all this info now goes in the study's meta.yaml instead (well this is still used for the studies that appear here)
    'katie' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/' + sample + '/' + sample + '.tsv',
        'extra-args' : {'all' : ['--seq-column', 'nucleotide']},
        'find-new-alleles' : False,  # reads are really short, so it'd be pretty meaningless
        'initial-gl-dirs' : {'igh' : 'germlines/imgt-and-intronic'},  # subdir of datascripts/meta/
    },
    'bf520-synth' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/vlad-processed-data-with-seed-seqs/' + sample + '.fasta',
        'aa-translation-fname-fcn' : lambda datadir, sample: datadir + '/processed_data/' + sample + '/04_igblast_out/' + sample + '.igblast.prod.scrub.clon.fasta',  # NOTE this is now also the input file
        'extra-args' : laura_extra_args,
    },
    'mg505-synth' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/vlad-processed-data-with-seed-seqs/' + sample + '.fasta',
        'aa-translation-fname-fcn' : lambda datadir, sample: datadir + '/processed_data/' + sample + '/04_igblast_out/' + sample + '.igblast.prod.scrub.clon.fasta',  # NOTE this is now also the input file
        'extra-args' : laura_extra_args,
    },
    'rubelt-heritable-influence' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/sequence_data/bcell/' + sample + '.fasta',
    },
    'jason-influenza' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/processed/patients/' + sample + '.csv',
        'extra-args' : {'all' : ['--seq-column', 'SEQUENCE_INPUT', '--timepoint-column', 'TIME_POINT']},  # both this and jason-mg use some nucleotide string as SEQUENCE_INPUT that isn't always unique, so I'm just letting partis use the line number
    },
    'cui-et-al' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/processed/' + sample + '.tsv',
        'extra-args' : {'all' : ['--name-column', 'SEQUENCE_ID', '--seq-column', 'SEQUENCE_INPUT']},
    },
    'jason-mg' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/processed/patients/' + sample + '.csv',
        'extra-args' : {'all' : ['--seq-column', 'SEQUENCE_INPUT']},  # both this and jason-mg use some nucleotide string as SEQUENCE_INPUT that isn't always unique, so I'm just letting partis use the line number
    },
    'chaim-donor-45' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/' + sample + '_noFrameshifts.fa',
    },
    'chaim-vrc01-i-think' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/' + sample + '.fa',
    },
    'adaptive-billion-read' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/' + sample + '/shuffled.csv',
        'find-new-alleles' : False,  # not enough V to make it worthwhile
        'extra-args' : {'all' : ['--name-column', 'unique_id', '--seq-column', 'seq']},
    },
    'vollmers' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/' + sample + '/' + sample + '_Lineages.fasta',
        'find-new-alleles' : False,  # not enough V to make it worthwhile
    },
    'three-finger' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/' + sample + '.fa',
    },
    'sheng-gssp' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/' + sample + '.fa',
    },
    'davide-gl-valid' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/naive/' + sample + '_N_atleast-2_reheader.fasta',
    },
    'crotty-fna' : {
        'infname-fcn' : lambda datadir, sample: datadir + '/processed_fastas/' + sample + '.fasta',
        # 'initial-gl-dirs' : {l : 'imgt-plus-ramesh' for l in ['igh', 'igk', 'igl']},  # don't need this any more, since we added the ramesh genes to the default set (at about the same time, we also added the igblast/christopher genes)
    },
}

# ----------------------------------------------------------------------------------------
def get_infname(study, sample, mfo):
    if study in studies:
        return mfo['infname-fcn'](get_datadir(study, 'raw'), sample)  # studies[study]['infname-fcn'](get_datadir(study, 'raw'), sample)
    else:
        return mfo['infname']

# ----------------------------------------------------------------------------------------
def get_datadir(study, dtype, extra_str=None):
    rootdir = '/fh/fast/matsen_e'

    if dtype == 'raw':
        return_str = rootdir + '/data/' + study
    elif dtype == 'processed':
        return_str = rootdir + '/processed-data/partis/' + study
    else:
        raise Exception('unkown dtype %s' % dtype)

    if extra_str is not None:
        return_str += '/' + extra_str

    return return_str

# ----------------------------------------------------------------------------------------
def get_parameter_dir(study, sample, extra_str=None):
    return get_datadir(study, 'processed', extra_str=extra_str) + '/' + sample

# ----------------------------------------------------------------------------------------
def full_sample(metafos, shorthand):  # <shorthand> can also be full sample
    if shorthand in metafos:  # full sample, not actually a shorthand
        return shorthand
    for sample in metafos:
        if shorthand == metafos[sample].get('shorthand', sample):  # .get() isn't actually because we want the value of 'sample' (it would've triggered the previous if clause), but just so we can have 'shorthand' be absent from the dict if it's not set (to keep the meta yaml files shorter))
            return sample
    raise Exception('couldn\'t find a sample for shorthand \'%s\'' % shorthand)

# ----------------------------------------------------------------------------------------
def full_dataset(metafos, shorthand):  # backwards compatibility
    return full_sample(metafos, shorthand)

# ----------------------------------------------------------------------------------------
def output_exists(args, outpath, read_existing_output):  # [originally] copied from partis/python/compareutils.py NOTE shit, there's also one of these in python/utils.py
    # designed such that if this fcn returns True, the calling function should return immediately

    def delete_output(fn):
        if os.path.isdir(fn):
            raise Exception('warning: output %s is a directory, comment this by hand a few times to make sure it\'s working properly' % fn)
            check_call(['rm', '-r', fn])
        else:
            os.remove(fn)

    if not os.path.exists(outpath):
        return False

    if args.overwrite:
        print '                      overwriting %s' % outpath
        delete_output(outpath)
        return False
    if args.rm:
        print '                      removing %s' % outpath
        delete_output(outpath)
        return True  # NOTE args.overwrite returns False -- args.rm is for deleting without immediately re-running

    # handle zero empty dirs/zero length files (by default we leave 'em be)
    if os.path.isdir(outpath):
        if len(os.listdir(outpath)) == 0:
            if args.rm_zero_length:
                print '                      removing empty dir %s' % outpath
                os.rmdir(outpath)
                return False
            else:
                print '                      leaving empty dir %s (set --rm-zero-length to delete)' % outpath
                return True
    else:
        if os.stat(outpath).st_size == 0:
            if args.rm_zero_length:
                print '                      deleting zero length %s' % outpath
                os.remove(outpath)
                return False
            else:
                print '                      leaving zero length %s (set --rm-zero-length to delete)' % outpath
                return True

    if read_existing_output:
        print '                      output exists, proceeding to read: %s' % outpath
        return False  # because it's not the output now, it's the input, see...

    print '                      output exists, skipping (%s)' % outpath
    return True

# ----------------------------------------------------------------------------------------
def annotate_to_get_light_chain_locus(l_uid, l_seq, debug=False):  # holy shit I would really rather not do it this way (but laura's on vacation, and it should be ok)
    def get_locus_mfreq_with_vsearch(locus):
        vs_info = utils.run_vsearch('search', {l_uid : l_seq}, '/tmp/%d' % random.randint(0, 999999), threshold=0.3, glfo=glutils.read_glfo('data/germlines/human', locus), get_annotations=True, expect_failure=True)
        if l_uid not in vs_info['annotations']:
            return 999.
        return vs_info['annotations'][l_uid]['v_mut_freq']
    # def get_locus_mfreq_with_partis(locus):  # works fine, but the vsearch way is way faster, and I only wrote this because I forgot the vsearch fcn can handle glfos
    #     infname = '/tmp/in.fa'
    #     outfname = '/tmp/out.yaml'
    #     with open(infname, 'w') as infile:
    #         infile.write('>%s\n%s\n' % (l_uid, l_seq))
    #     cmd = './bin/partis cache-parameters --infname %s --outfname %s --only-sm --locus %s --leave-default-germline --debug 1' % (infname, outfname, locus)
    #     utils.simplerun(cmd, swallow_stdout=True, debug=False)
    #     _, annotation_list, _ = utils.read_output(outfname)
    #     assert len(annotation_list) == 1 and annotation_list[0]['unique_ids'][0] == l_uid
    #     mfreq = annotation_list[0]['mut_freqs'][0]
    #     os.remove(infname)
    #     os.remove(outfname)
    #     return mfreq

    if debug:
        print '    %s: getting mfreq for igk/igl' % l_uid
    mfreqs = {}
    for locus in ['igk', 'igl']:
        mfreqs[locus] = get_locus_mfreq_with_vsearch(locus)
        if debug:
            print '        %s  %.3f' % (locus, mfreqs[locus])
    sorted_locis, sorted_mfreqs = zip(*sorted(mfreqs.items(), key=operator.itemgetter(1)))
    if sorted_mfreqs[1] / sorted_mfreqs[0] < 2.5:
        raise Exception('not enough separation between mfreqs for igk/igl')
    if sorted_mfreqs[0] < 0. or sorted_mfreqs[0] > 0.4:
        raise Exception('weird minimum mfreq %f' % sorted_mfreqs[0])
    if debug:
        print '    return: %s' % sorted_locis[0]
    return sorted_locis[0]

# ----------------------------------------------------------------------------------------
def get_uids_and_loci(seedfos, subject, line, dbg_str):
    h_uid, l_uid, h_locus, l_locus = None, None, 'igh', None

    h_seq = line['heavy_nt'].upper().replace('-', '').replace('V', 'N')
    l_seq = line['light_nt'].upper().replace('-', '').replace('V', 'N')

    extra_info_strs = []
    if ', ' in line['ab_name']:  # if it's a list of the heavy, light uids, try to guess the locus info
        try:
            def remove_crap(tmpstr):
                if 'or' in tmpstr:
                    tmpstr, orstr = [s.strip() for s in tmpstr.split('or')]
                    # extra_info_strs.append(orstr)
                return tmpstr

            ab_name_str = line['ab_name']
            if '(' in ab_name_str and ')' in ab_name_str:
                leftstr, parenstr, rightstr = [t.strip() for s in ab_name_str.split('(') for t in s.split(')')]
                ab_name_str = leftstr + rightstr
                # extra_info_strs.append(parenstr)
            h_uid_str, l_uid_str = [remove_crap(ustr) for ustr in ab_name_str.split(', ')]

            isostr = h_uid_str.split('.')[1]
            assert isostr in ['G', 'G1', 'G2', 'G3']
            h_uid = h_uid_str

            isostr = l_uid_str.split('.')[1]
            assert isostr in ['K', 'L']
            l_locus = 'ig%s' % isostr[0].lower()
            l_uid = l_uid_str

        except:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
            print utils.pad_lines(''.join(lines))
            print '   %s failed to guess for \'%s\'' % (utils.color('red', 'error'), line['ab_name'])
            return None, None, None, None, None, None, dbg_str
    else:  # otherwise see if it's in the existing seed info
        for locus in seedfos[subject]:
            seed_uid = line['ab_name']
            if seed_uid not in seedfos[subject][locus]:  # different conventions for different seeds
                seed_uid = seed_uid + '-' + locus
            if seed_uid not in seedfos[subject][locus]:
                continue
            if locus in ['igh']:
                assert h_uid is None  # shouldn't be able to happen twice
                h_uid = seed_uid
            elif locus in ['igk', 'igl']:
                l_uid = seed_uid
                l_locus = locus
            else:
                assert False

    if h_uid is None:
        h_uid = line['ab_name'] + '-' + h_locus
    annotated_for_l_locus = False
    if l_uid is None:
        l_uid = line['ab_name'] + '-' + h_locus
        l_locus = annotate_to_get_light_chain_locus(l_uid, l_seq)
        annotated_for_l_locus = True

    dbg_str += '%-20s %-20s  %s %s        %25s   %s' % (utils.color('purple', 'x', width=20, padside='right') if h_uid is None else h_uid,
                                                        utils.color('purple', 'x', width=20, padside='right') if l_uid is None else l_uid,
                                                        h_locus,
                                                        utils.color('green' if annotated_for_l_locus else None, l_locus),
                                                        utils.color('blue', line['ab_name'], width=25), ', '.join(extra_info_strs))

    return h_uid, l_uid, h_locus, l_locus, h_seq, l_seq, dbg_str

# ----------------------------------------------------------------------------------------
def parse_laura_neut_csv(study, seedfos):  # adds neutralization info from laura's csv file to <seedfos> (actually I think the csv files are something I parsed out of an excel file she sent)
    def add_info(subject, locus, uid, seq, line, neutfo, dbg_str):
        dbg_str += '   %s:' % locus
        if uid in seedfos[subject][locus]:
            if seq != seedfos[subject][locus][uid]['seq']:
                a, b = utils.color_mutants(seq, seedfos[subject][locus][uid]['seq'], align=True, return_ref=True)
                print '            %s different sequences for %s\n            %s\n            %s' % (utils.color('red', 'warning'), uid, a, b)
            dbg_str += ' %s in seedfos' % utils.color('green', 'already')
        else:
            seedfos[subject][locus][uid] = {'uid' : uid, 'seq' : seq}
            seedfos[subject][locus][uid].update({k : line[k] for k in ['donor_age', 'donor_sex', 'timepoint_postinfection']})
            dbg_str += ' %s to seedfos' % utils.color('yellow', 'adding')

        existing_ic50_vals = set()
        if 'ic50_vals' in seedfos[subject][locus][uid]:
            existing_ic50_vals |= set(seedfos[subject][locus][uid]['ic50_vals'])
        else:
            seedfos[subject][locus][uid]['ic50_vals'] = {}
        new_vals = set(neutfo) - existing_ic50_vals
        if len(new_vals) > 0:
            seedfos[subject][locus][uid]['ic50_vals'].update({k : neutfo[k] for k in new_vals})
            dbg_str += ' (added %d ic50 vals)' % len(new_vals)

        return dbg_str

    tmp_subject_list = ['bf520', 'mg505', 'bg505']
    neut_csv_fnames = glob.glob(metadir + '/laura-neut-database*.csv')  # arg, not the best way to combine them
    print '  reading neut info from %d files: %s    (%s locus means we ran annotation to figure out the light chain locus)' % (len(neut_csv_fnames), ' '.join(neut_csv_fnames), utils.color('green', 'green'))
    for ncfn in neut_csv_fnames:
        print 'starting file %d: %s' % (neut_csv_fnames.index(ncfn), ncfn)
        with open(ncfn) as nfile:
            reader = csv.DictReader(nfile)
            last_subject = None
            for line in reader:
                if line['ab_name'] == '':
                    continue
                subject = line['donor']
                if subject in tmp_subject_list:
                    subject = subject.upper()
                if subject not in seedfos:
                    continue

                dbg_str = '  %-7s' % (subject if subject != last_subject else '')
                last_subject = subject

                h_uid, l_uid, h_locus, l_locus, h_seq, l_seq, dbg_str = get_uids_and_loci(seedfos, subject, line, dbg_str)

                # NOTE the neutralization info is of course a property of this particular pairing, but at the moment we only handle single sequences, so for sequences that occur with more than one pairing, it's kind of whackadoodle
                neutfo = {}
                for key in [k for k in line if k.split('_')[0] == 'ic50']:
                    virus_name = key[key.find('_') + 1 : ]
                    for stmp in tmp_subject_list:
                        if stmp in virus_name:
                            virus_name = virus_name.replace(stmp, stmp.upper())
                    nvalue = line[key]
                    if nvalue == 'na' or nvalue == 'ND' or nvalue == '':
                        continue
                    nvalue = float(nvalue)
                    # if nvalue >= 20:  # equivalent to negative results
                    #     # maybe in future
                    neutfo[virus_name] = nvalue

                if h_uid is not None and h_seq != '':
                    dbg_str = add_info(subject, h_locus, h_uid, h_seq, line, neutfo, dbg_str)
                if l_uid is not None and l_seq != '':
                    dbg_str = add_info(subject, l_locus, l_uid, l_seq, line, neutfo, dbg_str)

                print dbg_str

# # ----------------------------------------------------------------------------------------
# # transfers neut info in seed file to affinity info in input meta file (this is old: the better/current way is to just write an input meta file initially and use that, but I'm keeping this for backwards compatibility with e.g. pc64, wu-focused, and some of laura's data)
# def write_neut_info(outfname, study, seedfos, locus=None, restrict_to_autologous=False, val_choice_type='min', debug=False):
#     def get_autolog_ic50s(subject, vals):
#         return {n : v for n, v in vals.items() if subject.lower() in n.lower()}
#     def remove_zeros(vals):
#         non_zero_vals = {n : v for n, v in vals.items() if v > 0}
#         if len(non_zero_vals) < len(vals):
#             print '    removed %d zero-valued ic50s (no idea what these mean)' % (len(vals) - len(non_zero_vals))
#         return non_zero_vals

#     if debug:
#         print '  reading ic50 neut info from seed file'
#     outfo = {}
#     for subject in seedfos:
#         if debug:
#             print '   %-12s           ic50' % ''
#             print '   %-12s       min     max' % subject
#         for locus in ['igh', 'igk', 'igl'] if locus is None else [locus]:
#             if debug:
#                 print '    %s' % locus
#             for uid, sfo in seedfos[subject][locus].items():
#                 if 'ic50_vals' not in sfo or len(sfo['ic50_vals']) == 0:
#                     continue
#                 ic50s = sfo['ic50_vals']
#                 if restrict_to_autologous:
#                     ic50s = get_autolog_ic50s(subject, ic50s)
#                 ic50s = remove_zeros(ic50s)
#                 if debug:
#                     print '    %15s  %.3f  %.3f    %s' % (uid, min(ic50s.values()), max(ic50s.values()), ' '.join(n for n, v in ic50s.items() if v < 20.))
#                 if val_choice_type == 'min':
#                     combo_ic50 = min(ic50s.values())
#                 elif val_choice_type == 'mean':
#                     combo_ic50 = numpy.mean(ic50s.values())
#                 else:
#                     assert False
#                 outfo[uid] = {'affinity' : 1. / combo_ic50}
#     if not os.path.exists(os.path.dirname(outfname)):
#         os.makedirs(os.path.dirname(outfname))
#     if debug:
#         print '   writing to %s' % outfname
#     with open(outfname, 'w') as outfile:
#         yaml.dump(outfo, outfile)

# ----------------------------------------------------------------------------------------
def write_cluster_summary_tables(outfname, seedfos, args, n_max_clusters=500):  # for search: print_seed_summary_table print seed summary
    # ----------------------------------------------------------------------------------------
    def write_file(summary_fname, clusterfos):
        if len(clusterfos) == 0:
            print '  nothing to write to %s' % summary_fname
            return
        print '  writing to %s' % summary_fname
        with open(summary_fname, 'w') as sfile:
            writer = csv.DictWriter(sfile, clusterfos.values()[0][0].keys())  # NOTE this uses the last <outfo> from the loop above
            writer.writeheader()
            for clusterstr, outfos in sorted(clusterfos.items(), key=lambda x: x[0].count(':'), reverse=True):
                for iseed, outfo in enumerate(outfos):
                    if len(outfos) > 1 and 'clonal-seeds' in outfo:
                        outfo['clonal-seeds'] = [ofo['uid'] for ofo in outfos]
                    writer.writerow(outfo)
    # ----------------------------------------------------------------------------------------
    def add_to_info(clusterfos, cluster, seedid=None):
        keystr = ':'.join(cluster)
        line = annotations[keystr]
        family_size_rank = sorted_cluster_sizes.index(len(line['unique_ids'])) + 1  # this finds the *first* instance of this len, which (properly) gives ties the same rank
        outfo = [
            ('family-size', len(line['unique_ids'])),
            ('family-size-fraction', len(line['unique_ids']) / float(repertoire_size)),
            ('family-size-rank', family_size_rank),
            ('v_gene', line['v_gene']),
            ('d_gene', line['d_gene']),
            ('j_gene', line['j_gene']),
            ('cdr3_length', line['cdr3_length']),
        ]
        if seedid is None:
            outfo += [
                ('mean_percent_shm_nt', '%.1f' % (100 * numpy.mean(line['mut_freqs']))),
                ('median_percent_shm_nt', '%.1f' % (100 * numpy.median(line['mut_freqs']))),
            ]
        else:
            indelfo = utils.per_seq_val(line, 'indelfos', seedid)
            outfo = [('uid', seedid)] + outfo
            outfo += [
                ('seed_percent_shm_nt', '%.1f' % (100 * utils.per_seq_val(line, 'mut_freqs', seedid))),
                ('family_percent_shm_nt', '%.1f' % (100 * numpy.mean(line['mut_freqs']))),
                ('family_percent_shm_nt_rank', shm_sorted_clusters.index(line['unique_ids']) + 1),
                ('indels', indelfo if indelutils.has_indels(indelfo) else None),
                ('clonal-seeds', []),
            ]
        outfo = collections.OrderedDict(outfo)
        if keystr not in clusterfos:  # we make a separate list for the seeds in each cluster, so in the output file we can have clonal ones together
            clusterfos[keystr] = []
        clusterfos[keystr].append(outfo)
    # ----------------------------------------------------------------------------------------
    def write_seed_cluster_info():
        clusterfos = {}  # just a way to sort the seeds so that clonal ones appear next to each other
        for sfo in seedfos.values():
            clusters = [c for c in cpath.partitions[cpath.i_best] if sfo['uid'] in c]
            if len(clusters) == 0:
                if sfo['uid'] in outfname and not (args.paired_loci and '+'+sfo['uid'] in outfname):  # hackey way to see if we expect this seed to be in this file (if --paired is set, we only read the igh output file, so expect to have all the light chain uids missing)
                    print 'couldn\'t find %s in %s' % (sfo['uid'], outfname)
                continue
            elif len(clusters) > 1:
                raise Exception('found %d clusters (rather than 1) for \'%s\'' % (len(clusters), sfo['uid']))  # I think this shouldn't happen, even for seed partitioning, if we're looking at the best partition

            add_to_info(clusterfos, clusters[0], seedid=sfo['uid'])

        write_file(utils.replace_suffix(outfname, '-seed-summary.csv'), clusterfos)
    # ----------------------------------------------------------------------------------------
    def write_partition_info():
        clusterfos = {}
        for cluster in size_sorted_clusters:
            add_to_info(clusterfos, cluster)
            if len(clusterfos) >= n_max_clusters:
                print '    stopping with %d (of %d) clusters (skipped: %s)' % (len(clusterfos), len(size_sorted_clusters), utils.cluster_size_str(size_sorted_clusters))
                break

        write_file(utils.replace_suffix(outfname, '-partition-summary.csv'), clusterfos)
    # ----------------------------------------------------------------------------------------
    glfo, annotation_list, cpath = utils.read_output(outfname)
    if len(annotation_list) == 0 or len(cpath.partitions) == 0:
        print '    empty partitions or annotation list in %s' % outfname
        return
    annotations = utils.get_annotation_dict(annotation_list)
    size_sorted_clusters = sorted(cpath.partitions[cpath.i_best], key=len, reverse=True)
    sorted_cluster_sizes = [len(c) for c in size_sorted_clusters]
    repertoire_size = sum(len(c) for c in cpath.partitions[cpath.i_best])
    shm_sorted_clusters, _ = zip(*sorted([(c, numpy.mean(annotations[':'.join(c)]['mut_freqs'])) for c in cpath.partitions[cpath.i_best] if ':'.join(c) in annotations], key=operator.itemgetter(1), reverse=True))  # not really sure why they're sometimes missing, but I think it's just the occasional failed annotation

    write_seed_cluster_info()
    write_partition_info()

# ----------------------------------------------------------------------------------------
def get_subst_data(args, locus, region='v', debug=False):
    fname = '%s/data/substitution-profiles/GSSPs_for_V%s_genes_with_at_least_300_lineages.txt' % (args.partis_dir, locus[2].upper())
    all_aas = sorted(utils.get_all_amino_acids(no_stop=True))
    twd = 4
    with open(fname) as dfile:
        reader = csv.DictReader(dfile, delimiter='\t')
        subst_info = {}
        last_gene = None
        for line in reader:
            gene = line['%sgene' % region.upper()]
            assert '*' not in gene  # they just write it for all alleles together i guess
            gene += '*x'
            utils.split_gene(gene)  # just to check it's a valid gene

            if gene not in subst_info:
                # if len(subst_info) > 0:
                #     sys.exit()
                subst_info[gene] = []
                if debug:
                    print '    %s' % utils.color_gene(gene)
                    print '       ipos  germlines   freq     %s' % '  '.join(('%'+str(twd)+'s')%aa for aa in all_aas)
            else:
                assert last_gene == gene  # make sure a gene isn't in there twice

            sfos = subst_info[gene]

            newfo = {}
            newfo['ipos'] = int(line['pos']) - 1
            assert len(sfos) == 0 or newfo['ipos'] == sfos[-1]['ipos'] + 1

            newfo['germlines'] = line['germ'].split(',')  # there's multiple choices, presumably at positions where the alleles differ
            assert len(set(newfo['germlines']) & set(all_aas)) == len(newfo['germlines'])  # make sure everything in the germlines columns is actually an AA
            newfo['mut_freq'] = float(line['freq'])  # I think this is observed mutation frequency? I don't really have any idea

            # for tmp_aa in [aa for aa in utils.all_amino_acids if aa != '*']:
            newfo['aa_freqs'] = {}
            for tmp_aa in [k for k in line if k in all_aas]:
                newfo['aa_freqs'][tmp_aa] = float(line[tmp_aa])
            total = sum(newfo['aa_freqs'].values())
            if total > 0:
                if not utils.is_normed(total, this_eps=0.05):  # some of them are off by quite a bit (i saw 1.02), presumably from rounding
                    print '  %s total not normed: %.6f      %s' % (utils.color('red', 'warning'), total, newfo['aa_freqs'].values())

            sfos.append(newfo)

            if debug:
                def fstr(aa):
                     if aa in newfo['aa_freqs']:
                         return utils.color('blue', '0', width=twd, padside='right') if newfo['aa_freqs'][aa]==0. else (('%'+str(twd)+'.2f') % newfo['aa_freqs'][aa])
                     else:
                         return utils.color('red', '-', width=twd, padside='right')
                print '      %3d     %-6s     %.2f       %s' % (newfo['ipos'], ' '.join(newfo['germlines']), newfo['mut_freq'], '  '.join(fstr(a) for a in all_aas))

            last_gene = gene

    return subst_info

# ----------------------------------------------------------------------------------------
# if <args.inferred_intermediate_fname> is None we only write them for the seed seq
def write_subst_prob_table(outfname, parameter_dir, locus, seedfos, args, debug=False):
    glfo, annotation_list, cpath = utils.read_output(outfname, glfo_dir=parameter_dir + '/hmm/germline-sets', locus=locus)
    if len(annotation_list) == 0 or len(cpath.partitions) == 0:
        print '    empty partitions or annotation list in %s' % outfname
        return
    annotations = {':'.join(adict['unique_ids']) : adict for adict in annotation_list}  # collect the annotations in a dictionary so they're easier to access

    extra_seqfos = []
    if args.inferred_intermediate_fname is not None:
        extra_seqfos = utils.read_fastx(args.inferred_intermediate_fname)

    subst_info = get_subst_data(args, locus, debug=False)

    assert len(seedfos) == 1
    for sfo in seedfos.values():
        clusters = [c for c in cpath.partitions[cpath.i_best] if sfo['uid'] in c]
        if len(clusters) == 0:
            if sfo['uid'] in outfname:  # hackey way to see if we expect this seed to be in this file
                print 'couldn\'t find %s in %s' % (sfo['uid'], outfname)
            continue
        elif len(clusters) > 1:
            raise Exception('found %d clusters (rather than 1) for \'%s\'' % (len(clusters), sfo['uid']))  # I think this shouldn't happen, even for seed partitioning, if we're looking at the best partition
        keystr = ':'.join(clusters[0])
        line = annotations[keystr]

        # get annotations for the inferred intermediates
        if len(extra_seqfos) > 0:
            utils.add_seqs_to_line(line, extra_seqfos, glfo)

        utils.add_naive_seq_aa(line)
        utils.add_seqs_aa(line)
        fv_len = len(line['fv_insertion'])
        if fv_len % 3 != 0:
            raise Exception('at least for now we need the fv insertion to be of length a multiple of three so the translation works')
        if line['v_5p_del'] > 0:
            raise Exception('would need to handle this case below (just some simple padding but i don\'t want to do it now)')
        aa_start, aa_stop = fv_len / 3, len(line['naive_seq_aa']) - int(len(line['jf_insertion']) / 3.)  # TODO not sure if i need the jf part
        naive_seq_aa = line['naive_seq_aa'][aa_start : aa_stop]
        output_infos = []
        for tuid in [sfo['uid']] + [s['name'] for s in extra_seqfos]:
            seq_aa = utils.per_seq_val(line, 'seqs_aa', tuid)[aa_start : aa_stop]
            if debug:
                print '      %s' % utils.color('blue', tuid)
                utils.color_mutants(naive_seq_aa, seq_aa, print_result=True, amino_acid=True, extra_str='           ')

            sfogenes = [g for g in subst_info if utils.are_alleles(line['v_gene'], g)]
            if len(sfogenes) != 1:
                raise Exception('couldn\'t find exactly one gene to match %s from among %s' % (utils.color_gene(line['v_gene']), utils.color_genes(subst_info.keys())))
            if debug:
                print '         using info from %s for gene from annotation %s' % (utils.color_gene(sfogenes[0]), utils.color_gene(line['v_gene']))
            sinfos = subst_info[sfogenes[0]]
            outfo = {'name' : tuid, 'naive_seq_aa' : [], 'mature_seq_aa' : [],
                     'overall_freq' : [],  # fraction of observed sequences in which this position had a non-synonymous mutation (so might make more sense to normalize to 1, since as-is they're scaled to the shm rate of whatever repertoire they were calculated on)
                     'per_base_freq' : []  # given that it's mutated (i.e. the previous line), propensity for each position to mutate to the aa that we observed in the mature sequence
            }
            if debug:
                print '           ipos  ipos+1  naive  mature  overall  per-base'
            for sfo in sinfos:
                naive_aa = naive_seq_aa[sfo['ipos']]
                if naive_aa not in sfo['germlines']:
                    print '    %s germline from annotation %s at %d not among germlines in subst_info %s' % (utils.color('red', 'warning'), naive_aa, sfo['ipos'], sfo['germlines'])
                mature_aa = seq_aa[sfo['ipos']]
                outfo['naive_seq_aa'].append(naive_aa)
                outfo['mature_seq_aa'].append(mature_aa)
                outfo['overall_freq'].append(sfo['mut_freq'])
                if mature_aa in utils.ambiguous_amino_acids or mature_aa == naive_aa:  # position not mutated
                    outfo['per_base_freq'].append(None)
                else:
                    outfo['per_base_freq'].append(sfo['aa_freqs'][mature_aa])
                    if debug:
                        print '           %3d %3d     %s --> %s        %.2f     %.2f' % (sfo['ipos'], sfo['ipos']+1, naive_aa, mature_aa if mature_aa!=naive_aa else '', sfo['mut_freq'], sfo['aa_freqs'][mature_aa])
            output_infos.append(outfo)

        if len(output_infos) == 0:
            print '  nothing to write to %s' % output_infos
            continue
        sprob_fname = utils.replace_suffix(outfname, '-subst-probs.csv')
        print '  writing to %s' % sprob_fname
        with open(sprob_fname, 'w') as outfile:
            for ifo, outfo in enumerate(output_infos):
                if ifo == 0:
                    outfile.write('0-based index,,%s\n' % ','.join(str(i) for i in range(len(outfo['naive_seq_aa']))))
                    outfile.write('1-based index,,%s\n' % ','.join(str(i+1) for i in range(len(outfo['naive_seq_aa']))))
                    outfile.write('naive,,%s\n' % ','.join(outfo['naive_seq_aa']))
                def fstr(val): return '%.2f'%val if val is not None else ''
                outfile.write('%s,,%s\n' % (outfo['name'], ','.join(outfo['mature_seq_aa'])))
                outfile.write('%s,overall freq,%s\n' % (outfo['name'], ','.join(fstr(f) for f in outfo['overall_freq'])))
                outfile.write('%s,per-base freq,%s\n' % (outfo['name'], ','.join(fstr(f) for f in outfo['per_base_freq'])))
