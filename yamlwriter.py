import sys
import copy
import os
import yaml
import subprocess

import heads
import preprocess

# ----------------------------------------------------------------------------------------
class YamlWriter(object):
    def __init__(self, args, yamldir):
        self.yamlpath = yamldir + '/info.yaml'
        if args.overwrite or not os.path.exists(self.yamlpath):
            print '  init yaml at %s' % self.yamlpath
            self.init_file(args)  # loops over all samples
        else:
            print '  read existing yaml at %s' % self.yamlpath
            self.read_file()

    # ----------------------------------------------------------------------------------------
    def init_file(self, args):
        self.yamlfo = {}

        # first get partis version info
        version_str = subprocess.check_output([args.partis_dir + '/bin/partis', 'version'])
        version_str_lines = version_str.strip().split('\n')
        if len(version_str_lines) != 2:
            raise Exception('expected two lines from ./bin/partis version but got %d:\n"%s"' % (len(version_str_lines), version_str))
        tmpstr, commit = version_str_lines[0].split()
        if tmpstr != 'commit:':
            raise Exception('unexpected commit line: %s' % version_str_lines[0])
        tmpstr, tag = version_str_lines[1].split()[:2]
        if tmpstr != 'tag:':
            raise Exception('unexpected tag line: %s' % version_str_lines[1])
        self.yamlfo['partis-version'] = {'tag' : tag, 'commit' : commit}
        self.yamlfo['id'] = args.study
        if args.extra_str is not None:
            self.yamlfo['id'] += '-' + args.extra_str

        self.yamlfo['subjects'] = {subject : {} for subject in set([mfo['subject'] for mfo in args.metafo.values()])}

        # then initialize the main sample dict, and add each sample's metafo
        self.yamlfo['samples'] = {
            sample : {
                'meta' : {k : v for k, v in mfo.items() if k not in ['infname-fcn', 'aa-translation-fname-fcn']},  # yaml crashes when it tries to read these, since it's got a lambda fcn in it
                'seeds' : {seedstr : {} for seedstr in heads.subset_seed_info(args.seedfos, mfo['subject'], mfo['locus'])},  # file looks nicer if we convert from OrderedDict
            } for sample, mfo in args.metafo.items()}

        print '  initializing %s' % self.yamlpath
        self.write_yaml()

    # ----------------------------------------------------------------------------------------
    def read_file(self):
        with open(self.yamlpath) as yamlfile:
            self.yamlfo = yaml.safe_load(yamlfile)

    # ----------------------------------------------------------------------------------------
    def write_yaml(self):
        with open(self.yamlpath, 'w') as yamlfile:
            yaml.dump(self.yamlfo, yamlfile, width=150)

    # ----------------------------------------------------------------------------------------
    def edit(self, args, outpath, sample, seedstr=None, extra_logstr=None):
        print ''

        if args.action == 'cache-parameters':
            assert outpath.endswith('/hmm/hmms')
            self.modify_yamlfo(sample, outpath, ['parameter-dir'], outpath[:-len('/hmm/hmms')])
            self.modify_yamlfo(sample, outpath, ['glfo-dir'], outpath[:-len('/hmms')] + '/germline-sets')
            sw_cache_path = outpath[:-len('/hmm/hmms')] + '/sw-cache.yaml'
            if os.path.exists(sw_cache_path.replace('.yaml', '.csv')):
                sw_cache_path = sw_cache_path.replace('.yaml', '.csv')
            self.modify_yamlfo(sample, outpath, ['sw-cache'], sw_cache_path)
            per_seq_metafname = preprocess.get_per_seq_metafo_fname(heads.get_infname(args.study, sample, args.metafo[sample]))
            if os.path.exists(per_seq_metafname):
                self.modify_yamlfo(sample, outpath, ['per-sequence-meta-file'], per_seq_metafname)
        elif 'partition' in args.action:
            if 'glfo-dir' not in self.yamlfo['samples'][sample]:
                raise Exception('need to write yaml info for \'cache-parameters\' before writing for \'partition\'')
            keylist = []
            val = outpath
            if args.action == 'seed-partition':
                keylist += ['seeds', seedstr]
            if args.logstr != '':
                keylist += ['other-partitions', heads.logstr_str(args, extra_logstr)]
            self.modify_yamlfo(sample, outpath, keylist + ['partition-file'], val)
            if '.csv' in outpath:
                self.modify_yamlfo(sample, outpath, keylist + ['cluster-annotation-file'], val.replace('.csv', '-cluster-annotations.csv'))
        else:
            assert False

    # ----------------------------------------------------------------------------------------
    def addval(self, samplefo, keylist, val):
        modified = False
        keylist = copy.deepcopy(keylist)
        subdict = samplefo
        while len(keylist) > 0:
            name = keylist[0]
            keylist = keylist[1:]
            if len(keylist) == 0:  # made it to the actual key -- so add the val
                if name in subdict:
                    if subdict[name] != val:
                        raise Exception('subdict[%s] existing val doesn\'t match new val:\n   %s\n   %s' % (name, subdict[name], val))
                    print '    already in yaml'
                else:
                    print '    adding %s' % val
                    subdict[name] = val
                    modified = True
            else:
                if name not in subdict:
                    subdict[name] = {}
                subdict = subdict[name]  # drill down one more layer

        return modified

    # ----------------------------------------------------------------------------------------
    def modify_yamlfo(self, sample, outpath, keylist, val):
        if os.path.exists(outpath):
            modified = self.addval(self.yamlfo['samples'][sample], keylist, val)
            if modified:
                print '      rewriting yaml'
                self.write_yaml()
        else:
            print '     missing output %s' % outpath
