#!/usr/bin/env python3
"""Compare v/d/j identity columns from AIRR/changeo TSVs to partis mute-freqs plotdirs.

Reads {v,d,j}_identity columns from TSV files, converts to mutation frequency
(1 - identity/100), builds histograms matching the partis mute-freqs/overall
format, writes them to an output plotdir, then runs compare-plotdirs.py to make
comparison plots.

Example:
  bin/cf-identity-vs-mfreqs.py --tsv /fh/fast/matsen_e/data/vanwinkle-170/raw/master_IGH_changeo.tsv \
      --mfreq-plotdirs /path/to/hmm/mute-freqs/overall:/path/to/multi-hmm/mute-freqs/overall \
      --names changeo:hmm:multi-hmm \
      --outdir _output/identity-vs-mfreqs
"""
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import argparse
import csv
import os
import sys
import colored_traceback.always

import partis.utils as utils
from partis.hist import Hist

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    description='Compare AIRR/changeo identity columns to partis mute-freq histograms.'
)
parser.add_argument('--tsv', required=True, help='AIRR/changeo TSV file (e.g. master_IGH_changeo.tsv)')
parser.add_argument('--mfreq-plotdirs', required=True, help='colon-separated partis mute-freqs/overall plotdir(s) with {v,d,j}_mean-freq.csv files')
parser.add_argument('--outdir', required=True, help='Output directory for comparison plots')
parser.add_argument('--names', default='changeo:partis', help='Colon-separated labels: first for the tsv, then one per --mfreq-plotdirs entry')
parser.add_argument('--n-bins', type=int, default=40, help='Number of bins for mutation frequency histograms')
parser.add_argument('--xmax', type=float, default=0.4, help='Max mutation frequency for histogram x-axis')
parser.add_argument('--normalize', action='store_true', default=True, help='Normalize histograms (default True)')
args = parser.parse_args()
args.mfreq_plotdirs = utils.get_arg_list(args.mfreq_plotdirs)
args.names = utils.get_arg_list(args.names)
if len(args.names) != len(args.mfreq_plotdirs) + 1:
    raise Exception('--names should have %d entries (1 for tsv + %d for --mfreq-plotdirs) but got %d' % (len(args.mfreq_plotdirs) + 1, len(args.mfreq_plotdirs), len(args.names)))

# ----------------------------------------------------------------------------------------
def read_identity_values(tsv_path):
    """Read v/d/j_identity from a single TSV, return {region: [mute_freq_values]}."""
    identity_cols = {'v': 'v_identity', 'd': 'd_identity', 'j': 'j_identity'}
    mfreqs = {r: [] for r in identity_cols}
    print('  reading %s' % tsv_path)
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row.get('v_identity') == 'v_identity':  # skip duplicate header lines
                continue
            for region, col in identity_cols.items():
                val = row.get(col, '')
                if val in ('', 'NA', 'None', None):
                    continue
                identity = float(val)
                mfreqs[region].append(1.0 - identity / 100.0)
    for region in mfreqs:
        print('    %s: %d values' % (region, len(mfreqs[region])))
    return mfreqs

# ----------------------------------------------------------------------------------------
def make_mfreq_hists(mfreq_values, n_bins, xmax):
    """Build {region: Hist} from mutation frequency value lists."""
    hists = {}
    for region, values in mfreq_values.items():
        if len(values) == 0:
            continue
        hist = Hist(n_bins, 0.0, xmax, value_list=values)
        if args.normalize:
            hist.normalize()
        hists[region] = hist
    return hists

# ----------------------------------------------------------------------------------------
def write_plotdir(hists, outdir):
    """Write hists to CSVs matching the partis mute-freqs/overall naming convention."""
    utils.prep_dir(outdir, wildlings=['*.csv'])
    for region, hist in hists.items():
        hist.write(os.path.join(outdir, '%s_mean-freq.csv' % region))
    print('  wrote identity-derived mfreq hists to %s' % outdir)

# ----------------------------------------------------------------------------------------
# read identity values and build histograms
mfreq_values = read_identity_values(args.tsv)
identity_hists = make_mfreq_hists(mfreq_values, args.n_bins, args.xmax)

# write to a temporary plotdir
identity_plotdir = os.path.join(args.outdir, 'identity-mfreqs')
write_plotdir(identity_hists, identity_plotdir)

# run compare-plotdirs.py
all_plotdirs = [identity_plotdir] + args.mfreq_plotdirs
cmd = 'bin/compare-plotdirs.py --plotdirs %s --names %s --outdir %s/comparison --normalize --extra-stats mean --translegend=-0.25:0' % (
    ':'.join(all_plotdirs),
    ':'.join(args.names),
    args.outdir,
)
print('  running: %s' % cmd)
utils.simplerun(cmd)
print('  comparison plots in %s/comparison' % args.outdir)
