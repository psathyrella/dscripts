#!/usr/bin/env python3
import os

""" Split germline seqs from different regions (V/D/J) from the joint/merged fasta downloaded from VDJbase (https://vdjbase.org/genetable)
    into separate per-region fasta files (e.g. ighv.fasta, ighd.fasta, ighj.fasta) for each locus.
"""

gdir = os.path.dirname(os.path.abspath(__file__))
for locus in ['igh', 'igk', 'igl']:
    inpath = os.path.join(gdir, locus, '%s-all.fasta' % locus)
    regions = ['v', 'd', 'j'] if locus == 'igh' else ['v', 'j']
    outfiles = {r: open(os.path.join(gdir, locus, '%s%s.fasta' % (locus, r)), 'w') for r in regions}
    with open(inpath) as infile:
        current_region = None
        for line in infile:
            if line.startswith('>'):
                current_region = line[4].lower()  # e.g. >IGHV... -> 'v'
                if current_region not in regions:
                    current_region = None  # skip non-V/D/J (e.g. constant region)
                    continue
            if current_region is not None:
                outfiles[current_region].write(line)
    for f in outfiles.values():
        f.close()
    print('%s: %s' % (locus, ', '.join('%s%s.fasta' % (locus, r) for r in regions)))
