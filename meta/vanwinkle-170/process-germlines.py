#!/usr/bin/env python3
""" Merge the newly-downloaded germline set with VDJbase alleles downloaded for the vanwinkle-170 cohort,
    remove duplicate sequences, and write the combined set to germlines/processed/. Also compares the two input sets.
    The idea is to get every gene that could possibly be in the changeo annotations for vanwinkle (somehow we still
    end up missing some).
"""
import partis.utils as utils
import partis.glutils as glutils
import colored_traceback.always

bdir = '/fh/fast/matsen_e/data/vanwinkle-170/germlines'
for locus in utils.sub_loci('ig'):
    print('%s' % utils.color('red', locus))
    processed_ref_glfo = glutils.read_glfo('data/germlines/human', locus)
    new_imgt_glfo = glutils.read_glfo('data/germlines/imgt-download/human', locus, template_glfo=processed_ref_glfo, remove_duplicates=True, skip_orfs=False, debug=True)
    print(new_imgt_glfo['seqs']['d'].keys())
    vdj_glfo = glutils.read_glfo(bdir+'/raw', locus, template_glfo=new_imgt_glfo, remove_duplicates=True, debug=True)
    print(utils.color('blue', 'merge'))
    merged_glfo, _ = glutils.get_merged_glfo(new_imgt_glfo, vdj_glfo, debug=True)
    print(utils.color('blue', 'write'))
    glutils.write_glfo(bdir+'/processed', merged_glfo, debug=True)
    print(utils.color('blue', 'read'))
    merged_glfo = glutils.read_glfo(bdir+'/processed', locus, debug=True)
    print(utils.color('blue', 'comparing %s' % locus))
    glutils.compare_glfos([new_imgt_glfo, vdj_glfo], ['default', 'vdjbase'], locus)
