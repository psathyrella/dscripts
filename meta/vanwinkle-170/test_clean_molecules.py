"""Tests for the delete/collapse/uid-rewrite pass.

The fixtures are hand-built so a wrong answer is visible: the collapse test
uses CONSCOUNTs that cannot be confused with either input, and the collision
test uses one UMI on two sequences, which is the case collapse must *not*
merge.

Run with `pytest meta/vanwinkle-170/test_clean_molecules.py` from the
datascripts checkout.
"""

import collections
import csv
import hashlib
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from clean_molecules import (  # noqa: E402
    UID_HASH_LEN,
    clean_donor_tsv,
    load_reference_genes,
    read_exclusions,
    resolve_ambi_call,
    split_master_tsv,
)

COLUMNS = ['sample_id', 'sequence_id', 'sequence', 'v_call']


def digest(sequence):
    return hashlib.md5(sequence.encode()).hexdigest()


def sequence_id(umi, conscount, dupcount):
    return '%s|C_CALL=IGHM|CONSCOUNT=%d|DUPCOUNT=%d' % (umi, conscount, dupcount)


def write_tsv(path, donor, rows):
    """rows: list of (umi, sequence, conscount, dupcount)."""
    with open(path, 'w') as tfile:
        tfile.write('\t'.join(COLUMNS) + '\n')
        for umi, sequence, conscount, dupcount in rows:
            tfile.write('\t'.join([donor, sequence_id(umi, conscount, dupcount), sequence, 'IGHV1-2*02']) + '\n')


def read_tsv(path):
    with open(path) as tfile:
        return list(csv.DictReader(tfile, delimiter='\t'))


def uids(rows):
    return [r['sequence_id'].split('|')[0] for r in rows]


def counts(row, key):
    return int([f for f in row['sequence_id'].split('|') if f.startswith(key + '=')][0].split('=')[1])


def test_deletes_every_copy_of_a_shared_molecule(tmp_path):
    """Both rows go, not just one: a surviving copy is still foreign material."""
    path = str(tmp_path / 'd.tsv')
    write_tsv(path, 'D', [
        ('AAAACCCCGGGG', 'ACGTAAAA', 6, 1),   # shared, held twice
        ('AAAACCCCGGGG', 'ACGTAAAA', 7, 2),
        ('TTTTCCCCGGGG', 'ACGTTTTT', 5, 1),   # kept
    ])
    clean_donor_tsv(path, {('AAAACCCCGGGG', digest('ACGTAAAA'))})
    rows = read_tsv(path)
    assert [r['sequence'] for r in rows] == ['ACGTTTTT']


def test_collapse_sums_conscount_and_dupcount(tmp_path):
    """One molecule's read support is split across the two concatenated sets."""
    path = str(tmp_path / 'd.tsv')
    write_tsv(path, 'D', [
        ('AAAACCCCGGGG', 'ACGTAAAA', 6, 1),
        ('AAAACCCCGGGG', 'ACGTAAAA', 7, 2),
    ])
    clean_donor_tsv(path, set())
    rows = read_tsv(path)
    assert len(rows) == 1
    assert counts(rows[0], 'CONSCOUNT') == 13, 'kept one row instead of adding the counts'
    assert counts(rows[0], 'DUPCOUNT') == 3


def test_one_umi_on_two_sequences_stays_two_rows(tmp_path):
    """Two different molecules that drew the same UMI are not one molecule.

    Collapse must not merge them, and the digest is what keeps their uids apart
    -- which is the whole reason the uid is not the bare UMI.
    """
    path = str(tmp_path / 'd.tsv')
    write_tsv(path, 'D', [
        ('AAAACCCCGGGG', 'ACGTAAAA', 6, 1),
        ('AAAACCCCGGGG', 'ACGTCCCC', 7, 1),
    ])
    clean_donor_tsv(path, set())
    rows = read_tsv(path)
    assert len(rows) == 2
    assert len(set(uids(rows))) == 2, 'the two molecules got the same uid, so partis would rename one'
    assert all(u.startswith('AAAACCCCGGGG.') for u in uids(rows))
    assert counts(rows[0], 'CONSCOUNT') == 6 and counts(rows[1], 'CONSCOUNT') == 7


def test_uid_is_the_molecule_key_across_donors(tmp_path):
    """The same molecule in two donors gets the same uid.

    That is what lets the sharing gate assert `no uid under two donor_ids`
    instead of recovering a UMI from a partis uid.
    """
    shared = ('AAAACCCCGGGG', 'ACGTAAAA')
    paths = []
    for donor in ('D1', 'D2'):
        path = str(tmp_path / ('%s.tsv' % donor))
        write_tsv(path, donor, [shared + (6, 1), ('TTTTCCCCGGGG', 'ACGT' + donor, 5, 1)])
        clean_donor_tsv(path, set())
        paths.append(path)
    first, second = [uids(read_tsv(p)) for p in paths]
    assert first[0] == second[0] == 'AAAACCCCGGGG.' + digest('ACGTAAAA')[:UID_HASH_LEN]
    assert first[1] != second[1], 'different molecules must not share a uid'


def test_other_columns_survive_unchanged(tmp_path):
    path = str(tmp_path / 'd.tsv')
    write_tsv(path, 'D', [('AAAACCCCGGGG', 'ACGTAAAA', 6, 1)])
    clean_donor_tsv(path, set())
    row = read_tsv(path)[0]
    assert (row['sample_id'], row['sequence'], row['v_call']) == ('D', 'ACGTAAAA', 'IGHV1-2*02')


def write_exclusion_files(tmp_path, table_rows, shared_rows):
    table = str(tmp_path / 'sharing.csv')
    shared = str(tmp_path / 'shared.csv')
    with open(table, 'w') as tfile:
        writer = csv.DictWriter(tfile, fieldnames=['locus', 'donor', 'decision'], lineterminator='\n')
        writer.writeheader()
        writer.writerows(table_rows)
    with open(shared, 'w') as sfile:
        writer = csv.DictWriter(sfile, fieldnames=['locus', 'donor', 'umi', 'sequence_md5', 'conscount'], lineterminator='\n')
        writer.writeheader()
        writer.writerows(shared_rows)
    return table, shared


def test_read_exclusions_gives_kept_donors_an_empty_set(tmp_path):
    """A kept donor with nothing to delete must be a key, not a missing lookup."""
    table, shared = write_exclusion_files(
        tmp_path,
        [{'locus': 'igh', 'donor': 'A', 'decision': 'keep'},
         {'locus': 'igh', 'donor': 'B', 'decision': 'drop'}],
        [])
    dropped, delete_keys = read_exclusions(table, shared)
    assert dropped == {('igh', 'B')}
    assert delete_keys == {('igh', 'A'): set()}


def test_read_exclusions_rejects_deletions_for_a_dropped_donor(tmp_path):
    """A dropped donor's whole file goes, so naming its rows means the two
    artifacts disagree about who was kept."""
    table, shared = write_exclusion_files(
        tmp_path,
        [{'locus': 'igh', 'donor': 'B', 'decision': 'drop'}],
        [{'locus': 'igh', 'donor': 'B', 'umi': 'AAAACCCCGGGG', 'sequence_md5': digest('ACGT'), 'conscount': 6}])
    with pytest.raises(Exception, match='not a kept donor'):
        read_exclusions(table, shared)


def test_read_exclusions_rejects_an_unknown_decision(tmp_path):
    table, shared = write_exclusion_files(
        tmp_path, [{'locus': 'igh', 'donor': 'A', 'decision': 'maybe'}], [])
    with pytest.raises(Exception, match='unexpected decision'):
        read_exclusions(table, shared)


MASTER_COLUMNS = ['sample_id', 'sequence_id', 'sequence', 'v_call']

# the genes the tests' reference germline set contains. IGLJ2-3 and IGKV1-13
# are deliberately absent: that is what makes them the interesting cases.
REFERENCE = {'IGHV1-2', 'IGHV1-69', 'IGHV3-23', 'IGHV3-30', 'IGHD4-4', 'IGHJ4',
             'IGKV1-33', 'IGKV1-37', 'IGKV1-39', 'IGKV2-28', 'IGKV2-40', 'IGKV1D-13',
             'IGLJ1', 'IGLJ2', 'IGLJ3', 'IGLJ3-1', 'IGLJ3-2', 'IGLJ3-4'}


def write_master(path, rows, header_every=None):
    """rows: list of (sample_id, umi, sequence, v_call).

    `header_every` repeats the header line every N data rows, which the real
    master files do.
    """
    headline = '\t'.join(MASTER_COLUMNS)
    lines = [headline]
    for i, (sample, umi, sequence, v_call) in enumerate(rows):
        if header_every is not None and i > 0 and i % header_every == 0:
            lines.append(headline)
        lines.append('\t'.join([sample, sequence_id(umi, 6, 1), sequence, v_call]))
    with open(path, 'w') as mfile:
        mfile.write('\n'.join(lines) + '\n')


def out_paths_for(tmp_path, sids):
    return collections.OrderedDict((s, str(tmp_path / ('%s.tsv' % s))) for s in sids)


def test_split_routes_rows_by_sample_id(tmp_path):
    master = str(tmp_path / 'master.tsv')
    write_master(master, [
        ('A', 'AAAACCCCGGGG', 'ACGTAAAA', 'IGHV1-2*02'),
        ('B', 'TTTTCCCCGGGG', 'ACGTTTTT', 'IGHV1-2*02'),
        ('A', 'GGGGCCCCGGGG', 'ACGTGGGG', 'IGHV1-2*02'),
    ])
    out = out_paths_for(tmp_path, ['A', 'B'])
    assert split_master_tsv(master, out, REFERENCE) == {'A': 2, 'B': 1}
    assert [r['sequence'] for r in read_tsv(out['A'])] == ['ACGTAAAA', 'ACGTGGGG']
    assert [r['sequence'] for r in read_tsv(out['B'])] == ['ACGTTTTT']


def test_split_does_not_give_a_subject_another_subjects_rows(tmp_path):
    """`RG113` must not collect `RG1130`'s rows.

    The grep this replaces was `grep '^<sid>'`, unanchored on the right, so it
    would have.  No id in this cohort is a prefix of another, so nothing was
    wrong -- but routing on the sample_id column cannot go wrong at all.
    """
    master = str(tmp_path / 'master.tsv')
    write_master(master, [
        ('RG113', 'AAAACCCCGGGG', 'ACGTAAAA', 'IGHV1-2*02'),
        ('RG1130', 'TTTTCCCCGGGG', 'ACGTTTTT', 'IGHV1-2*02'),
        ('RG1130', 'GGGGCCCCGGGG', 'ACGTGGGG', 'IGHV1-2*02'),
    ])
    out = out_paths_for(tmp_path, ['RG113', 'RG1130'])
    assert split_master_tsv(master, out, REFERENCE) == {'RG113': 1, 'RG1130': 2}
    assert [r['sequence'] for r in read_tsv(out['RG113'])] == ['ACGTAAAA']


def test_split_skips_the_repeated_header_lines(tmp_path):
    """The master files carry their header line many times over."""
    master = str(tmp_path / 'master.tsv')
    rows = [('A', 'AAAACCCCGG%02d' % i, 'ACGT%04d' % i, 'IGHV1-2*02') for i in range(9)]
    write_master(master, rows, header_every=3)
    out = out_paths_for(tmp_path, ['A'])
    assert split_master_tsv(master, out, REFERENCE) == {'A': 9}
    written = read_tsv(out['A'])
    assert len(written) == 9
    assert all(r['sample_id'] == 'A' for r in written), 'a header line was written as data'


def test_split_writes_only_the_wanted_subjects(tmp_path):
    master = str(tmp_path / 'master.tsv')
    write_master(master, [
        ('A', 'AAAACCCCGGGG', 'ACGTAAAA', 'IGHV1-2*02'),
        ('DROPPED', 'TTTTCCCCGGGG', 'ACGTTTTT', 'IGHV1-2*02'),
    ])
    out = out_paths_for(tmp_path, ['A'])
    assert split_master_tsv(master, out, REFERENCE) == {'A': 1}
    assert not os.path.exists(str(tmp_path / 'DROPPED.tsv'))


def test_split_resolves_ambi_gene_names(tmp_path):
    """Folded into the split, replacing a separate `sed -i` read-and-rewrite."""
    master = str(tmp_path / 'master.tsv')
    write_master(master, [
        ('A', 'AAAACCCCGGGG', 'ACGTAAAA', 'IGHV1-69ambi*01'),
        ('A', 'TTTTCCCCGGGG', 'ACGTTTTT', 'IGLJ2-3ambi*01'),
        ('A', 'GGGGCCCCGGGG', 'ACGTGGGG', 'IGHV1-69ambi*01,IGHV1-69*02'),
    ])
    out = out_paths_for(tmp_path, ['A'])
    split_master_tsv(master, out, REFERENCE)
    assert [r['v_call'] for r in read_tsv(out['A'])] == [
        'IGHV1-69*01', 'IGLJ2*01', 'IGHV1-69*01,IGHV1-69*02']


def test_split_leaves_rows_without_ambi_untouched(tmp_path):
    master = str(tmp_path / 'master.tsv')
    write_master(master, [('A', 'AAAACCCCGGGG', 'ACGTAAAA', 'IGHV1-69*01')])
    out = out_paths_for(tmp_path, ['A'])
    split_master_tsv(master, out, REFERENCE)
    assert read_tsv(out['A'])[0]['v_call'] == 'IGHV1-69*01'


def test_split_does_not_match_ambi_across_a_column_boundary(tmp_path):
    """The `sed -i 's/-[^-]*ambi//g'` this replaces ran over the whole line.

    Its `[^-]*` crosses tabs, so a hyphen in one column could pair with an
    `ambi` in another and delete everything between them, taking the columns'
    contents with it.
    """
    master = str(tmp_path / 'master.tsv')
    write_master(master, [('A', 'AAAACCCCGGGG', 'ACGTAAAA', 'IGLJ2-3ambi*01')])
    out = out_paths_for(tmp_path, ['A'])
    split_master_tsv(master, out, REFERENCE)
    written = read_tsv(out['A'])[0]
    assert written['sequence'] == 'ACGTAAAA', 'the sequence column was damaged'
    assert written['sequence_id'].startswith('AAAACCCCGGGG|')


AMBI_NAMES_IN_THE_COHORT = [
    # every `ambi` gene name in all three master files, with what it resolves to
    ('IGHV3-30ambi*01', 'IGHV3-30*01'),
    ('IGHV3-23ambi*01', 'IGHV3-23*01'),
    ('IGHV1-69ambi*01', 'IGHV1-69*01'),
    ('IGHD4-4ambi*01', 'IGHD4-4*01'),
    ('IGHV3-30ambi2*01', 'IGHV3-30*01'),      # the marker carries a suffix
    ('IGHV3-30ambi3*01', 'IGHV3-30*01'),
    ('IGHV3-30ambi-42*01', 'IGHV3-30*01'),    # ... which can contain a hyphen
    ('IGKV1-39ambi*01', 'IGKV1-39*01'),
    ('IGKV2-28ambi*01', 'IGKV2-28*01'),
    ('IGKV1-33ambi*01', 'IGKV1-33*01'),
    ('IGKV2-40ambi*01', 'IGKV2-40*01'),
    ('IGKV1-37ambi*01', 'IGKV1-37*01'),
    ('IGKV1-39ambi*01_c162a', 'IGKV1-39*01_c162a'),  # partis novel-allele suffix survives
    ('IGLJ2-3ambi*01', 'IGLJ2*01'),           # IGLJ2-3 is not a gene: IGLJ2 or IGLJ3
]


@pytest.mark.parametrize('call,expected', AMBI_NAMES_IN_THE_COHORT)
def test_resolve_every_ambi_name_in_the_cohort(call, expected):
    """The checked-in sed turned IGHV1-69ambi*01 into IGHV1*01, not a gene."""
    unresolved = collections.Counter()
    assert resolve_ambi_call(call, REFERENCE, unresolved) == expected
    assert unresolved == {}


def test_a_gene_the_reference_lacks_is_kept_and_counted():
    """IGKV1-13 is real; this reference carries IGKV1D-13 instead.

    The master calls it 850 times without the marker, so the marked and
    unmarked calls have to come out the same. Dropping the trailing -13 would
    give IGKV1, which is not a gene -- so the name is kept and reported.
    """
    unresolved = collections.Counter()
    assert resolve_ambi_call('IGKV1-13ambi*01', REFERENCE, unresolved) == 'IGKV1-13*01'
    assert unresolved == {'IGKV1-13': 1}


def test_load_reference_genes_reads_a_germline_set(tmp_path):
    gdir = tmp_path / 'glfo'
    (gdir / 'igh').mkdir(parents=True)
    (gdir / 'igh' / 'ighv.fasta').write_text('>IGHV1-69*01\nACGT\n>IGHV1-69*02\nACGA\n')
    (gdir / 'igh' / 'ighj.fasta').write_text('>IGHJ4*02\nACGT\n')
    assert load_reference_genes(str(gdir), 'igh') == {'IGHV1-69', 'IGHJ4'}


def test_load_reference_genes_raises_on_an_empty_set(tmp_path):
    (tmp_path / 'igh').mkdir()
    with pytest.raises(Exception, match='no germline genes'):
        load_reference_genes(str(tmp_path), 'igh')


def test_split_raises_when_a_wanted_subject_has_no_rows(tmp_path):
    """Silently writing a header-only file would shrink the cohort without saying so."""
    master = str(tmp_path / 'master.tsv')
    write_master(master, [('A', 'AAAACCCCGGGG', 'ACGTAAAA', 'IGHV1-2*02')])
    out = out_paths_for(tmp_path, ['A', 'MISSING'])
    with pytest.raises(Exception, match='no rows in'):
        split_master_tsv(master, out, REFERENCE)


def test_an_unmarked_call_beside_a_marked_one_is_untouched(tmp_path):
    """A call column can hold several calls with only some of them marked.

    Resolution must not touch the unmarked ones: dropping a trailing -N off a
    gene the reference happens to lack would silently rename it.
    """
    master = str(tmp_path / 'master.tsv')
    write_master(master, [('A', 'AAAACCCCGGGG', 'ACGTAAAA', 'IGHV1-69ambi*01,IGKV1-13*01')])
    out = out_paths_for(tmp_path, ['A'])
    split_master_tsv(master, out, REFERENCE)
    assert read_tsv(out['A'])[0]['v_call'] == 'IGHV1-69*01,IGKV1-13*01'
