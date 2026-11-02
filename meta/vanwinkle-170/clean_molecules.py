"""Split the master changeo TSV per donor and clean each one, for `process.py`.

Three transformations, in this order, all keyed on the molecule -- the
(UMI, sequence) pair:

1. **delete** every row of any molecule the donor shares with another donor
2. **collapse** rows that are the same molecule, summing CONSCOUNT and DUPCOUNT
3. **rewrite the uid** as `<umi>.<digest[:8]>`

Order matters: a molecule's two rows carry the same digest, so rewriting the
uid before collapsing would leave them colliding.

Why step 3 exists.  partis truncates a fasta header at the first `|`
(`partis/utils.py`, `read_fastx`), so the uid is the bare UMI, and on a
duplicate it renames the second to `<uid>-2` (`choose_non_dup_id`).  In a
properly collapsed set a UMI appears once, so this should never fire -- but
about a dozen donors hold two separately-collapsed sets concatenated together,
and in those the same UMI can carry two genuinely *different* molecules
(measured: 773 igh, 1,190 igk, 263 igl UMIs).  Step 2 cannot merge those, since
they are not one molecule.  Appending the digest makes the uid unique in every
donor without discarding anything, so partis has nothing to rename and
`--crash-on-duplicate-uids` becomes a real guard.

It also makes the uid *identical across donors* for the same molecule, so the
cross-donor sharing gate reads the molecule key straight off the id rather than
recovering a UMI from a partis uid.
"""
import collections
import csv
import hashlib
import os
import re

UID_HASH_LEN = 8  # hex chars of the sequence digest appended to the uid

# `ambi` marks an ambiguous gene call, appended to the gene name and sometimes
# with a suffix of its own: IGHV1-69ambi*01, IGHV3-30ambi2*01, IGHV3-30ambi-42*01.
AMBI_RE = re.compile(r'ambi[^*]*')

CALL_COLUMNS = ('v_call', 'd_call', 'j_call')

SAMPLE_COLUMN = 'sample_id'


# ----------------------------------------------------------------------------------------
def load_reference_genes(glfo_dir, locus):
    """Gene names, alleles stripped, from a germline set directory."""
    genes = set()
    for region in 'vdj':
        fname = '%s/%s/%s%s.fasta' % (glfo_dir, locus, locus, region)
        if not os.path.exists(fname):  # light chains have no d
            continue
        with open(fname) as gfile:
            genes |= set(l[1:].split('*')[0].strip() for l in gfile if l.startswith('>'))
    if len(genes) == 0:
        raise Exception('no germline genes under %s/%s' % (glfo_dir, locus))
    return genes


# ----------------------------------------------------------------------------------------
def resolve_ambi_call(call, reference_genes, unresolved):
    """Turn one `ambi`-marked gene call into the gene it names.

    The marker is appended to the gene name, and gene names contain hyphens, so
    two shapes look alike and only the reference can tell them apart:

      IGHV1-69ambi*01   -- IGHV1-69 is a gene; `ambi` is a marker on it
      IGLJ2-3ambi*01    -- IGLJ2-3 is *not* a gene; it means IGLJ2 or IGLJ3

    So: strip the marker; if that names a gene, keep it; otherwise drop the
    trailing -N and keep the first candidate.  If neither is in the reference
    the stripped name is kept and counted in <unresolved> -- IGKV1-13 is a real
    gene this reference lacks (it carries IGKV1D-13), and the master calls it
    850 times without the marker, so the marked and unmarked calls have to end
    up as the same thing.

    Do NOT do this with a whole-line regex.  The `sed -i 's/-[^-]*ambi//g'` this
    replaces ran over every column of the tsv, matched across tab boundaries,
    and -- by eating from the last hyphen -- turned IGHV1-69ambi*01 into
    IGHV1*01, which is not a gene and so vanished from every germline set built
    downstream.
    """
    stripped = AMBI_RE.sub('', call)
    gene = stripped.split('*')[0]
    if gene in reference_genes:
        return stripped
    trimmed = gene.rsplit('-', 1)[0]
    if trimmed in reference_genes:
        return stripped.replace(gene, trimmed, 1)
    unresolved[gene] += 1
    return stripped


# ----------------------------------------------------------------------------------------
def resolve_ambi_row(parts, call_indices, reference_genes, unresolved):
    """Rewrite a row's gene-call columns in place. Only those columns."""
    for icol in call_indices:
        if 'ambi' not in parts[icol]:
            continue
        parts[icol] = ','.join(resolve_ambi_call(c, reference_genes, unresolved) if 'ambi' in c else c
                               for c in parts[icol].split(','))  # a column can hold several calls, only some of them marked, and an unmarked one must come through untouched



# ----------------------------------------------------------------------------------------
def split_master_tsv(master_fn, out_paths, reference_genes):
    """One pass over the master changeo TSV, writing each wanted subject's rows.

    This replaces a `head -n1` plus `grep '^<sid>'` per subject, which re-read
    the whole file once per subject: 180 subjects against an 11-17 GB master is
    ~6 TB of reads per locus per rebuild, measured at ~88 s per subject, so
    about 13 h for the cohort.  One pass is 34 GB and a few minutes.

    It also removes the prefix-matching hazard in that grep, which is
    unanchored on the right: it would have handed a subject its own rows *plus*
    those of any subject whose id extends it.  No id in this cohort is a prefix
    of another, so nothing was wrong -- but routing on the sample_id column
    cannot go wrong at all.

    Ambiguous gene calls are resolved here too, per column rather than by a
    regex over the whole line -- see resolve_ambi_call.

    Args:
        master_fn: the master_<LOCUS>_changeo.tsv.
        out_paths: subject id -> path to write. Only these subjects are written.
        reference_genes: gene names from load_reference_genes, for resolving
            `ambi` calls.

    Returns:
        subject id -> number of rows written.
    """
    counts = collections.OrderedDict((sid, 0) for sid in out_paths)
    unresolved = collections.Counter()
    with open(master_fn) as mfile:
        headline = mfile.readline()
        columns = headline.rstrip('\n').split('\t')
        i_sample = columns.index(SAMPLE_COLUMN)
        call_indices = [columns.index(c) for c in CALL_COLUMNS if c in columns]
        handles = {}
        try:
            for sid, path in out_paths.items():
                handles[sid] = open(path, 'w')
                handles[sid].write(headline)
            for line in mfile:
                sid = line.split('\t', i_sample + 1)[i_sample]
                if sid not in handles:  # a subject we don't want, or one of the header lines the master file repeats throughout (their sample_id column reads 'sample_id', which is nobody)
                    continue
                if 'ambi' in line:  # cheap substring test, so untouched rows never get split into columns
                    parts = line.rstrip('\n').split('\t')
                    resolve_ambi_row(parts, call_indices, reference_genes, unresolved)
                    line = '\t'.join(parts) + '\n'
                handles[sid].write(line)
                counts[sid] += 1
        finally:
            for handle in handles.values():
                handle.close()
    empty = [sid for sid, n in counts.items() if n == 0]
    if len(empty) > 0:
        raise Exception('no rows in %s for %d subject(s) we meant to write: %s' % (master_fn, len(empty), ' '.join(empty)))
    print('    split %s into %d subject files, %d rows' % (os.path.basename(master_fn), len(counts), sum(counts.values())))
    for gene, n in sorted(unresolved.items(), key=lambda kv: -kv[1]):
        print('      note: %d calls of \'%s\', which is not in the reference germline set, so partis will drop them (as it does that gene\'s unmarked calls)' % (n, gene))
    return counts


# ----------------------------------------------------------------------------------------
def molecule_digest(sequence):
    """Hex digest identifying a sequence, shared with data-central's sharing table."""
    return hashlib.md5(sequence.encode()).hexdigest()


# ----------------------------------------------------------------------------------------
def read_exclusions(sharing_fn, shared_fn):
    """Read the two artifacts written by data-central scripts/vanwinkle_sharing.py.

    Returns:
        (dropped, delete_keys) where <dropped> is the set of (locus, donor)
        pairs to skip entirely, and <delete_keys> maps a kept (locus, donor) to
        the set of (umi, digest) molecules to delete from it.
    """
    dropped, delete_keys = set(), {}
    with open(sharing_fn) as sfile:
        for row in csv.DictReader(sfile):
            if row['decision'] == 'drop':
                dropped.add((row['locus'], row['donor']))
            elif row['decision'] == 'keep':
                delete_keys[(row['locus'], row['donor'])] = set()  # so a kept donor missing from <shared_fn> reads as 'nothing to delete' rather than as a missing key
            else:
                raise Exception('unexpected decision \'%s\' in %s' % (row['decision'], sharing_fn))
    with open(shared_fn) as sfile:
        for row in csv.DictReader(sfile):
            key = (row['locus'], row['donor'])
            if key not in delete_keys:
                raise Exception('%s has rows to delete in %s but is not a kept donor in %s' % (str(key), shared_fn, sharing_fn))
            delete_keys[key].add((row['umi'], row['sequence_md5']))
    print('    read %d dropped (locus, donor) pairs and %d molecules to delete from %d kept ones' % (len(dropped), sum(len(v) for v in delete_keys.values()), len(delete_keys)))
    return dropped, delete_keys


# ----------------------------------------------------------------------------------------
def count_index(fields, key):
    """Index of the single <key>=N field in a split changeo sequence_id."""
    indices = [i for i, f in enumerate(fields) if f.startswith(key + '=')]
    if len(indices) != 1:
        raise Exception('found %d %s fields in sequence id \'%s\'' % (len(indices), key, '|'.join(fields)))
    return indices[0]


# ----------------------------------------------------------------------------------------
def add_counts(fields, other):
    """Add <other>'s CONSCOUNT and DUPCOUNT into <fields>, in place.

    Two rows of one molecule are that molecule's read support split across the
    two separately-collapsed sets concatenated into the file, so the counts add.
    """
    for key in ('CONSCOUNT', 'DUPCOUNT'):
        idx, odx = count_index(fields, key), count_index(other, key)
        total = int(fields[idx].split('=')[1]) + int(other[odx].split('=')[1])
        fields[idx] = '%s=%d' % (key, total)


# ----------------------------------------------------------------------------------------
def clean_donor_tsv(sfn, delete_keys):
    """Rewrite <sfn> in place: delete, then collapse, then make uids unique.

    In that order, since a molecule's two rows carry the same digest and so
    would still collide if the uid were rewritten first.

      1. delete every row whose molecule (umi, digest) is in <delete_keys>
         -- the wide rule, so a molecule whose only partner was a dropped donor
         goes too
      2. collapse rows sharing a molecule to one, summing CONSCOUNT/DUPCOUNT
      3. replace the uid (first |-separated field of sequence_id, which is all
         partis keeps of the header) with <umi>.<digest[:8]>

    After this the uid is unique within the donor, so partis has nothing to
    rename, and identical across donors for the same molecule, so the split's
    ids carry the molecule key directly.
    """
    with open(sfn) as sfile:
        headline = sfile.readline()
        cols = headline.rstrip('\n').split('\t')
        i_sid, i_seq = cols.index('sequence_id'), cols.index('sequence')
        kept = collections.OrderedDict()
        n_read, n_deleted, n_collapsed = 0, 0, 0
        for line in sfile:
            parts = line.rstrip('\n').split('\t')
            fields = parts[i_sid].split('|')
            mol = (fields[0], molecule_digest(parts[i_seq]))
            n_read += 1
            if mol in delete_keys:
                n_deleted += 1
            elif mol in kept:
                add_counts(kept[mol][1], fields)
                n_collapsed += 1
            else:
                kept[mol] = (parts, fields)
    with open(sfn, 'w') as sfile:
        sfile.write(headline)
        for (umi, digest), (parts, fields) in kept.items():
            fields[0] = '%s.%s' % (umi, digest[:UID_HASH_LEN])
            parts[i_sid] = '|'.join(fields)
            sfile.write('\t'.join(parts) + '\n')
    print('        %s: %d rows in, %d deleted, %d collapsed, %d out' % (os.path.basename(sfn), n_read, n_deleted, n_collapsed, len(kept)))


# ----------------------------------------------------------------------------------------
