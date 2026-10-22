import sys
import pandas as pd


blast_file = sys.argv[1]
output = sys.argv[2]
PARALOG_RATIO = float(sys.argv[3]) if len(sys.argv) > 3 else 0.90
MULTI_RATIO = float(sys.argv[4]) if len(sys.argv) > 4 else 0.98
MERGE_DIST = int(sys.argv[5]) if len(sys.argv) > 5 else 1000
assert PARALOG_RATIO <= MULTI_RATIO, 'PARALOG_RATIO must be <= MULTI_RATIO'

ALT_COLS = ['sseqid', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'qseq', 'sseq']

df = pd.read_csv(blast_file, sep='\t')
df.columns = df.columns.str.strip()
for col in ['pident', 'bitscore', 'length', 'evalue', 'sstart', 'send']:
    df[col] = pd.to_numeric(df[col], errors='coerce')
df['slow'] = df[['sstart', 'send']].min(axis=1)


def loci_of(group):
    """One row per locus (its best HSP), sorted best first."""
    group = group.sort_values(['sseqid', 'slow'])
    new_locus = (group['sseqid'] != group['sseqid'].shift()) | \
                (group['slow'] - group['slow'].shift() > MERGE_DIST)
    group = group.assign(locus=new_locus.cumsum())
    group = group.sort_values(['bitscore', 'pident', 'evalue'],
                              ascending=[False, False, True])
    return group.drop_duplicates('locus')


rows = []
for q, group in df.groupby('qseqid', sort=False):
    loci = loci_of(group)
    best = loci['bitscore'].iloc[0]
    second = loci['bitscore'].iloc[1] if len(loci) > 1 else 0.0
    ratio = round(second / best, 3) if best > 0 else 0.0

    if ratio >= MULTI_RATIO:
        sel = loci[loci['bitscore'] >= MULTI_RATIO * best].copy()
        sel['Status'] = 'Multi-mapped'
    else:
        sel = loci.iloc[[0]].copy()
        if ratio >= PARALOG_RATIO:
            sel['Status'] = 'Unique_paralog'
            alt = loci.iloc[1]
            for c in ALT_COLS:
                if c in alt.index:
                    sel['alt_' + c] = alt[c]
        else:
            sel['Status'] = 'Unique'
    sel['second_best_ratio'] = ratio
    rows.append(sel)

best_hits = pd.concat(rows, ignore_index=True).drop(columns=['slow', 'locus'])
best_hits.to_csv(output, index=False, sep='\t')

print(best_hits.drop_duplicates('qseqid')['Status'].value_counts().to_string())
print(f'PARALOG_RATIO={PARALOG_RATIO}, MULTI_RATIO={MULTI_RATIO}, MERGE_DIST={MERGE_DIST}')
