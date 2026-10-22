import sys
import pandas as pd


def walk_alignment(row, marker_q):
    step = 1 if row['sstart'] < row['send'] else -1
    qi, si = int(row['qstart']), int(row['sstart'])
    for qc, sc in zip(str(row['qseq']), str(row['sseq'])):
        if qc != '-' and qi == marker_q:
            return si if sc != '-' else None
        if qc != '-':
            qi += 1
        if sc != '-':
            si += step
    return None


def calc_position(row, name=None):
    offset = marker_length.get(name if name is not None else row['Name'], 50)
    marker_q = offset + 1
    if not (int(row['qstart']) <= marker_q <= int(row['qend'])):
        return None, 'Marker_outside_alignment'
    if has_aln:
        pos = walk_alignment(row, marker_q)
        return (pos, None) if pos is not None else (None, 'Marker_in_gap')
    dist = offset - (int(row['qstart']) - 1)
    if row['sstart'] < row['send']:
        return int(row['sstart']) + dist, None
    return int(row['sstart']) - dist, None


pos_file = sys.argv[1]
len_data = sys.argv[2]
out = sys.argv[3]

df = pd.read_csv(pos_file, sep='\t')
length_df = pd.read_csv(len_data, sep=r'\s+')
marker_length = dict(zip(length_df['Marker'], length_df['Length']))
has_aln = {'qseq', 'sseq'}.issubset(df.columns)
if not has_aln:
    print('Warning: no qseq/sseq columns - indels near the marker are not corrected')

df = df.rename(columns={'qseqid': 'Name', 'sseqid': 'Chr'})
res = [calc_position(r) for _, r in df.iterrows()]
df['Position'] = [str(int(p)) if f is None else 'n.a.' for p, f in res]
flag = pd.Series([f is not None for _, f in res], index=df.index)
df.loc[flag, 'Status'] = [f for _, f in res if f is not None]
df.loc[flag, 'Chr'] = 'n.a.'
if flag.any():
    print(f'{flag.sum()} hits without a usable marker position:')
    print(df.loc[flag, 'Status'].value_counts().to_string())


def alt_position(row):
    if row.get('Status') != 'Unique_paralog' or pd.isna(row.get('alt_sseqid')):
        return ''
    alt = {k: row['alt_' + k] for k in ['qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq']
           if 'alt_' + k in row.index}
    pos, fail = calc_position(alt, name=row['Name'])
    return f"{row['alt_sseqid']}:{int(pos)}" if fail is None else f"{row['alt_sseqid']}:n.a."


cols = ['Name', 'Chr', 'Position', 'Status']
if 'second_best_ratio' in df.columns:
    cols.append('second_best_ratio')
if 'alt_sseqid' in df.columns:
    df['alt_locus'] = [alt_position(r) for _, r in df.iterrows()]
    cols.append('alt_locus')
df[cols].to_csv(out, sep='\t', index=False)
print(f'Number of markers in file {len(df)}')
