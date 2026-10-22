import sys
import pandas as pd

blast_file = sys.argv[1]
output = sys.argv[2]
failed = sys.argv[3]
flank = sys.argv[4]

MIN_PIDENT = float(sys.argv[5]) if len(sys.argv) > 5 else 90.0
MIN_COVERAGE = float(sys.argv[6]) if len(sys.argv) > 6 else 0.80

colnames = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
]

with open(blast_file) as fh:
    first = fh.readline().split()
if len(first) == 14:
    colnames += ['qseq', 'sseq']

blast_df = pd.read_csv(blast_file, sep=r'\s+', names=colnames)

numeric_cols = ['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                'sstart', 'send', 'evalue', 'bitscore']
for col in numeric_cols:
    blast_df[col] = pd.to_numeric(blast_df[col], errors='coerce')
    
expected_length = 2*int(flank)

passed_filter = (
    (blast_df['pident'] >= MIN_PIDENT) &
    (blast_df['length'] >= expected_length * MIN_COVERAGE)
)
filtered_df = blast_df[passed_filter]

all_qseqid = set(blast_df['qseqid'].unique())
filtered_qseqid = set(filtered_df['qseqid'].unique())
qseqid_not_in_filtered = sorted(all_qseqid - filtered_qseqid)

sorted_df = filtered_df.sort_values(by=['pident', 'bitscore'], ascending=[False, False])

sorted_df.to_csv(output, index=False, sep="\t")

with open(failed, 'w') as f:
    for q in qseqid_not_in_filtered:
        f.write(f'{q}\tn.a.\tn.a.\tFailed_Filter\n')

print(f"The initial number of markers: {len(blast_df.drop_duplicates(subset=['qseqid']))}")
print(f"Number of markers after filtration {len(sorted_df.drop_duplicates(subset=['qseqid']))} ")

print("Unique qseqid that failed the filter:")
print(qseqid_not_in_filtered)