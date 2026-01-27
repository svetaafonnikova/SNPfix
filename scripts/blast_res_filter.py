import sys
import pandas as pd

# Input files
blast_file = sys.argv[1]
output = sys.argv[2]
failed = sys.argv[3]



# Define column names as per your BLAST tabular format
colnames = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
]

# Load the file (no header in file)
blast_df = pd.read_csv(blast_file, sep=r'\s+', names=colnames)

# Convert relevant columns to numeric
numeric_cols = ['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                'sstart', 'send', 'evalue', 'bitscore']
for col in numeric_cols:
    blast_df[col] = pd.to_numeric(blast_df[col], errors='coerce')

# Apply filters
passed_filter = (
    (blast_df['mismatch'] <= 1) &
    (blast_df['evalue'] < 1e-20) &
    (blast_df['length'] > 51)
)
filtered_df = blast_df[passed_filter]

# Find qseqid present in original but not in filtered
all_qseqid = set(blast_df['qseqid'].unique())
filtered_qseqid = set(filtered_df['qseqid'].unique())
qseqid_not_in_filtered = sorted(all_qseqid - filtered_qseqid)

# Sort by pident and bitscore
sorted_df = filtered_df.sort_values(by=['pident', 'bitscore'], ascending=[False, False])

# Save filtered results
sorted_df.to_csv(output, index=False, sep="\t")

# Save markers failed to pass
with open(failed, 'w') as f:
    #f.write(f'Name\tChr\Position\n')
    for q in qseqid_not_in_filtered:
        f.write(f'{q}\tn.a\n')

# Output
print(f"The initial number of markers: {len(blast_df.drop_duplicates(subset=['qseqid']))}")
print(f"Number of markers after filtration {len(sorted_df.drop_duplicates(subset=['qseqid']))} ")

print("Unique qseqid that failed the filter:")
print(qseqid_not_in_filtered)