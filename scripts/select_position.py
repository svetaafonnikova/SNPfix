import sys
import pandas as pd

# Input files
blast_file = sys.argv[1]
chip = sys.argv[2]
output = sys.argv[3]



# Create dictionary
chip_dict = {}
with open(chip, encoding='utf-8') as f:
    next(f)  # skip header
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) < 3:
            continue
        name, _, pmap = fields[:3]
        if pmap.lower() != 'n.a.':
            chip_dict[name] = pmap



# Load the CSV file
df = pd.read_csv(blast_file, sep='\t')
df.columns = df.columns.str.strip()

# Convert columns to numeric for correct sorting/comparison
df['pident'] = pd.to_numeric(df['pident'], errors='coerce')
df['evalue'] = pd.to_numeric(df['evalue'], errors='coerce')


def select_best(group):
    print(group)
    # 1. Highest pident
    max_pident = group['pident'].max()
    group = group[group['pident'] == max_pident]
    # 2. Lowest evalue
    min_evalue = group['evalue'].min()
    group = group[group['evalue'] == min_evalue]
    # 3. If still tied, use preference dict
    if len(group) > 1:
        qseqid = group['qseqid'].iloc[0]
        preferred = chip_dict.get(qseqid)
        if preferred:
            preferred_group = group[group['sseqid'] == preferred]
            if not preferred_group.empty:
                return preferred_group.iloc[0]
     #Otherwise, just take the first
    return group.iloc[0]

# Apply the selection for each qseqid
#best_hits = (
#    df.groupby('qseqid', as_index=False)
#    .apply(lambda g: select_best(g.reset_index(drop=True)))
#    .reset_index(drop=True)
#)
best_hits = (
    df.groupby('qseqid', as_index=False)
    .apply(lambda g: select_best(g.reset_index(drop=True)), include_groups=True)
    .reset_index(drop=True)
)


# Save or display the result
best_hits.to_csv(output, index=False, sep='\t')

