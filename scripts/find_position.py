import sys
import pandas as pd



# Function to calculate position
def calc_position(row):
    marker_id = row['Name']
    size = marker_length.get(marker_id, 50)

    # Forward strand
    if row['sstart'] < row['send']:
        return row['sstart'] + size
    # Reverse strand
    else:
        return row['sstart'] - size

# Import data
pos = sys.argv[1]
len_data = sys.argv[2]
out = sys.argv[3]

# Read data
df = pd.read_csv(pos, sep='\\s+')
length_df = pd.read_csv(len_data, sep='\\s+')
marker_length = dict(zip(length_df['Marker'], length_df['Length']))


# Rename columns before applying function
df = df.rename(columns={'qseqid': 'Name', 'sseqid': 'Chr'})

# Calculate position for each row
df['Position'] = df.apply(calc_position, axis=1)

df['Position'] = df['Position'].astype(int).astype(str)


# Save only the requested columns as tab-separated
df[['Name', 'Chr', 'Position']].to_csv(out, sep='\t', index=False)

print(f'Number of markers in file {len(df)}')
