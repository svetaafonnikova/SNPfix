import sys
import pandas as pd



outpos = sys.argv[1]  
markers_na = sys.argv[2]
markers_failed_blast = sys.argv[3]
failed_filt = sys.argv[4]
output_file = sys.argv[5]

header = pd.read_csv(outpos, sep='\t', nrows=0).columns.tolist()

df1 = pd.read_csv(outpos, sep='\t')

df2 = pd.read_csv(markers_na, sep='\t', header=None, names=header)
df3 = pd.read_csv(markers_failed_blast, sep='\t', header=None, names=header)
df4 = pd.read_csv(failed_filt, sep='\t', header=None, names=header)


merged_df = pd.concat([df1, df2, df3, df4], ignore_index=True)


merged_df["Position"] = pd.to_numeric(merged_df["Position"], errors='coerce').astype('Int64')

print(merged_df.shape)
print(merged_df.columns.tolist())
print(merged_df.dtypes)
print(merged_df.head())

# Sort by Chr, then Position
final_df = merged_df.sort_values(by=["Chr", "Position"])


final_df.to_csv(output_file, sep='\t', index=False)
print("Save completed")
