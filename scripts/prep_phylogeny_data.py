import pandas as pd
import os

# 1. Load SV data
sv_df = pd.read_csv("results/clustered_svs.bed", sep='\t', names=['chrom', 'start', 'end', 'type', 'sampleID', 'clusterID'])

# 2. Load Nextclade data
clade_file = "results/phylogeny/nextclade_assignments.tsv"
if not os.path.exists(clade_file):
    print("Nextclade file not found!")
    exit()

clades = pd.read_csv(clade_file, sep='\t')

# DEBUG: Print first 3 names to see the format
print(f"SV sample example: '{sv_df['sampleID'].iloc[0]}'")
print(f"Nextclade sample example: '{clades['seqName'].iloc[0]}'")

# 3. Aggressive Cleaning Function
def clean(name):
    # Remove extension, take first part before space/pipe/underscore
    name = str(name).replace(".fasta", "").replace(".fna", "")
    name = name.split()[0].split('|')[0].split('_')[0]
    return name.strip()

sv_df['match_id'] = sv_df['sampleID'].apply(clean)
clades['match_id'] = clades['seqName'].apply(clean)

# 4. Perform the Merge
final_df = sv_df.merge(clades[['match_id', 'clade']], on='match_id', how='inner')

# 5. Summary
if len(final_df) == 0:
    print("\nSTILL NO MATCH. Check the examples above carefully.")
else:
    print(f"\nSuccess! Matched {final_df['match_id'].nunique()} unique samples.")
    # Create the binary presence/absence table
    pivot = final_df.groupby(['match_id', 'clade']).apply(lambda x: 1 if (x['clusterID'] == 4).any() else 0).reset_index()
    pivot.columns = ['sample', 'clade', 'Has_Cluster4']
    
    print("\n--- Structural Variant Distribution ---")
    summary = pivot.groupby('clade')['Has_Cluster4'].mean() * 100
    print("Percentage of each Clade carrying Cluster 4 GAP:")
    print(summary)

    final_df.to_csv("results/phylogeny/final_tree_metadata.csv", index=False)
