import pandas as pd
import os
import matplotlib.pyplot as plt

# 1. Load Data
sv_df = pd.read_csv("results/clustered_svs.bed", sep='\t', names=['chrom', 'start', 'end', 'type', 'sampleID', 'clusterID'])
clade_file = "results/phylogeny/nextclade_assignments.tsv"

if not os.path.exists(clade_file):
    print("Error: Nextclade file not found!")
    exit()

clades = pd.read_csv(clade_file, sep='\t')

# 2. Advanced Cleaning Function (Handles Accession.Version vs Accession)
def clean(name):
    name = str(name).replace(".fasta", "").replace(".fna", "")
    # Splits by space, pipe, underscore, and dot, then takes the first part
    return name.split()[0].split('|')[0].split('_')[0].split('.')[0].strip()

sv_df['match_id'] = sv_df['sampleID'].apply(clean)
clades['match_id'] = clades['seqName'].apply(clean)

# 3. Create a Presence/Absence Matrix
# We want to know for EVERY unique sample: What is its clade and does it have Cluster 4?
all_ids = clades[['match_id', 'clade']].drop_duplicates()
samples_with_c4 = sv_df[sv_df['clusterID'] == 4]['match_id'].unique()

all_ids['Has_Cluster4'] = all_ids['match_id'].apply(lambda x: 1 if x in samples_with_c4 else 0)

# 4. Final Merge Check
matched_count = all_ids['match_id'].nunique()
print(f"Successfully matched {matched_count} samples!")

if matched_count > 0:
    # 5. Calculate Statistics
    summary = all_ids.groupby('clade')['Has_Cluster4'].mean() * 100
    print("\n--- Structural Variant Distribution ---")
    print(summary)

    # 6. Plotting
    plt.figure(figsize=(10, 6))
    summary.plot(kind='bar', color=['#2c3e50', '#e74c3c', '#3498db', '#27ae60'])
    plt.title("Prevalence of Structural Hotspot Cluster 4 across Mpox Clades", fontsize=14)
    plt.ylabel("Percentage of Genomes carrying the GAP (%)", fontsize=12)
    plt.xlabel("Viral Clade", fontsize=12)
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    plt.savefig("results/phylogeny/SV_Clade_Distribution.png")
    all_ids.to_csv("results/phylogeny/final_correlation_data.csv", index=False)
    print("\nPlot saved to results/phylogeny/SV_Clade_Distribution.png")
else:
    print("Check your input files; no IDs matched between SV results and Nextclade.")
