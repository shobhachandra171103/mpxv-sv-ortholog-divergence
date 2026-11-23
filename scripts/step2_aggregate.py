import pandas as pd
import glob
import os
import subprocess

# 1. Parse .diff files
all_svs = []
diff_files = glob.glob("results/individual_svs/*.diff")

print(f"Parsing {len(diff_files)} files...")

for f in diff_files:
    sample_id = os.path.basename(f).replace(".diff", "")
    with open(f, 'r') as file:
        for line in file:
            if line.startswith("/") or line.startswith("NUCMER") or line.startswith("[SEQ]") or not line.strip():
                continue
            cols = line.split()
            if cols[0].startswith(">"):
                continue

            try:
                chrom = cols[0]
                sv_type = cols[1]
                # Standardize: Start must be less than End
                s1, e1 = int(cols[2]), int(cols[3])
                start, end = min(s1, e1), max(s1, e1)
                all_svs.append([chrom, start, end, sv_type, sample_id])
            except (ValueError, IndexError):
                continue

if not all_svs:
    print("No data parsed.")
    exit()

# 2. Save Unsorted BED
df = pd.DataFrame(all_svs)
unsorted_bed = "results/temp_unsorted.bed"
df.to_csv(unsorted_bed, sep='\t', index=False, header=False)

# 3. Use System Sort (This is what BEDTools wants)
sorted_bed = "results/MPXV_SVs_all_genomes.bed"
print("Sorting BED file...")
subprocess.run(f"sort -k1,1 -k2,2n {unsorted_bed} > {sorted_bed}", shell=True, check=True)

# 4. Cluster SVs
clustered_path = "results/clustered_svs.bed"
print("Clustering SVs...")
subprocess.run(f"bedtools cluster -i {sorted_bed} -d 200 > {clustered_path}", shell=True, check=True)

# 5. Frequency Analysis
c_df = pd.read_csv(clustered_path, sep='\t', names=['chrom', 'start', 'end', 'type', 'sampleID', 'clusterID'])
hotspots = c_df.groupby('clusterID').agg({
    'sampleID': 'nunique',
    'type': 'first',
    'chrom': 'first',
    'start': 'min',
    'end': 'max'
}).reset_index().rename(columns={'sampleID': 'genome_count'})

# Filter for hotspots (at least 12 genomes)
final_hotspots = hotspots[hotspots['genome_count'] >= 12].sort_values(by='genome_count', ascending=False)
final_hotspots.to_csv("results/MPXV_SV_hotspots.tab", sep='\t', index=False)

print(f"\n--- Success! ---")
print(f"Total SVs: {len(df)}")
print(f"Hotspots found: {len(final_hotspots)}")
if not final_hotspots.empty:
    print("\nTop Hotspots:")
    print(final_hotspots.head(10).to_string(index=False))

# Cleanup
if os.path.exists(unsorted_bed):
    os.remove(unsorted_bed)
