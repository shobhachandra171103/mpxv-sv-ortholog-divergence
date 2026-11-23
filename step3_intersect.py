import subprocess
import pandas as pd
import os

# Files
HOTSPOTS = "results/MPXV_SV_hotspots.tab"
REPEATS = "NC_063383_repeats.bed"
OUTPUT = "results/SVs_overlapping_repeats.tab"

if not os.path.exists(REPEATS):
    print(f"Error: {REPEATS} not found. Run trf and the conversion script first.")
    exit(1)

print("Intersecting Hotspots with TRF Repeats...")

# awk command to rearrange columns for bedtools:
# Original tab: clusterID(1), count(2), type(3), chrom(4), start(5), end(6)
# New order for stdin: chrom(4), start(5), end(6), clusterID(1), count(2), type(3)
cmd = (
    f"tail -n +2 {HOTSPOTS} | "
    f"awk 'BEGIN {{OFS=\"\\t\"}} {{print $4, $5, $6, $1, $2, $3}}' | "
    f"bedtools intersect -a stdin -b {REPEATS} -wa -wb > {OUTPUT}"
)

try:
    subprocess.run(cmd, shell=True, check=True)
except subprocess.CalledProcessError:
    print("Intersection failed. Check if column 4, 5, 6 in hotspots.tab are valid.")
    exit(1)

# Analysis
if os.path.exists(OUTPUT) and os.path.getsize(OUTPUT) > 0:
    # Col 0-5: SV info | Col 6-9: Repeat info
    res = pd.read_csv(OUTPUT, sep='\t', header=None)
    print("\n--- Structural Variants in Repeat Regions ---")
    
    # Track which clusters we've already printed to keep it clean
    seen_clusters = set()
    for _, row in res.iterrows():
        cluster_id = row[3]
        if cluster_id not in seen_clusters:
            print(f"Hotspot Cluster {cluster_id} (Freq: {row[4]}) overlaps with {row[9]} at {row[7]}-{row[8]}")
            seen_clusters.add(cluster_id)
else:
    print("\nNo overlaps found. The identified SVs are likely in non-repetitive regions.")