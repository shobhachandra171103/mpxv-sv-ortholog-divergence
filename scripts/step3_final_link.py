import subprocess
import pandas as pd
import os

HOTSPOTS = "results/MPXV_SV_hotspots.tab"
DOMAINS = "results/mobile_domains.bed"
OUTPUT = "results/hotspots_vs_domains.tab"

# Ensure the output directory exists
os.makedirs("results", exist_ok=True)

# Rearrange Hotspots to BED and intersect
cmd = (
    f"tail -n +2 {HOTSPOTS} | "
    f"awk 'BEGIN {{OFS=\"\\t\"}} {{print $4, $5, $6, $1, $2, $3}}' | "
    f"bedtools intersect -a stdin -b {DOMAINS} -wa -wb > {OUTPUT}"
)

print("Linking Hotspots to Functional Domains...")
subprocess.run(cmd, shell=True, check=True)

if os.path.exists(OUTPUT) and os.path.getsize(OUTPUT) > 0:
    res = pd.read_csv(OUTPUT, sep='\t', header=None)
    print("\n--- Structural Variants Linked to Functional Repeat Domains ---")
    
    # Using 'seen' to prevent duplicate prints for the same overlap
    seen = set()
    for _, row in res.iterrows():
        cluster_id = row[3]
        attr = str(row[9]) # Ensure it's a string
        
        # Cleaner logic to find the product name
        product = "Unknown Protein"
        if "product=" in attr:
            product = attr.split("product=")[-1].split(";")[0]
        
        link_id = f"{cluster_id}_{product}"
        if link_id not in seen:
            print(f"Hotspot Cluster {cluster_id} (Freq: {row[4]}) overlaps with: {product}")
            seen.add(link_id)
else:
    print("\nNo direct overlaps. The SVs may be in intergenic regions or non-repeat genes.")
