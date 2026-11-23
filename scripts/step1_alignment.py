import os
import subprocess

# Configuration based on your actual filenames
REF = "NC_063383.fasta"
QUERY_DIR = "1236_genomes/"
OUT_DIR = "results/individual_svs/"

# Get list of all fasta files
genomes = [f for f in os.listdir(QUERY_DIR) if f.endswith(('.fasta', '.fa'))]

print(f"Found {len(genomes)} genomes. Starting alignment...")

for i, g_file in enumerate(genomes):
    sample_id = os.path.splitext(g_file)[0]
    query_path = os.path.join(QUERY_DIR, g_file)
    prefix = os.path.join(OUT_DIR, sample_id)
    
    # 1. Nucmer Alignment
    # --maxmatch is used for high sensitivity in viral genomes
    subprocess.run(["nucmer", "--maxmatch", "-p", prefix, REF, query_path], check=True)

    # 2. Delta-filter 
    # -q keeps the best alignment for each query segment
    delta_file = f"{prefix}.delta"
    filtered_delta = f"{prefix}.filter.delta"
    with open(filtered_delta, "w") as f_out:
        subprocess.run(["delta-filter", "-q", delta_file], stdout=f_out, check=True)

    # 3. Call SVs using show-diff
    diff_file = f"{prefix}.diff"
    with open(diff_file, "w") as f_out:
        subprocess.run(["show-diff", filtered_delta], stdout=f_out, check=True)

    if (i + 1) % 10 == 0:
        print(f"Completed {i+1}/{len(genomes)}...")

print("Alignment and SV calling complete.")
