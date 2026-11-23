import pandas as pd
from pycirclize import Circos
import matplotlib.pyplot as plt

# --- 1. SETUP ---
genome_size = 197209
sectors = {"MPXV": genome_size}
# Increased space and figure size to accommodate the legend on the right
circos = Circos(sectors, space=5)

# --- 2. LOAD DATA ---
sv_df = pd.read_csv("results/clean_circos.bed", sep="\t", header=None,
                    names=["chrom", "start", "end", "type", "freq"])

for col in ['start', 'end', 'freq']:
    sv_df[col] = pd.to_numeric(sv_df[col], errors='coerce')
sv_df = sv_df.dropna()

repeats_df = pd.read_csv("NC_063383_repeats.bed", sep="\t", header=None)
repeats_df = repeats_df[[0, 1, 2]]
repeats_df.columns = ["chrom", "start", "end"]

# --- 3. CREATE PLOT ---
# Color Palette
color_gap = "#D62728"      # Sophisticated Crimson Red
color_brk = "#1F77B4"      # Deep Steel Blue
color_repeat = "#2CA02C"   # Forest Green (Better visibility for repeats)

for sector in circos.sectors:
    # Track 1: Outer Coordinates
    track1 = sector.add_track((95, 100))
    track1.axis(fc="#34495E", ec="none") 
    # Removed label_color to prevent the TypeError
    track1.xticks_by_interval(interval=20000, label_size=7, 
                              label_formatter=lambda v: f"{v/1000:.0f}k")

    # Track 2: Tandem Repeats
    # Placed slightly outer to show their relationship with SVs
    for _, row in repeats_df.iterrows():
        sector.rect(start=float(row["start"]), end=float(row["end"]), 
                   r_lim=(85, 92), fc=color_repeat, ec="none", alpha=0.5)
    
    # Track 3: SV Frequency
    max_f = sv_df["freq"].max()
    for _, row in sv_df.iterrows():
        color = color_gap if row["type"] == "GAP" else color_brk
        h_ratio = float(row["freq"]) / max_f
        r_top = 50 + (h_ratio * 30)
        
        sector.rect(start=float(row["start"]), end=float(row["end"]), 
                   r_lim=(50, r_top), fc=color, ec="none", alpha=0.8)
# --- 4. EXPORT & LEGEND ---
fig = circos.plotfig()

# Defining the legend handles
handles = [
    plt.Rectangle((0,0),1,1, fc=color_gap, label="Genomic GAP (Deletion)"),
    plt.Rectangle((0,0),1,1, fc=color_brk, label="Breakpoint (BRK)"),
    plt.Rectangle((0,0),1,1, fc=color_repeat, label="Tandem Repeats (ITR/STR)")
]

# Moving Legend to the Right
fig.legend(handles=handles, 
           loc="center left", 
           bbox_to_anchor=(1.05, 0.5), 
           frameon=True, 
           fontsize=10, 
           title="MPXV Pangenome Key")

# Fix: Use pad_inches instead of bbox_inches per the error message
circos.savefig("results/MPXV_Circos_Final_v2.png", pad_inches=0.5)
print("\nSuccess! Final presentation map saved to results/MPXV_Circos_Final_v2.png")