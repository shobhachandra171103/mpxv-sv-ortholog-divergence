# MPXV Structural Variants & Ortholog Divergence (`mpxv-structural-variants-orthologs`)

> **Objective 3** of the B.Sc. (Life Sciences) 4th-Year Dissertation: *"Multi-Scale Genomic Architecture of Monkeypox Virus Human Adaptation: Mutation Hotspots, Structural Variants, and Ortholog Divergence in Sustained-Transmission Lineages"* (Department of Zoology, Deshbandhu College, University of Delhi, 2026).

## Overview

This repository documents the workflow, scripts, and analytical outputs for investigating **structural variants (SVs)** and **ortholog divergence** across Monkeypox virus (MPXV) and *Orthopoxvirus* genomes, complementing genome-wide mutation hotspot mapping (Objective 1) and pangenome/phylogenetic analysis (Objective 2).

## Methodology & Pipeline

1. **Sequence Comparison & Structural Variant Detection:**
   * Compare assemblies/genomes to detect structural rearrangements, insertions/deletions, or copy/coordinate differences across lineages.
2. **Ortholog Identification & Divergence Analysis:**
   * Identify orthologous genes/proteins and quantify sequence divergence across MPXV and *Orthopoxvirus* comparisons.
3. **Summary & Visualization:**
   * Process tables and generate figures in Python (`pandas`, `Biopython`, `matplotlib`/`seaborn`) and R (`tidyverse`, `ggplot2`, `circlize`/`ComplexHeatmap`).

## Repository Structure

```text
mpxv-structural-variants-orthologs/
├── README.md
├── scripts/          # Bash, Python, and R scripts for structural variant & ortholog analyses
├── results/          # Summary tables and generated figures
└── docs/             # Workflow notes and parameter summaries
