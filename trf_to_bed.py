import sys

# Process the .dat file created by TRF
input_file = "NC_063383.fasta.2.7.7.80.10.50.500.dat"
output_file = "NC_063383_repeats.bed"

with open(input_file, 'r') as f, open(output_file, 'w') as out:
    for line in f:
        # TRF data lines start with an integer (the start position)
        parts = line.split()
        if len(parts) > 13 and parts[0].isdigit():
            start = parts[0]
            end = parts[1]
            period = parts[2]
            copy_num = parts[3]
            # Format: chrom, start, end, name (period_size)
            out.write(f"NC_063383.1\t{start}\t{end}\tTRF_period_{period}\n")

print(f"Created {output_file}")