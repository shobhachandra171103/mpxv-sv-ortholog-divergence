import os

# Keywords to look for in your existing annotation files
keywords = ['transposase', 'integrase', 'mobile element', 'insertion sequence']

# Check your annotation results (assuming you have .gff files from Objective 1/2)
# If you don't have GFFs for all 1236, we can search the protein sequences
print("Searching for Mobile Element signatures in annotations...")
# (Logic to parse GFF or Protein FASTA)