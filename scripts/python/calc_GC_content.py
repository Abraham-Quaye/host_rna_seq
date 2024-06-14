#!/usr/bin/env python3

import pysam
import glob
import os
import csv

# Path to the subdirectory containing BAM files
subdirectory_path = "results/hisat2/*.bam"

# Function to calculate GC content
def calculate_gc_content(bam_file):
    bam = pysam.AlignmentFile(bam_file, 'rb')
    total_bases = 0
    gc_bases = 0
    for read in bam.fetch():
        seq = read.query_sequence
        total_bases += len(seq)
        gc_bases += seq.count('G') + seq.count('C')
    bam.close()
    return (gc_bases / total_bases) * 100 if total_bases > 0 else 0

# Use glob to find all BAM files in the subdirectory and sort them alphabetically
bam_files = sorted(glob.glob(subdirectory_path))

# Create 'results' subdirectory if it doesn't exist
results_dir = 'results' 
os.makedirs(results_dir, exist_ok=True)

# Path to the output CSV file within the 'results' subdirectory
output_csv_path = os.path.join(results_dir, 'gc_content_results.csv')

# Open a CSV file to write the results
with open(output_csv_path, 'w', newline='') as csvfile:
    fieldnames = ['bamfile', 'gc_content']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    # Write the header
    writer.writeheader()

    # Iterate over the sorted list of BAM files and write the results
    for bam_file in bam_files:
        gc_content = calculate_gc_content(bam_file)
        # Write the results to the CSV file
        writer.writerow({'bamfile': os.path.basename(bam_file), 'gc_content': f'{gc_content:.2f}'})

print(f"GC content calculation is complete. Results have been saved to '{output_csv_path}'.")
