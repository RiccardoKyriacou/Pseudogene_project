#!/usr/bin/env python3
import os
import glob
import argparse

# Author: Riccardo Kyriacou
# Date: 06/06/2024
# Usage: python3 count_genes.py --path /path/to/gff/files/

def count_genes(gff_file_content, species, gene_count_file):
    count = 0
    for line in gff_file_content:
        if "gene" in line:
            count += 1
    print(species, count)
    gene_count_file.write(f"{species}\t{count}\n")

def count_pseudo(gff_file_content, species, pseudo_count_file):
    count = 0
    for line in gff_file_content:
        if "pseudogene" in line:
            count += 1
    print(species, count)
    pseudo_count_file.write(f"{species}\t{count}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", type=str, help="Path to GFF files", required=True)
    args = parser.parse_args()

    # Check for valid path
    if not args.path.endswith('/'):
        args.path += '/'

    for gff in glob.glob(args.path + '*.gff3'):
        # Get species name from file list
        sp_name = os.path.basename(gff).split("_HiC")[0]
        # Read gff file contents
        with open(gff, 'r', encoding='utf-8') as f:
            gff_content = f.readlines()
        # Output files in append mode
        with open("gene_counts.tsv", "a", encoding='utf-8') as outf, open("pseudo_counts.tsv", "a", encoding='utf-8') as outf1:
            print(f"Counting genes for {sp_name}...")
            count_genes(gff_content, sp_name, outf)
            count_pseudo(gff_content, sp_name, outf1)
            
if __name__ == "__main__":
    main()
