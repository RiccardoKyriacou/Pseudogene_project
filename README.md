# Pseudogene_project
Collection of scripts used to investigate Pseudogenes and stop-codon-readthrough in mammalian genomes

# Contents 
## get_gff3.py
- Script to download gff3/proteins/transcrpts from DNAzoo.com
- Usage: python3 get_gff.py --input [species text file] --type [gff3/proteins/transcripts]

## count_genes.py
- Script to count genes and pseudogenes for GFF files

## pseudogene_search.py
- Work in progress script that 1) masks GFF elements (genes) from a complete genome 2) Creates a list of genes from genome 3) BLASTn genes against masked genome to find pseudogenes  
- Usage: python pseudogene_search.py -g [GENOME FILE TO MASK] -a [GFF3 file]
- Usage with slurm: sbatch --export=GENOME=genome.fasta,GFF=gff.fasta pseudogene_search.slurm
