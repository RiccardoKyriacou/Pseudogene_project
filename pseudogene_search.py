import os
import argparse
from Bio import SeqIO
from subprocess import call as unix
from subprocess import run, PIPE

'''
Usage: python pseudogene_search.py -g pseudogene_search_test/genome_test.fasta -a pseudogene_search_test/gff_test.gff3

Usage with slurm: sbatch --export=GENOME=genome.fasta,GFF=gff.fasta pseudogene_search.slurm

sbatch --export=GENOME=GCF_000001405.38_GRCh38.p12_genomic.fna,GFF=genomic.gff pseudogene_search.slurm
'''

# Function to get masked genome and list of protein coding sequences that have been masked out 
def get_masked_genome(genome_fasta, gff_file, genome_name):
    # Read genome sequences from FASTA file
    genome_records = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    # Dictionary comprehension to store scaffold sequences as lists of nucleotides
    genome_seqs = {scaffold.id: list(scaffold.seq) for scaffold in genome_records.values()}

    # List to store CDS nucleotide sequences
    cds_sequence_lst = []
    
    # Parse GFF file and mask regions
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            columns = line.strip().split('\t')
            scaffold = columns[0]
            feature_type = columns[2]
            # Here, we want to get the start position for the feature, but python lists start indexing at 0: 
            # 0 1 2 3 4 5 6 7 8 9 ...
            # A G C T A G C T A G ...
            start = int(columns[3]) - 1  # Convert to 0-based index 
            end = int(columns[4])  # End is exclusive
            # Check cds and scaffold id match 
            if feature_type == "CDS" and scaffold in genome_seqs:
                
                # Extract CDS sequence from genome and join into a single string
                cds_sequence = "".join(genome_seqs[scaffold][start:end])
                cds_sequence_lst.append(cds_sequence)
         
                # Mask each CDS in genome_seq using "N"
                for i in range(start, end):
                    genome_seqs[scaffold][i] = 'N'
    
    # Write masked genome to output file
    with open(f"masked_genome_{genome_name}", "w") as f:
        for scaffold_id, seq_list in genome_seqs.items():
            masked_genome_seq = "".join(seq_list)
            f.write(f">{scaffold_id}\n{masked_genome_seq}\n")

    # Create a temporary FASTA file for the CDS sequences
    cds_fasta = "cds.fasta"
    with open(cds_fasta, "w") as f:
        for i, cds in enumerate(cds_sequence_lst, start=1):
            # NOTE IDK why but this seems to remove loads of CDS reagions that are just Ns 
            # Gets us to about 24000 CDS for humans which seems right
            if "NNNNNN" not in cds:
                f.write(f">cds{i}\n{cds}\n")
    
    # Return list of sequences that have been masked and the masked genome sequence
    print(f"Total number of masked CDS {len(cds_sequence_lst)}")
    

# Fucntion to BLAST cds sequences against masked genomes 
def blast_masked_genome(masked_genome_file):
    # Run BLAST for each CDS with evalue cutoff 1e-5
    blast_cmd = f"tblastn -query cds.fasta -subject {masked_genome_file} -out BLAST_outfile -outfmt '6 qseq sseq pident length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 1e-5"

    unix(blast_cmd, shell=True)

    # Clean up temporary files
    #os.remove("cds.fasta")

# Main 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mask CDS regions in a genome based on GFF annotations.")
    parser.add_argument("-g", "--genome", type=str, help="Path to genome FASTA file", required=True)
    parser.add_argument("-a", "--annotation", type=str, help="Path to GFF file", required=True)
    args = parser.parse_args()

    # Get genome name
    genome_name = args.genome.split("/")[-1]
    # Get masked genome and list of cds regions masked out 
    #print(f"Masking genome...")
    #get_masked_genome(args.genome, args.annotation, genome_name)
    # tblastn masked genome
    print(f"Performing tblastn on masked genome...")
    masked_genome_file = f"masked_genome_{genome_name}"
    blast_masked_genome(masked_genome_file)
    
    #Quick BLAST count
    with open("BLAST_outfile", "r") as f:
        count = 0
        for line in f:
            count += 1
        print(f"Number of BLASTn hits = {count}")


#NOTE : 1) Mask out all exons from genome
#       2) Make list of all genes (62,000)
#       3) BLAST list of genes against exon masked genome 



# NOTE use this is previous fucntion doesnt work, BLASTDM made
# def blast_masked_genome_withdb(cds_sequence_lst, masked_genome_file):
#     # Create a temporary FASTA file for the CDS sequences
#     cds_fasta = "cds.fasta"
#     with open(cds_fasta, "w") as f:
#         for i, cds in enumerate(cds_sequence_lst, start=1):
#             f.write(f">cds{i}\n{cds}\n")
    
#     # Create BLAST database from masked genome
#     db_name = "masked_genome_db"
#     makeblastdb_cmd = f"makeblastdb -in {masked_genome_file} -dbtype nucl -parse_seqids -out {db_name}"
#     unix(makeblastdb_cmd, shell=True)
    
#     # Run BLAST for each CDS with evalue cutoff 1e-5
#     blast_results = []
#     blastn_cmd = f"tblastn -query cds.fasta -db {db_name} -outfmt '6 qseq sseq pident length mismatch gapopen qstart qend sstart send evalue bitscore' -evalue 1e-5"

#     result = run(blastn_cmd, capture_output=True, text=True, check=True)
#     blast_results.append(result.stdout)
    
#     # Clean up temporary files
#     os.remove(cds_fasta)
#     for ext in [".nhr", ".nin", ".nsq"]:
#         os.remove(f"{db_name}{ext}")
    
#     return blast_results

