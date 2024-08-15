import os
import argparse
from Bio import SeqIO
from subprocess import call as unix

'''
sbatch --export=GENOME=GCF_000001405.38_GRCh38.p12_genomic.fna,GFF=genomic.gff pseudogene_search.slurm
'''

def get_outgroup_blastdb(outgroup_cds_1, outgroup_cds_2, outgroup_cds_3, outgroup_cds_4, outgroup_cds_5):
    # Define the output file for combined sequences
    combined_file = "combined_outgroup_cds.fasta"
    
    # List of outgroup CDS files
    outgroup_cds_files = [outgroup_cds_1, outgroup_cds_2, outgroup_cds_3, outgroup_cds_4, outgroup_cds_5]
    
    # Create dict of sacc : sp_name
    species_dict = {}
    for outgroup_file in outgroup_cds_files:
        sp_name = os.path.basename(outgroup_file).split(".")[0]
        with open(outgroup_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                header = record.id
                # Header in form: 
                # >lcl|NC_036902.1_cds_XP_024208912.1_77757 [gene=LOC112206852] ...
                sacc = header.split("|")[1].split(" ")[0].strip()   
                species_dict[sacc] = sp_name
    
    for i, (key, value) in enumerate(species_dict.items()):
        if i >= 50:
            break
        print(f"{key}: {value}")

    # Combine all input files into one
    with open(combined_file, 'w') as outf:
        for outgroup_file in outgroup_cds_files:
            with open(outgroup_file, 'r') as f:
                outf.write(f.read())
    
    # Define the name of the BLAST database
    db_name = "mammals_blast_database"
    
    # Create the BLAST database from the combined file
    makeblastdb_cmd = f"makeblastdb -in {combined_file} -dbtype nucl -parse_seqids -out {db_name}"
    unix(makeblastdb_cmd, shell=True)

    return species_dict
  
# Function to get masked human genome 
def get_masked_genome(genome_fasta, gff_file, genome_name):
    # Read genome sequences from FASTA file
    genome_records = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    # Dictionary comprehension to store scaffold sequences as lists of nucleotides
    genome_seqs = {scaffold.id: list(scaffold.seq) for scaffold in genome_records.values()}
    
    # Parse GFF file and mask regions
    with open(gff_file, 'r') as gff:
        mask_count = 0 
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

            # Mask genome exons
            if feature_type == "exon" and scaffold in genome_seqs:
                # Mask each CDS in genome_seq using "N"
                for i in range(start, end):
                    genome_seqs[scaffold][i] = 'N'
                    mask_count += 1 # Count to check if all exons being masked
    
    # Write masked genome to output file
    with open(f"masked_genome_{genome_name}", "w") as f:
        for scaffold_id, seq_list in genome_seqs.items():
            masked_genome_seq = "".join(seq_list)
            f.write(f">{scaffold_id}\n{masked_genome_seq}\n")

# Fucntion to BLAST masked genome against mammal blastDB with evalue and %identity cutoff 
def blast_masked_genome(blast_db_file, masked_genome, perc_identity, sp_dict):
    # Define the output file name for BLAST results
    blast_output = "BLAST_outfile"
    
    # Construct the BLAST command
    blast_cmd = (
        f"blastn -query {masked_genome} -db {blast_db_file} "
        f"-out {blast_output} "
        f"-outfmt '6 qseq sseq pident length mismatch gapopen qstart qend sstart send evalue bitscore sacc' "
        f"-evalue 1e-5 -perc_identity {perc_identity}"
    )

    unix(blast_cmd, shell=True)

    # Post-process the BLAST output to include species names
    with open(blast_output, "r") as f, open("BLAST_outfile_with_species", "w") as outf:
        for line in f:
            fields = line.strip().split('\t')
            sacc = fields[12]
            species_name = sp_dict.get(sacc, "Unknown")
            outf.write(f"{line.strip()}\t{species_name}\n")


# Main 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mask CDS regions in a genome based on GFF annotations.")
    parser.add_argument("-g", "--genome", type=str, help="Path to genome FASTA file", required=True)
    parser.add_argument("-a", "--annotation", type=str, help="Path to GFF file", required=True)
    args = parser.parse_args()

    # Get outgroup genomes 
    outgroup_cds_1 = "mouse_cds.fna"
    outgroup_cds_2 = "goat_cds.fna"
    outgroup_cds_3 = "cow_cds.fna"
    outgroup_cds_4 = "chimp_cds.fna"
    outgroup_cds_5 = "beluga_cds.fna"

    #Get blastdb file
    species_dict = get_outgroup_blastdb(outgroup_cds_1, outgroup_cds_2, outgroup_cds_3, outgroup_cds_4, outgroup_cds_5)

    # Get genome name
    genome_name = args.genome.split("/")[-1]

    # Get masked genome and list of cds regions masked out 
    print(f"Masking genome...")
    get_masked_genome(args.genome, args.annotation, genome_name)
    
    # BLAST masked genome
    print(f"Performing BLAST on masked genome...")
    blast_db_file = "mammals_blast_database"
    masked_genome_file = f"masked_genome_{genome_name}"
    perc_identity = 95
    blast_masked_genome(blast_db_file, masked_genome_file, perc_identity, species_dict)
    
    #Quick BLAST count
    with open("BLAST_outfile", "r") as f:
        count = 0
        for line in f:
            count += 1
        print(f"Number of BLASTn hits = {count}")

    # Cleanup 
    if not os.path.exists("blast_DB"):
        os.makedirs("blast_DB")
    unix("mv mammals_blast_database* blast_DB/", shell=True)
    
    if not os.path.exists("data"):
        os.makedirs("data")
    unix("mv *.fna *.gff *.fasta data/", shell=True)

        
'''
Notes:
1) Look at distribution of % identity again 
2) Run again with species hits in the BLAST file
3) Look at overlap of hits
'''