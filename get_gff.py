#!/usr/bin/env python3
import os
import json
import requests
import argparse
from subprocess import call as unix

# Author: Riccardo Kyriacou
# Date: 01/05/2024
# Usage: python3 get_gff.py -i ../data/DNAZoo_species_list.txt

# Function to get species list from input file 
def get_species_lst(species_list_file):
    sp_lst = []
    with open(species_list_file) as f:
        for line in f:
            species_name = line.strip()  # Remove leading/trailing whitespace
            sp_lst.append(species_name)

    return sp_lst

# Code to make borders for error messages, makes std inpt easier to read
def print_border_1(message):
    border = "-" * (len(message) + 2)
    print(f"{border}\n{message}\n{border}\n")

def print_border_2(message):
    border = "+" * (len(message) + 2)
    print(f"{border}\n{message}\n{border}\n")

# Function to download genomes
def genome_download(sp, failed_downloads_file, file_type):
    
    # Construct the URL for the README.json file specific to each species
    url = f"https://dnazoo.s3.wasabisys.com/{sp}/README.json"

    try:
        # Retrieve the README.json file from DNAzoo
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for non-2xx status codes

        # Parse JSON data
        data = response.json()

        # Get the name of the file
        chromlength_assembly_name = data["chromlengthAssembly"]["name"]

        # Check if file already downloaded or if file size > 0 
        if (os.path.exists(f"{chromlength_assembly_name}.fasta_v2.functional.{file_type}.gz")) or ( os.path.getsize(f"{chromlength_assembly_name}.fasta_v2.functional.{file_type}.gz") < 0):
            print(f"{chromlength_assembly_name}.fasta_v2.functional.{file_type}.gz already exists. Skipping download.\n")
        elif (os.path.exists(f"{chromlength_assembly_name}.functional_v1.{file_type}gz")) or (( os.path.getsize(f"{chromlength_assembly_name}.functional_v1.{file_type}gz") < 0)):
            print(f"{chromlength_assembly_name}.functional_v1.{file_type}.gz already exists. Skipping download.\n")
        # If file not downloaded, use wget
        else:
            # Use wget to download either .fasta_v2.functional.gff3.gz or .functional_v1.gff3.gz
            # wget || wget - This will execute second command only if first failed
            wget = f"wget https://dnazoo.s3.wasabisys.com/{sp}/{chromlength_assembly_name}.fasta_v2.functional.{file_type}.gz || wget https://dnazoo.s3.wasabisys.com/{sp}/{chromlength_assembly_name}.functional_v1.{file_type}.gz"
            unix(wget, shell=True)

            # Check if wget worked by checking if the file exists
            if os.path.exists(f"{chromlength_assembly_name}.fasta_v2.functional.{file_type}.gz"):
                # TODO: Include these lines??
                print_border_1(f"Downloaded GFF file for {sp}")
            elif os.path.exists(f"{chromlength_assembly_name}.functional_v1.{file_type}.gz"):
                print_border_1(f"Downloaded GFF file for {sp}")
            else:
                print_border_2(f"{sp} download unsuccessful. File not found.")
                # Write failed species to file
                failed_downloads_file.write(f"{sp}\n")  

    # Except clauses for if the .JSON file doesn't exist or cannot be parsed
    except requests.exceptions.RequestException as e:
        print_border_2(f"Error downloading {url}: {e}")
        failed_downloads_file.write(f"{sp}\n")  # Write failed species to file
    except KeyError as e:
        print_border_2(f"Error parsing JSON data for {sp}: {e}")
        failed_downloads_file.write(f"{sp}\n")  # Write failed species to file
        
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="Species list file", required=True)
    parser.add_argument("-t", "--type", type=str, help="type of file (gff3, proteins, transcripts)", required=False, default = "gff3")
    args = parser.parse_args()

    # Code to check the file type flag
    if args.type == "gff3":
        pass 
    elif args.type == "proteins":
        pass
    elif args.type == "transcripts":
        pass
    else:
        print("Error in --type flag, options: gff3, proteins, transcripts")

    # Make gff dir
    os.makedirs("gff", exist_ok=True)

    sp_lst = get_species_lst(args.input)

    # Change to gff dir for download
    os.chdir("gff")

    # Open the file for writing failed downloads
    with open("failed_downloads.tsv", "w") as ouf:
        for species in sp_lst:
            # Download species
            # NOTE: Parallel(n_jobs=threads)(delayed(genome_download)(sp) for sp in sp_lst) if joblib allowed
            genome_download(species, ouf, args.type)

    print_border_1(f"Unzipping files")
    unix(f"for file in *.gz; do gzip -d '$file'; done", shell=True)


if __name__ == "__main__":
    main()
