#!/usr/bin/env python3

"""
Script Function: Standardize chromosome names in a FASTA file by assigning species prefixes 
and sequential numbers to chromosome headers.
"""

import re
import sys
import os

def extract_species_prefix(species_name):
    """
    Determines the species prefix based on the file name.

    Parameters:
        species_name (str): The name of the species derived from the file name.

    Returns:
        str: A two-letter species prefix (e.g., 'AT' for Arabidopsis).
    """
    if 'arabidopsis' in species_name.lower():
        return 'AT'
    elif 'chlamydomonas' in species_name.lower():
        return 'CR'
    elif 'cuscuta' in species_name.lower():
        return 'CU'
    else:
        return 'RO'  # Default prefix for unknown species

def extract_chrom_number(chrom_name):
    """
    Extracts a numeric identifier from the chromosome name.

    Parameters:
        chrom_name (str): The chromosome name from the FASTA header.

    Returns:
        int or None: The last 5 digits of the last number in the chromosome name, 
        or None if no numbers are found.
    """
    chrom_base = chrom_name.split('.')[0]  # Remove suffix after a dot
    numbers = re.findall(r'\d+', chrom_base)  # Extract all numbers
    if numbers:
        last_number = numbers[-1]
        return int(last_number[-5:].lstrip('0') or '0')  # Convert to integer, remove leading zeros
    return None

def build_chromosome_mapping(fasta_file, species_prefix):
    """
    Creates a mapping of original chromosome names to standardized names.

    Parameters:
        fasta_file (str): Path to the input FASTA file.
        species_prefix (str): Prefix to add to chromosome names.

    Returns:
        dict: Mapping of original chromosome names to new standardized names.
    """
    chrom_numbers = {}

    # Extract chromosome names and numbers from the FASTA headers
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line[1:].strip()
                chrom_name = header.split()[0]  # Take only the first word
                chrom_number = extract_chrom_number(chrom_name)
                if chrom_number is not None:
                    chrom_numbers[chrom_name] = chrom_number

    # Sort chromosomes by the extracted number and assign new names
    sorted_chroms = sorted(chrom_numbers.items(), key=lambda x: x[1])
    chrom_name_mapping = {chrom_name: f"{species_prefix}{idx}" for idx, (chrom_name, _) in enumerate(sorted_chroms, start=1)}

    return chrom_name_mapping

def process_fasta_file(fasta_file, output_file, chrom_name_mapping, species_prefix):
    """
    Processes a FASTA file to replace chromosome names with standardized names.

    Parameters:
        fasta_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file with updated headers.
        chrom_name_mapping (dict): Mapping of original to new chromosome names.
        species_prefix (str): Prefix to add to chromosome names.

    Writes:
        A new FASTA file with standardized chromosome names.
    """
    with open(fasta_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                header = line[1:].strip()
                chrom_name = header.split()[0]  # Take only the first word
                chrom_name_with_prefix = species_prefix + chrom_name  # Add species prefix
                new_chrom_name = chrom_name_mapping.get(chrom_name, chrom_name_with_prefix)
                fout.write(f">{new_chrom_name}\n")
            else:
                fout.write(line)

    print(f"Processed FASTA file: {output_file}")

if __name__ == "__main__":
    """
    Main script execution:
        - Reads a FASTA file.
        - Standardizes chromosome names based on species prefix and numeric identifiers.
        - Outputs a new FASTA file with updated headers.
    """
    if len(sys.argv) != 2:
        print("Usage: python rename_chromosomes_fasta.py <absolute_path_to_fasta_file>")
        sys.exit(1)

    fasta_path = sys.argv[1]
    if not os.path.isfile(fasta_path):
        print(f"Error: File '{fasta_path}' does not exist.")
        sys.exit(1)

    species_name = os.path.basename(fasta_path).split('.')[0]
    species_prefix = extract_species_prefix(species_name)
    output_file = fasta_path.replace('.fna', '_renamed.fna')
    chrom_name_mapping = build_chromosome_mapping(fasta_path, species_prefix)
    process_fasta_file(fasta_path, output_file, chrom_name_mapping, species_prefix)
