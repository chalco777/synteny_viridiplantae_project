#!/usr/bin/env python3

"""
Script Function: Process a GFF file to standardize chromosome names by extracting species prefixes
and assigning sequential numbers to chromosomes.
"""

import re

def extract_species_prefix(chrom_name):
    """
    Extracts the species prefix from a chromosome name.

    Parameters:
        chrom_name (str): Chromosome name from the GFF file.

    Returns:
        str: The first two characters of the chromosome name.
    """
    return chrom_name[:2]

def extract_chrom_number(chrom_name):
    """
    Extracts a numeric identifier from the chromosome name.

    Parameters:
        chrom_name (str): Chromosome name from the GFF file.

    Returns:
        int or None: The last 5 digits of the last number in the name, or None if no numbers are found.
    """
    chrom_base = chrom_name.split('.')[0]  # Remove suffix after a dot
    numbers = re.findall(r'\d+', chrom_base)  # Extract all numbers
    if numbers:
        last_number = numbers[-1]
        return int(last_number[-5:].lstrip('0') or '0')  # Convert to integer, remove leading zeros
    return None

def build_chromosome_mapping(gff_file):
    """
    Creates a mapping of original chromosome names to standardized names.

    Parameters:
        gff_file (str): Path to the input GFF file.

    Returns:
        dict: Mapping of original chromosome names to new names.
    """
    species_chrom_numbers = {}

    # Extract chromosome and species information
    with open(gff_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split('\t')
            chrom_name = fields[0]
            species_prefix = extract_species_prefix(chrom_name)
            chrom_number = extract_chrom_number(chrom_name)
            if chrom_number is None:
                continue
            if species_prefix not in species_chrom_numbers:
                species_chrom_numbers[species_prefix] = {}
            species_chrom_numbers[species_prefix][chrom_name] = chrom_number

    # Create mapping for standardized chromosome names
    chrom_name_mapping = {}
    for species_prefix, chrom_dict in species_chrom_numbers.items():
        sorted_chroms = sorted(chrom_dict.items(), key=lambda x: x[1])  # Sort by extracted number
        for idx, (chrom_name, _) in enumerate(sorted_chroms, start=1):
            new_chrom_name = f"{species_prefix}{idx}"
            chrom_name_mapping[chrom_name] = new_chrom_name

    return chrom_name_mapping

def replace_chrom_names_in_gff(gff_file, output_file, chrom_name_mapping):
    """
    Replaces chromosome names in a GFF file with standardized names.

    Parameters:
        gff_file (str): Path to the input GFF file.
        output_file (str): Path to the output GFF file.
        chrom_name_mapping (dict): Mapping of original to new chromosome names.
    """
    with open(gff_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if not line.strip():
                fout.write(line)
                continue
            fields = line.strip().split('\t')
            chrom_name = fields[0]
            new_chrom_name = chrom_name_mapping.get(chrom_name, chrom_name)
            fields[0] = new_chrom_name
            fout.write('\t'.join(fields) + '\n')

if __name__ == "__main__":
    """
    Main script execution:
        - Processes the input GFF file.
        - Standardizes chromosome names based on species prefix and numeric identifiers.
        - Outputs a new GFF file with updated chromosome names.
    """
    gff_file = '/media/crowfoot2/DATOS/CHALCO/eproceso/mcscanx_analysis_relaxed/all.gff'
    output_file = '/media/crowfoot2/DATOS/CHALCO/eproceso/mcscanx_analysis_relaxed/all_new.gff'

    # Build chromosome name mapping
    chrom_name_mapping = build_chromosome_mapping(gff_file)

    # Replace chromosome names in the GFF file
    replace_chrom_names_in_gff(gff_file, output_file, chrom_name_mapping)

    print("The file 'all_new.gff' has been generated with updated chromosome names.")
