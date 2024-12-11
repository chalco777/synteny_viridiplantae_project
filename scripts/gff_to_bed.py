#!/usr/bin/env python3

"""
Script Function: Extract mRNAs that have corresponding CDS features and convert them into BED format.
"""

import sys

def gff_to_bed(gff_file, bed_file, species_prefix):
    """
    Converts a GFF file to a BED file containing mRNAs with corresponding CDS entries.

    Parameters:
        gff_file (str): Path to the input GFF file.
        bed_file (str): Path to the output BED file.
        species_prefix (str): Prefix to add to chromosome names in the output.

    Functionality:
        - Extracts mRNA and CDS information from the GFF file.
        - Links mRNAs to their corresponding CDS using the 'Parent' and 'protein_id' attributes.
        - Writes the BED file with chromosome, start, end, and protein ID for valid mRNAs.
    """
    mrna_info = {}  # Stores mRNA information: mRNA_ID -> {'chrom', 'start', 'end'}
    mrna_to_protein = {}  # Maps mRNA_ID to protein_id

    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, feature_type, start, end, strand, attributes = fields[0], fields[2], int(fields[3]) - 1, fields[4], fields[6], fields[8]

            # Process attributes into a dictionary
            attr_dict = dict(attr.split('=') for attr in attributes.strip(';').split(';') if '=' in attr)

            if feature_type == 'mRNA':
                mrna_id = attr_dict.get('ID', '')
                if mrna_id:
                    mrna_info[mrna_id] = {'chrom': chrom, 'start': start, 'end': int(end)}
            elif feature_type == 'CDS':
                parent_ids = attr_dict.get('Parent', '').split(',')
                protein_id = attr_dict.get('protein_id', '')
                if protein_id:
                    for parent_id in parent_ids:
                        mrna_to_protein[parent_id] = protein_id

    with open(bed_file, 'w') as bed:
        for mrna_id, info in mrna_info.items():
            if mrna_id in mrna_to_protein:
                chrom = species_prefix + info['chrom']
                start = info['start']
                end = info['end']
                protein_id = mrna_to_protein[mrna_id]
                bed.write(f"{chrom}\t{start}\t{end}\t{protein_id}\n")


def gorgonias(gff_file, bed_file, species_prefix):
    """
    Converts a GFF file to a BED file containing all mRNA features.

    Parameters:
        gff_file (str): Path to the input GFF file.
        bed_file (str): Path to the output BED file.
        species_prefix (str): Prefix to add to chromosome names in the output.

    Functionality:
        - Extracts mRNA features and their positions from the GFF file.
        - Writes the BED file with chromosome, start, end, and mRNA ID.
    """
    mrna_info = {}  # Stores mRNA information: mRNA_ID -> {'chrom', 'start', 'end'}

    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, feature_type, start, end, attributes = fields[0], fields[2], int(fields[3]) - 1, fields[4], fields[8]

            # Process attributes into a dictionary
            attr_dict = dict(attr.split('=') for attr in attributes.strip(';').split(';') if '=' in attr)

            if feature_type == 'mRNA':
                mrna_id = attr_dict.get('ID', '')
                if mrna_id:
                    mrna_info[mrna_id] = {'chrom': chrom, 'start': start, 'end': int(end)}

    with open(bed_file, 'w') as bed:
        for mrna_id, info in mrna_info.items():
            chrom = species_prefix + info['chrom']
            start = info['start']
            end = info['end']
            bed.write(f"{chrom}\t{start}\t{end}\t{mrna_id}\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 gff_to_bed.py <input.gff> <output.bed> <species_prefix>")
        sys.exit(1)

    # Call the desired function (e.g., gorgonias or gff_to_bed) based on your use case
    gorgonias(sys.argv[1], sys.argv[2], sys.argv[3])
