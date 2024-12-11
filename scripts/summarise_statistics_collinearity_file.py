#!/usr/bin/env python3

"""
Script Function:
This script processes collinearity and BED files to analyze syntenic relationships among genes. It generates two outputs:
1. A TSV file summarizing the count of syntenic genes for each chromosome, arranged from the maximum to the lowest count value
2. A TSV file detailing interactions between chromosome pairs, including syntenic block counts and total genes involved in each chromosome pair and arranged according from higher to lower according to the syntenic block counts column.

Input:
- A collinearity file describing syntenic alignments.
- A BED file mapping genes to chromosomes.

Output:
- A syntenic gene count per chromosome.
- Chromosomal interaction statistics.
"""

import os
from collections import defaultdict

def process_collinearity(file_path, bed_file, output_file, interaction_output_file):
    """
    Processes collinearity and BED files to extract syntenic relationships and chromosomal interactions.

    Parameters:
        file_path (str): Path to the collinearity file.
        bed_file (str): Path to the BED file mapping genes to chromosomes.
        output_file (str): Path for the syntenic gene count output.
        interaction_output_file (str): Path for the chromosomal interaction output.

    Functionality:
        - Reads the BED file to map genes to chromosomes.
        - Analyzes the collinearity file to identify syntenic blocks and gene interactions.
        - Outputs syntenic gene counts per chromosome.
        - Outputs chromosomal interaction statistics, including block counts and total genes involved.
    """
    # Map genes to their respective chromosomes using the BED file
    gene_to_chromosome = {}
    with open(bed_file, 'r') as bed:
        for line in bed:
            fields = line.strip().split('\t')
            gene = fields[1]
            chromosome = fields[0]
            gene_to_chromosome[gene] = chromosome

    # Initialize data structures for syntenic gene and interaction analysis
    syntenic_genes = set()
    chromosome_counts = defaultdict(int)
    interaction_counts = defaultdict(int)
    interaction_genes = defaultdict(set)

    with open(file_path, 'r') as collinearity:
        current_alignment = None
        current_pair = None

        for line in collinearity:
            # Skip non-relevant lines
            if (line.startswith("#") and "Alignment" not in line) or not line.strip():
                continue

            # Identify syntenic block information from alignment headers
            if line.startswith("#"):
                fields = line.strip().split()
                for field in fields:
                    if '&' in field:
                        pair = field
                        if field.split('&')[0][:2] != field.split('&')[1][:2]:
                            current_pair = pair
                            interaction_counts[current_pair] += 1
                            current_alignment = current_pair
                        break
                continue

            # Process gene pairs within syntenic blocks
            fields = line.strip().split("\t")
            if len(fields) < 3:
                continue

            gene1, gene2 = fields[1], fields[2]
            if gene1[:2] == gene2[:2]:
                continue

            # Count unique syntenic genes and update chromosome statistics
            if gene1 in gene_to_chromosome and gene1 not in syntenic_genes:
                chromosome_counts[gene_to_chromosome[gene1]] += 1
                syntenic_genes.add(gene1)

            if gene2 in gene_to_chromosome and gene2 not in syntenic_genes:
                chromosome_counts[gene_to_chromosome[gene2]] += 1
                syntenic_genes.add(gene2)

            # Track genes involved in chromosomal interactions
            if current_alignment:
                interaction_genes[current_alignment].update([gene1, gene2])

    # Write syntenic gene counts to the output file
    with open(output_file, 'w') as out:
        out.write("Chromosome\tSyntenic_Gene_Count\n")
        for chromosome, count in sorted(chromosome_counts.items(), key=lambda x: x[1], reverse=True):
            out.write(f"{chromosome}\t{count}\n")

    # Write chromosomal interaction statistics to the output file
    with open(interaction_output_file, 'w') as out_interaction:
        out_interaction.write("Chromosome_Pair\tSyntenic_Block_Count\tTotal_Genes_Involved\n")
        sorted_interactions = sorted(interaction_counts.items(), key=lambda x: x[1], reverse=True)
        for pair, block_count in sorted_interactions:
            total_genes = len(interaction_genes[pair])
            out_interaction.write(f"{pair}\t{block_count}\t{total_genes}\n")

    print(f"Gene count file generated: {output_file}")
    print(f"Interaction statistics file generated: {interaction_output_file}")

if __name__ == "__main__":
    """
    Main Execution:
        - Validates the existence of input files.
        - Processes the collinearity and BED files.
        - Outputs syntenic gene statistics and chromosomal interactions.
    """
    collinearity_path = "/media/crowfoot2/DATOS/CHALCO/eproceso/mcscanx_analysis_stringent/all_new.collinearity"
    bed_path = "/media/crowfoot2/DATOS/CHALCO/eproceso/mcscanx_analysis_stringent/all_new.gff"
    output_path = "/media/crowfoot2/DATOS/CHALCO/eproceso/mcscanx_analysis_stringent/syntenic_gene_counts.tsv"
    interaction_output_path = "/media/crowfoot2/DATOS/CHALCO/eproceso/mcscanx_analysis_stringent/chromosome_interactions.tsv"

    if not os.path.exists(collinearity_path) or not os.path.exists(bed_path):
        print("Error: Ensure the input files exist.")
    else:
        process_collinearity(collinearity_path, bed_path, output_path, interaction_output_path)
