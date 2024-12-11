#!/usr/bin/env python3

"""
Script Function:
This script processes EggNOG annotation files and a GFF file to analyze gene-to-GO term associations,
summarize gene counts per GO term, and report chromosome-level distributions of specific GO terms.

Outputs:
1. A list of genes associated with specific GO terms.
2. A summary of gene counts for each GO term across species.
3. A detailed summary of GO term distributions across chromosomes.
"""

import os
from collections import defaultdict

# Define specific GO terms of interest
specific_gos = {
    "GO:0019253", "GO:0010110", "GO:0019685", "GO:0015979", "GO:0006796",
    "GO:0040011", "GO:1902019", "GO:0009399", "GO:1902025"
}

# Species prefixes mapped to their EggNOG filenames
species_prefix = {
    "arabidopsis_eggnog.tsv": "AT",
    "chlamydomonas_eggnog.tsv": "CR",
    "cuscuta_eggnog.tsv": "CU",
    "gorgonias_eggnog.tsv": "RO",
}

# File paths and directories
input_dir = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog"
gff_file = "/media/crowfoot2/DATOS/CHALCO/eproceso/mcscanx_analysis_relaxaed_names/all_new.gff"
# Output files
output_file = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/go_genes.tsv"
summary_file = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/go_summary.tsv"
chromosome_file = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/go_chromosome_summary.tsv"

# Data structures for storing results
go_to_genes = defaultdict(lambda: defaultdict(list))
go_summary = defaultdict(lambda: defaultdict(int))
chromosome_summary = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
gene_to_chromosome = {}

# Map genes to chromosomes from the GFF file
with open(gff_file, "r") as gff:
    for line in gff:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) > 3:
            gene_name = fields[1]
            if gene_name:
                gene_to_chromosome[gene_name] = fields[0]  # Chromosome is in column 1

# Process EggNOG annotation files
for filename in os.listdir(input_dir):
    if filename.endswith("_eggnog.tsv"):
        prefix = species_prefix.get(filename, "")
        filepath = os.path.join(input_dir, filename)

        with open(filepath, "r") as file:
            for line in file:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) > 9:  # Column 10 contains GO terms
                    genes = fields[0]
                    if fields[9] != '-':
                        go_terms = fields[9].split(",")  # GO terms are comma-separated
                        for go in go_terms:
                            if go in specific_gos:
                                gene_id = f"{prefix}{genes}"
                                go_to_genes[go][prefix].append(gene_id)
                                go_summary[go][prefix] += 1
                                # Map genes to chromosomes for GO terms
                                chrom = gene_to_chromosome.get(gene_id, "Unknown")
                                chromosome_summary[go][prefix][chrom] += 1

# Write genes associated with GO terms to file
with open(output_file, "w") as out:
    for go, genes_dict in go_to_genes.items():
        line = [go]
        for prefix, genes in genes_dict.items():
            line.extend(genes)
        out.write("\t".join(line) + "\n")

# Write summary of gene counts per GO term
with open(summary_file, "w") as summary:
    summary.write("GO\tarabidopsis\tchlamydomonas\tcuscuta\troridula\tTotal\n")
    for go, counts_dict in go_summary.items():
        ar = counts_dict.get("AT", 0)
        cr = counts_dict.get("CR", 0)
        cu = counts_dict.get("CU", 0)
        ro = counts_dict.get("RO", 0)
        total = ar + cr + cu + ro
        summary.write(f"{go}\t{ar}\t{cr}\t{cu}\t{ro}\t{total}\n")

# Write chromosome-level summary of GO term distributions
with open(chromosome_file, "w") as chrom_file:
    chrom_file.write("GO\tSpecies\tChromosome\tGene_Count\n")
    for go, species_dict in chromosome_summary.items():
        for prefix, chrom_dict in species_dict.items():
            for chrom, count in chrom_dict.items():
                chrom_file.write(f"{go}\t{prefix}\t{chrom}\t{count}\n")

print(f"Files generated: {output_file}, {summary_file}, {chromosome_file}")
