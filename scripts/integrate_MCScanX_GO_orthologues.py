#!/usr/bin/env python3

import os
import requests
from collections import defaultdict

# Absolute paths to the files (update these paths as needed)
N0_TSV_PATH = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/N0.tsv"
EGGNOG_DIR = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/"
COLLINEARITY_PATH = "/media/crowfoot2/DATOS/CHALCO/eproceso/mcscanx_analysis_stringent/all_new.collinearity"
GO_OBO_URL = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
GFF_DIR = "/media/crowfoot2/DATOS/CHALCO/eproceso/gff/"

# Output files
OUTPUT_GO_PATH = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/syntenic_blocks_go_stringent.tsv"
OUTPUT_OG_PATH = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/syntenic_blocks_og_stringent.tsv"
OUTPUT_GO_FILTERED_PATH = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/syntenic_go_filtered_stringent.tsv"
OUTPUT_KEGG_PATH = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/syntenic_kegg_stringent.tsv"
# New output file for protein IDs and locus_tags
OUTPUT_PROTEIN_LOCUS_PATH = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/syntenic_blocks_protein_locus.tsv"
# New output files for interspecies blocks
OUTPUT_GO_FILTERED_INTERSPECIES_PATH = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/syntenic_go_filtered_interspecies_stringent.tsv"
OUTPUT_KEGG_INTERSPECIES_PATH = "/media/crowfoot2/DATOS/CHALCO/eproceso/eggnog/syntenic_kegg_interspecies_stringent.tsv"

def download_go_ontology(go_obo_url):
    """
    Downloads the GO ontology file and parses it into a dictionary.
    Also identifies obsolete terms and stores them.
    """
    print("Downloading GO ontology file...")
    response = requests.get(go_obo_url)
    if response.status_code != 200:
        raise Exception(f"Failed to download GO ontology file: {response.status_code}")
    go_obo_content = response.text

    # Parse the GO ontology content
    print("Parsing GO ontology file...")
    go_dict = {}
    obsolete_terms = set()
    in_term = False
    go_id = None
    name = None
    is_obsolete = False
    for line in go_obo_content.split('\n'):
        line = line.strip()
        if line == '[Term]':
            in_term = True
            go_id = None
            name = None
            is_obsolete = False
        elif line == '':
            if in_term and go_id and name:
                go_dict[go_id] = {'name': name, 'is_obsolete': is_obsolete}
                if is_obsolete:
                    obsolete_terms.add(go_id)
            in_term = False
        elif in_term:
            if line.startswith('id: GO:'):
                go_id = line[4:]  # Get everything after 'id: '
            elif line.startswith('name:'):
                name = line[6:]  # Get everything after 'name: '
            elif line.startswith('is_obsolete: true'):
                is_obsolete = True
    print(f"Total GO terms parsed: {len(go_dict)}")
    print(f"Total obsolete GO terms: {len(obsolete_terms)}")
    return go_dict, obsolete_terms

def download_kegg_pathways():
    """
    Downloads the list of KEGG pathways and parses it into a dictionary.
    """
    print("Downloading KEGG pathways list...")
    response = requests.get('http://rest.kegg.jp/list/pathway')
    if response.status_code != 200:
        raise Exception(f"Failed to download KEGG pathways list: {response.status_code}")
    data = response.text

    # Parse the data into a dictionary
    kegg_pathway_dict = {}
    for line in data.strip().split('\n'):
        entry, description = line.split('\t')
        pathway_id = entry.replace('path:', '')
        kegg_pathway_dict[pathway_id] = description
    print(f"Total KEGG pathways parsed: {len(kegg_pathway_dict)}")
    return kegg_pathway_dict
def parse_gff_files(gff_dir):
    """
    Parses the GFF files to build mappings from protein IDs to locus_tags for each species.
    Returns a dictionary where keys are species prefixes and values are dictionaries mapping protein IDs to locus_tags.
    """
    print("Parsing GFF files to build protein ID to locus_tag mappings...")
    species_gff_files = {
        'AT': 'arabidopsis.gff',
        'CU': 'cuscuta.gff',
        'CR': 'chlamydomonas.gff',
        'RO': 'roridula.gff'
    }
    species_protein_to_locus = {}
    for prefix, filename in species_gff_files.items():
        filepath = os.path.join(gff_dir, filename)
        protein_to_locus = {}
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                attributes = fields[8]
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                if 'protein_id' in attr_dict and 'locus_tag' in attr_dict:
                    protein_id = attr_dict['protein_id']
                    locus_tag = attr_dict['locus_tag']
                    protein_to_locus[protein_id] = locus_tag
        species_protein_to_locus[prefix] = protein_to_locus
        print(f"Species {prefix}: {len(protein_to_locus)} protein ID to locus_tag mappings")
    return species_protein_to_locus

def parse_n0_tsv(n0_tsv_path):
    """
    Parses the N0.tsv file to build a mapping from gene IDs to OGs, adding species prefixes to gene IDs.
    """
    print("Parsing N0.tsv to build gene to OG mapping...")
    gene_to_og = {}
    with open(n0_tsv_path, 'r') as f:
        header = f.readline().strip().split('\t')
        species_columns = header[3:]  # Skip first three columns
        species_prefixes = {}
        # Define prefixes based on column names
        for idx, col_name in enumerate(species_columns):
            if 'arabidopsis' in col_name.lower():
                species_prefixes[idx] = 'AT'
            elif 'cuscuta' in col_name.lower():
                species_prefixes[idx] = 'CU'
            elif 'clamydomonas' in col_name.lower() or 'chlamydomonas' in col_name.lower():
                species_prefixes[idx] = 'CR'
            elif 'gorgonias' in col_name.lower():
                species_prefixes[idx] = 'RO'
            else:
                species_prefixes[idx] = ''
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split('\t')
            og_id = fields[1]
            for idx, genes_str in enumerate(fields[3:]):
                if genes_str:
                    genes = [gene.strip() for gene in genes_str.split(',') if gene.strip()]
                    prefix = species_prefixes[idx]
                    for gene in genes:
                        gene_with_prefix = prefix + gene if prefix else gene
                        gene_to_og[gene_with_prefix] = og_id
    print(f"Total genes mapped to OGs: {len(gene_to_og)}")
    return gene_to_og

def parse_eggnog_files(eggnog_dir):
    """
    Parses the eggnog files to build mappings from gene IDs to GO terms and KEGG pathways,
    adding species prefixes to gene IDs.
    """
    print("Parsing eggnog files to build gene to GO and KEGG pathway mappings...")
    gene_to_go = defaultdict(set)
    gene_to_kegg_pathways = defaultdict(set)
    # Mapping from filenames to species prefixes
    eggnog_files = {
        "arabidopsis_eggnog.tsv": "AT",
        "cuscuta_eggnog.tsv": "CU",
        "chlamydomonas_eggnog.tsv": "CR",
        "gorgonias_eggnog.tsv": "RO"
    }
    for filename, prefix in eggnog_files.items():
        filepath = os.path.join(eggnog_dir, filename)
        with open(filepath, 'r') as f:
            headers = None
            for line in f:
                if line.startswith('#'):
                    headers = line.strip().split('\t')
                    query_idx = headers.index('#query')
                    gos_idx = headers.index('GOs')
                    kegg_pathway_idx = headers.index('KEGG_Pathway')  # Columna 13
                    continue
                fields = line.strip().split('\t')
                if len(fields) <= max(gos_idx, kegg_pathway_idx):
                    continue
                gene_id = fields[query_idx]
                gene_id_with_prefix = prefix + gene_id
                # Process GOs
                gos = fields[gos_idx]
                if gos:
                    go_terms = [go.strip() for go in gos.split(',') if go.strip()]
                    gene_to_go[gene_id_with_prefix].update(go_terms)
                # Process KEGG pathways
                kegg_pathways = fields[kegg_pathway_idx]
                if kegg_pathways:
                    pathways = [pw.strip() for pw in kegg_pathways.split(',') if pw.strip()]
                    # Filter to include only pathways starting with 'map'
                    pathways = [pw for pw in pathways if pw.startswith('map')]
                    gene_to_kegg_pathways[gene_id_with_prefix].update(pathways)
    print(f"Total genes mapped to GOs: {len(gene_to_go)}")
    print(f"Total genes mapped to KEGG pathways: {len(gene_to_kegg_pathways)}")
    return gene_to_go, gene_to_kegg_pathways

def parse_collinearity_file(collinearity_path):
    """
    Parses the collinearity file to extract syntenic blocks and their genes, including the comparison names.
    Also extracts species information for each block.
    Returns a dictionary of blocks with their genes, total gene counts, and species information.
    """
    print("Parsing collinearity file to get syntenic blocks...")
    blocks = {}
    block_species = {}
    with open(collinearity_path, 'r') as f:
        current_block = None
        block_gene_count = {}
        for line in f:
            line = line.strip()
            if line.startswith('## Alignment'):
                # Extract block name including comparison
                parts = line.split()
                block_id = parts[2].replace(":", "")  # e.g., Alignment 0:
                comparison = " ".join(parts[5:])  # e.g., AT1&AT1 minus
                full_block_name = f"Block{block_id} {comparison}"
                current_block = full_block_name
                blocks[current_block] = set()
                block_gene_count[current_block] = 0

                # Extract species from comparison
                # Assuming the comparison is like 'AT1&AT3 plus'
                chroms = parts[6]  # e.g., 'AT1&AT3'
                print(chroms)
                print(parts)
                species1 = chroms.split('&')[0][:2]  # Get first two letters as species code
                species2 = chroms.split('&')[1][:2]
                block_species[current_block] = (species1, species2)
            elif line.startswith('#') or line.startswith('##') or not line:
                continue
            else:
                # Parse gene pairs
                fields = line.split("\t")
                if len(fields) >= 4:
                    gene1 = fields[1]
                    gene2 = fields[2]
                    # Keep prefixes as is
                    blocks[current_block].update([gene1, gene2])
                    block_gene_count[current_block] += 2  # Counting both genes
    print(f"Total syntenic blocks found: {len(blocks)}")
    return blocks, block_gene_count, block_species

def process_blocks(blocks, block_gene_count, block_species, gene_to_go, gene_to_kegg_pathways, gene_to_og, go_dict, kegg_pathway_dict, obsolete_terms, species_protein_to_locus):
    """
    Processes each block, collects GO terms, KEGG pathways, and OGs, sorts the blocks and terms as required,
    and writes the outputs.
    Also generates additional files for interspecies blocks.
    """
    print("Processing blocks to collect GO terms, KEGG pathways, OGs, and descriptions...")
    # Sort blocks by total gene count in descending order
    sorted_blocks = sorted(blocks.items(), key=lambda x: len(x[1]), reverse=True)

    # Define general GO terms to exclude
    general_go_terms = {
        'GO:0003674',  # molecular_function
        'GO:0008150',  # biological_process
        'GO:0005575',  # cellular_component
    }

    with open(OUTPUT_GO_PATH, 'w') as go_out, \
         open(OUTPUT_OG_PATH, 'w') as og_out, \
         open(OUTPUT_GO_FILTERED_PATH, 'w') as go_filtered_out, \
         open(OUTPUT_KEGG_PATH, 'w') as kegg_out, \
         open(OUTPUT_GO_FILTERED_INTERSPECIES_PATH, 'w') as go_filtered_interspecies_out, \
         open(OUTPUT_KEGG_INTERSPECIES_PATH, 'w') as kegg_interspecies_out, \
         open(OUTPUT_PROTEIN_LOCUS_PATH, 'w') as protein_locus_out:
        # Write headers
        go_out.write("Block\tGO_Term\tGO_Description\tGO_Count\tGenes\n")
        og_out.write("Block\tOG\tOG_Count\tGenes\n")
        go_filtered_out.write("Block\tGO_Term\tGO_Description\tGO_Count\tGenes\n")
        kegg_out.write("Block\tKEGG_Pathway_ID\tKEGG_Pathway_Description\tPathway_Count\tGenes\n")
        go_filtered_interspecies_out.write("Block\tGO_Term\tGO_Description\tGO_Count\tGenes\n")
        kegg_interspecies_out.write("Block\tKEGG_Pathway_ID\tKEGG_Pathway_Description\tPathway_Count\tGenes\n")
        # Header for the new output file
        protein_locus_out.write("Block\tGene_ID\tProtein_ID\tLocus_Tag\n")

        for block_id, genes in sorted_blocks:
            # Determine if block is interspecies
            species1, species2 = block_species[block_id]
            is_interspecies = species1 != species2

            # Collect GO terms
            go_counts = defaultdict(list)
            # Collect KEGG pathways
            kegg_counts = defaultdict(list)
            og_counts = defaultdict(list)

            # For the new output
            processed_genes = set()

            for gene in genes:
                # Get species prefix (first two letters)
                species_prefix = gene[:2]
                # Remove the prefix to get the protein_id
                protein_id = gene[2:]
                # Get the species mapping
                protein_to_locus = species_protein_to_locus.get(species_prefix)
                if protein_to_locus:
                    locus_tag = protein_to_locus.get(protein_id, 'NA')
                else:
                    locus_tag = 'NA'
                # Write to the new output file
                protein_locus_out.write(f"{block_id}\t{gene}\t{protein_id}\t{locus_tag}\n")

                # Avoid processing the same gene multiple times
                if gene in processed_genes:
                    continue
                processed_genes.add(gene)

                # Get GO terms
                if gene in gene_to_go:
                    for go_term in gene_to_go[gene]:
                        go_counts[go_term].append(gene)
                # Get KEGG pathways
                if gene in gene_to_kegg_pathways:
                    for kegg_pathway in gene_to_kegg_pathways[gene]:
                        kegg_counts[kegg_pathway].append(gene)
                # Get OG
                if gene in gene_to_og:
                    og_id = gene_to_og[gene]
                    og_counts[og_id].append(gene)

            # Sort GO terms by count in descending order
            sorted_go_counts = sorted(go_counts.items(), key=lambda x: len(x[1]), reverse=True)
            # Sort KEGG pathways by count in descending order
            sorted_kegg_counts = sorted(kegg_counts.items(), key=lambda x: len(x[1]), reverse=True)
            # Write GO terms with descriptions
            for go_term, gene_list in sorted_go_counts:
                go_info = go_dict.get(go_term)
                if go_info:
                    go_description = go_info['name']
                    is_obsolete = go_info['is_obsolete']
                else:
                    go_description = "Description not found"
                    is_obsolete = False
                # Write to the main GO output file
                go_out.write(f"{block_id}\t{go_term}\t{go_description}\t{len(gene_list)}\t{','.join(gene_list)}\n")
                # Check if the term should be included in the filtered file
                if (
                    not is_obsolete and
                    go_term not in general_go_terms and
                    'obsolete' not in go_description.lower() and
                    go_description != "Description not found"
                ):
                    # Write to the filtered GO output file
                    go_filtered_out.write(f"{block_id}\t{go_term}\t{go_description}\t{len(gene_list)}\t{','.join(gene_list)}\n")
                    if is_interspecies:
                        # Write to the interspecies GO filtered output file
                        go_filtered_interspecies_out.write(f"{block_id}\t{go_term}\t{go_description}\t{len(gene_list)}\t{','.join(gene_list)}\n")
            # Write KEGG pathways with descriptions
            for kegg_id, gene_list in sorted_kegg_counts:
                kegg_description = kegg_pathway_dict.get(kegg_id, "Description not found")
                kegg_out.write(f"{block_id}\t{kegg_id}\t{kegg_description}\t{len(gene_list)}\t{','.join(gene_list)}\n")
                if is_interspecies:
                    # Write to the interspecies KEGG output file
                    kegg_interspecies_out.write(f"{block_id}\t{kegg_id}\t{kegg_description}\t{len(gene_list)}\t{','.join(gene_list)}\n")
            # Write OGs
            for og_id, gene_list in og_counts.items():
                og_out.write(f"{block_id}\t{og_id}\t{len(gene_list)}\t{','.join(gene_list)}\n")

        print("Processing complete.")
        print(f"GO term output written to: {OUTPUT_GO_PATH}")
        print(f"Filtered GO term output written to: {OUTPUT_GO_FILTERED_PATH}")
        print(f"KEGG pathway output written to: {OUTPUT_KEGG_PATH}")
        print(f"OG output written to: {OUTPUT_OG_PATH}")
        print(f"Interspecies filtered GO term output written to: {OUTPUT_GO_FILTERED_INTERSPECIES_PATH}")
        print(f"Interspecies KEGG pathway output written to: {OUTPUT_KEGG_INTERSPECIES_PATH}")
        print(f"Protein ID to locus_tag output written to: {OUTPUT_PROTEIN_LOCUS_PATH}")

def main():
    # Step 1: Download and parse the GO ontology to build go_dict
    go_dict, obsolete_terms = download_go_ontology(GO_OBO_URL)

    # Step 1b: Download and parse KEGG pathways to build kegg_pathway_dict
    kegg_pathway_dict = download_kegg_pathways()

    # Step 2: Parse N0.tsv to build gene_to_og mapping
    gene_to_og = parse_n0_tsv(N0_TSV_PATH)

    # Step 3: Parse eggnog files to build gene_to_go and gene_to_kegg_pathways mappings
    gene_to_go, gene_to_kegg_pathways = parse_eggnog_files(EGGNOG_DIR)

    # Step 4: Parse collinearity file to get syntenic blocks and gene counts
    blocks, block_gene_count, block_species = parse_collinearity_file(COLLINEARITY_PATH)

    # Step 5: Parse GFF files to build protein_id to locus_tag mappings
    species_protein_to_locus = parse_gff_files(GFF_DIR)

    # Step 6: Process blocks and write outputs
    process_blocks(
        blocks, block_gene_count, block_species,
        gene_to_go, gene_to_kegg_pathways, gene_to_og,
        go_dict, kegg_pathway_dict, obsolete_terms,
        species_protein_to_locus
    )

if __name__ == "__main__":
    main()