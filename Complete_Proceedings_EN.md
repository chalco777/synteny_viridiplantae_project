# Analysis of the synteny and conservation of metabolic pathways in photosynthetic organisms: from unicellular algae to specialized plants

**Authors:**

- Armas Torres Valery Leonor  
- Cahuana Mamani, Nataly Caterine  
- Chalco González Adrián Alexandro  
- Alvarao Caballero Johan Isaac

**Teachers:**

- Claudia Machicado  
- Felipe One

LIMA - PERU

2024

**Index:**

- [Premise](#premise)
- [1. Introduction](#1-introduction)
- [2. Goals](#2-goals)
- [3. Materials and methods](#3-materials-and-methods)
  - [3.1. Sequence selection parameters](#31-sequence-selection-parameters)
  - [3.2. Evaluation of continuity and integrity of assemblies](#32-evaluation-of-continuity-and-integrity-of-assemblies)
  - [3.3. Assign functions to genes or protein sequences](#33-assign-functions-to-genes-or-protein-sequences)
  - [3.4. Identification of conserved orthologous gene groups: OrthoFinder](#34-identification-of-conserved-orthologous-gene-groups-orthofinder)
  - [3.5. Preparing GFF Annotation Files for MCScanX](#35-preparing-gff-annotation-files-for-mcscanx)
  - [3.6. Preparing Fasta Protein Sequences for MCScanX in Bash](#36-preparing-fasta-protein-sequences-for-mcscanx-in-bash)
  - [3.7. Execution in MCScanX](#37-execution-in-mcscanx)
  - [3.8. MCScanX data processing: Obtaining statistics of syntenic blocks and analysis of differential and enriched GO terms](#38-mcscanx-data-processing-obtaining-statistics-of-syntenic-blocks-and-analysis-of-differential-and-enriched-go-terms)
- [4. Results](#4-results)
- [5. Discussion](#5-discussion)
- [6. Conclusions](#6-conclusions)
- [7. Literature](#7-literature)
- [8. Attachments](#8-attachments)

# Premise

Throughout evolution, chlorophyll-based photosynthetic organisms have colonized different environments and adapted to different habitat conditions, ranging from aquatic algae to parasitic and carnivorous plants. The aim is to identify synteny and its rearrangements among representative species of algae, photosynthetic, parasitic, and carnivorous plants, and thus identify within the conserved blocks which gene clusters related to metabolic pathways remain conserved over time and evolution. It is suggested that a representative species from each group, with a small genome, be selected and that the analysis of gene clusters be limited to functional groups to encompass processes conserved over time.

## 1. Introduction

The evolutionary history of terrestrial plants was born in the sea, within the group of green algae charophytes (Jianchao Ma, 2022). This evolutionary process involved not only the gain and co-option of genes, but also the loss of others. In addition, there are studies on gene families among photosynthetic species, such as *Arabidopsis thaliana* and *Chlamydomonas reinhardtii*. However, until now it has not been compared with *Cuscuta australis* (parasitic plant) and *Roridula gorgonias* (carnivorous plant). Gene families conserved in terrestrial plants are often related to stress response, ion transport, and metabolites (Jianchao Ma, 2022).

Therefore, in this work we will evaluate the shared syntenic blocks in these four species: three terrestrial plants (*Arabidopsis thaliana*, *Cuscuta australis* and *Roridula gorgonias*) and a unicellular alga (*Chlamydomonas reinhardtii*).

The objective is to determine the synteny of conserved metabolic pathways among these photosynthetic species. It is expected to identify the conservation of some key metabolic pathways among the analyzed organisms, such as: ion transport (such as K, Fe, Cu, etc.), cell cycle regulation, abiotic stress response mechanisms (particularly relevant in terrestrial plants), photosynthesis (possibly partially in *Cuscuta australis*), metabolism of secondary compounds for biotic adaptation and embryonic development (Guo, 2012; Sun et al., 2018; Jianchao Ma, 2022; Fleck et al., 2023).

To achieve these objectives, a series of bioinformatics tools are used. Genomic data were downloaded from public databases such as NCBI and Dryad. The quality of the genomes will be assessed using the tools **QUAST** and **BUSCO**, available in Galaxy. Programming languages such as Python, Bash and R will be used for file editing, data processing, and analysis. Annotation and identification of conserved orthologous genes will be performed using **eggNOG-mapper** and **OrthoFinder**, respectively. R will be used to graphically represent the results, complemented by Excel in specific cases, such as the generation of bar charts. Syntenic blocks will be determined using **MCScanX**, and their visualization will be managed with **SynVisio**. Finally, the enrichment of GO terms in the syntenic blocks will be analyzed using the **DAVID** and **PANTHER** platforms.

## 2. Goals

**Main objective:** To determine the synteny of conserved metabolic pathways among photosynthetic species: a unicellular alga, a model plant, a parasitic plant, and a carnivorous plant.

**Secondary objectives:**

- Genome exclusion criteria  
- Identify conserved orthologous genes among the selected species, using tools such as **eggNOG-mapper** and **OrthoFinder**.  
- Analyze the synteny of orthologous genes of interest in specific metabolic pathways using **MCScanX** and **SynVisio**.

## 3. Materials and methods

**Important note:** *From this point on, this work mentions several scripts in Python or R, of which only screenshots are taken as a general sample. The reader can access each of these scripts by name in the "Scripts" folder.*

Guide flowchart  
The numbers correspond to the section where the table is discussed. QUAST and BUSCO are not considered.  
![][image2]

### 3.1. Sequence selection parameters

For the sequence selection process, the following criteria were considered: photosynthetic organism, availability of the genome, proteins and annotations in databases such as NCBI or Dryad, and genome size between 100 Mb to 300 Mb (small), for efficient computational management given our limited resources.

The NCBI search was filtered using known photosynthetic organisms. Priority was given to those with assembled genomes, accompanied by genome sequence files (FASTA), protein files (FASTA), and complete annotations (GFF). Organisms that met the selection criteria were eligible for downloading their sequences in the following formats:

- Genome and proteins in format **FASTA**  
- Annotations in format **GFF**

Table 1. General display of the selected organisms

|   | *Arabidopsis thaliana* | *Chlamydomonas reinhardtii* | *Cuscuta australis* | *Utricularia gibba* |
|---|------------------------|-----------------------------|---------------------|---------------------|
| Identifier | GCF_0000001735.4 | GCF_000002595.2 | GCA_003260385.1 | GCA_002189035.1 |
| Source | RefSeq | RefSeq | GenBank | GenBank |
| Downloaded documents | Genome sequences (FASTA) → fna; Annotation Functions → GFF; Proteins (FASTA)→ faa | Genome sequences (FASTA) → fna; Annotation Functions → GFF; Proteins (FASTA)→ faa | Genome sequences (FASTA) → fna; Annotation Functions → GFF; Proteins (FASTA)→ faa | Genome sequences (FASTA) → fna; Sequence and annotations → GBFF; There are no proteins |

Initially, it was considered to use the carnivorous plant sequence *Utricularia gibba*. However, as it lacked protein sequence and annotations, it was decided to use *Roridula gorgonias*. Although annotated genomes were searched for *Utricularia gibba* in multiple sources and supplementary file pages of articles, the necessary genome annotation was not found.  
So, Table 1 with the genomes used looks like this:

|   | *Arabidopsis thaliana* | *Chlamydomonas reinhardtii* | *Cuscuta australis* | *Roridula gorgonias* |
|---|---------------------|---------------------------|-------------------|--------------------|
| Identifier | GCF_0000001735.4 | GCF_000002595.2 | GCA_003260385.1 | _______________ |
| Source | RefSeq | RefSeq | GenBank | Dryad |
| Downloaded documents | Genome sequences (FASTA) → fna; Annotation Functions → GFF; Proteins (FASTA)→ faa | Genome sequences (FASTA) → fna; Annotation Functions → GFF; Proteins (FASTA)→ faa | Genome sequences (FASTA) → fna; Annotation Functions → GFF; Proteins (FASTA)→ faa | Genome sequences (FASTA) → fna; Annotation Functions → GFF; Proteins (FASTA)→ faa |

### 3.2. Evaluation of continuity and integrity of assemblies

The quality of the assemblies was assessed using the **QUAST** and **BUSCO** tools within Galaxy Europe.

**QUAST**

1. Selection of input sequences in fna (fasta) format - Input: assembled genomes.  
2. Selecting the type of organisms → choose organisms **Eukaryotes**  
3. Selecting output files HTML and PDF report.  
4. The other parameters are by default.

**Assembly Integrity Analysis Workflow with BUSCO**

Protein analysis in BUSCO. It was divided into two operations: the aquatic and terrestrial species lineages, the Chlorophyta lineage (*Chlamydomonas reinhardtii*) and Eudicots (*Arabidopsis thaliana*, *Cuscuta australis*, *Roridula gorgonias*). The protein sequences were uploaded, in the first operation this configuration was used:

- Lineage data source: Download lineage data  
- Mode: Annotated gene sets (protein)  
- Auto-detect or select lineage: Select lineage  
- Lineage: Chlorophyta

and in the following one this configuration was used:

- Lineage data source: Download lineage data  
- Mode: Annotated gene sets (protein)  
- Auto-detect or select lineage: Select lineage  
- Lineage: Eudicots

It returns several results, a table is highlighted to see the fragmented, complete, and missing genes, and a bar graph is compiled in Excel.

### 3.3. Assign functions to genes or protein sequences

To assign functions to the protein sequences of each species, the tool **EggNOG Mapper** (Functional sequence annotation by orthology) contained in Galaxy Europe was used. Protein files were uploaded here in faa format.

1. Working parameters  
   1. Input: The protein sequence in faa format.  
   2. Sequence type: proteins  
   3. Scoring and gap cost matrix: BLOSUM62  
   4. Base for scoring: Diamond  
   5. Diamond Sensitivity Mode: Sensitive  
   6. The other parameters are by default.

2. Output files  
   1. seed_orthologs (not used)  
   2. Annotations: The GO term identifiers and KEGG pathways mentioned below were extracted from this file.

### 3.4. Identification of conserved orthologous gene groups: OrthoFinder

**OrthoFinder** is a tool for identifying orthologous (and paralogous) genes between different genomes. Its GitHub is at [https://github.com/davidemms/OrthoFinder](https://github.com/davidemms/OrthoFinder)

**Installing OrthoFinder:** Clone the repository from GitHub and install the dependencies.

![][image3]  
![][image4]  
![][image5]  
![][image6]

**Data preparation:**

- Gather protein sequence files in format **FASTA** for each genome.  
- Place the files in a folder and make sure they are named clearly.

![][image7]

**Running OrthoFinder:** Run OrthoFinder pointing to the directory where you have the FASTA files.

![][image8]

**Review of the results:** OrthoFinder creates a results folder (for example, `Results_date`) containing key files:

- **Orthogroups.csv:** Lists the orthogroups with the genes of each genome.  
- **Orthologs:** Relationship of orthologs between genomes.  
- **Gene_Trees:** Phylogenetic trees of genes.  
- **Viewing results:** Using the script **figures.Rmd** in R for the results placed in the figures X

![][image9]

### 3.5. Preparing GFF Annotation Files for MCScanX

(**to get all_new.gff**)

To run MCScanX, two files are required: first, the tabular result of performing an all vs all blast on the protein sequences (all against all). Second, a single gff file (which is really a bed) corresponding to all these sequences. Each protein sequence corresponds to a single mRNA in the genome and gff file of each species. It should be noted that a gene can have multiple mRNAs, and therefore, multiple protein sequences. We have worked from mRNAs.

1. First, the gff_to_bed script below was created and executed (the script name has been changed, it was previously edit_gff), which extracts only the mRNAs from the gff and puts them into my new bed file. It also adds a species prefix (AT, CR, CU or RO) to the identifier of each chromosome in the chromosomes column. This is so that it can be displayed correctly when performing synteny analysis, and to identify which species each chromosome comes from. The script below consists of two functions (each one starts with def), one for *Roridula gorgonias*, and the other three for the rest of the species. The reason is that the structure of the gff of *Roridula* was different from those of the other three species because it didn't come from the NCBI Refseq. The script has three inputs: the original input or gff, the output in the bed file, and finally, the prefix, which is two capital letters that vary between each species.

![][image10]

```bash
# Executing the script in the terminal  
python3 gff_to_bed.py arabidopsis.gff arabidopsis.bed AT  
python3 gff_to_bed.py chlamydomonas.gff chlamydomonas.bed CR  
python3 gff_to_bed.py cuscuta.gff cuscuta.bed CU  
python3 gff_to_bed.py roridula.gff roridula.bed RG
```
2. Using the Linux commands shown, it was confirmed that the number of features in each bed of each species was equal to the number of proteins in the fasta protein file of each species.

![][image11]  
![][image12]

3. All the beds of the species were concatenated (or joined) into a single bed called all.bed  
4. In order to improve gene distinction by understanding the future output of MCScanX and to better manage SynVisio, the following bash script was created to edit the all.bed file. What was done was to add the species prefixes AT, CR, CU and RO to the identifiers of each protein sequence in the BED, so that also the proteins, for example in *Arabidopsis*, have AT at the beginning, and not just the chromosomes.

```bash
# Creating the script  
cat << 'EOF' > add_prefixes_bed.awk  
BEGIN { OFS = "\t" }  
{ prefix = substr($1, 1, 2); $4 = prefix $4; print $0 }  
EOF  

# Script execution  
awk -f add_prefixes_bed.awk all.bed > all_prefixed.bed  
# After adding the prefixes to the genes, I reorder the columns  
awk '{print $1, $4, $2, $3}' OFS='\t' all_prefixed.bed > all.gff
```

5. This would terminate the process, but when running MCScanX, it's better to shorten the name of each chromosome, for example, from ATNC0009.1 to AT1. To do this, I used the Python script simplify_chromosome_name_gff.py to rename all the chromosomes in all.gff and simplify them.

![][image13]

**EXTRA STEP:**

- Chromosome names in the whole-genome fasta of each species were also simplified with simplify_chromosome_name_fasta.py

![][image14]

- Afterwards, the script get_length_chromosomes_fasta.py was used to calculate the lengths of each chromosome for each species and visualize the largest genomes with SynVisio (since most are small chromosomes). He also saw where the largest synteny blocks might be.

![][image15]

- The script output looks like this, where the chromosomes are ordered from largest to smallest:

![][image16]

### 3.6. Preparing Fasta Protein Sequences for MCScanX in Bash

1. First the other words were removed from the header and I kept only the first word (header accession)

```bash
sed 's/^\(\>[^ ]*\).*/\1/' arabidopsis_protein.faa > arabidopsis_protein_edited.faa  
sed 's/^\(\>[^ ]*\).*/\1/' chlamydomonas_protein.faa > chlamydomonas_protein_edited.faa  
sed 's/^\(\>[^ ]*\).*/\1/' cuscuta_protein.faa > cuscuta_protein_edited.faa  
sed 's/^\(\>[^ ]*\).*/\1/' roridula_protein.faa > roridula_protein_edited.faa
```

2. Added the corresponding prefix to the accession of each header (of each protein sequence) of the fasta protein file

```bash
sed -i 's/^>\(.*\)/>AT\1/' arabidopsis_protein_edited.faa  
sed -i 's/^>\(.*\)/>CR\1/' chlamydomonas_protein_edited.faa  
sed -i 's/^>\(.*\)/>CU\1/' cuscuta_protein_edited.faa  
sed -i 's/^>\(.*\)/>RO\1/' roridula_protein_edited.faa
```

Example of the change in the Arabidopsis thaliana header, before and after the execution of the previous commands (1 and 2).  
![][image17]![][image18]

3. Next, the database for each species was created to perform the multiple intra- and interspecies BLASTP (as mentioned, MCScanX requires the result of a BLAST as input):

```bash
makeblastdb -in arabidopsis_protein_edited.faa -dbtype prot -out arabidopsis_db  
makeblastdb -in chlamydomonas_protein_edited.faa -dbtype prot -out chlamydomonas_db  
makeblastdb -in cuscuta_protein_edited.faa -dbtype prot -out cuscuta_db  
makeblastdb -in roridula_protein_edited.faa -dbtype prot -out roridula_db
```

4. Performing intra- and interspecific BLASTP (between sequences of the same species and between sequences of different species)

In the commands below, the -outfmt 6 parameter indicates that the output will be in tabular format, which is required by MCScanX. This format contains predefined columns such as identities, e-values, and start/stop positions of the alignments. The -max_target_seqs 5 parameter limits the maximum number of target sequences (hits) reported per query, showing only the top 5 matches based on the score. This optimizes the output by focusing on the most relevant hits and is also recommended when running MCScanX.

```bash
# Arabidopsis vs Arabidopsis  
blastp -query arabidopsis_protein_edited.faa -db arabidopsis_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out arabidopsis_vs_arabidopsis.blast  
# Arabidopsis vs Chlamydomonas  
blastp -query arabidopsis_protein_edited.faa -db chlamydomonas_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out arabidopsis_vs_chlamydomonas.blast  
# Arabidopsis vs Cuscuta  
blastp -query arabidopsis_protein_edited.faa -db cuscuta_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out arabidopsis_vs_cuscuta.blast  
# Arabidopsis vs Roridula  
blastp -query arabidopsis_protein_edited.faa -db roridula_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out arabidopsis_vs_roridula.blast  
# Chlamydomonas vs Arabidopsis  
blastp -query chlamydomonas_protein_edited.faa -db arabidopsis_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out chlamydomonas_vs_arabidopsis.blast  
# Chlamydomonas vs Chlamydomonas  
blastp -query chlamydomonas_protein_edited.faa -db chlamydomonas_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out chlamydomonas_vs_chlamydomonas.blast  
# Chlamydomonas vs Cuscuta  
blastp -query chlamydomonas_protein_edited.faa -db cuscuta_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out chlamydomonas_vs_cuscuta.blast  
# Chlamydomonas vs Roridula  
blastp -query chlamydomonas_protein_edited.faa -db roridula_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out chlamydomonas_vs_roridula.blast  
# Cuscuta vs Arabidopsis  
blastp -query cuscuta_protein_edited.faa -db arabidopsis_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out cuscuta_vs_arabidopsis.blast  
# Cuscuta vs Chlamydomonas  
blastp -query cuscuta_protein_edited.faa -db chlamydomonas_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out cuscuta_vs_chlamydomonas.blast  
# Cuscuta vs Cuscuta  
blastp -query cuscuta_protein_edited.faa -db cuscuta_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out cuscuta_vs_cuscuta.blast  
# Cuscuta vs Roridula  
blastp -query cuscuta_protein_edited.faa -db roridula_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out cuscuta_vs_roridula.blast  
# Roridula vs Arabidopsis  
blastp -query roridula_protein_edited.faa -db arabidopsis_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out roridula_vs_arabidopsis.blast  
# Roridula vs Chlamydomonas  
blastp -query roridula_protein_edited.faa -db chlamydomonas_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out roridula_vs_chlamydomonas.blast  
# Roridula vs Cuscuta  
blastp -query roridula_protein_edited.faa -db cuscuta_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out roridula_vs_cuscuta.blast  
# Roridula vs Roridula  
blastp -query roridula_protein_edited.faa -db roridula_db -evalue 1e-5 -outfmt 6 -num_threads 28 -max_target_seqs 5 -out roridula_vs_roridula.blast
```

5. The results of all BLASTPs were concatenated.

```bash
cat *.blast > all_vs_all.blast
```

### 3.7. Execution in MCScanX

The MCSCanX program first detects homologous genes with the given blastp as input, and then detects homologous genes that are close together, within a specified range, to construct collinear syntenic blocks.

Input: all.blast and all.gff

1. First it was run with these default parameters

```bash
./MCScanX /media/crowfoot2/DATOS/CHALCO/eproceso/mcscanx_analysis/all -s 5 -m 25 -w 5
```

2. Then it was run with the parameters with more relaxed values, to allow obtaining more collinear syntenic blocks between the species, because they are distant.

```bash
./MCScanX /media/crowfoot2/DATOS/CHALCO/eproceso/mcscanx_analysis_relaxed/all_new -s 3 -m 50 -w 10
```

**Explanation of the parameters**  
-w OVERLAP_WINDOW: Maximum distance (in number of genes) between matching BLAST genes to be collapsed as part of the same block. In the “second” relaxed run (step 2), this parameter was increased to allow generating syntenic blocks even between relatively distant BLAST homologous genes (up to 10 genes distant).  
-m MAX_GAPS: Maximum number of gaps (missing or non-homologous genes) allowed within a syntenic block. In the relaxed run in step 2, this parameter was increased to 50 to allow the syntenic block to be maintained even if there are up to 50 non-homologous genes.  
-s MATCH_SIZE: Minimum number of genes required to consider a block syntenic. The value of this parameter was reduced to 3 in the second run, to allow calling smaller blocks (up to a minimum of 3 genes).

It is important to highlight that although both forms were worked with, both with stricter parameters by default (step 1), and with relaxed ones (step 2), after the inspection of the data carried out in “MCScanX data processing”, it was decided to use the strict parameters by default, because consistent syntenic blocks with similar functions (of GO terms) were obtained.

![][image19]

**Use:** *Here, the use of two steps does not illustrate that the program was run only twice. MCScanX was actually run multiple times, as the gff input and blast input were fine-tuned. There were also similar issues to those mentioned in the following "[Issue](https://github.com/wyp1125/MCScanX/issues/53)" on its GitHub page (with the order of the incoming bed or gff columns):*

*Overall, the MCScanX input documentation needs improvement.*

3. MCScanX output visualization, all. collinearity file:

![][image20]  
Note that each line “\#\# Alignment X” defines a syntenic block. The first blocks are between the same species (note the prefixes of the genes compared in both columns, both are AT, so the genes compared in these syntenic blocks are from *Arabidopsis*)

* This all. collinearity file, along with the gff accepted by MCScanX, served as input for SynVisio. The SynVisio web address is: [https://SynVisio. github. io/\#/](https://SynVisio. github. io/#/). 

![][image21]

* The chromosomes were used to visualize:

![][image22]

* And thus figures 6, 11 and 12 were produced.

### 3.8. MCScanX data processing: Obtaining statistics of syntenic blocks and analysis of differential and enriched GO terms

1. First, we sought to obtain general statistics on syntenic blocks. Since we are interested in blocks between species (not within a single species), we obtained the number of genes participating in interspecies blocks on each chromosome to determine which chromosomes contribute most to interspecies conservation. We also obtained the number of interspecies chromosome pairs with the largest number of blocks. syntenicand genes. All of the above was obtained by processing the results of thefile all. collinearity with the script en python “summarise_statistics_collinearity_file. py”:

![][image23]

2. Number of genes participating in interspecies collinearity blocks on each chromosome, ordered from chromosomes with the most to the fewest genes:

![][image24]

3. Number of chromosome pairs, the count of syntenic blocks between each pair, as well as the total number of genes in each, ordered from pairs with the most syntenic blocks to the fewest syntenic blocks:

![][image25]

4. After obtaining general statistics for the syntenic blocks, we considered how to evaluate the enrichment of GO terms among genes from different species and analyze differential and conserved GO terms. We first analyzed **differentialy** represented GO terms, computing their count across different species. To do this, we used the Python script “count_go_per_eggnogfile. py.” This script analyzes the annotations output of each species eggNOG file to (i) calculate the number of times a given GO appears in the annotated file for each species, (ii) provide all the genes associated with each GO, and (iii) generate a tsv file with the count of each GO on each chromosome for each species.

![][image26]

5. The main result was a summary of the number of occurrences of each GO in each species. This was performed with GOs associated with photosynthesis, nitrate import, and cilia, all using the R script. ** figures. Rmd **(see STEP 5 in the script) and the following was obtained:

![][image27]

From this table Table 3 and Figure 7 were produced.

6. The next step was to analyze the GO terms** preserved**between species, which is the main objective of the project. To this end, the largest Python script ever created, with multiple functions, called "integrate_MCScanX_GO_orthologs. py" was generated. This script takes as input the Eggnog annotation results for each species, the names associated with each GO term and KEGG pathway from the APIs of each of these databases, the ortholog groups detected by OrthoFinder and the collinearity of MCScanX, the original gff annotation files for each species, and generates what is described in the following steps.

![][image28]

7. The probably most important result of the “integrate_MCScanX_GO_orthologs. py” script is the tsv file with the count of GO terms present in each syntenic block, which basically allows me to know which GOs are enriched in that syntenic block conserved between species, how many proteins correspond to that GO and what IDs of those proteins:

![][image29]

8. The same was obtained for the KEGG metabolic pathways, although in the end they did not contribute as much to the analysis (GO terms were enough).

![][image30]

9. Another relevant result was obtaining the locus_tag for each protein in each syntenic block. In the image below, for example, we see for Block431 (as the reader will notice, this is an intraspecies syntenic block) each of the proteins in *A. thaliana *corresponding in the second column, and their locus_tag in the fourth column (the third column is the same as the second, it refers to the same protein, just without the two-letter species prefix at the beginning of its code). And why do I want to get the locus_tag of each protein in each syntenic block? We'll see later.

![][image31]

10. Next, the first result (step 7) of the GO count in each interspecies syntenic block was used to easily locate some large syntenic blocks in SynVisio, visualize them, and be able to verify that they are enriched in certain GO terms, as seen in the image below. Although initially the enrichment inspection was visual, only seeing the description of the most frequent GO terms in each syntenic block (see image from step 7), then the tag locus of the proteins of my syntenic block of interest (from step 9) was used and introduced into the [DAVID](https://davidbioinformatics.nih.gov/tools.jsp) platform, which resulted in the enrichment of GO terms in the blocksyntenicpecific, which is seen in Figure 11.

11. On the other hand, comparisons were made between all the blocks. syntenicbetween each pair of species, but, as mentioned above, only between* Arabidopsis * and * Cuscuta * A good number of conserved blocks were obtained, we decided to focus on these species. Then, all syntenic blocks were chosen between* Arabidopsis * and * Cuscuta * , I extracted all of their locus_tags (from step 9) and entered them into Panther via the [gene ontology platform](https://geneontology.org/), DAVID was not used because it did not accept such a large number of locus_tags. The result was a table, which was downloaded and processed with the R script **figures. Rmd** (see STEP 11 in the script), to filter out those with the highest Fold_Enrichment and p-value adjusted for multiple testing. The output was the following HTML file.

![][image32]

From this table the most enriched GO terms were examined when considering all collinear syntenic blocks between * Arabidopsis* and * Cuscuta *and figures 8, 9 and 10 were obtained.

## 4. Results

Figure 1. Continuity results of the assembled genomes

![][image33]

Figure 2. Integrity results of the assembled genomes

![][image34]

Figure 3. Gene duplication events

![][image35]

Table 2. Matrix of common orthologs between species.  
![][image36]

Each entry (i, j) is the number of proteins in j that have orthologs in species i.*Chlamydomonas reinhardtii* has few orthologous sequences in common with the other species. While *Cuscuta* and *Roridula* have a high number of orthologous genes in common with *Arabidopsis* and among themselves.

Figure 4. Percentage of sequences from each species assigned to singletons (orthologous groups of a single sequence), orthologous groups of two sequences, and orthologous groups with more than three sequences.![][image37]  
Figure 5. Percentage of genes of each species assigned to singletons, orthologous groups only detected in that species (Species specific), and orthologous groups shared with at least one species (Shared orthogroups).

![][image38]

Figure 6. Collinearity of three terrestrial photosynthetic species

| ![][image39] |  ![][image40]  Share of Collinear Genesindicates that of the total genes analyzed, \~16% are found in collinear blocks, that is, conserved genomic regions. |
| :---- | :---- |

The fragmentscolored lines represent chromosomes, while gray lines indicate common syntenic blocks.

Table 3. GO term count for each species  
![][image41]  
Each of the GOs is shown with their respective names, in addition to the number of genes contained in each GO according to the species.

Figure 7. Bar chart of counts by GO terms (number of genes) and species

![][image42]

Figure 8. GO terms enriched in syntenic blocks between* Arabidopsis * and * Cuscuta *.  

![][image43]

GO terms associated with the regulation of intracellular signals (kinase activity) and cell cycle (actin nucleation, replication, cyclins, cytokinesis) are observed.

Figure 9. GO terms enriched in syntenic blocks between* Arabidopsis * and * Cuscuta *. Ionic transport.

![][image44]

GO terms associated with ionic transport are observed

Figure 10. GO terms enriched in syntenic blocks between* Arabidopsis * and * Cuscuta *. 

![][image45]

GO terms associated with embryonic development and pattern establishment are observed.

Figure 11. Largest syntenic block

![][image46]

Figure 12. Global panorama of synteny between Arabidopsis and Cuscuta

![][image47]

## 5. Discussion

From the analysis carried out with QUAST (see figure 1), the assembly with the greatest continuity (higher N50 and lower L50) corresponds to *Arabidopsis thaliana*, followed by *Chlamydomonas reinhardtii*. This is because both species are widely studied models in molecular and evolutionary biology. Furthermore, these genomes were obtained from the NCBI RefSeq Assembly database, ensuring their quality. Furthermore, the genome size is relatively small, where *Roridula gorgonias* has the largest genome among the four organisms, however, it is the most fragmented. This genome was downloaded from Dryad, a platform where research data from other articles are available and can be reused (Dryad, 2024). In this case, the genome sequence of *Roridula gorgonias* obtained from the article titled “Annotated genome sequences of the carnivorous plant *Roridula gorgonias* and a non-carnivorous relative, *Clethra arborea*". Within this, the limitations of the data set are mentioned, since it is a draft of the genome, therefore, some genes located at the end of the scaffolds may be incomplete, and repeated regions may be missing or poorly assembled (Hartmann, S. et al., 2020). This is supported by the QUAST results seen at the beginning. Furthermore, the fragmented genome agrees with the BUSCO result (Figure 2), where *Roridula* appears to be the species with the most fragmented single-copy orthologs.	

In Figure 3 (Gene duplication events), we see the clustering of the species. The scale of 0.45 substitutions per site indicates that the total genome of the clade of *Roridula gorgonias* and *Cuscuta* has diverged in nucleotide distance almost as much from that of *Arabidopsis*, like that of *Arabidopsis* has diverged from that of *Chlamydomonas*, which is a result of the great evolutionary distance of *Chlamydomonas* with respect to other species, and the large number of changes suffered by *Roridula* and *Cuscuta* with respect to *Arabidopsis*, probably due to genome reduction. Values such as 5527 in “chlamydomonas_5527” represent the number of gene duplications with at least 50% support occurring on that branch. Gene family expansion by duplication in *Chlamydomonas* has already been reported in other articles, due to its evolutionary distance from other terrestrial plants (Guo, 2012), specific adaptations to aquatic environments and characteristics of the adaptation of its genes to unicellular life. Likewise, the clade that leads to *Arabidopsis* has **24,751 duplications**, which does not share with *Roridula gorgonias* and *Cuscuta*. *Cuscuta* and *Roridula* have suffered reductions in their genome, so the duplications of *Arabidopsis* reflect the large number of duplications and expansions of gene families in angiosperms (Guo, 2012; Cannon et al., 2004). As *Roridula* and *Cuscuta* are the only angiosperm species against which the OrthoFinder program can compare, a good portion of the number of orthologs being obtained only in *Arabidopsis* actually corresponds to most dicotyledons, but is assigned only to *Arabidopsis* because they are no longer in *Roridula* and *Cuscuta*. 

Also, a good number of the reported duplications may be due to the rounds of duplication, gene expansion and rearrangements detected in the specific Brassicales order of *Arabidopsis* (Cannon et al., 2004). This is also consistent with the reduced number of proteins in both carnivorous and parasitic species when compared to *Arabidopsis* as a reference (Sun et al., 2018; Fleck & Jobson, 2023). Likewise, it cannot be ruled out that the high number of duplications and orthologs in *Arabidopsis* may also occur due to a more detailed annotation in its genome, since it is a model plant.

In Table 2 (matrix of orthologs in common between species) we are shown the orthologous proteins common between species. Where *Chlamydomonas reinhardtii* has few orthologous sequences in common with the other species (approximately only 25% of its proteins). For its part, *Cuscuta* and *Roridula* have a high number of orthologous genes in common with *Arabidopsis* and among themselves. In fact, 72% of the sequences of *Cuscuta* (13042/18157) and 78% of those of *Roridula* , are orthologous with those of *Arabidopsis* (17778/22655). It is also interesting to note that a similar number (70%) are shared among themselves. These orthologous proteins could be part of an essential proteome, with basic cellular functions, shared among all these Angiosperms.

In Figure 4 we see that *Chlamydomonas* is the species with the highest number of singletons or genes without any homologues (orphan genes). They also have the highest number of sequences participating in orthologous groups of two genes, which are probably also orthologous groups specific to this species. Furthermore, Figure 5 shows us that the overwhelming majority of their proteins belong to species-specific orthologous groups, *Chlamydomonas*, in addition to singletons (which by definition are only present in the species in which they are detected). The percentage of proteins in shared orthologous groups is 25%, which was also verified in Table 2 (Matrix of orthologs in common between species). These results are consistent with the report by Guo (2013). Likewise, *Arabidopsis* has a high number of orthologous groups reported as unique to this species, for reasons already mentioned such as the gene expansion suffered in the lineage of Brassicales and angiosperms, the loss of many of these genes during the reduction of the genomes of parasitic and carnivorous plants, as well as the better annotation of the genome of *Arabidopsis*.

According to the result of MCScanX (see figure 6) the syntenic blocks were found in common between *Chlamydomonas* and other species, which is expected given the low percentage of orthologous genes found in common with other species and the high number of singletons or orphan genes (Figure 5 and Table 2). *Chlamydomonas* has experienced a very significant evolutionary divergence of approximately 725 MA (Guo 2012, see Figure 3. Gene duplication events). Moreover, in the collinearity image, it can be seen that the genome of *Roridula* is highly fragmented compared to *Cuscuta* and *Arabidopsis*. Although *Roridula* has the largest genome among the four organisms, from the beginning there were a number of limitations in the data set, given that it is a draft genome without assembly at the scafold level. This prevented the recovery of large syntenic blocks and for this reason, the identification of syntenic blocks focused exclusively on *A. thaliana* and *C. australis*. 

Before analyzing the syntenic blocks, we will describe some differential functions detected based on the GO terms. Table 3 and Figure 7 show the GO term count (number of genes with that term) and species. *Arabidopsis thaliana* has the largest number of genes in each GO term, except for the GO for cilia-dependent cell motility regulation (which as expected is unique to *Chlamydomonas*, a single-celled organism that moves with cilia). As explained before, this may be due to greater annotation and gene duplication or even, to a lesser extent, horizontal gene transfer (HGT). Although the role of HGT in multicellular eukaryotic organisms is much less clear and is sometimes considered anecdotal (Jianchao Ma, 2022). However, we do not rule out horizontal gene transfer (HGT) because two significant HGT events occurred in the evolution of land plants. Furthermore, the vast majority of genes acquired in both events have been conserved in descendant groups and are involved in stress responses, ion and metabolite transport, as well as growth and development; and specialized metabolism as in the case of carnivorous and parasitic plants (Jianchao Ma, 2022). Likewise, it is observed that *Cuscuta* has a small number of genes involved in most pathways, including photosynthesis, nitrate import, and the dark reaction of photosynthesis. This may be due to gene loss events resulting from its parasitic lifestyle. This hypothesis is supported by previous studies, which report that *Cuscuta* has lost 11.7% of the orthologs normally conserved in autotrophic plants (Guiling Sun, 2018), including structural genes involved in the cycle of electrons obtained in photosystem I (the genes etc) (Guiling Sun, 2018). The case of *Roridula* is very similar to that of *Cuscuta*, and expected, because it is also known that carnivorous species have a reduced genome and suffer changes in obtaining nutrients (Fleck & Jobson, 2023). In both species, a decline in genes for the dark phase of photosynthesis and its regulation is observed, due to the fact that they obtain energy and carbon from other sources.

Finally, in the main result, the GO terms enriched in syntenic blocks between *Arabidopsis* and *Cuscuta* are related to the regulation of intracellular signals (kinase activity), and cell cycle (actin nucleation, replication, cyclins, cytokinesis) (see figure 8), ionic transport (see figure 9), and embryonic development and pattern establishment (see figure 10). These results correspond to what was mentioned in the article by Jianchao Ma (2022), who points out some genes involved in transport of metabolites and ions, as well as in growth and development were acquired during two major HGT events and have often accumulated, duplicated, or functionally differentiated into descendant groups, thus contributing to the diversification and long-term evolution of land plants (Jianchao Ma, 2022). In addition, we decided to examine which genes or metabolic pathways were present in the largest syntenic block between *Arabidopsis* and *Cuscuta*. It turned out that the largest syntenic block (see Figure 11) contains conserved genes for DNA regulatory transcription factors. These transcription factors are grouped into 47 collinear genes, which contain master regulators of response to abiotic stress, light stimuli, radiation, response to chemicals, bacterial defense, and regulation of macromolecule production (see Appendices). Therefore, we understand which are a critical part of the adaptation to terrestrial life between these two plants, to deal with biotic and abiotic stresses that they have in common. Likewise, it makes sense that this result was not obtained with *Chlamydomonas*, an aquatic organism. It is also possible that *Roridula* had similar results to those found in *Arabidopsis* and *Cuscuta*, however, to corroborate this, a less fragmented genome with a greater number of annotations is needed.

## 6. Conclusions

* Most of the proteins in *Chlamydomonas* are unique to this species, lack GO terms, and did not find orthologous groups or syntenic blocks when comparing them with other species, probably because they must face particular stresses of the aquatic environment.  
* *Cuscuta* and *Roridula* have a reduced number of proteins compared to *Arabidopsis*, including reduction in proteins associated with photosynthesis.  
* The highest number of syntenic blocks was found between *Arabidopsis* and *Cuscuta*, and although both species have numerous orthologous genes with *Roridula*, no syntenic blocks were formed due to the fragmentation of the latter's genome.  
* The conserved metabolic pathways between *Arabidopsis* and *Cuscuta* include genes regulating the cell cycle, ion transport, embryonic development, intracellular signaling cascades, and transcription factors regulating the response to biotic and abiotic stress.

## 7. Literature

● Cannon, S. B., Mitra, A., Baumgarten, A., Young, N. D., & May, G. (2004). The roles of segmental and tandem gene duplication in the evolution of large gene families in *Arabidopsis thaliana*. BMC plant biology, 4, 1-21.  
● Dryad. (2024). Who we are. https://datadryad.org/stash/about  
● Fleck, S. J., & Jobson, R. W. (2023). Molecular Phylogenomics Reveals the Deep Evolutionary History of Carnivory across Land Plants. Plants, 12(19), 3356. https://doi.org/10.3390/plants12193356  
● Guiling Sun, Yuxing Xu, Hui Liu, Ting Sun, Jingxiong Zhang, Christian Hettenhausen, Guojing Shen, Jinfeng Qi, Yan Qin, Jing Li, Lei Wang, Wei Chang, Zhenhua Guo, Ian T. Baldwin and Jianqiang Wu. (2018). Large-scale gene losses underlie the genome evolution of the parasitic plant *Cuscuta australis*. Nature Communications. Vol 9. https://www.nature.com/articles/s41467-018-04721-8  
● Guo Y. L. (2012). Gene family evolution in green plants with emphasis on the origination and evolution of *Arabidopsis thaliana* genes. The Plant Journal. https://onlinelibrary.wiley.com/doi/10.1111/tpj.12089   
● Hartmann, S., Preick, M., Abelt, S., Scheffel, A. and Hofreiter, M. (2020). Annotated genome sequences of the carnivorous plant *Roridula gorgonias* and a non-carnivorous relative, *Clethra arborea*. BMC Research Notes. https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-020-05254-4  
● Jianchao Ma, Shuanghua Wang, Xiaojing Zhu, Guiling Sun, Guanxiao Chang, Linhong Li, Xiangyang Hu, Shouzhou Zhang, Yun Zhou, Chun-Peng Song, and Jinling Huang. (2022). Major episodes of horizontal gene transfer drove the evolution of land plants. Molecular Plant. Volume 15, pag. 857-871. https://www.sciencedirect.com/science/article/pii/S1674205222000491?via%3Dihub  
● Sun, G., Xu, Y., Liu, H. et al. (2018). Large-scale gene losses underlie genome evolution of parasitic plant *Cuscuta australis*. NatCommun 9, 2683. https://www.nature.com/article/s41467-018-04721-8#citeas 

## 8. Attachments

Visualization of the GO terms present in the largest syntenic block in common between *Arabidopsis* and *Cuscuta*. You can see their description in the third column, and the number of genes in which they appear in the fourth column.

![][image48]  
![][image49]  
![][image50]  
![][image51]

[image2]: docs/image2.png
[image3]: docs/image3.png
[image4]: docs/image4.png
[image5]: docs/image5.png
[image6]: docs/image6.png
[image7]: docs/image7.png
[image8]: docs/image8.png
[image9]: docs/image9.png
[image10]: docs/image10.png
[image11]: docs/image11.png
[image12]: docs/image12.png
[image13]: docs/image13.png
[image14]: docs/image14.png
[image15]: docs/image15.png
[image16]: docs/image16.png
[image17]: docs/image17.png
[image18]: docs/image18.png
[image19]: docs/image19.png
[image20]: docs/image20.png
[image21]: docs/image21.png
[image22]: docs/image22.png
[image23]: docs/image23.png
[image24]: docs/image24.png
[image25]: docs/image25.png
[image26]: docs/image26.png
[image27]: docs/image27.png
[image28]: docs/image28.png
[image29]: docs/image29.png
[image30]: docs/image30.png
[image31]: docs/image31.png
[image32]: docs/image32.png
[image33]: docs/image33.png
[image34]: docs/image34.png
[image35]: docs/image35.png
[image36]: docs/image36.png
[image37]: docs/image37.png
[image38]: docs/image38.png
[image39]: docs/image39.png
[image40]: docs/image40.png
[image41]: docs/image41.png
[image42]: docs/image42.png
[image43]: docs/image43.png
[image44]: docs/image44.png
[image45]: docs/image45.png
[image46]: docs/image46.png
[image47]: docs/image47.png
[image48]: docs/image48.png
[image49]: docs/image49.png
[image50]: docs/image50.png
[image51]: docs/image51.png
