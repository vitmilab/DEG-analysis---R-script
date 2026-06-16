# Comparative Transcriptomic Analysis Identifies Conserved Abiotic Stress-Responsive Hub Genes in Sorghum bicolor

This repository contains the R scripts used for the differential gene expression analysis performed in the study:
“Integrative Transcriptomic Analysis Reveals Conserved Abiotic Stress-Responsive Hub Genes in Sorghum bicolor”

The study integrates publicly available microarray and RNA-seq datasets to identify conserved abiotic stress-responsive hub genes associated with heat–drought and saline–alkali stress conditions in sorghum.

## Files Included
GSE48205.R - R script for microarray differential expression analysis of dataset GSE48205 using GEOquery and limma | 
GSE140928.R	- R script for RNA-seq differential expression analysis of dataset GSE140928 using DESeq2
External Transcriptomic Validation
GSE80699_external_validation.R
Extraction of candidate hub gene expression values from GSE80699, calculation of log2 fold changes using published RPKM values, generation of Table S4, and creation of Supplementary Figure S2.

## Data Availability
The datasets analysed in this study are publicly available in the NCBI Gene Expression Omnibus (GEO):
GSE48205
GSE140928
GSE80699

## Analysis Workflow

--> The analysis includes:
1. Data preprocessing
2. Normalization
3. Differential gene expression analysis
4. DEG filtering using:
5. Adjusted p-value < 0.05
6. |log2 Fold Change| > 1
7. Volcano plot generation

--> The overlapping differentially expressed genes (DEGs) obtained from both datasets were used for downstream analyses including:
1. Protein–protein interaction (PPI) network construction
2. Hub gene identification
3. Gene Ontology enrichment analysis
4. Software Requirements

--> External Validation

The four prioritized candidate hub genes were evaluated using the independent drought-stress RNA-seq dataset GSE80699.

Validation included:

1. Extraction of published RPKM expression values
2. Calculation of log2 fold changes between drought and control samples
3. Generation of validation summary tables
4. Visualization of expression patterns across drought-sensitive and drought-tolerant sorghum genotypes

--> The scripts were executed in R using the following packages:
1. GEOquery
2. limma
3. DESeq2
4. ggplot2
5. umap
