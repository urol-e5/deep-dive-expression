# F-Ptuh Code Directory

## Overview

This directory contains the complete bioinformatics pipeline for analyzing *Porites tuahiniensis* multi-omics data, including RNA-seq, small RNA-seq, and whole genome bisulfite sequencing (WGBS) data. The analysis focuses on understanding gene expression regulation through miRNA-mRNA-lncRNA interactions, DNA methylation patterns, and functional annotation in this coral species.

## Workflow Schematic

```mermaid
graph TD
    A[Raw Data] --> B[Quality Control]
    B --> C[Data Processing]
    C --> D[Alignment & Quantification]
    D --> E[Functional Analysis]
    E --> F[Network Analysis]
    
    A1[RNA-seq FastQ] --> B1[FastQC/MultiQC]
    A2[sRNA-seq FastQ] --> B2[FastQC/MultiQC]
    A3[WGBS FastQ] --> B3[FastQC/MultiQC]
    
    B1 --> C1[RNA Trimming]
    B2 --> C2[sRNA Trimming]
    B3 --> C3[WGBS Trimming]
    
    C1 --> D1[Hisat2 Alignment]
    C2 --> D2[ShortStack Analysis]
    C3 --> D3[Bismark Alignment]
    
    D1 --> E1[Gene Expression Analysis]
    D2 --> E2[miRNA Discovery]
    D3 --> E3[Methylation Analysis]
    
    E1 --> F1[mRNA-miRNA Interactions]
    E2 --> F2[miRNA-lncRNA Interactions]
    E3 --> F3[Correlation Networks]
    
    F1 --> G[Integrated Network Analysis]
    F2 --> G
    F3 --> G
```

## Key Results

- **Gene Expression**: Identified differentially expressed genes across *P. tuahiniensis* samples
- **miRNA Discovery**: Discovered and annotated miRNAs using ShortStack 4.1.0
- **DNA Methylation**: Analyzed genome-wide methylation patterns using Bismark
- **Regulatory Networks**: Constructed miRNA-mRNA-lncRNA interaction networks
- **Functional Annotation**: Performed GO enrichment analysis and functional annotation
- **Species-Specific Analysis**: Comparative analysis with other coral species

## File Descriptions

### Data Quality Control and Processing

- **`00.00-F-Ptuh-WGBS-reads-FastQC-MultiQC.Rmd`**: Downloads and performs quality control on WGBS sequencing reads using FastQC and MultiQC
- **`01-Ptuh-RNA-trimming-FastQC.Rmd`**: Trims RNA-seq reads and performs quality control on both raw and trimmed reads
- **`01.00-F-Ptuh-WGBS-trimming-cutadapt-FastQC-MultiQC.Rmd`**: Trims WGBS reads using cutadapt and performs quality control
- **`01.00-F-Ptuh-WGBS-trimming-fastp-FastQC-MultiQC.Rmd`**: Alternative WGBS trimming using fastp

### Reference and Annotation

- **`00.00-genome-GFF-formatting.Rmd`**: Formats and prepares genome GFF files for analysis
- **`02-Ptuh-reference-annotation.Rmd`**: Annotates the *P. tuahiniensis* reference genome with functional information and GO terms
- **`04-Ptuh-genome-explore.Rmd`**: Explores and analyzes the *P. tuahiniensis* genome characteristics
- **`16-Ptuh-annotate-UTRs.Rmd`**: Annotates 5' and 3' UTR regions in the genome

### RNA-seq Analysis

- **`03-Ptuh-RNA-summary.Rmd`**: Summarizes RNA-seq gene expression data and performs differential expression analysis
- **`03.1-Ptuh-sRNA-summary.Rmd`**: Summarizes small RNA-seq data and identifies different RNA classes
- **`03.2-Ptuh-lncRNA-summary.Rmd`**: Identifies and summarizes long non-coding RNAs
- **`06-Ptuh-Hisat.qmd`**: RNA-seq alignment using Hisat2
- **`06.2-Ptuh-Hisat.qmd`**: Additional Hisat2 analysis

### Small RNA Analysis

- **`05-Ptuh-sRNA-ShortStack_4.1.0.Rmd`**: Discovers and annotates small RNAs using ShortStack 4.1.0
- **`05.1-Ptuh-sRNA-ShortStack-rename.Rmd`**: Renames ShortStack output files for consistency

### DNA Methylation Analysis

- **`12-Ptuh-WGBS.Rmd`**: Performs whole genome bisulfite sequencing analysis using Bismark
- **`12-Ptuh-methylation.qmd`**: Analyzes DNA methylation patterns and differential methylation

### Interaction Analysis

- **`11-Ptuh-mRNA-miRNA-interactions.Rmd`**: Identifies and analyzes miRNA-mRNA interactions
- **`11.01-Ptuh-mRNA-miRNA-interactions-CDS_5UTR.Rmd`**: Analyzes miRNA interactions with CDS and 5' UTR regions
- **`11.1-Ptuh-mRNA-miRNA-interactions-functional-enrichment.Rmd`**: Performs functional enrichment analysis on miRNA-mRNA interactions
- **`11.2-Ptuh-known-miRNA-interactions-functional-enrichment.Rmd`**: Functional enrichment analysis for known miRNA interactions
- **`11.11-Ptuh-mRNA-miRNA-interactions-FE-CDS.Rmd`**: Functional enrichment analysis for CDS regions
- **`11.12-Ptuh-mRNA-miRNA-interactions-FE-5UTR.Rmd`**: Functional enrichment analysis for 5' UTR regions
- **`11.13-Ptuh-mRNA-miRNA-interactions-FE-3UTR.Rmd`**: Functional enrichment analysis for 3' UTR regions
- **`11.14-Ptuh-mRNA-miRNA-interactions-FE-pooled.Rmd`**: Pooled functional enrichment analysis
- **`14-Ptuh-miRNA-lncRNA-BLASTs-miRanda.Rmd`**: Analyzes miRNA-lncRNA interactions using BLAST and miRanda
- **`15-Ptuh-miRNA-lncRNA-PCC.Rmd`**: Calculates Pearson correlation coefficients between miRNAs and lncRNAs

### Network Analysis

- **`09-Ptuh-mixOmics.Rmd`**: Performs multi-omics integration analysis using mixOmics
- **`10-Ptuh-mRNA-lncRNA-correlation-networks.Rmd`**: Constructs correlation networks between mRNAs and lncRNAs
- **`10.1-Ptuh-mRNA-lncRNA-correlation-PCC.Rmd`**: Calculates Pearson correlation coefficients between mRNAs and lncRNAs
- **`31-Ptuh-miRNA-lncRNA-network.Rmd`**: Integrates all interaction data into a comprehensive network

### Functional Analysis

- **`13-Ptuh-mRNA-GO-enrichment.Rmd`**: GO enrichment analysis for mRNAs
- **`30.00-Ptua-transcriptome-GOslims.Rmd`**: GO-slim analysis of the transcriptome

### Data Management

- **`08-Ptuh-lncRNA-matrix.qmd`**: Creates lncRNA count matrix
- **`17-Ptuh-lncRNA.Rmd`**: Comprehensive lncRNA analysis
- **`18-Ptuh-lncRNA-matrix.qmd`**: Updated lncRNA matrix creation

### Utility Scripts

- **`references.bib`**: Bibliography file for citations
- **`.gitignore`**: Git ignore file for version control

## Data Flow

1. **Raw Data Processing**: Quality control and trimming of sequencing reads
2. **Alignment**: Mapping reads to reference genome using appropriate tools
3. **Quantification**: Counting reads and estimating expression levels
4. **Discovery**: Identifying novel RNAs (miRNAs, lncRNAs) and methylation sites
5. **Interaction Analysis**: Predicting and validating molecular interactions
6. **Network Construction**: Building regulatory networks from interaction data
7. **Functional Analysis**: Enrichment analysis and functional annotation
8. **Integration**: Combining all analyses into comprehensive networks

## Key Software Tools

- **FastQC/MultiQC**: Quality control
- **Hisat2**: RNA-seq alignment
- **ShortStack**: Small RNA analysis
- **Bismark**: DNA methylation analysis
- **miRanda**: miRNA target prediction
- **mixOmics**: Multi-omics integration
- **DESeq2**: Differential expression analysis
- **R/Bioconductor**: Statistical analysis and visualization

## Output Structure

All outputs are organized in the `../output/` directory with subdirectories corresponding to each analysis step. Key outputs include:

- Quality control reports (FastQC, MultiQC)
- Alignment files (BAM, SAM)
- Count matrices (gene, miRNA, lncRNA)
- Interaction predictions (miRNA-mRNA, miRNA-lncRNA)
- Network files (CSV, GML)
- Functional enrichment results
- Visualization plots (PDF, PNG)

## Citation

Please cite the original papers for the software tools used in this analysis, as well as any relevant publications from the E5 coral research group.
