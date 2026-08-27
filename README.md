# deep-dive-expression

[![DOI](https://img.shields.io/badge/DOI-pending-blue)]()

A comprehensive multi-omics bioinformatics pipeline for analyzing gene expression, ncRNA co-expression, and DNA methylation interactions in three species of stony corals. This work builds on the ncRNA landscape analysis from the [deep-dive](https://github.com/urol-e5/deep-dive) project, and is currently available as a preprint on bioRxiv: 

[**Cross-talk among miRNAs, lncRNAs, and DNA methylation in three coral species reveal conserved epigenetic regulatory architecture**](https://doi.org/10.64898/2026.07.19.739451)
Kathleen M. Durkin, Jill Ashey, Javier A. Rodriguez-Casariego, Zoe Dellaert, Zachary Bengtsson, Samuel J. White, Jose M. Eirin-Lopez, Hollie M. Putnam, Steven B. Roberts

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Prerequisites](#prerequisites)
- [Citation](#citation)
- [Related Projects](#related-projects)
- [Contact](#contact)

## Overview

This repository contains a complete bioinformatics pipeline for analyzing matched whole-genome bisulfite sequencing (WGBS), RNA-sequencing, and small RNA-sequencing data from three stony coral species (*Acropora pulchra*, *Porites evermanni*, and *Pocillopora tuahiniensis*). Analyses include:

- **Gene Expression Analysis**: RNA-seq data processing
- **Small RNA Analysis**: miRNA discovery and annotation, and target prediction
- **Long Non-coding RNA**: lncRNA identification and characterization
- **DNA Methylation**: Whole genome bisulfite sequencing (WGBS) analysis and methylation pattern characterization
- **Regulatory Networks**: Integration of mRNA, miRNA, lncRNA, and DNA methylation data, including identification of candidate epi-miRNAs and ceRNA triads
- **Functional Annotation**: GO enrichment analysis and functional characterization
- **Comparative Analysis**: Cross-species comparisons and orthology analysis

Here are some resources to get your started!

- [DIVE genome browser](https://urol-e5.github.io/deep-dive-genome-browser/) - View ncRNA features across all three species in our interactive genome browser.
- [CENE network viewer](https://urol-e5.github.io/CENE/) - explore mRNA-miRNA-lncRNA interaction networks, including epi-miRNAs and ceRNA networks, across the three species.
- [OSF Project Link](https://osf.io/aw53f/) - Access large genomic reference files in our Open Science Framework (OSF) storage directory. 
- [![SRA](https://img.shields.io/badge/SRA-PRJNA1201098-blue)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1201098/) - Access raw RNA sequences stored on NCBI. 
- [![SRA](https://img.shields.io/badge/SRA-PRJNA1236666-blue)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1236666/) - Access raw small RNA sequences stored on NCBI. 

- [**Sample metadata**](M-multi-species/data/e5_deep_dive_WGBS_metadata.csv) - per-sample colony IDs, species, collection date/site (Mo'orea,
    2020-03-05), extraction & library-prep details. Colony ID prefixes: `ACR-` (*A. pulchra*), `POR-` (*P. evermanni*), `POC-` (*P. tuahiniensis*).

Full genomic resources, sequencing data files, and count matrices can be found in the [project wiki](https://github.com/urol-e5/deep-dive-expression/wiki).

## Repository Structure

The repository is organized into species-specific directories and a multi-species comparison directory:

```
deep-dive-expression/
├── data/                # Reference files used throughout the below directories
├── D-Apul/              # Acropora pulchra analyses
│   ├── code/               # Analysis scripts (*.Rmd, *.qmd)
│   ├── data/               # Input data and references
│   └── output/             # Analysis outputs (organized by script)
├── E-Peve/              # Porites evermanni analyses
│   ├── code/
│   ├── data/
│   └── output/
├── F-Ptuh/              # Pocillopora tuahiniensis analyses
│   ├── code/
│   ├── data/
│   └── output/
├── M-multi-species/     # Cross-species comparative analyses
│   ├── code/
│   ├── data/
│   └── output/
├── manuscript_files/    # Figures and tables included in manuscript
├── references.bib       # Shared bibliography
└── README.md            # This file
```

Each species directory (`D-Apul`, `E-Peve`, `F-Ptuh`) follows the same structure:

- **`code/`**: Contains all analysis scripts with numbered prefixes (e.g., `01-*.Rmd`, `02-*.Rmd`)
  - See species-specific README files for detailed descriptions of each script
- **`data/`**: Contains input data, reference genomes, and annotations
  - Large data files (FASTQ, BAM) are linked via the [project wiki](https://github.com/urol-e5/deep-dive-expression/wiki)
- **`output/`**: Contains subdirectories matching script names (e.g., output from `01-script.Rmd` goes in `output/01-script/`)


## Prerequisites

### System Requirements
- Linux/Unix environment (Ubuntu 20.04+ recommended)
- Minimum 32 GB RAM (64 GB+ recommended for large genomes)
- 500 GB+ available disk space
- Multi-core processor (8+ cores recommended)

### Software Dependencies

The following software tools are required:

#### Core Bioinformatics Tools
- **FastQC** (v0.11.9+): Quality control
- **MultiQC** (v1.12+): Aggregated QC reports
- **fastp** (v0.23+) / **cutadapt** (v3.5+): Adapter trimming
- **HISAT2** (v2.2.1+): RNA-seq alignment
- **StringTie** (v2.2+): Transcript assembly and quantification
- **Bismark** (v0.23+): WGBS alignment and methylation calling
- **ShortStack** (v4.1.0): Small RNA analysis
- **Bowtie** / **Bowtie2**: Short read alignment
- **SAMtools** (v1.14+): BAM file manipulation
- **BEDtools** (v2.30+): Genomic interval operations

#### Target Prediction and Interaction Analysis
- **miRanda** (v3.3a+): miRNA target prediction

#### Statistical Analysis and Visualization
- **R** (v4.1+): Statistical computing
- **RStudio** (recommended): IDE for R

#### R/Bioconductor Packages
```r
# Data manipulation and visualization
install.packages(c("tidyverse", "ggplot2", "pheatmap", "RColorBrewer"))

# Bioconductor packages
BiocManager::install(c(
  "DESeq2",           # Differential expression
  "edgeR",            # RNA-seq analysis
  "WGCNA",            # Co-expression networks
  "methylKit",        # DNA methylation analysis
  "GenomicRanges",    # Genomic intervals
  "rtracklayer",      # Genome annotations
  "clusterProfiler",  # GO enrichment
  "biomaRt"           # Annotation databases
))
```

#### Other Tools
- **CPC2**: Coding potential calculator
- **GffCompare**: Transcript comparison
- **OrthoFinder**: Orthology analysis
- **BLAST+**: Sequence similarity searches

Refer to individual script headers for specific version requirements and conda environment configurations.

#### Reporting Issues

Found a bug or have a suggestion? Please [open an issue](https://github.com/urol-e5/deep-dive-expression/issues) with:
- Clear description of the problem or suggestion
- Steps to reproduce (for bugs)
- Expected vs. actual behavior
- System information (OS, R version, package versions)

#### Code of Conduct

Please be respectful and constructive in all interactions. We aim to maintain a welcoming and inclusive environment for all contributors.

## Citation

If you use this pipeline or data in your research, please cite:

```bibtex
@misc{deepdive-expression,
  title = {deep-dive-expression: Multi-omics analysis of coral gene expression and regulatory networks},
  author = {E5 Coral Research Group},
  year = {2026},
  url = {https://github.com/urol-e5/deep-dive-expression},
  note = {GitHub repository}
}
```

## Related Projects

- **[deep-dive](https://github.com/urol-e5/deep-dive)**: Predecessor project focusing on ncRNA landscape characterization
- **[E5 Coral Project](https://e5coral.org/)**: Broader coral resilience and epigenetics research initiative

## Contact

For questions, issues, or collaborations:

- **Issues**: [GitHub Issues](https://github.com/urol-e5/deep-dive-expression/issues)
- **Wiki**: [Project Wiki](https://github.com/urol-e5/deep-dive-expression/wiki)
- **E5 Coral Project**: Visit the main E5 project page for team information

---

**Last Updated**: August 2026

**Repository maintained by**: E5 Coral Research Group
