# deep-dive-expression

[![DOI](https://img.shields.io/badge/DOI-pending-blue)]()

A comprehensive multi-omics bioinformatics pipeline for analyzing gene expression, ncRNA co-expression, and DNA methylation interactions in three species of stony corals. This work builds on the ncRNA landscape analysis from the [deep-dive](https://github.com/urol-e5/deep-dive) project.

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Species Coverage](#species-coverage)
- [Analysis Workflows](#analysis-workflows)
- [Datasets](#datasets)
- [Prerequisites](#prerequisites)
- [Getting Started](#getting-started)
- [Usage Examples](#usage-examples)
- [File Naming Conventions](#file-naming-conventions)
- [Software Requirements](#software-requirements)
- [Data Access](#data-access)
- [Contributing](#contributing)
- [Citation](#citation)
- [Related Projects](#related-projects)
- [Contact](#contact)

## Overview

This repository contains a complete bioinformatics pipeline for analyzing multi-omics data from three coral species, focusing on:

- **Gene Expression Analysis**: RNA-seq data processing, differential expression analysis, and co-expression networks
- **Small RNA Analysis**: miRNA discovery, annotation, and target prediction using ShortStack 4.1.0
- **Long Non-coding RNA**: lncRNA identification, characterization, and interaction analysis
- **DNA Methylation**: Whole genome bisulfite sequencing (WGBS) analysis and methylation pattern characterization
- **Regulatory Networks**: Integration of mRNA-miRNA-lncRNA interactions with methylation data
- **Functional Annotation**: GO enrichment analysis and functional characterization
- **Comparative Analysis**: Cross-species comparisons and orthology analysis

### Key Features

- **Reproducible workflows**: All analyses documented in R Markdown/Quarto format
- **Modular structure**: Species-specific and multi-species analysis modules
- **Comprehensive documentation**: Each analysis script includes detailed methods and citations
- **Quality control**: FastQC/MultiQC reports at every processing step
- **Standardized outputs**: Consistent output organization across all analyses

## Repository Structure

The repository is organized into species-specific directories and a multi-species comparison directory:

```
deep-dive-expression/
├── D-Apul/              # Acropora pulchra analyses
│   ├── code/            # Analysis scripts (*.Rmd, *.qmd)
│   ├── data/            # Input data and references
│   └── output/          # Analysis outputs (organized by script)
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
├── references.bib       # Shared bibliography
└── README.md           # This file
```

### Directory Organization

Each species directory (`D-Apul`, `E-Peve`, `F-Ptuh`) follows the same structure:

- **`code/`**: Contains all analysis scripts with numbered prefixes (e.g., `01-*.Rmd`, `02-*.Rmd`)
  - See species-specific README files for detailed descriptions of each script
- **`data/`**: Contains input data, reference genomes, and annotations
  - Large data files (FASTQ, BAM) are linked via the [project wiki](https://github.com/urol-e5/deep-dive-expression/wiki)
- **`output/`**: Contains subdirectories matching script names (e.g., output from `01-script.Rmd` goes in `output/01-script/`)

## Species Coverage

This project focuses on three ecologically important Pacific coral species:

| Code | Species | Common Name | Genome |
|------|---------|-------------|---------|
| **D-Apul** | *Acropora pulchra* | Table coral | [NCBI](https://github.com/urol-e5/deep-dive-expression/wiki) |
| **E-Peve** | *Porites evermanni* | Evermann's coral | [NCBI](https://github.com/urol-e5/deep-dive-expression/wiki) |
| **F-Ptuh** | *Pocillopora tuahiniensis* | Cauliflower coral | [NCBI](https://github.com/urol-e5/deep-dive-expression/wiki) |

Detailed genome information and resources can be found on our [species descriptions and genomic resources wiki page](https://github.com/urol-e5/deep-dive-expression/wiki).

## Analysis Workflows

The pipeline includes the following major analysis components:

### 1. Data Quality Control and Processing
- FastQC quality assessment of raw sequencing reads
- MultiQC summary reports
- Adapter trimming and quality filtering (fastp, cutadapt)
- Post-trimming quality validation

### 2. Reference Genome Annotation
- Functional annotation using BLAST and InterProScan
- Gene Ontology (GO) term assignment
- Repeat element identification (RepeatMasker)
- UTR annotation

### 3. RNA-seq Analysis
- Read alignment with HISAT2
- Gene expression quantification with StringTie
- Differential expression analysis with DESeq2
- Co-expression network construction (WGCNA)

### 4. Small RNA Analysis
- Small RNA discovery and annotation (ShortStack 4.1.0)
- miRNA identification and naming
- siRNA cluster detection
- Small RNA expression profiling

### 5. Long Non-coding RNA Analysis
- lncRNA identification pipeline
- Coding potential assessment (CPC2)
- lncRNA expression quantification
- lncRNA count matrix generation

### 6. DNA Methylation Analysis
- WGBS alignment with Bismark
- Methylation calling and quantification
- Differential methylation analysis (methylKit)
- Context-specific methylation patterns (CpG, CHG, CHH)

### 7. Interaction Analysis
- miRNA-mRNA target prediction (miRanda, RNAhybrid)
- miRNA-lncRNA interaction analysis
- siRNA-mRNA interaction mapping
- Correlation analysis (Pearson correlation coefficients)

### 8. Network Analysis
- Regulatory network construction
- Network visualization and export (GML format)
- Hub gene identification
- Module detection

### 9. Functional Enrichment
- GO term enrichment analysis
- Pathway analysis
- Functional annotation of interaction targets

### 10. Comparative Analysis (M-multi-species)
- Cross-species miRNA comparison
- Orthology analysis (OrthoFinder)
- Conserved target identification
- Comparative expression analysis

## Datasets

The project analyzes three types of high-throughput sequencing data:

### RNA-seq (mRNA)
- **Purpose**: Gene expression profiling and differential expression analysis
- **Platform**: Illumina paired-end sequencing
- **Processing**: Quality filtering → HISAT2 alignment → StringTie quantification

### Small RNA-seq
- **Purpose**: miRNA and siRNA discovery and quantification
- **Platform**: Illumina single-end sequencing
- **Processing**: Quality filtering → ShortStack analysis → Target prediction

### Whole Genome Bisulfite Sequencing (WGBS)
- **Purpose**: DNA methylation profiling
- **Platform**: Illumina paired-end sequencing
- **Processing**: Quality filtering → Bismark alignment → Methylation calling

Access to raw FASTQ files and aligned BAM files is provided through the [project wiki](https://github.com/urol-e5/deep-dive-expression/wiki).

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
- **RNAhybrid** (v2.1+): RNA-RNA interaction prediction

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

### Installation

Most tools can be installed using conda/mamba:

```bash
# Create a conda environment
conda create -n deep-dive python=3.9

# Install core tools
conda install -c bioconda fastqc multiqc hisat2 stringtie samtools bedtools bismark shortstack

# Install target prediction tools
conda install -c bioconda miranda rnahybrid

# Install R and Bioconductor in a separate environment
conda create -n deep-dive-r r-base=4.2 r-tidyverse bioconductor-deseq2
```

Refer to individual script headers for specific version requirements and conda environment configurations.

## Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/urol-e5/deep-dive-expression.git
cd deep-dive-expression
```

### 2. Set Up Your Environment

Install required software (see [Prerequisites](#prerequisites)) or use the provided conda environments specified in individual analysis scripts.

### 3. Download Reference Data

Reference genomes and annotation files are available through the [project wiki](https://github.com/urol-e5/deep-dive-expression/wiki). Download the appropriate files for your species of interest:

```bash
# Example for D-Apul
cd D-Apul/data
# Follow wiki instructions to download genome and annotation files
```

### 4. Access Sequencing Data

Raw FASTQ files and aligned BAM files are linked in the [project wiki](https://github.com/urol-e5/deep-dive-expression/wiki). Large data files are not stored in the repository to keep it lightweight.

### 5. Run Analyses

Navigate to the species directory and run scripts in numerical order:

```bash
cd D-Apul/code
# Open and run scripts in RStudio or from command line
```

## Usage Examples

### Example 1: Quality Control of RNA-seq Data

```bash
# Navigate to code directory
cd D-Apul/code

# Run FastQC analysis (open in RStudio or run with R)
R -e "rmarkdown::render('01-Apul-RNA-trimming-FastQC.Rmd')"

# View MultiQC report in output directory
firefox ../output/01-Apul-RNA-trimming-FastQC/multiqc_report.html
```

### Example 2: Small RNA Discovery with ShortStack

```bash
# Set up ShortStack conda environment (see script for details)
conda activate ShortStack-4.1.0_env

# Navigate to code directory
cd D-Apul/code

# Run ShortStack analysis
R -e "rmarkdown::render('11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome.Rmd')"

# Results will be in ../output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/
```

### Example 3: miRNA Target Prediction

```bash
cd D-Apul/code

# Run miRNA-mRNA interaction analysis
R -e "rmarkdown::render('16-Apul-RNAhybrid.Rmd')"

# View interaction predictions
head ../output/16-Apul-RNAhybrid/Apul-RNAhybrid-mRNA-results.txt
```

### Example 4: Cross-Species Comparison

```bash
cd M-multi-species/code

# Compare miRNAs across species
R -e "rmarkdown::render('04-miRNA-comparison.Rmd')"

# View comparison results
ls -l ../output/04-miRNA-comparison/
```

## File Naming Conventions

To maintain consistency and organization:

### Code Files
- **Format**: `XX-SpeciesCode-Description.Rmd` or `.qmd`
- **Numbering**: Two-digit prefix (e.g., `01`, `02`, `11`)
  - `00-0X`: Quality control and preprocessing
  - `01-09`: Core analyses (alignment, quantification)
  - `10-19`: Advanced analyses (interactions, networks)
  - `20+`: Specialized analyses

**Examples**:
- `01-Apul-RNA-trimming-FastQC.Rmd`
- `11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome.Rmd`
- `16-Apul-RNAhybrid.Rmd`

### Output Directories
- **Format**: Match the code file name without extension
- **Location**: `output/XX-SpeciesCode-Description/`

**Example**: Output from `01-Apul-RNA-trimming-FastQC.Rmd` → `output/01-Apul-RNA-trimming-FastQC/`

### Best Practices
- Use **relative paths** in all scripts (e.g., `../data/`, `../output/`)
- Commit and push changes frequently
- Document any deviations from standard naming
- Keep script numbers sequential when adding new analyses

## Software Requirements

### Conda Environments

Many scripts use conda/mamba environments. Example setup for ShortStack:

```bash
# Create environment
mamba create -n ShortStack-4.1.0_env -c bioconda shortstack=4.1.0

# Activate environment
conda activate ShortStack-4.1.0_env

# Find conda path (needed for R scripts)
which conda
# Example output: /home/user/mambaforge/condabin/conda
```

Update the conda path in the R script before running:

```r
# In the R script
shortstack_conda_env_name <- c("ShortStack-4.1.0_env")
shortstack_cond_path <- c("/home/user/mambaforge/condabin/conda")
```

### Version Information

Each script includes specific version requirements in its header. Refer to individual files for detailed software versions and parameters used.

### Common Issues

- **Memory errors**: Increase available RAM or reduce number of parallel processes
- **Path issues**: Ensure all paths are relative and reference genomes are downloaded
- **Conda environment**: Activate the correct environment before running scripts
- **R package versions**: Use BiocManager for Bioconductor packages to ensure compatibility

## Data Access

### Large Data Files

Due to size constraints, large data files are not stored in this repository:

- **Raw sequencing data** (FASTQ files): Linked in [project wiki](https://github.com/urol-e5/deep-dive-expression/wiki)
- **Aligned reads** (BAM files): Linked in [project wiki](https://github.com/urol-e5/deep-dive-expression/wiki)
- **Reference genomes**: Download instructions in [wiki](https://github.com/urol-e5/deep-dive-expression/wiki)

### `.gitignore` Configuration

The following file types are excluded from version control:
- Sequencing data: `*.fastq`, `*.fastq.gz`, `*.fq`, `*.fq.gz`
- Alignments: `*.bam`, `*.sam`
- Genome files: `*.fa`, `*.fasta`
- Intermediate files: Build directories, cache files

See `.gitignore` for the complete list.

## Contributing

We welcome contributions to improve analyses, add new features, or fix bugs!

### How to Contribute

1. **Fork the repository**
   ```bash
   # Click 'Fork' button on GitHub
   ```

2. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

3. **Make your changes**
   - Follow existing code style and naming conventions
   - Add documentation to your scripts
   - Test your changes thoroughly

4. **Commit your changes**
   ```bash
   git add .
   git commit -m "Description of changes"
   ```

5. **Push to your fork**
   ```bash
   git push origin feature/your-feature-name
   ```

6. **Create a Pull Request**
   - Go to the original repository
   - Click "New Pull Request"
   - Describe your changes and their purpose

### Contribution Guidelines

- **Code style**: Follow existing R/Rmd formatting conventions
- **Documentation**: Include comments and markdown explanations
- **File naming**: Use the established two-digit prefix system
- **Relative paths**: Always use relative paths (e.g., `../data/`, `../output/`)
- **Testing**: Verify that your code runs without errors
- **Commit messages**: Write clear, descriptive commit messages
- **Pull requests**: Reference relevant issues and provide detailed descriptions

### Reporting Issues

Found a bug or have a suggestion? Please [open an issue](https://github.com/urol-e5/deep-dive-expression/issues) with:
- Clear description of the problem or suggestion
- Steps to reproduce (for bugs)
- Expected vs. actual behavior
- System information (OS, R version, package versions)

### Code of Conduct

Please be respectful and constructive in all interactions. We aim to maintain a welcoming and inclusive environment for all contributors.

## Citation

If you use this pipeline or data in your research, please cite:

```bibtex
@misc{deepdive-expression,
  title = {deep-dive-expression: Multi-omics analysis of coral gene expression and regulatory networks},
  author = {E5 Coral Research Group},
  year = {2024},
  url = {https://github.com/urol-e5/deep-dive-expression},
  note = {GitHub repository}
}
```

### Key Software Citations

Please also cite the key software tools used:

- **ShortStack**: Axtell, M.J. (2013). ShortStack: comprehensive annotation and quantification of small RNA genes. RNA, 19(6), 740-751.
- **HISAT2**: Kim, D., et al. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature Biotechnology, 37(8), 907-915.
- **Bismark**: Krueger, F., & Andrews, S.R. (2011). Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. Bioinformatics, 27(11), 1571-1572.
- **DESeq2**: Love, M.I., et al. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550.

See `references.bib` for additional citations.

## Related Projects

- **[deep-dive](https://github.com/urol-e5/deep-dive)**: Predecessor project focusing on ncRNA landscape characterization
- **E5 Coral Project**: Broader coral resilience and epigenetics research initiative

## Contact

For questions, issues, or collaborations:

- **Issues**: [GitHub Issues](https://github.com/urol-e5/deep-dive-expression/issues)
- **Wiki**: [Project Wiki](https://github.com/urol-e5/deep-dive-expression/wiki)
- **E5 Coral Project**: Visit the main E5 project page for team information

---

**Last Updated**: October 2024

**Repository maintained by**: E5 Coral Research Group
