# miRNA Integration and Network Analysis

This directory contains R scripts for comprehensive miRNA analysis across three coral species: *Acropora pulchra* (D-Apul), *Porites evermanni* (E-Peve), and *Pocillopora tuahiniensis* (F-Ptuh).

## Scripts

### 1. `99-miRNA_integration_analysis.R`
**Purpose**: Integrates miRNA expression and interaction data across all three species and generates comprehensive visualizations.

**Key Features**:
- Loads miRNA expression data from ShortStack results
- Integrates miRNA-mRNA and miRNA-lncRNA interaction data
- Performs CPM normalization for expression data
- Generates expression distribution plots, heatmaps, and correlation analyses
- Creates species comparison visualizations
- Saves processed data for further analysis

**Outputs**: Saved to `../output/99-miRNA-integration-analysis/`
- Expression plots and heatmaps
- Interaction statistics
- Correlation distributions
- Species comparisons
- Processed data files (CSV format)

### 2. `99-miRNA_network_analysis.R`
**Purpose**: Performs advanced network analysis on miRNA interaction data.

**Key Features**:
- Creates interaction networks for each species
- Calculates network metrics (density, clustering, centrality)
- Identifies hub miRNAs based on multiple centrality measures
- Performs community detection
- Generates network visualizations
- Prepares data for interactive visualization tools

**Outputs**: Saved to `../output/99-miRNA-integration-analysis/`
- Network plots for each species
- Network comparison metrics
- Hub miRNA identification
- Interactive network data files

### 3. `99-install_packages.R`
**Purpose**: Installs all required R packages for the analysis.

**Required Packages**:
- **CRAN**: dplyr, tidyr, ggplot2, stringr, reshape2, corrplot, igraph, ggraph, tidygraph, network, sna, intergraph, patchwork, viridis, RColorBrewer, readr
- **Bioconductor**: edgeR, limma

## Data Sources

The scripts expect the following data structure:

```
D-Apul/output/09-Apul-mRNA-miRNA-interactions/miranda_PCC_miRNA_mRNA.csv
D-Apul/output/28-Apul-miRNA-lncRNA-interactions/miranda_PCC_miRNA_lncRNA.csv
D-Apul/output/03.1-Apul-sRNA-summary/Apul_miRNA_ShortStack_counts_formatted.txt

E-Peve/output/10-Peve-mRNA-miRNA-interactions/Peve-miranda_PCC_miRNA_mRNA.csv
E-Peve/output/15-Peve-miRNA-lncRNA-PCC/miranda_PCC_miRNA_lncRNA.csv
E-Peve/output/03.1-Peve-sRNA-summary/Peve_miRNA_ShortStack_counts_formatted.txt

F-Ptuh/output/11-Ptuh-mRNA-miRNA-interactions/three_prime_interaction/Ptuh-miranda_PCC_miRNA_mRNA.csv
F-Ptuh/output/15-Ptuh-miRNA-lncRNA-PCC/miranda_PCC_miRNA_lncRNA.csv
F-Ptuh/output/03.1-Ptuh-sRNA-summary/Ptuh_miRNA_ShortStack_counts_formatted.txt
```

## Installation

1. **Install R packages**:
   ```r
   source("install_packages.R")
   ```

2. **Verify data files exist** in the expected locations

## Usage

### Run Integration Analysis
```r
source("99-miRNA_integration_analysis.R")
```

### Run Network Analysis
```r
source("99-miRNA_network_analysis.R")
```

## Output Structure

All outputs are saved to `../output/99-miRNA-integration-analysis/`:

- **Visualizations**: PNG plots for expression, interactions, networks
- **Data Files**: CSV files with processed data and statistics
- **Network Data**: Files ready for interactive visualization tools

## Key Findings

### Expression Analysis
- **Total miRNAs analyzed**: 121 across all species
- **Expression range**: 5.27 to 642,916.9 CPM
- **Species-specific patterns**: Each species shows unique miRNA expression profiles

### Network Analysis
- **A. pulchra**: 274 nodes, 253 edges, density 0.007
- **P. evermanni**: 202 nodes, 182 edges, density 0.009
- **P. tuahiniensis**: 205 nodes, 194 edges, density 0.009

### Hub miRNAs
Each species has 3-4 key hub miRNAs with high connectivity:
- **A. pulchra**: Cluster_5981 (24 connections)
- **P. evermanni**: Cluster_9149 (28 connections)
- **P. tuahiniensis**: Cluster_2973 & Cluster_1938 (30 connections each)

## Advanced Analysis Concepts

### Network Metrics
- **Density**: Measures network connectivity (0-1 scale)
- **Clustering**: Indicates functional module formation
- **Centrality**: Identifies key regulatory nodes
- **Community Detection**: Reveals functional clusters

### Hub Identification
Hubs are identified using multiple centrality measures:
- **Degree**: Number of direct connections
- **Betweenness**: Control over information flow
- **Closeness**: Average distance to other nodes
- **Eigenvector**: Influence based on neighbors' importance

## Troubleshooting

### Common Issues
1. **Package conflicts**: Use `install_packages.R` to ensure correct versions
2. **File path errors**: Verify data files exist in expected locations
3. **Memory issues**: Large datasets may require increased memory allocation

### Error Messages
- **"cannot open file"**: Check file paths and permissions
- **"object not found"**: Ensure data loading completed successfully
- **Package conflicts**: Restart R session after package installation

## Future Enhancements

- **Interactive visualizations**: Integration with Cytoscape or Gephi
- **Functional enrichment**: GO term analysis for target genes
- **Time series analysis**: Expression changes across conditions
- **Machine learning**: Predictive modeling of miRNA-target interactions

## Citation

If you use these scripts in your research, please cite:
- R packages: igraph, tidygraph, ggplot2
- Analysis methods: Network analysis, Centrality measures
- Data sources: ShortStack miRNA prediction, miRanda target prediction

## Contact

For questions or issues with these scripts, please refer to the project documentation or contact the development team.

---

**Note**: This analysis integrates data from multiple coral species to provide comparative insights into miRNA regulatory networks. Ensure all data files are properly formatted and accessible before running the scripts.

