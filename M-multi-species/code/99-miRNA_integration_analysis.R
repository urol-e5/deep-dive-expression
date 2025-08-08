# miRNA Integration Analysis for Deep-Dive Expression Project
# This script integrates miRNA results across three coral species:
# D-Apul (Acropora pulchra)
# E-Peve (Porites evermanni) 
# F-Ptuh (Pocillopora tuahiniensis)

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(edgeR)
library(reshape2)
library(corrplot)
library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(readr)

# Set working directory and create output directory
setwd(".")
# Create output directory
output_dir <- "M-multi-species/output/99-miRNA-integration-analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to normalize counts using CPM
normalize_cpm <- function(counts_matrix) {
  dge <- DGEList(counts = counts_matrix)
  dge <- calcNormFactors(dge)
  cpm_matrix <- cpm(dge)
  return(cpm_matrix)
}

# Function to calculate correlation statistics
calc_correlation_stats <- function(x, y) {
  if (length(x) != length(y)) {
    return(c(correlation = NA, p_value = NA))
  }
  result <- cor.test(x, y, method = "pearson", use = "complete.obs")
  return(c(correlation = result$estimate, p_value = result$p.value))
}

# 1. LOAD AND PROCESS MIRNA EXPRESSION DATA
cat("Loading miRNA expression data...\n")

# Load miRNA count matrices for each species
apul_miRNA_counts <- read.table("D-Apul/output/03.1-Apul-sRNA-summary/Apul_miRNA_ShortStack_counts_formatted.txt", 
                               header = TRUE, row.names = 1)
peve_miRNA_counts <- read.table("E-Peve/output/03.1-Peve-sRNA-summary/Peve_miRNA_ShortStack_counts_formatted.txt", 
                               header = TRUE, row.names = 1)
ptuh_miRNA_counts <- read.table("F-Ptuh/output/03.1-Ptuh-sRNA-summary/Ptuh_miRNA_ShortStack_counts_formatted.txt", 
                               header = TRUE, row.names = 1)

# Normalize counts using CPM
apul_miRNA_cpm <- normalize_cpm(apul_miRNA_counts)
peve_miRNA_cpm <- normalize_cpm(peve_miRNA_counts)
ptuh_miRNA_cpm <- normalize_cpm(ptuh_miRNA_counts)

# Calculate mean expression for each miRNA
apul_miRNA_summary <- data.frame(
  miRNA = rownames(apul_miRNA_cpm),
  mean_expression = rowMeans(apul_miRNA_cpm),
  species = "A. pulchra",
  stringsAsFactors = FALSE
)

peve_miRNA_summary <- data.frame(
  miRNA = rownames(peve_miRNA_cpm),
  mean_expression = rowMeans(peve_miRNA_cpm),
  species = "P. evermanni",
  stringsAsFactors = FALSE
)

ptuh_miRNA_summary <- data.frame(
  miRNA = rownames(ptuh_miRNA_cpm),
  mean_expression = rowMeans(ptuh_miRNA_cpm),
  species = "P. tuahiniensis",
  stringsAsFactors = FALSE
)

# Combine all species data
all_miRNA_expression <- rbind(apul_miRNA_summary, peve_miRNA_summary, ptuh_miRNA_summary)

# 2. LOAD MIRNA-MRNA INTERACTION DATA
cat("Loading miRNA-mRNA interaction data...\n")

# Load miRNA-mRNA interactions for each species (using correct file names)
apul_mirna_mrna <- read.csv("D-Apul/output/09-Apul-mRNA-miRNA-interactions/miranda_PCC_miRNA_mRNA.csv")
peve_mirna_mrna <- read.csv("E-Peve/output/10-Peve-mRNA-miRNA-interactions/Peve-miranda_PCC_miRNA_mRNA.csv")
ptuh_mirna_mrna <- read.csv("F-Ptuh/output/11-Ptuh-mRNA-miRNA-interactions/three_prime_interaction/Ptuh-miranda_PCC_miRNA_mRNA.csv")

# Add species information
apul_mirna_mrna$species <- "A. pulchra"
peve_mirna_mrna$species <- "P. evermanni"
ptuh_mirna_mrna$species <- "P. tuahiniensis"

# Combine interaction data
all_mirna_mrna_interactions <- rbind(apul_mirna_mrna, peve_mirna_mrna, ptuh_mirna_mrna)

# 3. LOAD MIRNA-LNCRNA INTERACTION DATA
cat("Loading miRNA-lncRNA interaction data...\n")

# Load miRNA-lncRNA interactions for each species
apul_mirna_lncrna <- read.csv("D-Apul/output/28-Apul-miRNA-lncRNA-interactions/miranda_PCC_miRNA_lncRNA.csv")
peve_mirna_lncrna <- read.csv("E-Peve/output/15-Peve-miRNA-lncRNA-PCC/miranda_PCC_miRNA_lncRNA.csv")
ptuh_mirna_lncrna <- read.csv("F-Ptuh/output/15-Ptuh-miRNA-lncRNA-PCC/miranda_PCC_miRNA_lncRNA.csv")

# Add species information
apul_mirna_lncrna$species <- "A. pulchra"
peve_mirna_lncrna$species <- "P. evermanni"
ptuh_mirna_lncrna$species <- "P. tuahiniensis"

# Combine interaction data
all_mirna_lncrna_interactions <- rbind(apul_mirna_lncrna, peve_mirna_lncrna, ptuh_mirna_lncrna)

# 4. GENERATE VISUALIZATIONS

# 4.1 miRNA Expression Distribution Across Species
cat("Generating miRNA expression distribution plot...\n")

expression_plot <- ggplot(all_miRNA_expression, aes(x = species, y = log2(mean_expression + 1), fill = species)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(title = "miRNA Expression Distribution Across Coral Species",
       x = "Species",
       y = "Log2(Mean Expression + 1)",
       fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "miRNA_expression_distribution.png"), 
       expression_plot, width = 10, height = 6, dpi = 300)

# 4.2 miRNA Expression Heatmap
cat("Generating miRNA expression heatmap...\n")

# Prepare data for heatmap
expression_matrix <- all_miRNA_expression %>%
  pivot_wider(names_from = species, values_from = mean_expression, values_fill = 0)

# Convert to matrix with miRNA names as row names
expression_matrix_mat <- as.matrix(expression_matrix[, -1])
rownames(expression_matrix_mat) <- expression_matrix$miRNA

# Log transform for better visualization
expression_matrix_log <- log2(expression_matrix_mat + 1)

# Create heatmap
heatmap_plot <- ggplot(reshape2::melt(expression_matrix_log), 
                       aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_viridis(option = "magma", name = "Log2(Expression + 1)") +
  labs(title = "miRNA Expression Heatmap Across Species",
       x = "Species",
       y = "miRNA") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "miRNA_expression_heatmap.png"), 
       heatmap_plot, width = 8, height = 12, dpi = 300)

# 4.3 miRNA-mRNA Interaction Network Analysis
cat("Analyzing miRNA-mRNA interaction networks...\n")

# Calculate interaction statistics for each species
interaction_stats <- all_mirna_mrna_interactions %>%
  group_by(species) %>%
  summarise(
    total_interactions = n(),
    unique_mirnas = n_distinct(miRNA),
    unique_mrnas = n_distinct(mRNA),
    mean_correlation = mean(PCC.cor, na.rm = TRUE),
    significant_interactions = sum(p_value < 0.05, na.rm = TRUE),
    .groups = 'drop'
  )

# Plot interaction statistics
interaction_stats_plot <- ggplot(interaction_stats, aes(x = species, y = total_interactions, fill = species)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_viridis(discrete = TRUE, option = "C") +
  labs(title = "miRNA-mRNA Interaction Statistics by Species",
       x = "Species",
       y = "Total Interactions",
       fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "miRNA_mRNA_interaction_stats.png"), 
       interaction_stats_plot, width = 10, height = 6, dpi = 300)

# 4.4 Correlation Distribution Analysis
cat("Analyzing correlation distributions...\n")

# Filter significant interactions
significant_interactions <- all_mirna_mrna_interactions %>%
  filter(p_value < 0.05)

correlation_plot <- ggplot(significant_interactions, aes(x = PCC.cor, fill = species)) +
  geom_density(alpha = 0.6) +
  scale_fill_viridis(discrete = TRUE, option = "B") +
  labs(title = "Distribution of Significant miRNA-mRNA Correlations",
       x = "Pearson Correlation Coefficient",
       y = "Density",
       fill = "Species") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "correlation_distribution.png"), 
       correlation_plot, width = 10, height = 6, dpi = 300)

# 4.5 Top miRNA Analysis
cat("Analyzing top miRNAs by expression and interactions...\n")

# Calculate top miRNAs by expression
top_mirnas_expression <- all_miRNA_expression %>%
  group_by(species) %>%
  top_n(10, mean_expression) %>%
  arrange(species, desc(mean_expression))

# Calculate top miRNAs by number of interactions
top_mirnas_interactions <- all_mirna_mrna_interactions %>%
  group_by(species, miRNA) %>%
  summarise(
    interaction_count = n(),
    mean_correlation = mean(PCC.cor, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  group_by(species) %>%
  top_n(10, interaction_count) %>%
  arrange(species, desc(interaction_count))

# Plot top miRNAs by expression
top_expression_plot <- ggplot(top_mirnas_expression, 
                             aes(x = reorder(miRNA, mean_expression), y = mean_expression, fill = species)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  facet_wrap(~species, scales = "free_y") +
  coord_flip() +
  labs(title = "Top 10 miRNAs by Expression Level",
       x = "miRNA",
       y = "Mean Expression (CPM)",
       fill = "Species") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 8))

ggsave(file.path(output_dir, "top_mirnas_expression.png"), 
       top_expression_plot, width = 12, height = 8, dpi = 300)

# 4.6 Network Analysis
cat("Performing network analysis...\n")

# Create a simplified network for visualization
# Focus on top miRNAs and their interactions
top_network_data <- all_mirna_mrna_interactions %>%
  filter(p_value < 0.05) %>%
  group_by(miRNA) %>%
  filter(n() >= 5) %>%  # Only miRNAs with at least 5 interactions
  ungroup()

# Check if we have data for network analysis
if (nrow(top_network_data) > 0) {
  cat("Creating network with", nrow(top_network_data), "interactions\n")
  
  # Create edges using base R to avoid namespace conflicts
  edges <- data.frame(
    from = top_network_data$miRNA,
    to = top_network_data$mRNA,
    weight = abs(top_network_data$PCC.cor),
    species = top_network_data$species,
    stringsAsFactors = FALSE
  )
  
  # Create nodes
  all_nodes <- unique(c(edges$from, edges$to))
  nodes <- data.frame(
    id = all_nodes,
    type = ifelse(all_nodes %in% edges$from, "miRNA", "mRNA"),
    stringsAsFactors = FALSE
  )
  
  # Create graph
  graph <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  
  # Calculate network metrics
  network_metrics <- data.frame(
    metric = c("Number of nodes", "Number of edges", "Average degree", "Density"),
    value = c(
      vcount(graph),
      ecount(graph),
      mean(degree(graph)),
      edge_density(graph)
    )
  )
  
  # Save network metrics
  write.csv(network_metrics, file.path(output_dir, "network_metrics.csv"), row.names = FALSE)
  
  cat("Network analysis completed with", vcount(graph), "nodes and", ecount(graph), "edges\n")
} else {
  cat("No significant interactions found for network analysis\n")
}

# 4.7 Species Comparison Analysis
cat("Performing species comparison analysis...\n")

# Compare miRNA expression patterns between species
species_comparison <- all_miRNA_expression %>%
  group_by(species) %>%
  summarise(
    total_mirnas = n(),
    mean_expression = mean(mean_expression),
    median_expression = median(mean_expression),
    expression_variance = var(mean_expression),
    .groups = 'drop'
  )

# Plot species comparison
species_comparison_plot <- ggplot(species_comparison, aes(x = species, y = mean_expression, fill = species)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_viridis(discrete = TRUE, option = "A") +
  labs(title = "Average miRNA Expression by Species",
       x = "Species",
       y = "Mean Expression (CPM)",
       fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "species_comparison.png"), 
       species_comparison_plot, width = 10, height = 6, dpi = 300)

# 4.8 Functional Enrichment Analysis
cat("Performing functional enrichment analysis...\n")

# Analyze miRNA targets for functional patterns
target_analysis <- all_mirna_mrna_interactions %>%
  filter(p_value < 0.05) %>%
  group_by(miRNA, species) %>%
  summarise(
    target_count = n(),
    mean_correlation = mean(PCC.cor, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  group_by(species) %>%
  summarise(
    avg_targets_per_mirna = mean(target_count),
    total_targets = sum(target_count),
    .groups = 'drop'
  )

# Plot target analysis
target_analysis_plot <- ggplot(target_analysis, aes(x = species, y = avg_targets_per_mirna, fill = species)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_viridis(discrete = TRUE, option = "E") +
  labs(title = "Average Number of Targets per miRNA by Species",
       x = "Species",
       y = "Average Targets per miRNA",
       fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "target_analysis.png"), 
       target_analysis_plot, width = 10, height = 6, dpi = 300)

# 5. GENERATE SUMMARY REPORT
cat("Generating summary report...\n")

# Create summary statistics
summary_stats <- list(
  total_mirnas = nrow(all_miRNA_expression),
  total_interactions = nrow(all_mirna_mrna_interactions),
  significant_interactions = nrow(significant_interactions),
  species_count = length(unique(all_miRNA_expression$species)),
  expression_range_min = min(all_miRNA_expression$mean_expression),
  expression_range_max = max(all_miRNA_expression$mean_expression),
  correlation_range_min = min(all_mirna_mrna_interactions$PCC.cor, na.rm = TRUE),
  correlation_range_max = max(all_mirna_mrna_interactions$PCC.cor, na.rm = TRUE)
)

# Save summary statistics
summary_df <- data.frame(
  metric = names(summary_stats),
  value = unlist(summary_stats),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)

# Create comprehensive summary plot
summary_plot <- (expression_plot + correlation_plot) / 
                (interaction_stats_plot + species_comparison_plot) +
  plot_annotation(
    title = "Comprehensive miRNA Integration Analysis",
    subtitle = "Analysis of miRNA expression and interactions across three coral species",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                 plot.subtitle = element_text(size = 12, hjust = 0.5))
  )

ggsave(file.path(output_dir, "comprehensive_summary.png"), 
       summary_plot, width = 16, height = 12, dpi = 300)

# 6. SAVE PROCESSED DATA
cat("Saving processed data...\n")

# Save processed data for further analysis
write.csv(all_miRNA_expression, file.path(output_dir, "all_miRNA_expression.csv"), row.names = FALSE)
write.csv(all_mirna_mrna_interactions, file.path(output_dir, "all_mirna_mrna_interactions.csv"), row.names = FALSE)
write.csv(all_mirna_lncrna_interactions, file.path(output_dir, "all_mirna_lncrna_interactions.csv"), row.names = FALSE)
write.csv(significant_interactions, file.path(output_dir, "significant_interactions.csv"), row.names = FALSE)

cat("Analysis complete! Results saved in 'M-multi-species/output/99-miRNA-integration-analysis' directory.\n")

# 7. PRINT KEY FINDINGS
cat("\n=== KEY FINDINGS ===\n")
cat("1. Total miRNAs analyzed:", summary_stats$total_mirnas, "\n")
cat("2. Total interactions found:", summary_stats$total_interactions, "\n")
cat("3. Significant interactions (p < 0.05):", summary_stats$significant_interactions, "\n")
cat("4. Species analyzed:", summary_stats$species_count, "\n")
cat("5. Expression range (CPM):", round(summary_stats$expression_range_min, 2), "to", round(summary_stats$expression_range_max, 2), "\n")
cat("6. Correlation range:", round(summary_stats$correlation_range_min, 3), "to", round(summary_stats$correlation_range_max, 3), "\n")

# Print top miRNAs by expression
cat("\n=== TOP MIRNAS BY EXPRESSION ===\n")
for (species in unique(all_miRNA_expression$species)) {
  cat("\n", species, ":\n")
  top_species <- all_miRNA_expression %>%
    filter(species == !!species) %>%
    arrange(desc(mean_expression)) %>%
    head(5)
  print(top_species[, c("miRNA", "mean_expression")])
}

cat("\nAnalysis completed successfully!\n")
