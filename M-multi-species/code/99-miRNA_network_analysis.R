# Advanced miRNA Network Analysis
# This script performs detailed network analysis of miRNA regulatory networks
# across three coral species

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(tidygraph)
library(network)
library(sna)
library(intergraph)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(readr)

# Set working directory
setwd(".")

# Create output directory
output_dir <- "M-multi-species/output/99-miRNA-integration-analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Function to create and analyze networks
analyze_network <- function(interaction_data, species_name, min_interactions = 3) {
  
  # Filter for significant interactions and minimum interaction count
  filtered_data <- interaction_data %>%
    filter(p_value < 0.05) %>%
    group_by(miRNA) %>%
    filter(n() >= min_interactions) %>%
    ungroup()
  
  if (nrow(filtered_data) == 0) {
    cat("No significant interactions found for", species_name, "\n")
    return(NULL)
  }
  
  # Create edges using base R to avoid namespace conflicts
  edges <- data.frame(
    from = filtered_data$miRNA,
    to = filtered_data$mRNA,
    weight = abs(filtered_data$PCC.cor),
    correlation = filtered_data$PCC.cor,
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
  
  # Calculate network metrics using igraph functions
  metrics <- list(
    nodes = vcount(graph),
    edges = ecount(graph),
    density = edge_density(graph),
    avg_degree = mean(igraph::degree(graph)),
    avg_clustering = transitivity(graph, type = "global"),
    diameter = diameter(graph),
    avg_path_length = mean_distance(graph)
  )
  
  # Calculate centrality measures using igraph functions
  centrality <- data.frame(
    id = V(graph)$name,
    degree = igraph::degree(graph),
    betweenness = igraph::betweenness(graph),
    closeness = igraph::closeness(graph),
    eigenvector = eigen_centrality(graph)$vector,
    stringsAsFactors = FALSE
  )
  
  return(list(
    graph = graph,
    metrics = metrics,
    centrality = centrality,
    edges = edges,
    nodes = nodes
  ))
}

# Function to create network visualization
create_network_plot <- function(network_data, species_name) {
  if (is.null(network_data)) return(NULL)
  
  # Create tidygraph object
  tidy_graph <- as_tbl_graph(network_data$graph)
  
  # Add node attributes
  tidy_graph <- tidy_graph %>%
    activate(nodes) %>%
    mutate(
      type = ifelse(name %in% network_data$edges$from, "miRNA", "mRNA"),
      degree = centrality_degree(),
      betweenness = centrality_betweenness(),
      community = group_louvain()
    )
  
  # Create plot
  network_plot <- ggraph(tidy_graph, layout = "fr") +
    geom_edge_link(aes(alpha = weight), color = "gray50") +
    geom_node_point(aes(size = degree, color = type, shape = type)) +
    scale_color_manual(values = c("miRNA" = "#E64B35", "mRNA" = "#4DBBD5")) +
    scale_shape_manual(values = c("miRNA" = 16, "mRNA" = 17)) +
    scale_size_continuous(range = c(2, 8)) +
    labs(title = paste("miRNA-mRNA Network:", species_name),
         subtitle = paste("Nodes:", network_data$metrics$nodes, 
                         "Edges:", network_data$metrics$edges,
                         "Density:", round(network_data$metrics$density, 3))) +
    theme_graph() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  return(network_plot)
}

# 1. LOAD DATA
cat("Loading interaction data...\n")

# Load miRNA-mRNA interactions (using correct file names)
apul_interactions <- read.csv("D-Apul/output/09-Apul-mRNA-miRNA-interactions/miranda_PCC_miRNA_mRNA.csv")
peve_interactions <- read.csv("E-Peve/output/10-Peve-mRNA-miRNA-interactions/Peve-miranda_PCC_miRNA_mRNA.csv")
ptuh_interactions <- read.csv("F-Ptuh/output/11-Ptuh-mRNA-miRNA-interactions/three_prime_interaction/Ptuh-miranda_PCC_miRNA_mRNA.csv")

# Load miRNA-lncRNA interactions
apul_lncrna <- read.csv("D-Apul/output/28-Apul-miRNA-lncRNA-interactions/miranda_PCC_miRNA_lncRNA.csv")
peve_lncrna <- read.csv("E-Peve/output/15-Peve-miRNA-lncRNA-PCC/miranda_PCC_miRNA_lncRNA.csv")
ptuh_lncrna <- read.csv("F-Ptuh/output/15-Ptuh-miRNA-lncRNA-PCC/miranda_PCC_miRNA_lncRNA.csv")

# 2. NETWORK ANALYSIS
cat("Performing network analysis...\n")

# Analyze networks for each species
apul_network <- analyze_network(apul_interactions, "A. pulchra")
peve_network <- analyze_network(peve_interactions, "P. evermanni")
ptuh_network <- analyze_network(ptuh_interactions, "P. tuahiniensis")

# 3. CREATE NETWORK VISUALIZATIONS
cat("Creating network visualizations...\n")

# Create network plots
apul_plot <- create_network_plot(apul_network, "A. pulchra")
peve_plot <- create_network_plot(peve_network, "P. evermanni")
ptuh_plot <- create_network_plot(ptuh_network, "P. tuahiniensis")

# Save individual network plots
if (!is.null(apul_plot)) {
  ggsave(file.path(output_dir, "apul_network.png"), apul_plot, width = 12, height = 10, dpi = 300)
}
if (!is.null(peve_plot)) {
  ggsave(file.path(output_dir, "peve_network.png"), peve_plot, width = 12, height = 10, dpi = 300)
}
if (!is.null(ptuh_plot)) {
  ggsave(file.path(output_dir, "ptuh_network.png"), ptuh_plot, width = 12, height = 10, dpi = 300)
}

# 4. COMPARATIVE NETWORK ANALYSIS
cat("Performing comparative network analysis...\n")

# Combine network metrics for comparison
network_comparison <- data.frame(
  species = c("A. pulchra", "P. evermanni", "P. tuahiniensis"),
  nodes = c(
    ifelse(is.null(apul_network), 0, apul_network$metrics$nodes),
    ifelse(is.null(peve_network), 0, peve_network$metrics$nodes),
    ifelse(is.null(ptuh_network), 0, ptuh_network$metrics$nodes)
  ),
  edges = c(
    ifelse(is.null(apul_network), 0, apul_network$metrics$edges),
    ifelse(is.null(peve_network), 0, peve_network$metrics$edges),
    ifelse(is.null(ptuh_network), 0, ptuh_network$metrics$edges)
  ),
  density = c(
    ifelse(is.null(apul_network), 0, apul_network$metrics$density),
    ifelse(is.null(peve_network), 0, peve_network$metrics$density),
    ifelse(is.null(ptuh_network), 0, ptuh_network$metrics$density)
  ),
  avg_degree = c(
    ifelse(is.null(apul_network), 0, apul_network$metrics$avg_degree),
    ifelse(is.null(peve_network), 0, peve_network$metrics$avg_degree),
    ifelse(is.null(ptuh_network), 0, ptuh_network$metrics$avg_degree)
  ),
  avg_clustering = c(
    ifelse(is.null(apul_network), 0, apul_network$metrics$avg_clustering),
    ifelse(is.null(peve_network), 0, peve_network$metrics$avg_clustering),
    ifelse(is.null(ptuh_network), 0, ptuh_network$metrics$avg_clustering)
  )
)

# Create comparison plots
comparison_plot <- ggplot(network_comparison, aes(x = species, y = nodes, fill = species)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  labs(title = "Network Size Comparison",
       x = "Species",
       y = "Number of Nodes",
       fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

density_plot <- ggplot(network_comparison, aes(x = species, y = density, fill = species)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_viridis(discrete = TRUE, option = "E") +
  labs(title = "Network Density Comparison",
       x = "Species",
       y = "Network Density",
       fill = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Save comparison plots
ggsave(file.path(output_dir, "network_comparison.png"), comparison_plot, width = 10, height = 6, dpi = 300)
ggsave(file.path(output_dir, "network_density.png"), density_plot, width = 10, height = 6, dpi = 300)

# 5. CENTRALITY ANALYSIS
cat("Performing centrality analysis...\n")

# Analyze centrality measures for each species
centrality_analysis <- function(network_data, species_name) {
  if (is.null(network_data)) return(NULL)
  
  # Get top nodes by different centrality measures
  top_degree <- network_data$centrality %>%
    arrange(desc(degree)) %>%
    head(10)
  
  top_betweenness <- network_data$centrality %>%
    arrange(desc(betweenness)) %>%
    head(10)
  
  top_closeness <- network_data$centrality %>%
    arrange(desc(closeness)) %>%
    head(10)
  
  return(list(
    top_degree = top_degree,
    top_betweenness = top_betweenness,
    top_closeness = top_closeness
  ))
}

# Perform centrality analysis for each species
apul_centrality <- centrality_analysis(apul_network, "A. pulchra")
peve_centrality <- centrality_analysis(peve_network, "P. evermanni")
ptuh_centrality <- centrality_analysis(ptuh_network, "P. tuahiniensis")

# 6. COMMUNITY DETECTION
cat("Performing community detection...\n")

detect_communities <- function(network_data, species_name) {
  if (is.null(network_data)) return(NULL)
  
  # Detect communities using Louvain method
  communities <- cluster_louvain(network_data$graph)
  
  # Get community information
  community_info <- data.frame(
    node = V(network_data$graph)$name,
    community = communities$membership,
    stringsAsFactors = FALSE
  )
  
  # Calculate community statistics
  community_stats <- community_info %>%
    group_by(community) %>%
    summarise(
      size = n(),
      .groups = 'drop'
    ) %>%
    arrange(desc(size))
  
  return(list(
    communities = communities,
    community_info = community_info,
    community_stats = community_stats
  ))
}

# Detect communities for each species
apul_communities <- detect_communities(apul_network, "A. pulchra")
peve_communities <- detect_communities(peve_network, "P. evermanni")
ptuh_communities <- detect_communities(ptuh_network, "P. tuahiniensis")

# 7. HUB MIRNA ANALYSIS
cat("Analyzing hub miRNAs...\n")

# Identify hub miRNAs (high degree centrality)
identify_hub_mirnas <- function(network_data, species_name) {
  if (is.null(network_data)) return(NULL)
  
  # Filter for miRNAs only
  mirna_centrality <- network_data$centrality %>%
    filter(id %in% network_data$edges$from) %>%
    arrange(desc(degree))
  
  # Identify hubs (top 10% by degree)
  hub_threshold <- quantile(mirna_centrality$degree, 0.9)
  hub_mirnas <- mirna_centrality %>%
    filter(degree >= hub_threshold)
  
  return(hub_mirnas)
}

# Identify hub miRNAs for each species
apul_hubs <- identify_hub_mirnas(apul_network, "A. pulchra")
peve_hubs <- identify_hub_mirnas(peve_network, "P. evermanni")
ptuh_hubs <- identify_hub_mirnas(ptuh_network, "P. tuahiniensis")

# 8. CREATE COMPREHENSIVE REPORT
cat("Creating comprehensive report...\n")

# Save network metrics
write.csv(network_comparison, file.path(output_dir, "network_comparison_metrics.csv"), row.names = FALSE)

# Save hub miRNAs
if (!is.null(apul_hubs)) {
  write.csv(apul_hubs, file.path(output_dir, "apul_hub_mirnas.csv"), row.names = FALSE)
}
if (!is.null(peve_hubs)) {
  write.csv(peve_hubs, file.path(output_dir, "peve_hub_mirnas.csv"), row.names = FALSE)
}
if (!is.null(ptuh_hubs)) {
  write.csv(ptuh_hubs, file.path(output_dir, "ptuh_hub_mirnas.csv"), row.names = FALSE)
}

# 9. CREATE INTERACTIVE NETWORK VISUALIZATION
cat("Creating interactive network visualization...\n")

# Function to create interactive network data
create_interactive_network <- function(network_data, species_name) {
  if (is.null(network_data)) return(NULL)
  
  # Prepare nodes data
  nodes_data <- network_data$nodes %>%
    mutate(
      group = type,
      size = ifelse(type == "miRNA", 20, 10),
      color = ifelse(type == "miRNA", "#E64B35", "#4DBBD5")
    )
  
  # Prepare edges data
  edges_data <- network_data$edges %>%
    mutate(
      width = weight * 5,  # Scale weight for visualization
      color = "#666666"
    )
  
  return(list(
    nodes = nodes_data,
    edges = edges_data
  ))
}

# Create interactive network data for each species
apul_interactive <- create_interactive_network(apul_network, "A. pulchra")
peve_interactive <- create_interactive_network(peve_network, "P. evermanni")
ptuh_interactive <- create_interactive_network(ptuh_network, "P. tuahiniensis")

# Save interactive network data
if (!is.null(apul_interactive)) {
  write.csv(apul_interactive$nodes, file.path(output_dir, "apul_interactive_nodes.csv"), row.names = FALSE)
  write.csv(apul_interactive$edges, file.path(output_dir, "apul_interactive_edges.csv"), row.names = FALSE)
}
if (!is.null(peve_interactive)) {
  write.csv(peve_interactive$nodes, file.path(output_dir, "peve_interactive_nodes.csv"), row.names = FALSE)
  write.csv(peve_interactive$edges, file.path(output_dir, "peve_interactive_edges.csv"), row.names = FALSE)
}
if (!is.null(ptuh_interactive)) {
  write.csv(ptuh_interactive$nodes, file.path(output_dir, "ptuh_interactive_nodes.csv"), row.names = FALSE)
  write.csv(ptuh_interactive$edges, file.path(output_dir, "ptuh_interactive_edges.csv"), row.names = FALSE)
}

# 10. PRINT SUMMARY
cat("\n=== NETWORK ANALYSIS SUMMARY ===\n")
cat("A. pulchra network:\n")
if (!is.null(apul_network)) {
  cat("  Nodes:", apul_network$metrics$nodes, "\n")
  cat("  Edges:", apul_network$metrics$edges, "\n")
  cat("  Density:", round(apul_network$metrics$density, 3), "\n")
  cat("  Hub miRNAs:", nrow(apul_hubs), "\n")
} else {
  cat("  No significant network found\n")
}

cat("\nP. evermanni network:\n")
if (!is.null(peve_network)) {
  cat("  Nodes:", peve_network$metrics$nodes, "\n")
  cat("  Edges:", peve_network$metrics$edges, "\n")
  cat("  Density:", round(peve_network$metrics$density, 3), "\n")
  cat("  Hub miRNAs:", nrow(peve_hubs), "\n")
} else {
  cat("  No significant network found\n")
}

cat("\nP. tuahiniensis network:\n")
if (!is.null(ptuh_network)) {
  cat("  Nodes:", ptuh_network$metrics$nodes, "\n")
  cat("  Edges:", ptuh_network$metrics$edges, "\n")
  cat("  Density:", round(ptuh_network$metrics$density, 3), "\n")
  cat("  Hub miRNAs:", nrow(ptuh_hubs), "\n")
} else {
  cat("  No significant network found\n")
}

cat("\nNetwork analysis completed successfully!\n")

