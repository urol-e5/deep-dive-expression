31-Peve-miRNA-mRNA-lncRNA-network
================
Kathleen Durkin
2025-07-16

- <a href="#1-all-interactions" id="toc-1-all-interactions">1 All
  interactions</a>
- <a href="#2-pval--005" id="toc-2-pval--005">2 pval &lt; 0.05</a>
- <a href="#3-pval--001" id="toc-3-pval--001">3 pval &lt; 0.01</a>
- <a href="#4-save" id="toc-4-save">4 Save</a>
- <a href="#5-using-cytoscape" id="toc-5-using-cytoscape">5 Using
  Cytoscape</a>
- <a href="#6-plotting-with-igraph" id="toc-6-plotting-with-igraph">6
  Plotting with igraph</a>

Reran 10/09/2025 to incorporate changes due to updated lncRNA counts
matrix

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
library(tidyr)
library(igraph)
```

    ## 
    ## Attaching package: 'igraph'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     crossing

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

``` r
library(tidygraph)
```

    ## 
    ## Attaching package: 'tidygraph'

    ## The following object is masked from 'package:igraph':
    ## 
    ##     groups

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(ggraph)
knitr::opts_chunk$set(
  echo = TRUE,         # Display code chunks
  eval = TRUE         # Evaluate code chunks
)
```

I want to generate an interaction network plot showing the sum of
miRNA-mRNA-lncRNA interactions.

This will include all putative miRNA-mRNA interactions (miRanda +
significant PCC) and all putative miRNA-lncRNA interactions (miRanda +
significant PCC). However, it will *NOT* include mRNA-lncRNA links. This
is for two reasons:

1.  First, that dataset is simply too large – there were tens of
    millions of mRNA-lncRNA pairs with significantly correlated
    expression.

2.  Second, the only method we have to investigate potential mRNA-lncRNA
    links is expression correlation, which doesn’t provide any
    indication of the direction of action and is generally inappropriate
    to use as a sole predictor of causative relationships. In contrast,
    the miRNA-mRNA and miRNA-lncRNA putative interactions are primarily
    supported by binding predictions (through miRanda), and only use
    expression correlation as a method of validating interactions and
    evaluating their “polarity” (positive or negative relationship).
    When making a summary plot, I think it’s only appropriate to include
    putative interactions which we can predict with similar degrees of
    confidence. Including mRNA-lncRNA interactions in the final plot may
    imply we are as confident in those interactions truly existing as we
    are in the miRNA-mRNA and miRNA-lncRNA interactions.

Because I’m excluding the mRNA-lncRNA links, I think this will
functionally end up being a putative ceRNA (competative endogenous)
network.

For both `igraph` and `Cytoscape` I need two inputs, an “Edges”
dataframe and a “Nodes” dataframe (both should be saved as `.csv` files
for export to `Cytoscape`.

The “Edges” file should associate each node with all other nodes it
connects to. It should also contain edge-specific metadata. For example:

| source   | target  | correlation | correlation magnitude | correlation direction | correlation pval | binding pval |
|:---------|:--------|:------------|:----------------------|:----------------------|:-----------------|:-------------|
| miR-100  | Peve001 | -0.9        | 0.9                   | -1                    | 0.001            | 0.02         |
| miR-100  | Peve002 | 0.85        | 0.85                  | 1                     | 0.02             | 0.03         |
| lncRNA01 | Peve001 | -0.95       | 0.95                  | -1                    | 0.01             | 0.01         |

Note that there may be duplicates in both the “source” and “target”
columns, and there may be duplicate source-target combinations if they
have different attributes (e.g. representing different predicted binding
sites of the same pair).However, the rows should be unique.

The “Nodes” file contains metadata for every node included in the plot.
Importantly, the set of nodes listed in the “Nodes” file should match
exactly the set of nodes included in the “Edges” document. For example:

| id       | type   |
|:---------|:-------|
| Peve001  | gene   |
| Peve002  | gene   |
| miR-100  | miRNA  |
| lncRNA01 | lncRNA |

I’ll need the following files to compile the Cytoscape inputs:

- miRNA-mRNA interaction files (contains binding and coexpression
  information for miRNA-gene pairs):

- miRNA-3UTR:
  `deep-dive-expression/E-Peve/output/10-Peve-mRNA-miRNA-interactions/miranda_PCC_miRNA_mRNA.csv`

- miRNA-CDS:
  `~/deep-dive-expression/E-Peve/output/10.01-Peve-mRNA-miRNA-interactions-CDS_5UTR/miRanda_PCC_miRNA_CDS.csv`

- miRNA-5UTR:
  `~/deep-dive-expression/E-Peve/output/10.01-Peve-mRNA-miRNA-interactions-CDS_5UTR/miRanda_PCC_miRNA_5UTR.csv`

- miRNA-lncRNA interaction file (contains binding and coexpression
  information for miRNA-lncRNA pairs):

- `~/deep-dive-expression/E-Peve/output/15-Peve-miRNA-lncRNA-PCC/miranda_PCC_miRNA-lncRNA.csv`

Load packages:

``` r
library(dplyr)
library(tidyr)
library(igraph)
```

Load and format files:

``` r
# Load the three miRNA-mRNA files, and format so they have the same column contents and names

miRNA_3UTR <- read.csv("../output/10-Peve-mRNA-miRNA-interactions/Peve-miranda_PCC_miRNA_mRNA.csv") %>% dplyr::select(-X.1, -X)
# Add label for binding region
miRNA_3UTR$region <- "3UTR"


miRNA_CDS <- read.csv("../output/10.01-Peve-mRNA-miRNA-interactions-CDS_5UTR/miRanda_PCC_miRNA_CDS.csv") %>% dplyr::select(-X)
colnames(miRNA_CDS) <- c("miRNA", "mRNA_coord", "score", "energy", "query_start", "query_end", "subject_start", "subject_end", "total_bp_shared", "query_similar", "subject_similar", "mRNA", "PCC.cor", "p_value", "adjusted_p_value")
miRNA_CDS$query_start_end <- paste0(miRNA_CDS$query_start, " ", miRNA_CDS$query_end)
miRNA_CDS$subject_start_end <- paste0(miRNA_CDS$subject_start, " ", miRNA_CDS$subject_end)
# Add label for binding region
miRNA_CDS$region <- "CDS"
# Match column order of miRNA_3UTR
miRNA_CDS <- miRNA_CDS %>% dplyr::select(colnames(miRNA_3UTR))


miRNA_5UTR <- read.csv("../output/10.01-Peve-mRNA-miRNA-interactions-CDS_5UTR/miRanda_PCC_miRNA_5UTR.csv") %>% dplyr::select(-X)
colnames(miRNA_5UTR) <- c("miRNA", "mRNA_coord", "score", "energy", "query_start", "query_end", "subject_start", "subject_end", "total_bp_shared", "query_similar", "subject_similar", "mRNA", "PCC.cor", "p_value", "adjusted_p_value")
miRNA_5UTR$query_start_end <- paste0(miRNA_5UTR$query_start, " ", miRNA_5UTR$query_end)
miRNA_5UTR$subject_start_end <- paste0(miRNA_5UTR$subject_start, " ", miRNA_5UTR$subject_end)
# Add label for binding region
miRNA_5UTR$region <- "5UTR"
# Match column order of miRNA_3UTR
miRNA_5UTR <- miRNA_5UTR %>% dplyr::select(colnames(miRNA_3UTR))


miRNA_gene <- rbind(miRNA_3UTR, miRNA_CDS, miRNA_5UTR)
# Remove NA rows
miRNA_gene <- miRNA_gene %>% filter(!is.na(miRNA))

# Load and format the miRNA-lncRNA file
miRNA_lncRNA <- read.csv("../output/15-Peve-miRNA-lncRNA-PCC/miranda_PCC_miRNA_lncRNA.csv") %>% dplyr::select(-X.1, -X)
# Add label for binding region
miRNA_lncRNA$region <- "lncRNA"
```

For miRNA, annotate with assigned names (e.g., Peve-miR-100,
Peve-novel-7)

``` r
Peve_names <- read.csv("../../E-Peve/output/05-Peve-sRNA-ShortStack_4.1.0/ShortStack_out/Peve_Results_mature_named_miRNAs.csv") %>% dplyr::select(Name, given_miRNA_name)

miRNA_gene <- left_join(miRNA_gene, Peve_names, by = c("miRNA" = "Name"))
miRNA_lncRNA <- left_join(miRNA_lncRNA, Peve_names, by = c("miRNA" = "Name"))
```

# 1 All interactions

Format “Edges” file:

``` r
# Add correlation magnitude and direction columns
miRNA_gene$PCC_magnitude <- abs(miRNA_gene$PCC.cor)
miRNA_gene$PCC_direction <- sign(miRNA_gene$PCC.cor)
miRNA_gene$Alignment <- paste0(miRNA_gene$mRNA, ";", miRNA_gene$query_start_end, ";", miRNA_gene$subject_start_end)
# Select columns I want to keep in Edges file
miRNA_gene_edges <- miRNA_gene %>% dplyr::select(given_miRNA_name, mRNA, region, Alignment, energy, total_bp_shared, query_similar, subject_similar, PCC.cor, PCC_magnitude, PCC_direction, p_value)
# rename columns
miRNA_gene_edges <- miRNA_gene_edges %>% rename(source = given_miRNA_name, target = mRNA)

# Add correlation magnitude and direction columns
miRNA_lncRNA$PCC_magnitude <- abs(miRNA_lncRNA$PCC.cor)
miRNA_lncRNA$PCC_direction <- sign(miRNA_lncRNA$PCC.cor)
miRNA_lncRNA$Alignment <- paste0(miRNA_lncRNA$lncRNA, ";", miRNA_lncRNA$query_start_end, ";", miRNA_lncRNA$subject_start_end)
# Select columns I want to keep in Edges file (ensure in same order as in the miRNA_gene_edges file)
miRNA_lncRNA_edges <- miRNA_lncRNA %>% dplyr::select(given_miRNA_name, lncRNA, region, Alignment, energy, total_bp_shared, query_similar, subject_similar, PCC.cor, PCC_magnitude, PCC_direction, p_value)
# rename columns
miRNA_lncRNA_edges <- miRNA_lncRNA_edges %>% rename(source = given_miRNA_name, target = lncRNA)

# Combine miRNA-gene edges and miRNA-lncRNA edges
edges <- rbind(miRNA_gene_edges, miRNA_lncRNA_edges)

# Ensure we have no duplicate rows
nrow(edges)
```

    ## [1] 26714

``` r
nrow(edges %>% distinct())
```

    ## [1] 26708

``` r
# Check formatting/contents
head(edges)
```

    ##              source        target region                  Alignment energy
    ## 1 peve-mir-novel-24 Peve_00009100   3UTR Peve_00009100;2 21;322 343 -20.39
    ## 2 peve-mir-novel-24 Peve_00009100   3UTR Peve_00009100;2 21;661 682 -20.39
    ## 3     peve-mir-2023 Peve_00009103   3UTR Peve_00009103;2 20;302 321 -21.38
    ## 4 peve-mir-novel-15 Peve_00018268   3UTR   Peve_00018268;2 16;36 60 -22.14
    ## 5  peve-mir-novel-2 Peve_00028953   3UTR Peve_00028953;2 13;168 189 -23.43
    ## 6 peve-mir-novel-39 Peve_00018488   3UTR Peve_00018488;2 21;241 264 -22.50
    ##   total_bp_shared query_similar subject_similar    PCC.cor PCC_magnitude
    ## 1              19        73.68%          84.21% -0.4772517     0.4772517
    ## 2              19        73.68%          84.21% -0.4772517     0.4772517
    ## 3              18        83.33%          88.89%  0.4539645     0.4539645
    ## 4              17        76.47%          76.47% -0.6860633     0.6860633
    ## 5              11        90.91%          90.91% -0.9581593     0.9581593
    ## 6              21        71.43%          80.95%  0.8679339     0.8679339
    ##   PCC_direction   p_value
    ## 1            -1 0.6832660
    ## 2            -1 0.6832660
    ## 3             1 0.7000186
    ## 4            -1 0.5186746
    ## 5            -1 0.1848079
    ## 6             1 0.3308953

Format Nodes file:

``` r
# Make a df that contains all miRNA, genes, and lncRNA listed in the `source` and `target` columns of `edges`
nodes <- data.frame(
  # The `unique` argument ensures we remove duplicates
  id = unique(unname(unlist(edges[, c("source", "target")])))
)

# Add column identifying the type of each node (miRNA, lncRNA, or gene)
nodes <- nodes %>%
  mutate(type = case_when(
    grepl("mir", id) ~ "miRNA",
    grepl("Peve", id) ~ "gene",
    grepl("lncRNA", id) ~ "lncRNA",
    TRUE ~ "other"
  ))

# Check formatting/contents
head(nodes)
```

    ##                  id  type
    ## 1 peve-mir-novel-24 miRNA
    ## 2     peve-mir-2023 miRNA
    ## 3 peve-mir-novel-15 miRNA
    ## 4  peve-mir-novel-2 miRNA
    ## 5 peve-mir-novel-39 miRNA
    ## 6 peve-mir-novel-40 miRNA

# 2 pval \< 0.05

Edges:

``` r
edges_pval_0.05 <- edges %>% filter(p_value < 0.05)
nrow(edges_pval_0.05)
```

    ## [1] 1442

Nodes:

``` r
# Make a df that contains all miRNA, genes, and lncRNA listed in the `source` and `target` columns of `edges`
nodes_pval_0.05 <- data.frame(
  # The `unique` argument ensures we remove duplicates
  id = unique(unname(unlist(edges_pval_0.05[, c("source", "target")])))
)

# Add column identifying the type of each node (miRNA, lncRNA, or gene)
nodes_pval_0.05 <- nodes_pval_0.05 %>%
  mutate(type = case_when(
    grepl("mir", id) ~ "miRNA",
    grepl("Peve", id) ~ "gene",
    grepl("lncRNA", id) ~ "lncRNA",
    TRUE ~ "other"
  ))

nrow(nodes_pval_0.05)
```

    ## [1] 1293

# 3 pval \< 0.01

Edges:

``` r
edges_pval_0.01 <- edges %>% filter(p_value < 0.01)
nrow(edges_pval_0.01)
```

    ## [1] 455

Nodes:

``` r
# Make a df that contains all miRNA, genes, and lncRNA listed in the `source` and `target` columns of `edges`
nodes_pval_0.01 <- data.frame(
  # The `unique` argument ensures we remove duplicates
  id = unique(unname(unlist(edges_pval_0.01[, c("source", "target")])))
)

# Add column identifying the type of each node (miRNA, lncRNA, or gene)
nodes_pval_0.01 <- nodes_pval_0.01 %>%
  mutate(type = case_when(
    grepl("mir", id) ~ "miRNA",
    grepl("Peve", id) ~ "gene",
    grepl("lncRNA", id) ~ "lncRNA",
    TRUE ~ "other"
  ))

nrow(nodes_pval_0.01)
```

    ## [1] 444

# 4 Save

Save files

``` r
write.csv(edges_pval_0.05, "../output/31-Peve-miRNA-mRNA-lncRNA-network/edges_miRNA_mRNA_lncRNA_network_p0.05.csv", quote = FALSE)
write.csv(nodes_pval_0.05, "../output/31-Peve-miRNA-mRNA-lncRNA-network/nodes_miRNA_mRNA_lncRNA_network_p0.05.csv", quote = FALSE)
write.csv(edges_pval_0.01, "../output/31-Peve-miRNA-mRNA-lncRNA-network/edges_miRNA_mRNA_lncRNA_network_p0.01.csv", quote = FALSE)
write.csv(nodes_pval_0.01, "../output/31-Peve-miRNA-mRNA-lncRNA-network/nodes_miRNA_mRNA_lncRNA_network_p0.01.csv", quote = FALSE)
```

# 5 Using Cytoscape

To load a network into Cytoscape:

1.  Open Cytoscape and select File \> Import \> Network from File…

2.  Select “Edges” file. Ensure The source and target columns are
    appropriately identified before loading the file.

3.  To load “Nodes” file, select File \> Import \> Table from File…

# 6 Plotting with igraph

``` r
# Rename columns for igraph
colnames(edges_pval_0.01)[1:2] <- c("from", "to")         # For igraph edge input
colnames(nodes_pval_0.01)[1] <- "name"                    # For igraph vertex input; must match nodes in edges

# Build graph
g <- graph_from_data_frame(d = edges_pval_0.01, vertices = nodes_pval_0.01, directed = FALSE)

# Edge attributes
E(g)$edge_color <- ifelse(E(g)$PCC_direction > 0, "positive", "negative")

# Convert to tbl_graph
g_tbl <- as_tbl_graph(g)


p <- ggraph(g_tbl, layout = "fr") +
  geom_edge_link(aes(edge_width = PCC_magnitude, color = edge_color), alpha = 0.6) +
  geom_node_point(aes(color = type), size = 2) +
  scale_edge_width(range = c(0.5, 3)) +
  
  # Split color scales
  scale_edge_color_manual(
    values = c("positive" = "green3", "negative" = "red"),
    name = "Correlation Direction"
  ) +
  scale_color_manual(
    values = c("miRNA" = "orange", "gene" = "lightblue", "lncRNA" = "steelblue4"),
    name = "type"
  ) +

  theme_graph() +
  labs(title = "miRNA-lncRNA-mRNA Interaction Network")

print(p)
```

![](31-Peve-miRNA-mRNA-lncRNA-network_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
