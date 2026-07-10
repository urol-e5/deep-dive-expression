15-miRNA-mRNA-lncRNA-network-ceRNA
================
Kathleen Durkin
2026-07-09

- [1 Helper function](#1-helper-function)
- [2 Apul](#2-apul)
- [3 Peve](#3-peve)
- [4 Ptuh](#4-ptuh)
- [5 Using Cytoscape](#5-using-cytoscape)
- [6 Exporting a web session with a
  legend](#6-exporting-a-web-session-with-a-legend)

For each species, the species-specific network scripts
(`33-Apul-miRNA-mRNA-lncRNA-network`,
`31-Peve-miRNA-mRNA-lncRNA-network`, and `31-Ptuh-miRNA-lncRNA-network`)
generate Cytoscape-compatible network tables containing miRNA-mRNA and
miRNA-lncRNA edges supported by both predicted binding (miRanda) and
significant coexpression. However, those networks do not include any
lncRNA-mRNA links.

This script modifies those networks to make three additions/changes:

1.  **ceRNA-associated lncRNA-mRNA edges.** For each species, the ceRNA
    network analysis (`33.1-Apul-ceRNA-network`,
    `31.1-Peve-ceRNA-network`, `31.1-Ptuh-ceRNA-network`) identified
    triads in which an lncRNA and mRNA share a common miRNA, the lncRNA
    is negatively coexpressed with that miRNA, and the lncRNA is
    positively coexpressed with the mRNA. Here I want to add lncRNA-mRNA
    edges for these ceRNA-associated pairs. Only instances where both
    the lncRNA and mRNA are part of the same ceRNA triad are included.
    Unlike the miRNA-target edges, there was no binding prediction
    performed between lncRNA and mRNA — the interaction is only
    supported by significant coexpression within a ceRNA context.

2.  **An `interaction_support` column in the Edges file** indicating the
    type of evidence supporting each edge:

    - `"predicted binding + coexpression"`: miRNA-mRNA and miRNA-lncRNA
      edges (miRanda binding prediction + significant PCC)
    - `"coexpression alone"`: ceRNA-associated lncRNA-mRNA edges
      (significant PCC only, no predicted binding)

3.  **A `special_status` column in the Nodes file** flagging nodes with
    a “special” designation (from `12-miRNA-epimachinery`):

    - `"epi-miRNA"`: miRNA nodes that putatively target epigenetic /
      ncRNA machinery
    - `"ceRNA"`: lncRNA nodes that are putative ceRNAs (i.e. appear in
      the ceRNA triads)
    - `"epi-machinery"`: mRNA nodes that encode epigenetic / ncRNA
      machinery proteins

I also want to save “unique” versions of these networks that remove
duplicate edges for unique node pairs (these have been used to indicate
instances of repetitive binding, in which there are multiple predicted
binding sites for a single miRNA-target pair).

For `Cytoscape` I need two inputs: an “Edges” dataframe and a “Nodes”
dataframe (both should be saved as `.csv` files for export to
`Cytoscape`).

The “Edges” file should associate each node with all other nodes it
connects to. It should also contain edge-specific metadata. For example:

| source | target | region | PCC.cor | PCC_magnitude | PCC_direction | p_value | interaction_support |
|:---|:---|:---|:---|:---|:---|:---|:---|
| miR-100 | FUN001 | 3UTR | -0.9 | 0.9 | -1 | 0.001 | predicted binding + coexpression |
| miR-100 | lncRNA01 | lncRNA | -0.95 | 0.95 | -1 | 0.01 | predicted binding + coexpression |
| lncRNA01 | FUN001 | lncRNA_mRNA | 0.95 | 0.95 | 1 | 0.01 | coexpression alone |

Note that there may be duplicates in both the “source” and “target”
columns, and there may be duplicate source-target combinations if they
have different attributes (e.g. representing different predicted binding
sites of the same pair). However, the rows should be unique.

The “Nodes” file contains metadata for every node included in the plot.
Importantly, the set of nodes listed in the “Nodes” file should match
exactly the set of nodes included in the “Edges” document. For example:

| id       | type   | special_status | notes                           |
|:---------|:-------|:---------------|:--------------------------------|
| FUN001   | gene   | epi-machinery  | USP19-205 (Ubiquitin signaling) |
| miR-100  | miRNA  | epi-miRNA      |                                 |
| lncRNA01 | lncRNA | ceRNA          |                                 |

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
library(tidyr)
knitr::opts_chunk$set(
  echo = TRUE,         # Display code chunks
  eval = TRUE         # Evaluate code chunks
)
```

Load the miRNA-epimachinery targeting results (used to identify
epi-miRNAs and epimachinery mRNAs across all three species), and build a
lookup table mapping each epimachinery mRNA to the name(s) of the
epimachinery gene(s) it matches

``` r
epi_miRNA_targets <- read.csv("../output/12-miRNA-epimachinery/miRNAtargets_mach.csv")

# Build a lookup table mapping each epimachinery mRNA to the epimachinery gene name(s) it matches, along with the functional category label
epi_machinery_notes <- epi_miRNA_targets %>%
  group_by(target) %>%
  summarise(
    gene_names = paste(unique(unlist(strsplit(gene, ";\\s*"))), collapse = "; "),
    category   = paste(unique(category), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(notes = paste0(gene_names, " (", category, ")")) %>%
  dplyr::select(target, notes)
```

# 1 Helper function

To avoid repeating the same network-building code for each species,
define a function. Should: - take the species-specific inputs - add the
ceRNA-associated lncRNA-mRNA edges - annotate an `interaction_support`
column on the edges - annotate `special_status` and `notes` columns on
the nodes - save the Cytoscape-compatible network tables at both p \<
0.05 and p \< 0.01 coexpression significance levels.

``` r
build_ceRNA_network <- function(edges_file, ceRNA_file, species_label, gene_pattern,
                                strip_gene_prefix, outdir, file_prefix) {

  # Create output directory (if it doesn't exist)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # Load the existing miRNA-target edges (from the species-specific network script)
  edges <- read.csv(edges_file) %>% dplyr::select(-X)

  # Load the ceRNA triads
  # Each row = one miRNA-lncRNA-mRNA set with significant positive lncRNA-mRNA coexpression
  ceRNA <- read.delim(ceRNA_file, header = TRUE, sep = "\t")

  # Identify epi-miRNAs and epimachinery mRNAs for this species
  epi_miRNAs <- epi_miRNA_targets %>%
    filter(species == species_label) %>%
    pull(given_miRNA_name) %>%
    unique()

  epi_machinery_mRNAs <- epi_miRNA_targets %>%
    filter(species == species_label) %>%
    pull(target) %>%
    unique()

  # Identify ceRNA lncRNAs
  ceRNA_lncRNAs <- ceRNA$lncRNA %>% unique()

  # Build lncRNA-mRNA edges from the ceRNA triads
  # Strip "gene-" prefix from mRNA names if needed (Peve/Ptuh ceRNA files use it
  # to match the lncRNA-mRNA PCC file, but the network uses unprefixed names)
  mRNA_target <- if (strip_gene_prefix) sub("^gene-", "", ceRNA$mRNA) else ceRNA$mRNA

  lncRNA_mRNA_edges <- ceRNA %>%
    mutate(
      source = lncRNA,
      target = mRNA_target,
      region = "lncRNA_mRNA",
      Alignment = NA_character_,
      energy = NA_real_,
      total_bp_shared = NA_integer_,
      query_similar = NA_character_,
      subject_similar = NA_character_,
      PCC.cor = PCC_lncRNA_mRNA,
      PCC_magnitude = abs(PCC_lncRNA_mRNA),
      PCC_direction = sign(PCC_lncRNA_mRNA),
      p_value = p_value_lncRNA_mRNA,
      interaction_support = "coexpression alone"
    ) %>%
    dplyr::select(source, target, region, Alignment, energy, total_bp_shared,
                  query_similar, subject_similar, PCC.cor, PCC_magnitude,
                  PCC_direction, p_value, interaction_support) %>%
    distinct()

  # Add interaction_support to the existing miRNA-target edges
  edges <- edges %>%
    mutate(interaction_support = "predicted binding + coexpression")

  # Combine miRNA-target edges with ceRNA-associated lncRNA-mRNA edges
  # The ceRNA triads were built from the p < 0.05 miRNA-target edges, so at that
  # threshold every node is connected to at least one miRNA.
  # NOTE: At stricter thresholds (e.g. p < 0.01) some nodes may lose their miRNA edge while
  # their ceRNA-associated lncRNA-mRNA edge survives, producing nodes that are only
  # connected via coexpression-alone edges
  edges <- rbind(edges, lncRNA_mRNA_edges)

  # Helper to build a nodes df from an edges df
  build_nodes <- function(edges_df) {
    nodes <- data.frame(
      # The `unique` argument ensures we remove duplicates
      id = unique(unname(unlist(edges_df[, c("source", "target")])))
    )
    nodes <- nodes %>%
      mutate(type = case_when(
        grepl("mir", id) ~ "miRNA",
        grepl(gene_pattern, id) ~ "gene",
        grepl("lncRNA", id) ~ "lncRNA",
        TRUE ~ "other"
      )) %>%
      # Annotate special status: epi-miRNA, ceRNA, epi-machinery, or NA
      mutate(special_status = case_when(
        id %in% epi_miRNAs ~ "epi-miRNA",
        id %in% ceRNA_lncRNAs ~ "ceRNA",
        id %in% epi_machinery_mRNAs ~ "epi-machinery",
        TRUE ~ NA_character_
      )) %>%
      # For epimachinery mRNAs, add a notes column with the epimachinery gene name(s)
      left_join(epi_machinery_notes, by = c("id" = "target")) %>%
      mutate(notes = ifelse(is.na(special_status) | special_status != "epi-machinery", NA_character_, notes))
    return(nodes)
  }

  # pval < 0.05
  edges_pval_0.05 <- edges %>% filter(p_value < 0.05)
  nodes_pval_0.05 <- build_nodes(edges_pval_0.05)

  # pval < 0.01
  edges_pval_0.01 <- edges %>% filter(p_value < 0.01)
  nodes_pval_0.01 <- build_nodes(edges_pval_0.01)

  # Save
  # Interaction-level
  write.csv(edges_pval_0.05, file.path(outdir, paste0(file_prefix, "_edges_miRNA_mRNA_lncRNA_ceRNA_network_p0.05.csv")), quote = FALSE)
  write.csv(nodes_pval_0.05, file.path(outdir, paste0(file_prefix, "_nodes_miRNA_mRNA_lncRNA_ceRNA_network_p0.05.csv")), quote = FALSE)
  write.csv(edges_pval_0.01, file.path(outdir, paste0(file_prefix, "_edges_miRNA_mRNA_lncRNA_ceRNA_network_p0.01.csv")), quote = FALSE)
  write.csv(nodes_pval_0.01, file.path(outdir, paste0(file_prefix, "_nodes_miRNA_mRNA_lncRNA_ceRNA_network_p0.01.csv")), quote = FALSE)

  # Unique-feature level
  edges_pval_0.05_unique <- edges_pval_0.05 %>% dplyr::select(-region, -Alignment, -energy, -total_bp_shared, -query_similar, -subject_similar) %>% distinct()
  edges_pval_0.01_unique <- edges_pval_0.01 %>% dplyr::select(-region, -Alignment, -energy, -total_bp_shared, -query_similar, -subject_similar) %>% distinct()

  write.csv(edges_pval_0.05_unique, file.path(outdir, paste0(file_prefix, "_edges_miRNA_mRNA_lncRNA_ceRNA_network_p0.05_unique.csv")), quote = FALSE)
  write.csv(edges_pval_0.01_unique, file.path(outdir, paste0(file_prefix, "_edges_miRNA_mRNA_lncRNA_ceRNA_network_p0.01_unique.csv")), quote = FALSE)

  # Print summary counts
  cat("Species:", species_label, "\n")
  cat("  epi-miRNAs:", length(epi_miRNAs), "\n")
  cat("  epi-machinery mRNAs:", length(epi_machinery_mRNAs), "\n")
  cat("  ceRNA lncRNAs:", length(ceRNA_lncRNAs), "\n")
  cat("  ceRNA-associated lncRNA-mRNA edges added:", nrow(lncRNA_mRNA_edges), "\n")
  cat("  Total edges (p < 0.05):", nrow(edges_pval_0.05), "\n")
  cat("  Total nodes (p < 0.05):", nrow(nodes_pval_0.05), "\n")
  cat("    miRNA:", sum(nodes_pval_0.05$type == "miRNA"),
      " gene:", sum(nodes_pval_0.05$type == "gene"),
      " lncRNA:", sum(nodes_pval_0.05$type == "lncRNA"), "\n")
  cat("    special_status (p < 0.05): epi-miRNA:", sum(nodes_pval_0.05$special_status == "epi-miRNA", na.rm = TRUE),
      " ceRNA:", sum(nodes_pval_0.05$special_status == "ceRNA", na.rm = TRUE),
      " epi-machinery:", sum(nodes_pval_0.05$special_status == "epi-machinery", na.rm = TRUE), "\n")
  cat("  Total edges (p < 0.01):", nrow(edges_pval_0.01), "\n")
  cat("  Total nodes (p < 0.01):", nrow(nodes_pval_0.01), "\n")
  cat("    miRNA:", sum(nodes_pval_0.01$type == "miRNA"),
      " gene:", sum(nodes_pval_0.01$type == "gene"),
      " lncRNA:", sum(nodes_pval_0.01$type == "lncRNA"), "\n")
  cat("    special_status (p < 0.01): epi-miRNA:", sum(nodes_pval_0.01$special_status == "epi-miRNA", na.rm = TRUE),
      " ceRNA:", sum(nodes_pval_0.01$special_status == "ceRNA", na.rm = TRUE),
      " epi-machinery:", sum(nodes_pval_0.01$special_status == "epi-machinery", na.rm = TRUE), "\n\n")
}
```

# 2 Apul

``` r
build_ceRNA_network(
  edges_file       = "../../D-Apul/output/33-Apul-miRNA-mRNA-lncRNA-network/edges_miRNA_mRNA_lncRNA_network_p0.05.csv",
  ceRNA_file       = "../../D-Apul/output/33.1-Apul-ceRNA-network/pos_lncRNA_mRNA_hits.txt",
  species_label    = "A. pulchra",
  gene_pattern     = "FUN",
  strip_gene_prefix = FALSE,
  outdir           = "../output/15-miRNA-mRNA-lncRNA-network-ceRNA",
  file_prefix      = "Apul"
)
```

    ## Species: A. pulchra 
    ##   epi-miRNAs: 24 
    ##   epi-machinery mRNAs: 34 
    ##   ceRNA lncRNAs: 117 
    ##   ceRNA-associated lncRNA-mRNA edges added: 3171 
    ##   Total edges (p < 0.05): 5957 
    ##   Total nodes (p < 0.05): 2302 
    ##     miRNA: 39  gene: 1821  lncRNA: 442 
    ##     special_status (p < 0.05): epi-miRNA: 24  ceRNA: 117  epi-machinery: 34 
    ##   Total edges (p < 0.01): 1978 
    ##   Total nodes (p < 0.01): 874 
    ##     miRNA: 39  gene: 641  lncRNA: 194 
    ##     special_status (p < 0.01): epi-miRNA: 24  ceRNA: 111  epi-machinery: 13

# 3 Peve

Note: the Peve ceRNA file uses a `gene-` prefix on mRNA names (to match
the lncRNA-mRNA PCC file), so `strip_gene_prefix = TRUE` is used to
match the network’s mRNA naming.

``` r
build_ceRNA_network(
  edges_file       = "../../E-Peve/output/31-Peve-miRNA-mRNA-lncRNA-network/edges_miRNA_mRNA_lncRNA_network_p0.05.csv",
  ceRNA_file       = "../../E-Peve/output/31.1-Peve-ceRNA-network/pos_lncRNA_mRNA_hits.txt",
  species_label    = "P. evermanni",
  gene_pattern     = "Peve",
  strip_gene_prefix = TRUE,
  outdir           = "../output/15-miRNA-mRNA-lncRNA-network-ceRNA",
  file_prefix      = "Peve"
)
```

    ## Species: P. evermanni 
    ##   epi-miRNAs: 8 
    ##   epi-machinery mRNAs: 14 
    ##   ceRNA lncRNAs: 41 
    ##   ceRNA-associated lncRNA-mRNA edges added: 450 
    ##   Total edges (p < 0.05): 1892 
    ##   Total nodes (p < 0.05): 1293 
    ##     miRNA: 42  gene: 1091  lncRNA: 160 
    ##     special_status (p < 0.05): epi-miRNA: 8  ceRNA: 41  epi-machinery: 14 
    ##   Total edges (p < 0.01): 648 
    ##   Total nodes (p < 0.01): 580 
    ##     miRNA: 31  gene: 490  lncRNA: 59 
    ##     special_status (p < 0.01): epi-miRNA: 7  ceRNA: 35  epi-machinery: 8

# 4 Ptuh

Note: the Ptuh ceRNA file uses a `gene-` prefix on mRNA names (to match
the lncRNA-mRNA PCC file), so `strip_gene_prefix = TRUE` is used to
match the network’s mRNA naming.

``` r
build_ceRNA_network(
  edges_file       = "../../F-Ptuh/output/31-Ptuh-miRNA-mRNA-lncRNA-network/edges_miRNA_mRNA_lncRNA_network_p0.05.csv",
  ceRNA_file       = "../../F-Ptuh/output/31.1-Ptuh-ceRNA-network/pos_lncRNA_mRNA_hits.txt",
  species_label    = "P. tuahiniensis",
  gene_pattern     = "Pocillopora",
  strip_gene_prefix = TRUE,
  outdir           = "../output/15-miRNA-mRNA-lncRNA-network-ceRNA",
  file_prefix      = "Ptuh"
)
```

    ## Species: P. tuahiniensis 
    ##   epi-miRNAs: 11 
    ##   epi-machinery mRNAs: 16 
    ##   ceRNA lncRNAs: 161 
    ##   ceRNA-associated lncRNA-mRNA edges added: 3245 
    ##   Total edges (p < 0.05): 4711 
    ##   Total nodes (p < 0.05): 1266 
    ##     miRNA: 37  gene: 836  lncRNA: 393 
    ##     special_status (p < 0.05): epi-miRNA: 11  ceRNA: 161  epi-machinery: 16 
    ##   Total edges (p < 0.01): 1685 
    ##   Total nodes (p < 0.01): 600 
    ##     miRNA: 34  gene: 378  lncRNA: 188 
    ##     special_status (p < 0.01): epi-miRNA: 11  ceRNA: 144  epi-machinery: 5

# 5 Using Cytoscape

To load a network into Cytoscape:

1.  Open Cytoscape and select File \> Import \> Network from File…

2.  Select the “Edges” file. Ensure the source and target columns are
    appropriately identified before loading the file.

3.  To load the “Nodes” file, select File \> Import \> Table from File…

# 6 Exporting a web session with a legend

After styling the networks in Cytoscape (using the “default” visual
style, which maps node fill color to `type`, node shape to
`special_status`, edge color to `PCC_direction`, and edge line style to
`interaction_support`; using the yFiles Organic Layout), the session was
exported as a standalone web app via File \> Export \> Network to Web
Session… . It can be launched by clicking web_session/index.html and
selecting “View in Web Browser”

The webapp does not include a legend of visual features, so I then used
OpenCode GLM 5.2 to manually add a legend to the exported Cytoscape
webapp files, as described below:

Cytoscape web session exports include a `scripts/custom.js` file that is
intentionally left as a placeholder for user customization. This file
was used to inject an HTML/CSS legend overlay into the web app, so that
the meaning of each visual encoding is visible when viewing the network
in a browser. No other files in the web session export were modified.

The legend is added entirely through jQuery DOM injection in `custom.js`
— CSS is injected into a `<style>` tag in the document head, and the
legend HTML is appended to `<body>` on document ready. A ☰ button in the
legend header toggles its visibility. The RGB values for node colors and
edge colors are taken directly from the Cytoscape visual style mappings
(parsed from the exported `data/styles.js` file).

The full contents of the manually added `custom.js` (located at
`../output/15-miRNA-mRNA-lncRNA-network-ceRNA/web_session/scripts/custom.js`)
are reproduced below:

``` js
$( document ).ready(function(){

  // ---- Legend overlay ----
  // Describes the visual style mappings used in the "default" visual style:
  //   Node fill color  -> feature type (miRNA, mRNA, lncRNA)
  //   Node shape       -> special_status (epi-miRNA, ceRNA, epi-machinery, none)
  //   Edge line color  -> PCC_direction (positive, negative)
  //   Edge line style  -> interaction_support (predicted binding + coexpression, coexpression alone)

  var legendCSS = [
    '#cytoscape-legend {',
    '  position: fixed;',
    '  bottom: 16px;',
    '  right: 16px;',
    '  z-index: 9999;',
    '  background: rgba(255,255,255,0.95);',
    '  border: 1px solid #aaa;',
    '  border-radius: 8px;',
    '  padding: 12px 16px;',
    '  font-family: sans-serif;',
    '  font-size: 12px;',
    '  color: #333;',
    '  width: 280px;',
    '  box-shadow: 0 2px 12px rgba(0,0,0,0.2);',
    '  user-select: none;',
    '  line-height: 1.4;',
    '}',
    '#cytoscape-legend .legend-header {',
    '  display: flex;',
    '  align-items: center;',
    '  justify-content: space-between;',
    '  margin-bottom: 8px;',
    '}',
    '#cytoscape-legend .legend-title {',
    '  font-size: 14px;',
    '  font-weight: bold;',
    '}',
    '#cytoscape-legend .toggle-btn {',
    '  cursor: pointer;',
    '  font-size: 16px;',
    '  color: #999;',
    '  border: none;',
    '  background: none;',
    '  padding: 0 2px;',
    '}',
    '#cytoscape-legend .toggle-btn:hover { color: #333; }',
    '#cytoscape-legend.collapsed .legend-content { display: none; }',
    '#cytoscape-legend.collapsed { width: auto; }',
    '',
    '#cytoscape-legend .legend-section {',
    '  margin-bottom: 10px;',
    '}',
    '#cytoscape-legend .legend-section:last-child { margin-bottom: 0; }',
    '#cytoscape-legend .section-title {',
    '  font-size: 11px;',
    '  font-weight: bold;',
    '  color: #666;',
    '  text-transform: uppercase;',
    '  letter-spacing: 0.5px;',
    '  margin-bottom: 5px;',
    '  padding-bottom: 3px;',
    '  border-bottom: 1px solid #ddd;',
    '}',
    '#cytoscape-legend .legend-item {',
    '  display: flex;',
    '  align-items: center;',
    '  margin: 4px 0;',
    '}',
    '#cytoscape-legend .legend-item .swatch {',
    '  flex: 0 0 28px;',
    '  margin-right: 8px;',
    '  display: inline-flex;',
    '  justify-content: center;',
    '  align-items: center;',
    '}',
    '#cytoscape-legend .legend-item .desc {',
    '  flex: 1;',
    '}',
    '#cytoscape-legend .legend-item .desc .desc-name {',
    '  font-weight: 600;',
    '}',
    '#cytoscape-legend .legend-item .desc .desc-detail {',
    '  color: #777;',
    '  font-size: 11px;',
    '}',
    '',
    '/* Node shape swatches */',
    '#cytoscape-legend .sw {',
    '  width: 20px;',
    '  height: 20px;',
    '  display: inline-block;',
    '}',
    '#cytoscape-legend .sw-ellipse { border-radius: 50%; }',
    '#cytoscape-legend .sw-rectangle { border-radius: 2px; }',
    '#cytoscape-legend .sw-diamond {',
    '  transform: rotate(45deg);',
    '  width: 14px;',
    '  height: 14px;',
    '  margin: 3px;',
    '}',
    '#cytoscape-legend .sw-triangle {',
    '  width: 0;',
    '  height: 0;',
    '  border-left: 10px solid transparent;',
    '  border-right: 10px solid transparent;',
    '  border-bottom: 17px solid currentColor;',
    '  background: none;',
    '}',
    '',
    '/* Edge line swatches */',
    '#cytoscape-legend .sw-edge {',
    '  width: 28px;',
    '  height: 0;',
    '  display: inline-block;',
    '}',
  ].join('\n');

  $('<style>' + legendCSS + '</style>').appendTo('head');

  var legendHTML = [
    '<div id="cytoscape-legend">',
    '  <div class="legend-header">',
    '    <span class="legend-title">Legend</span>',
    '    <button class="toggle-btn" title="Show/Hide legend">&#9776;</button>',
    '  </div>',
    '  <div class="legend-content">',
    '',
    '    <div class="legend-section">',
    '      <div class="section-title">Node fill color &rarr; feature type</div>',
    '      <div class="legend-item">',
    '        <span class="swatch"><span class="sw sw-ellipse" style="background:rgb(255,127,0)"></span></span>',
    '        <span class="desc"><span class="desc-name">Orange</span> <span class="desc-detail">&mdash; miRNA</span></span>',
    '      </div>',
    '      <div class="legend-item">',
    '        <span class="swatch"><span class="sw sw-ellipse" style="background:rgb(166,206,227)"></span></span>',
    '        <span class="desc"><span class="desc-name">Light blue</span> <span class="desc-detail">&mdash; mRNA</span></span>',
    '      </div>',
    '      <div class="legend-item">',
    '        <span class="swatch"><span class="sw sw-ellipse" style="background:rgb(31,120,180)"></span></span>',
    '        <span class="desc"><span class="desc-name">Dark blue</span> <span class="desc-detail">&mdash; lncRNA</span></span>',
    '      </div>',
    '    </div>',
    '',
    '    <div class="legend-section">',
    '      <div class="section-title">Node shape &rarr; special status</div>',
    '      <div class="legend-item">',
    '        <span class="swatch"><span class="sw sw-rectangle" style="background:#ddd;border:1px solid #999"></span></span>',
    '        <span class="desc"><span class="desc-name">Rectangle</span> <span class="desc-detail">&mdash; epi-miRNA (targets epigenetic machinery)</span></span>',
    '      </div>',
    '      <div class="legend-item">',
    '        <span class="swatch"><span class="sw sw-triangle" style="color:#bbb"></span></span>',
    '        <span class="desc"><span class="desc-name">Triangle</span> <span class="desc-detail">&mdash; ceRNA (putative competing endogenous RNA)</span></span>',
    '      </div>',
    '      <div class="legend-item">',
    '        <span class="swatch"><span class="sw sw-diamond" style="background:#ddd;border:1px solid #999"></span></span>',
    '        <span class="desc"><span class="desc-name">Diamond</span> <span class="desc-detail">&mdash; epi-machinery (encodes epigenetic/ncRNA machinery)</span></span>',
    '      </div>',
    '      <div class="legend-item">',
    '        <span class="swatch"><span class="sw sw-ellipse" style="background:#ddd;border:1px solid #999"></span></span>',
    '        <span class="desc"><span class="desc-name">Circle</span> <span class="desc-detail">&mdash; no special status</span></span>',
    '      </div>',
    '    </div>',
    '',
    '    <div class="legend-section">',
    '      <div class="section-title">Edge color &rarr; coexpression direction</div>',
    '      <div class="legend-item">',
    '        <span class="swatch"><span class="sw-edge" style="border-top:3px solid rgb(51,160,44)"></span></span>',
    '        <span class="desc"><span class="desc-name">Green</span> <span class="desc-detail">&mdash; positive coexpression</span></span>',
    '      </div>',
    '      <div class="legend-item">',
    '        <span class="swatch"><span class="sw-edge" style="border-top:3px solid rgb(227,26,28)"></span></span>',
    '        <span class="desc"><span class="desc-name">Red</span> <span class="desc-detail">&mdash; negative coexpression</span></span>',
    '      </div>',
    '    </div>',
    '',
    '    <div class="legend-section">',
    '      <div class="section-title">Edge style &rarr; interaction support</div>',
    '      <div class="legend-item">',
    '        <span class="swatch"><span class="sw-edge" style="border-top:3px solid #555"></span></span>',
    '        <span class="desc"><span class="desc-name">Solid</span> <span class="desc-detail">&mdash; predicted binding + coexpression</span></span>',
    '      </div>',
    '      <div class="legend-item">',
    '        <span class="swatch"><span class="sw-edge" style="border-top:3px dashed #555"></span></span>',
    '        <span class="desc"><span class="desc-name">Dashed</span> <span class="desc-detail">&mdash; coexpression alone (no binding prediction performed between lncRNA and mRNA)</span></span>',
    '      </div>',
    '    </div>',
    '',
    '  </div>',
    '</div>',
  ].join('\n');

  $(legendHTML).appendTo('body');

  // Toggle legend visibility
  $('#cytoscape-legend .toggle-btn').on('click', function() {
    $('#cytoscape-legend').toggleClass('collapsed');
  });

});
```
