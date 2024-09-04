03-Apul-RNA-summary
================
Kathleen Durkin
2024-08-20

- <a href="#001-install-and-load-packages"
  id="toc-001-install-and-load-packages">0.0.1 Install and load
  packages</a>
- <a href="#1-load-data" id="toc-1-load-data">1 Load data</a>
  - <a href="#11-load-count-data" id="toc-11-load-count-data">1.1 Load count
    data</a>
  - <a href="#12-count-data-munging" id="toc-12-count-data-munging">1.2
    Count data munging</a>
- <a href="#2-summary-stats-and-visualizations"
  id="toc-2-summary-stats-and-visualizations">2 Summary stats and
  visualizations</a>
  - <a href="#21-expression-levels" id="toc-21-expression-levels">2.1
    Expression levels</a>
  - <a href="#22-transcript-counts" id="toc-22-transcript-counts">2.2
    Transcript counts</a>
  - <a href="#23-most-common-biological-processes"
    id="toc-23-most-common-biological-processes">2.3 Most common biological
    processes</a>

Gene expression summary for *Acropora pulchra* RNA-seq data.

- trimmed reads generated in `deep-dive` project, trimming and QC
  details in `01-Apul-RNA-trimming-FastQC`

- Reads aligned to *Acropora millipora* transcriptome downloaded from
  [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_013753865.1/),
  stored
  [here](https://gannet.fish.washington.edu/acropora/E5-deep-dive/Transcripts/Apul_GCF_013753865.1_rna.fna)
  as a part of [deep-dive genomic
  resources](https://github.com/urol-e5/deep-dive/wiki/Species-Characteristics-and-Genomic-Resources#genomic-resources).

### 0.0.1 Install and load packages

``` r
library(tidyverse)
library(ggplot2)
library(reshape2)
```

# 1 Load data

## 1.1 Load count data

Load in the count matrix we generated after kallisto pseudoalignment
using the Trinity abundance_estimates_to_matrix.pl script. We also need
to slightly reformat the count matrix

``` r
# Read in counts data. This is a gene-level counts matrix generated from kallisto transcript abundances using Trinity
Apul_counts_data_OG <- read_delim("../../../deep-dive/D-Apul/output/14-Apul-RNAseq-kallisto/kallisto/kallisto.isoform.counts.matrix") 
head(Apul_counts_data_OG)
```

    # A tibble: 6 × 6
      ...1      kallisto_quant_sampl…¹ kallisto_quant_sampl…² kallisto_quant_sampl…³
      <chr>                      <dbl>                  <dbl>                  <dbl>
    1 XM_04431…                  108                     46                    63   
    2 XR_00639…                    0                      0                     0   
    3 XM_02934…                  467.                   337.                  362.  
    4 XM_04432…                   43.5                   20.7                   5.14
    5 XR_00639…                  108.                     0                     0   
    6 XM_04431…                   20.4                    0                    21.3 
    # ℹ abbreviated names: ¹​kallisto_quant_sample140, ²​kallisto_quant_sample145,
    #   ³​kallisto_quant_sample150
    # ℹ 2 more variables: kallisto_quant_sample173 <dbl>,
    #   kallisto_quant_sample178 <dbl>

``` r
# Read in ID mapping of transcripts and associated GO terms etc.
Apul_IDmapping <- read_delim("..//output/02-Apul-reference-annotation/Apul_GCF_013753865.1_rna-IDmapping-2024_08_21.tab") %>%
  select(-...1)
head(Apul_IDmapping)
```

    # A tibble: 6 × 7
      V1             V3           V13 Protein.names  Organism Gene.Ontology..biolo…¹
      <chr>          <chr>      <dbl> <chr>          <chr>    <chr>                 
    1 XM_029323402.2 Q9XTR8 3.59e- 25 Lipase ZK262.… Caenorh… lipid catabolic proce…
    2 XM_029323410.2 Q9NUQ6 7.10e- 39 SPATS2-like p… Homo sa… <NA>                  
    3 XM_029323412.2 Q96EZ8 2.78e-134 Microspherule… Homo sa… chromatin remodeling …
    4 XM_029323419.2 P35991 2.61e-177 Tyrosine-prot… Mus mus… apoptotic process [GO…
    5 XM_029323420.2 P35991 1.57e-175 Tyrosine-prot… Mus mus… apoptotic process [GO…
    6 XM_029323421.2 Q3ZCL5 1.27e- 87 Arfaptin-2 (A… Bos tau… intracellular protein…
    # ℹ abbreviated name: ¹​Gene.Ontology..biological.process.
    # ℹ 1 more variable: Gene.Ontology.IDs <chr>

## 1.2 Count data munging

``` r
# We need to modify this data frame so that the row names are actually row names, instead of comprising the first column
Apul_counts_data <- Apul_counts_data_OG %>%
  column_to_rownames(var = "...1")

# Additional formatting
# Round all estimated counts to integers
Apul_counts_data <- round(Apul_counts_data, digits = 0)

# Remove all transcripts with 5 or fewer counts in all samples
Apul_counts_data <- Apul_counts_data[!apply(Apul_counts_data, 1, function(row) all(row < 6)), ]

# Remove the "kallisto_quant_" portion of the column names, to leave just the sample names
colnames(Apul_counts_data) <- sub("kallisto_quant_", "", colnames(Apul_counts_data))

# Reorder the columns into alphabetical order (to make it easier to create an associated metadata spreadsheet)
Apul_counts_data <- Apul_counts_data[, order(colnames(Apul_counts_data))]

Apul_sample_names <- names(Apul_counts_data)

head(Apul_counts_data)
```

                   sample140 sample145 sample150 sample173 sample178
    XM_044313592.1       108        46        63       105       131
    XR_006395886.1         0         0         0         6       139
    XM_029343024.2       467       337       362       570      1234
    XM_044327972.1        43        21         5        12        41
    XR_006392795.1       108         0         0         0         0
    XM_044315346.1        20         0        21        18       232

``` r
Apul_sample_names
```

    [1] "sample140" "sample145" "sample150" "sample173" "sample178"

``` r
Apul_counts_GO <- Apul_counts_data %>%
  rownames_to_column(var = "transcript") %>%
  left_join(Apul_IDmapping, by = c("transcript" = "V1"))

head(Apul_counts_GO)
```

          transcript sample140 sample145 sample150 sample173 sample178     V3 V13
    1 XM_044313592.1       108        46        63       105       131   <NA>  NA
    2 XR_006395886.1         0         0         0         6       139 Q8WZA2   0
    3 XM_029343024.2       467       337       362       570      1234 P16331   0
    4 XM_044327972.1        43        21         5        12        41   <NA>  NA
    5 XR_006392795.1       108         0         0         0         0   <NA>  NA
    6 XM_044315346.1        20         0        21        18       232   <NA>  NA
                                                                                                                                                                                                             Protein.names
    1                                                                                                                                                                                                                 <NA>
    2 Rap guanine nucleotide exchange factor 4 (Exchange factor directly activated by cAMP 2) (Exchange protein directly activated by cAMP 2) (EPAC 2) (cAMP-regulated guanine nucleotide exchange factor II) (cAMP-GEFII)
    3                                                                                                                                               Phenylalanine-4-hydroxylase (PAH) (EC 1.14.16.1) (Phe-4-monooxygenase)
    4                                                                                                                                                                                                                 <NA>
    5                                                                                                                                                                                                                 <NA>
    6                                                                                                                                                                                                                 <NA>
                  Organism
    1                 <NA>
    2 Homo sapiens (Human)
    3 Mus musculus (Mouse)
    4                 <NA>
    5                 <NA>
    6                 <NA>
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               Gene.Ontology..biological.process.
    1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        <NA>
    2 adaptive immune response [GO:0002250]; adenylate cyclase-activating G protein-coupled receptor signaling pathway [GO:0007189]; calcium-ion regulated exocytosis [GO:0017156]; G protein-coupled receptor signaling pathway [GO:0007186]; insulin secretion [GO:0030073]; positive regulation of GTPase activity [GO:0043547]; positive regulation of insulin secretion [GO:0032024]; Ras protein signal transduction [GO:0007265]; regulation of exocytosis [GO:0017157]; regulation of synaptic vesicle cycle [GO:0098693]
    3                                                                                                                                                                                                                                                                   L-phenylalanine catabolic process [GO:0006559]; L-phenylalanine metabolic process [GO:0006558]; protein hydroxylation [GO:0018126]; tyrosine biosynthetic process [GO:0006571]; tyrosine biosynthetic process, by oxidation of phenylalanine [GO:0019293]
    4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        <NA>
    5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        <NA>
    6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        <NA>
                                                                                                                                                                                                                       Gene.Ontology.IDs
    1                                                                                                                                                                                                                               <NA>
    2 GO:0002250; GO:0005085; GO:0005096; GO:0005829; GO:0005886; GO:0007186; GO:0007189; GO:0007265; GO:0016020; GO:0017156; GO:0017157; GO:0030073; GO:0030552; GO:0030674; GO:0031267; GO:0032024; GO:0043547; GO:0098686; GO:0098693
    3                                                                                                                         GO:0004505; GO:0005506; GO:0006558; GO:0006559; GO:0006571; GO:0016597; GO:0018126; GO:0019293; GO:0042802
    4                                                                                                                                                                                                                               <NA>
    5                                                                                                                                                                                                                               <NA>
    6                                                                                                                                                                                                                               <NA>

# 2 Summary stats and visualizations

## 2.1 Expression levels

Plot histograms of the expression levels in each sample

``` r
# Melt the count matrix into long format
Apul_counts_melted <- melt(Apul_counts_data, variable.name = "sample", value.name = "counts")

# Plot the expression level histograms for each sample
ggplot(Apul_counts_melted, aes(x = counts)) +
  geom_histogram(binwidth = 1, fill = "#408EC6", color = "black") +
  scale_x_log10() +  # Optional: Log-transform the x-axis for better visualization
  facet_wrap(~sample, scales = "free_y") +
  labs(title = "Gene Expression Level Histogram for Each Sample",
       x = "Expression Level (Counts)",
       y = "Frequency") +
  theme_minimal()
```

![](03-Apul-RNA-summary_files/figure-gfm/expression-level-histograms-1.png)<!-- -->

## 2.2 Transcript counts

First let’s check the total number of transcripts in each sample – keep
in mind this expression data has *not* been normalized yet, so there may
be different totals for each sample

``` r
# Calculate the total number of transcripts for each sample
total_transcripts <- colSums(Apul_counts_data)

# Create a data frame for plotting
total_transcripts_df <- data.frame(sample = names(total_transcripts),
                                   totals = total_transcripts)

# Plot the total number of transcripts for each sample
ggplot(total_transcripts_df, aes(x = sample, y = totals)) +
  geom_bar(stat = "identity", fill = "#408EC6", color = "black") +
  labs(title = "Total Number of Transcripts per Sample",
       x = "Sample",
       y = "Total Transcripts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
```

![](03-Apul-RNA-summary_files/figure-gfm/transcript-counts-plot-1.png)<!-- -->

Now let’s check the number of unique transcripts in each sample – that
is, how many genes are expressed in each sample? This should be pretty
much the same across samples, even without normalization.

``` r
# Calculate the number of unique transcripts (non-zero counts) for each sample
unique_transcripts <- colSums(Apul_counts_data > 0)

# Create a data frame for plotting
unique_transcripts_df <- data.frame(sample = names(unique_transcripts),
                                    uniques = unique_transcripts)

# Plot the total number of unique transcripts for each sample
ggplot(unique_transcripts_df, aes(x = sample, y = uniques)) +
  geom_bar(stat = "identity", fill = "#408EC6", color = "black") +
  labs(title = "Total Number of Unique Expressed Transcripts per Sample",
       x = "Sample",
       y = "Unique Transcripts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
```

![](03-Apul-RNA-summary_files/figure-gfm/total-unique-transcripts-plot-1.png)<!-- -->

## 2.3 Most common biological processes

Similar to the plot generated in `02-Apul-reference-annotation`, let’s
plot the biological processes most represented in these samples’
expression

``` r
# Rename the `Gene.Ontology..biological.process.` column to `Biological_Process`
colnames(Apul_counts_GO)[colnames(Apul_counts_GO) == "Gene.Ontology..biological.process."] <- "Biological_Process"

# Separate the `Biological_Process` column into individual biological processes
data_separated <- unlist(strsplit(Apul_counts_GO$Biological_Process, split = ";"))

# Trim whitespace from the biological processes
data_separated <- gsub("^\\s+|\\s+$", "", data_separated)

# Count the occurrences of each biological process
process_counts <- table(data_separated)
process_counts <- data.frame(Biological_Process = names(process_counts), Count = as.integer(process_counts))
process_counts <- process_counts[order(-process_counts$Count), ]

# Select the 20 most predominant biological processes
top_20_processes <- process_counts[1:20, ]

# Create a color palette for the bars
bar_colors <- rainbow(nrow(top_20_processes))

# Create a staggered vertical bar plot with different colors for each bar
barplot(top_20_processes$Count, names.arg = rep("", nrow(top_20_processes)), col = bar_colors,
        ylim = c(0, max(top_20_processes$Count) * 1.25),
        main = "Occurrences of the 20 Most Predominant Biological Processes", xlab = "Biological Process", ylab = "Count")
```

![](03-Apul-RNA-summary_files/figure-gfm/most-common-processes-1.png)<!-- -->

``` r
# Create a separate plot for the legend
png("../output/03-Apul-RNA-summary/GOlegend.png", width = 800, height = 600)
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = top_20_processes$Biological_Process, fill = bar_colors, cex = 1, title = "Biological Processes")
dev.off()
```

    png 
      2 

``` r
knitr::include_graphics("../output/03-Apul-RNA-summary/GOlegend.png")
```

<img src="../output/03-Apul-RNA-summary/GOlegend.png" width="800" />

``` bash
rm ../output/03-Apul-RNA-summary/GOlegend.png
```
