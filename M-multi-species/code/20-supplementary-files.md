20-supplementary-files
================
Kathleen Durkin
2026-06-15

- [1 Alignment summary tables](#1-alignment-summary-tables)
  - [1.1 RNA-seq](#11-rna-seq)
    - [1.1.1 1.1 Raw reads](#111-11-raw-reads)
    - [1.1.2 1.2 Trimmed reads](#112-12-trimmed-reads)
    - [1.1.3 1.3 Retained after trimming
      (%)](#113-13-retained-after-trimming-)
    - [1.1.4 1.4 MISSING HISAT2 alignment stats (overall rate, unique,
      multi)](#114-14-missing-hisat2-alignment-stats-overall-rate-unique-multi)
    - [1.1.5 1.5 MISSING Mapped reads (from featureCounts
      .summary)](#115-15-missing-mapped-reads-from-featurecounts-summary)
    - [1.1.6 1.6 MISSING Mismatch rate
      (%)](#116-16-missing-mismatch-rate-)
    - [1.1.7 1.7 Number of lncRNAs
      identified](#117-17-number-of-lncrnas-identified)
  - [1.2 sRNA-seq alignment
    statistics](#12-srna-seq-alignment-statistics)
    - [1.2.1 2.1 Species metadata
      (genome)](#121-21-species-metadata-genome)
    - [1.2.2 2.1 Raw reads](#122-21-raw-reads)
    - [1.2.3 2.2 Trimmed reads](#123-22-trimmed-reads)
    - [1.2.4 2.3 Retained after trimming
      (%)](#124-23-retained-after-trimming-)
    - [1.2.5 2.4 ShortStack alignment
      stats](#125-24-shortstack-alignment-stats)
    - [1.2.6 2.5 Mapped miRNA reads](#126-25-mapped-mirna-reads)
    - [1.2.7 2.6 Number of miRNAs
      identified](#127-26-number-of-mirnas-identified)

Code for any compiling/formating of supplementary figures/tables for the
DDE manuscript.

# 1 Alignment summary tables

Compiling per-sample alignment and processing statistics for all three
species

## 1.1 RNA-seq

Set up the metadata

``` r
species_meta <- tibble(
  species = c("D-Apul", "E-Peve", "F-Ptuh"),
  species_full = c("Acropora pulchra", "Porites evermanni",
                   "Pocillopora tuahiniensis"),
  genome = c("Apulcra-genome.fa", "Porites_evermanni_v1.fa",
             "Pocillopora_meandrina_HIv1.assembly.fasta"),
  rna_prefix = c("../../D-Apul", "../../E-Peve", "../../F-Ptuh")
)
species_meta
```

    ## # A tibble: 3 × 4
    ##   species species_full             genome                             rna_prefix
    ##   <chr>   <chr>                    <chr>                              <chr>     
    ## 1 D-Apul  Acropora pulchra         Apulcra-genome.fa                  ../../D-A…
    ## 2 E-Peve  Porites evermanni        Porites_evermanni_v1.fa            ../../E-P…
    ## 3 F-Ptuh  Pocillopora tuahiniensis Pocillopora_meandrina_HIv1.assemb… ../../F-P…

### 1.1.1 1.1 Raw reads

from MultiQC FastQC on raw reads

``` r
raw_reads <- rbind(
  
read_tsv("../../D-Apul/output/01-Apul-RNA-trimming-FastQC/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("../../E-Peve/output/01-Peve-RNA-trimming-FastQC/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("../../F-Ptuh/output/01-Ptuh-RNA-trimming-FastQC/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences')
)
```

    ## Rows: 18 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 20 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 20 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
raw_reads
```

    ## # A tibble: 15 × 2
    ##    Sample                    `Total Sequences`
    ##    <chr>                                 <dbl>
    ##  1 RNA-ACR-140-S1-TP2_R1_001          48199794
    ##  2 RNA-ACR-145-S1-TP2_R1_001          43288869
    ##  3 RNA-ACR-150-S1-TP2_R1_001          44216948
    ##  4 RNA-ACR-173-S1-TP2_R1_001          47976843
    ##  5 RNA-ACR-178-S1-TP2_R1_001          43185440
    ##  6 RNA-POR-71-S1-TP2_R1_001           51248209
    ##  7 RNA-POR-73-S1-TP2_R1_001           51750431
    ##  8 RNA-POR-76-S1-TP2_R1_001           50199897
    ##  9 RNA-POR-79-S1-TP2_R1_001           50299323
    ## 10 RNA-POR-82-S1-TP2_R1_001           49279402
    ## 11 RNA-POC-47-S1-TP2_R1_001           54705696
    ## 12 RNA-POC-48-S1-TP2_R1_001           52036501
    ## 13 RNA-POC-50-S1-TP2_R1_001           55792440
    ## 14 RNA-POC-53-S1-TP2_R1_001           53888583
    ## 15 RNA-POC-57-S1-TP2_R1_001           42934866

### 1.1.2 1.2 Trimmed reads

from MultiQC FastQC on trimmed reads

``` r
trimmed_reads <- rbind(
  
read_tsv("../../D-Apul/output/01-Apul-RNA-trimming-FastQC/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("../../E-Peve/output/01-Peve-RNA-trimming-FastQC/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("../../F-Ptuh/output/01-Ptuh-RNA-trimming-FastQC/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences')
)
```

    ## Rows: 10 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
trimmed_reads
```

    ## # A tibble: 15 × 2
    ##    Sample                    `Total Sequences`
    ##    <chr>                                 <dbl>
    ##  1 RNA-ACR-140-S1-TP2_R1_001          47710408
    ##  2 RNA-ACR-145-S1-TP2_R1_001          42864294
    ##  3 RNA-ACR-150-S1-TP2_R1_001          43712298
    ##  4 RNA-ACR-173-S1-TP2_R1_001          47501524
    ##  5 RNA-ACR-178-S1-TP2_R1_001          42677752
    ##  6 RNA-POR-71-S1-TP2_R1_001           50831351
    ##  7 RNA-POR-73-S1-TP2_R1_001           51385213
    ##  8 RNA-POR-76-S1-TP2_R1_001           49828147
    ##  9 RNA-POR-79-S1-TP2_R1_001           49976568
    ## 10 RNA-POR-82-S1-TP2_R1_001           48908730
    ## 11 RNA-POC-47-S1-TP2_R1_001           54219233
    ## 12 RNA-POC-48-S1-TP2_R1_001           51579450
    ## 13 RNA-POC-50-S1-TP2_R1_001           55302459
    ## 14 RNA-POC-53-S1-TP2_R1_001           53353608
    ## 15 RNA-POC-57-S1-TP2_R1_001           42513488

### 1.1.3 1.3 Retained after trimming (%)

``` r
read_counts <- raw_reads %>%
  left_join(trimmed_reads, by = c("Sample")) %>%
  mutate(retained_pct = round(trimmed_reads$'Total Sequences' / raw_reads$'Total Sequences' * 100, 2))

read_counts
```

    ## # A tibble: 15 × 4
    ##    Sample                   `Total Sequences.x` `Total Sequences.y` retained_pct
    ##    <chr>                                  <dbl>               <dbl>        <dbl>
    ##  1 RNA-ACR-140-S1-TP2_R1_0…            48199794            47710408         99.0
    ##  2 RNA-ACR-145-S1-TP2_R1_0…            43288869            42864294         99.0
    ##  3 RNA-ACR-150-S1-TP2_R1_0…            44216948            43712298         98.9
    ##  4 RNA-ACR-173-S1-TP2_R1_0…            47976843            47501524         99.0
    ##  5 RNA-ACR-178-S1-TP2_R1_0…            43185440            42677752         98.8
    ##  6 RNA-POR-71-S1-TP2_R1_001            51248209            50831351         99.2
    ##  7 RNA-POR-73-S1-TP2_R1_001            51750431            51385213         99.3
    ##  8 RNA-POR-76-S1-TP2_R1_001            50199897            49828147         99.3
    ##  9 RNA-POR-79-S1-TP2_R1_001            50299323            49976568         99.4
    ## 10 RNA-POR-82-S1-TP2_R1_001            49279402            48908730         99.2
    ## 11 RNA-POC-47-S1-TP2_R1_001            54705696            54219233         99.1
    ## 12 RNA-POC-48-S1-TP2_R1_001            52036501            51579450         99.1
    ## 13 RNA-POC-50-S1-TP2_R1_001            55792440            55302459         99.1
    ## 14 RNA-POC-53-S1-TP2_R1_001            53888583            53353608         99.0
    ## 15 RNA-POC-57-S1-TP2_R1_001            42934866            42513488         99.0

### 1.1.4 1.4 MISSING HISAT2 alignment stats (overall rate, unique, multi)

D-Apul has a consolidated log; E-Peve and F-Ptuh have per-sample stderr
files.

``` r
# Parse a single HISAT2 stderr block.
# Returns a one-row tibble with reads, overall_rate, unique_pct, multi_pct.
parse_hisat2_block <- function(lines) {
  text <- paste(lines, collapse = "\n")

  reads <- as.numeric(str_match(text, "^(\\d+) reads")[2])
  overall <- as.numeric(str_match(text, "(\\S+)% overall alignment rate")[2])

  # Concordant exactly 1 time  -> uniquely mapped
  unique_line <- str_match(text,
    "aligned concordantly exactly 1 time\n")[1]
  unique_pct <- as.numeric(str_extract(unique_line, "\\d+\\.\\d+"))

  # Concordant >1 times -> multi-mapped
  multi_line <- str_match(text,
    "aligned concordantly >1 times\n")[1]
  multi_pct <- as.numeric(str_extract(multi_line, "\\d+\\.\\d+"))

  tibble(
    hisat2_reads   = reads,
    overall_align  = overall,
    unique_mapped  = unique_pct,
    multi_mapped   = multi_pct
  )
}

# D-Apul: consolidated log, 5 blocks separated by blank lines
parse_hisat2_consolidated <- function(path, samples) {
  lines <- readLines(path)
  # Split into blocks at "N reads; of these:"
  block_starts <- which(str_detect(lines, "^\\d+ reads; of these:"))
  results <- list()
  for (i in seq_along(block_starts)) {
    end <- if (i < length(block_starts)) block_starts[i + 1] - 1 else length(lines)
    block <- lines[block_starts[i]:end]
    results[[i]] <- parse_hisat2_block(block) %>%
      mutate(sample = samples[i])
  }
  bind_rows(results)
}

# E-Peve / F-Ptuh: one stderr file per sample
parse_hisat2_per_file <- function(dir_path, sample_pattern, sample_replacement) {
  files <- list.files(dir_path, pattern = sample_pattern, full.names = TRUE)
  results <- list()
  for (f in files) {
    lines <- readLines(f)
    # Extract sample name from filename: RNA-POR-71_hisat.stderr -> POR-71
    sname <- basename(f) %>%
      str_remove("_hisat.stderr") %>%
      str_remove("^RNA-")
    results[[sname]] <- parse_hisat2_block(lines) %>%
      mutate(sample = sname)
  }
  bind_rows(results)
}

hisat2_stats <- bind_rows(
  parse_hisat2_consolidated(
    "../../D-Apul/output/07-Apul-Hisat/hisat.out",
    samples = c("ACR-140", "ACR-145", "ACR-150", "ACR-173", "ACR-178")
  ),
  parse_hisat2_per_file(
    "../../E-Peve/output/06-Peve-Hisat",
    pattern = "_hisat.stderr$"
  ),
  parse_hisat2_per_file(
    "../../F-Ptuh/output/06-Ptuh-Hisat",
    pattern = "_hisat.stderr$"
  )
) %>%
  mutate(species = case_when(
    str_detect(sample, "^ACR") ~ "D-Apul",
    str_detect(sample, "^POR") ~ "E-Peve",
    str_detect(sample, "^POC") ~ "F-Ptuh"
  )) %>%
  select(species, sample, everything())

hisat2_stats
```

### 1.1.5 1.5 MISSING Mapped reads (from featureCounts .summary)

``` r
# featureCounts .summary: rows = Status (Assigned, Unassigned_Unmapped, etc.),
# columns = BAM paths. Mapped reads = total - Unassigned_Unmapped.
# Total reads per sample = sum of all rows.

read_featurecounts_summary <- function(path, species_code) {
  df <- read_tsv(path, show_col_types = FALSE)
  # Column names are BAM paths; extract sample names
  bam_cols <- colnames(df)[-1]  # first col is "Status"
  samples <- str_extract(bam_cols, "(ACR|POR|POC)-\\d+")

  assigned <- df %>% filter(Status == "Assigned") %>% select(-Status) %>% as.numeric()
  unmapped <- df %>% filter(Status == "Unassigned_Unmapped") %>% select(-Status) %>% as.numeric()

  # Total = sum of all status rows
  total <- colSums(df %>% select(-Status) %>% mutate(across(everything, as.numeric)))

  tibble(
    species = species_code,
    sample = samples,
    mapped_reads = total - unmapped,
    total_reads_fc = total
  )
}

mapped_reads <- bind_rows(
  read_featurecounts_summary("../../D-Apul/output/32-Apul-lncRNA-matrix/Apul-lncRNA-counts.txt.summary", "D-Apul"),
  read_featurecounts_summary("../../E-Peve/output/18-Peve-lncRNA-matrix/Peve-lncRNA-counts.txt.summary", "E-Peve"),
  read_featurecounts_summary("../../F-Ptuh/output/18-Ptuh-lncRNA-matrix/Ptuh-lncRNA-counts.txt.summary", "F-Ptuh")
)
mapped_reads
```

### 1.1.6 1.6 MISSING Mismatch rate (%)

``` r
# NOTE: Mismatch rate is NOT available in this repository.
# There are no samtools stats / flagstat / idxstats or Picard outputs.
# The sorted BAMs are gitignored (available via the project wiki).
#
# To compute mismatch rate, run the following on the BAMs:
#
#   samtools stats <sample.sorted.bam> | grep "^SN" > <sample>.stats.txt
#
# Then parse the line:
#   SN  mismatch rate:  X.XXe-XX
# and convert to percentage.
#
# Placeholder:
mismatch_rate <- tibble(
  species = character(),
  sample = character(),
  mismatch_rate_pct = numeric()
)
# mismatch_rate  # leave empty — populate after running samtools stats on BAMs
```

### 1.1.7 1.7 Number of lncRNAs identified

``` r
count_lncRNAs <- function(filepath) {
  
df <- read_tsv(filepath, col_names = TRUE)
sample_cols <- grep(".sorted.bam$", colnames(df), value = TRUE)
sample_names <- str_remove(sample_cols, ".*pipeline\\.RNA\\.") %>% str_remove("\\.sorted\\.bam$")

# Count non-zero entries per sample column
nonzero_counts <- sapply(sample_cols, function(col) {
  sum(df[[col]] > 0, na.rm = TRUE)
})

# Build a summary tibble
summary <- tibble(
  sample = sample_names,
  n_unique_lncRNA = nonzero_counts
)
print(summary)
}

lncRNA_counts <- rbind(
count_lncRNAs("../output/01.6-lncRNA-pipeline/Apul-lncRNA-counts-filtered.txt"),
count_lncRNAs("../output/01.6-lncRNA-pipeline/Peve-lncRNA-counts-filtered.txt"),
count_lncRNAs("../output/01.6-lncRNA-pipeline/Ptuh-lncRNA-counts-filtered.txt")
)
```

    ## Rows: 31490 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): Geneid, Chr, Strand
    ## dbl (8): Start, End, Length, X.home.shared.8TB_HDD_02.zbengt.github.deep.div...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 5 × 2
    ##   sample  n_unique_lncRNA
    ##   <chr>             <int>
    ## 1 ACR.140           27515
    ## 2 ACR.145           21807
    ## 3 ACR.150           28248
    ## 4 ACR.173           24585
    ## 5 ACR.178           29713

    ## Rows: 10085 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): Geneid, Chr, Strand
    ## dbl (8): Start, End, Length, X.home.shared.8TB_HDD_02.zbengt.github.deep.div...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 5 × 2
    ##   sample n_unique_lncRNA
    ##   <chr>            <int>
    ## 1 POR.71            8577
    ## 2 POR.73            7763
    ## 3 POR.76            9498
    ## 4 POR.79            8007
    ## 5 POR.82            8681

    ## Rows: 16152 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): Geneid, Chr, Strand
    ## dbl (8): Start, End, Length, ...output.01.6.Ptuh.lncRNA.pipeline.RNA.POC.47....
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 5 × 2
    ##   sample        n_unique_lncRNA
    ##   <chr>                   <int>
    ## 1 POC.47.S1.TP2           12017
    ## 2 POC.48.S1.TP2           13866
    ## 3 POC.50.S1.TP2           14294
    ## 4 POC.53.S1.TP2           14409
    ## 5 POC.57.S1.TP2           12197

``` r
lncRNA_counts
```

    ## # A tibble: 15 × 2
    ##    sample        n_unique_lncRNA
    ##    <chr>                   <int>
    ##  1 ACR.140                 27515
    ##  2 ACR.145                 21807
    ##  3 ACR.150                 28248
    ##  4 ACR.173                 24585
    ##  5 ACR.178                 29713
    ##  6 POR.71                   8577
    ##  7 POR.73                   7763
    ##  8 POR.76                   9498
    ##  9 POR.79                   8007
    ## 10 POR.82                   8681
    ## 11 POC.47.S1.TP2           12017
    ## 12 POC.48.S1.TP2           13866
    ## 13 POC.50.S1.TP2           14294
    ## 14 POC.53.S1.TP2           14409
    ## 15 POC.57.S1.TP2           12197

## 1.2 sRNA-seq alignment statistics

### 1.2.1 2.1 Species metadata (genome)

``` r
# Genome files are referenced in the ShortStack Rmd scripts:
srna_genomes <- tibble(
  species = c("D-Apul", "E-Peve", "F-Ptuh"),
  genome = c("Apulcra-genome.fa", "Porites_evermanni_v1.fa",
             "Pocillopora_meandrina_HIv1.assembly.fasta")
)
srna_genomes
```

    ## # A tibble: 3 × 2
    ##   species genome                                   
    ##   <chr>   <chr>                                    
    ## 1 D-Apul  Apulcra-genome.fa                        
    ## 2 E-Peve  Porites_evermanni_v1.fa                  
    ## 3 F-Ptuh  Pocillopora_meandrina_HIv1.assembly.fasta

### 1.2.2 2.1 Raw reads

``` r
s_raw_reads <- rbind(
  
read_tsv("~/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("~/deep-dive/E-Peve/output/06-Peve-sRNAseq-trimming/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("~/deep-dive/F-Pmea/output/08-Pmea-sRNAseq-trimming/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences')
)
```

    ## Rows: 10 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 6 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
s_raw_reads$Sample <- str_remove(s_raw_reads$Sample, "_R1_001")

s_raw_reads
```

    ## # A tibble: 13 × 2
    ##    Sample              `Total Sequences`
    ##    <chr>                           <dbl>
    ##  1 sRNA-ACR-140-S1-TP2          18506972
    ##  2 sRNA-ACR-145-S1-TP2          21295979
    ##  3 sRNA-ACR-150-S1-TP2          22586680
    ##  4 sRNA-ACR-173-S1-TP2          19686271
    ##  5 sRNA-ACR-178-S1-TP2          17687362
    ##  6 sRNA-POR-73-S1-TP2           22960890
    ##  7 sRNA-POR-79-S1-TP2           19762559
    ##  8 sRNA-POR-82-S1-TP2           20675792
    ##  9 sRNA-POC-47-S1-TP2           27579504
    ## 10 sRNA-POC-48-S1-TP2           24884997
    ## 11 sRNA-POC-50-S1-TP2           26073937
    ## 12 sRNA-POC-53-S1-TP2           25726329
    ## 13 sRNA-POC-57-S1-TP2           26459749

### 1.2.3 2.2 Trimmed reads

``` r
s_trimmed_reads <- rbind(
  
read_tsv("~/deep-dive/D-Apul/output/08.2-Apul-sRNAseq-trimming-31bp-fastp-merged/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "sRNA-")) %>% select(Sample, 'Total Sequences'), 

read_tsv("~/deep-dive/E-Peve/output/06.2-Peve-sRNAseq-trimming-31bp-fastp-merged/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% select(Sample, 'Total Sequences') %>% mutate(Sample = paste0("sRNA-", Sample)), 

read_tsv("~/deep-dive/F-Pmea/output/08.2-Pmea-sRNAseq-trimming-31bp-fastp-merged/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "sRNA-")) %>% select(Sample, 'Total Sequences')
)
```

    ## Rows: 5 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 3 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 5 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
s_trimmed_reads$Sample <- str_remove(s_trimmed_reads$Sample, "-fastp-adapters-polyG-31bp-merged")

s_trimmed_reads
```

    ## # A tibble: 13 × 2
    ##    Sample              `Total Sequences`
    ##    <chr>                           <dbl>
    ##  1 sRNA-ACR-140-S1-TP2          15076001
    ##  2 sRNA-ACR-145-S1-TP2          17265123
    ##  3 sRNA-ACR-150-S1-TP2          18395038
    ##  4 sRNA-ACR-173-S1-TP2          16332044
    ##  5 sRNA-ACR-178-S1-TP2          14295261
    ##  6 sRNA-POR-73-S1-TP2           12238952
    ##  7 sRNA-POR-79-S1-TP2           12259576
    ##  8 sRNA-POR-82-S1-TP2           13933812
    ##  9 sRNA-POC-47-S1-TP2           12066776
    ## 10 sRNA-POC-48-S1-TP2           14045042
    ## 11 sRNA-POC-50-S1-TP2           12141955
    ## 12 sRNA-POC-53-S1-TP2           15732665
    ## 13 sRNA-POC-57-S1-TP2           14136486

### 1.2.4 2.3 Retained after trimming (%)

``` r
s_read_counts <- s_raw_reads %>%
  left_join(s_trimmed_reads, by = c("Sample")) %>%
  mutate(retained_pct = round(s_trimmed_reads$'Total Sequences' / s_raw_reads$'Total Sequences' * 100, 2))

s_read_counts
```

    ## # A tibble: 13 × 4
    ##    Sample              `Total Sequences.x` `Total Sequences.y` retained_pct
    ##    <chr>                             <dbl>               <dbl>        <dbl>
    ##  1 sRNA-ACR-140-S1-TP2            18506972            15076001         81.5
    ##  2 sRNA-ACR-145-S1-TP2            21295979            17265123         81.1
    ##  3 sRNA-ACR-150-S1-TP2            22586680            18395038         81.4
    ##  4 sRNA-ACR-173-S1-TP2            19686271            16332044         83.0
    ##  5 sRNA-ACR-178-S1-TP2            17687362            14295261         80.8
    ##  6 sRNA-POR-73-S1-TP2             22960890            12238952         53.3
    ##  7 sRNA-POR-79-S1-TP2             19762559            12259576         62.0
    ##  8 sRNA-POR-82-S1-TP2             20675792            13933812         67.4
    ##  9 sRNA-POC-47-S1-TP2             27579504            12066776         43.8
    ## 10 sRNA-POC-48-S1-TP2             24884997            14045042         56.4
    ## 11 sRNA-POC-50-S1-TP2             26073937            12141955         46.6
    ## 12 sRNA-POC-53-S1-TP2             25726329            15732665         61.2
    ## 13 sRNA-POC-57-S1-TP2             26459749            14136486         53.4

### 1.2.5 2.4 ShortStack alignment stats

from alignment_details.tsv

NOTE: the Github repo README.md states that the H tag means “H: Very
highly multi-mapped read (\>=50 hits)”. However,the log from running
ShortStack 4.1.0 (e.g., the file
~/deep-dive-expression/D-Apul/output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/shortstack.log)
states that the tag represents “Very highly multi-mapped (\>=20
hits)(H)”. It looks like the \>=50 threshold was used in previous
versions of ShortStack (confirmed by looking at deep-dive ShortStack
logs), but the updated 4.1.0 version switched to a \>=20 threshold. The
Github repo’s README is out of date.

``` r
# alignment_details.tsv columns:
#   readfile  mapping_type  sequence_length  sequence_count  read_count
# mapping_type (from https://github.com/MikeAxtell/ShortStack)
# U: Uniquely mapped (not a multimapper).
# P: Multimapper placed using the method set by option --mmap.
# R: Multimapper placed at random.
# H: Very highly multi-mapped read (>=20 hits).
# N: Unmapped reads.

read_alignment_details <- function(path, species_code) {
  df <- read_tsv(path, show_col_types = FALSE)

  # Extract sample name from readfile path
  df <- df %>%
    mutate(
      sample = str_extract(readfile, "(ACR|POR|POC)-\\d+"),
      species = species_code
    )

  # Aggregate by sample × mapping_type, sum sequence_count
  agg <- df %>%
    group_by(species, sample, mapping_type) %>%
    summarise(seq_count = sum(read_count), .groups = "drop")

  # Pivot to wide format
  agg_wide <- agg %>%
    pivot_wider(names_from = mapping_type, values_from = seq_count,
                values_fill = 0)

  # totals and percentages
  # U=unique, P=guided, R=random, H=very-high, N=not mapped
  agg_wide %>%
    mutate(
      total_reads = U + P + R + H + N,
      mapped_reads = U + P + R + H,
      overall_align_pct = round((mapped_reads) / total_reads * 100, 2),
      unique_pct = round(U / total_reads * 100, 2),
      guided_pct = round(P / total_reads * 100, 2),
      random_pct = round(R / total_reads * 100, 2),
      very_high_pct = round(H / total_reads * 100, 2),
    ) %>%
    select(
      species, sample, total_reads, overall_align_pct, unique_pct, guided_pct, 
      random_pct, very_high_pct, mapped_reads
    ) %>%
    rename(
      `Overall alignment (%)` = overall_align_pct,
      `Uniquely mapped reads (%)` = unique_pct,
      `Multi-mapped reads placed with guidance (%)` = guided_pct,
      `Multi-mapped reads placed randomly (%)` = random_pct,
      `Very high multi-mapping reads (>=50 hits; %)` = very_high_pct,
      `Mapped sRNA reads` = mapped_reads
    )
}

srna_align <- bind_rows(
  read_alignment_details("../../D-Apul/output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out/alignment_details.tsv", "D-Apul"),
  read_alignment_details("../../E-Peve/output/05-Peve-sRNA-ShortStack_4.1.0/ShortStack_out/alignment_details.tsv", "E-Peve"),
  read_alignment_details("../../F-Ptuh/output/05-Ptuh-sRNA-ShortStack_4.1.0/ShortStack_out/alignment_details.tsv", "F-Ptuh")
)
srna_align
```

    ## # A tibble: 13 × 9
    ##    species sample  total_reads `Overall alignment (%)` Uniquely mapped reads (…¹
    ##    <chr>   <chr>         <dbl>                   <dbl>                     <dbl>
    ##  1 D-Apul  ACR-140    15076180                    80.0                      33.2
    ##  2 D-Apul  ACR-145    17265150                    78.8                      28.9
    ##  3 D-Apul  ACR-150    18395068                    77.8                      35.0
    ##  4 D-Apul  ACR-173    16332641                    78.1                      31.1
    ##  5 D-Apul  ACR-178    14295278                    80.1                      31.9
    ##  6 E-Peve  POR-73     12238972                    73.6                      30.7
    ##  7 E-Peve  POR-79     12259629                    79.3                      40.3
    ##  8 E-Peve  POR-82     13933829                    81.4                      43.7
    ##  9 F-Ptuh  POC-47     12068141                    62.4                      17.8
    ## 10 F-Ptuh  POC-48     14046678                    62.4                      12.9
    ## 11 F-Ptuh  POC-50     12142097                    55.9                      16.2
    ## 12 F-Ptuh  POC-53     15732765                    75.8                      18.3
    ## 13 F-Ptuh  POC-57     14137264                    83.2                      16.8
    ## # ℹ abbreviated name: ¹​`Uniquely mapped reads (%)`
    ## # ℹ 4 more variables: `Multi-mapped reads placed with guidance (%)` <dbl>,
    ## #   `Multi-mapped reads placed randomly (%)` <dbl>,
    ## #   `Very high multi-mapping reads (>=50 hits; %)` <dbl>,
    ## #   `Mapped sRNA reads` <dbl>

### 1.2.6 2.5 Mapped miRNA reads

from ShortStack Counts.txt

``` r
# Counts.txt is a TSV: Coords, Name, MIRNA, <sample1>_condensed, ...
# Filter rows where MIRNA == "Y", then sum each sample column.

read_mirna_counts <- function(path, species_code) {
  df <- read_tsv(path, show_col_types = FALSE)
  sample_cols <- colnames(df)[-(1:3)]  # skip Coords, Name, MIRNA

  # Extract clean sample names from column names
  # e.g. "sRNA-ACR-140-S1-TP2-..._condensed" -> "ACR-140"
  samples <- str_extract(sample_cols, "(ACR|POR|POC)-\\d+")

  mirna_only <- df %>% filter(MIRNA == "Y")

  # Sum each sample column across all miRNA rows
  mapped_mirna <- colSums(mirna_only[, sample_cols])

  tibble(
    species = species_code,
    sample = samples,
    mapped_mirna_reads = mapped_mirna
  )
}

srna_mirna_reads <- bind_rows(
  read_mirna_counts("../../D-Apul/output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out/Counts.txt", "D-Apul"),
  read_mirna_counts("../../E-Peve/output/05-Peve-sRNA-ShortStack_4.1.0/ShortStack_out/Counts.txt", "E-Peve"),
  read_mirna_counts("../../F-Ptuh/output/05-Ptuh-sRNA-ShortStack_4.1.0/ShortStack_out/Counts.txt", "F-Ptuh")
)
srna_mirna_reads
```

    ## # A tibble: 13 × 3
    ##    species sample  mapped_mirna_reads
    ##    <chr>   <chr>                <dbl>
    ##  1 D-Apul  ACR-140             226397
    ##  2 D-Apul  ACR-145             444883
    ##  3 D-Apul  ACR-150             668500
    ##  4 D-Apul  ACR-173             373547
    ##  5 D-Apul  ACR-178             455227
    ##  6 E-Peve  POR-73              373247
    ##  7 E-Peve  POR-79              279129
    ##  8 E-Peve  POR-82              399594
    ##  9 F-Ptuh  POC-47              623726
    ## 10 F-Ptuh  POC-48              368009
    ## 11 F-Ptuh  POC-50              512651
    ## 12 F-Ptuh  POC-53              430756
    ## 13 F-Ptuh  POC-57              178801

### 1.2.7 2.6 Number of miRNAs identified

``` r
count_expressed_mirnas <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)

  sample_cols <- grep("merged_condensed", colnames(df), value = TRUE)

  as.data.frame(t(df %>%
    filter(MIRNA == "Y") %>%
    summarise(across(all_of(sample_cols), ~ sum(.x > 0, na.rm = TRUE)))))
}

miRNA_counts <- rbind(
count_expressed_mirnas("../../D-Apul/output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out/Counts.txt"),
count_expressed_mirnas("../../E-Peve/output/05-Peve-sRNA-ShortStack_4.1.0/ShortStack_out/Counts.txt"),
count_expressed_mirnas("../../F-Ptuh/output/05-Ptuh-sRNA-ShortStack_4.1.0/ShortStack_out/Counts.txt")
)
colnames(miRNA_counts) <- c("miRNA_count")
miRNA_counts$Sample <- rownames(miRNA_counts)
miRNA_counts <- miRNA_counts %>% select(Sample, miRNA_count) %>% remove_rownames()
miRNA_counts$Sample <- miRNA_counts$Sample %>% str_remove("-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed") %>% str_remove("sRNA-")

miRNA_counts
```

    ##     Sample miRNA_count
    ## 1  ACR-140          38
    ## 2  ACR-145          39
    ## 3  ACR-150          39
    ## 4  ACR-173          39
    ## 5  ACR-178          38
    ## 6   POR-73          44
    ## 7   POR-79          42
    ## 8   POR-82          43
    ## 9   POC-47          37
    ## 10  POC-48          37
    ## 11  POC-50          37
    ## 12  POC-53          37
    ## 13  POC-57          37
