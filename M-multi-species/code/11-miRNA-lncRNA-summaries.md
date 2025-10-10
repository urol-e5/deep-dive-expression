11-miRNA-lncRNA-summaries
================
Kathleen Durkin
2025-10-09

- <a href="#01-apul" id="toc-01-apul">0.1 Apul</a>
- <a href="#02-peve" id="toc-02-peve">0.2 Peve</a>
- <a href="#03-ptuh" id="toc-03-ptuh">0.3 Ptuh</a>

Load packages

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
library(readr)
library(stringr)
library(tidyr)
```

## 0.1 Apul

Reminder that miRanda output files are formatted as such: “miRNA”,
“lncRNA”, “score”, “energy”, “query_start_end”, “subject_start_end”,
“total_bp_shared”, “query_similar”, “subject_similar”

Jill has previously suggested filtering miRNA-lncRNA binding results to
retain only instances in which at least 75% of the miRNA query is
complementary, so I will include this as a filter option in the below
summary stats (query_similar \> 75%). Note, however, that the query may
not be the entire miRNA, it is the portion of the miRNA that is
predicted to bind. So, for example, the query may be positions 2 to 10
on the miRNA (likely the seed) and have 100% query similarity. This
would indicate that those 8 nucleotides of the miRNA are perfectly
complementary to a set of 8 nucleotides in the lncRNA, NOT that the
entire \~22nt long miRNA is perfectly complementary. For now, though,
I’ll include filters for just “query_similar” because that’s what Jill
used. We shouldn’t need to use that filter in final results, because the
join of predicted binding + coexpression is already a sufficient
reduction of network size.

``` r
apul_miranda <- "../../D-Apul/output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict-parsed.txt"
apul_pcc_link <- "https://gannet.fish.washington.edu/kdurkin1/ravenbackups/deep-dive-expression/D-Apul/output/28-Apul-miRNA-lncRNA-interactions/Apul-PCC_miRNA_lncRNA.csv"
apul_pcc_loc <- "../../D-Apul/output/28-Apul-miRNA-lncRNA-interactions/Apul-PCC_miRNA_lncRNA.csv"
apul_miranda_pcc <- "../../D-Apul/output/28-Apul-miRNA-lncRNA-interactions/miranda_PCC_miRNA_lncRNA.csv"

# ---- Optionally PCC download file if missing ----
if (!file.exists(apul_pcc_loc)) {
  download.file(apul_pcc_link, apul_pcc_loc)
}

# ---- miRanda binding interactions ----
cat("\nNum. of miRNA-lncRNA binding interactions (miRanda):\n")
```

    ## 
    ## Num. of miRNA-lncRNA binding interactions (miRanda):

``` r
num_clusters <- sum(grepl("^>Cluster", readLines(apul_miranda)))
print(num_clusters)
```

    ## [1] 12390

``` r
cat("\nNum. of miRNA-lncRNA bindings >75% (miRanda):\n")
```

    ## 
    ## Num. of miRNA-lncRNA bindings >75% (miRanda):

``` r
miranda_df <- read.table(apul_miranda, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
miranda_df$V8 <- as.numeric(str_remove(miranda_df$V8, "%"))
num_75 <- sum(miranda_df$V8 >= 75, na.rm = TRUE)
print(num_75)
```

    ## [1] 6620

``` r
# ---- Significant PCC correlations ----
pcc_df <- read_csv(apul_pcc_loc, show_col_types = FALSE)
```

    ## New names:
    ## • `` -> `...1`

``` r
count_sig <- sum(pcc_df[[5]] < 0.05, na.rm = TRUE)
count_sig_pos <- sum(pcc_df[[5]] < 0.05 & pcc_df[[4]] > 0, na.rm = TRUE)
percent <- if (count_sig > 0) round((count_sig_pos / count_sig) * 100, 2) else 0

cat("\nSignificant PCC correlations + % positive\n")
```

    ## 
    ## Significant PCC correlations + % positive

``` r
cat("Significant PCC (p<0.05):       ", count_sig, "\n")
```

    ## Significant PCC (p<0.05):        62486

``` r
cat("Significant (p<0.05), pos. PCC: ", count_sig_pos, "\n")
```

    ## Significant (p<0.05), pos. PCC:  47187

``` r
cat("Percent positive:               ", percent, "%\n")
```

    ## Percent positive:                75.52 %

``` r
# ---- Overlap (binding AND sig. PCC) ----
miranda_pcc <- read_csv(apul_miranda_pcc, show_col_types = FALSE)
```

    ## New names:
    ## • `` -> `...1`

``` r
# Clean columns 9, 12, 13 (remove quotes/percent)
miranda_pcc <- miranda_pcc %>%
  mutate(
    across(c(9, 12, 13), ~ as.numeric(str_remove_all(as.character(.x), '["%]')))
  )

count_sig2 <- sum(miranda_pcc[[13]] < 0.05, na.rm = TRUE)
count_sig_pos2 <- sum(miranda_pcc[[13]] < 0.05 & miranda_pcc[[12]] > 0, na.rm = TRUE)
percent2 <- if (count_sig2 > 0) round((count_sig_pos2 / count_sig2) * 100, 2) else 0

cat("\nOverlap (binding AND sig. PCC) + % positive\n")
```

    ## 
    ## Overlap (binding AND sig. PCC) + % positive

``` r
cat("Significant PCC (p<0.05):       ", count_sig2, "\n")
```

    ## Significant PCC (p<0.05):        564

``` r
cat("Significant (p<0.05), pos. PCC: ", count_sig_pos2, "\n")
```

    ## Significant (p<0.05), pos. PCC:  433

``` r
cat("Percent positive:               ", percent2, "%\n")
```

    ## Percent positive:                76.77 %

``` r
# ---- Overlap >75% (binding >75% AND PCC) ----
count_sig3 <- sum(miranda_pcc[[13]] < 0.05 & miranda_pcc[[9]] > 75, na.rm = TRUE)
count_sig_pos3 <- sum(miranda_pcc[[13]] < 0.05 & miranda_pcc[[12]] > 0 & miranda_pcc[[9]] > 75, na.rm = TRUE)
percent3 <- if (count_sig3 > 0) round((count_sig_pos3 / count_sig3) * 100, 2) else 0

cat("\nOverlap >75% (binding >75% AND PCC) + % positive\n")
```

    ## 
    ## Overlap >75% (binding >75% AND PCC) + % positive

``` r
cat("Significant PCC (p<0.05):       ", count_sig3, "\n")
```

    ## Significant PCC (p<0.05):        268

``` r
cat("Significant (p<0.05), pos. PCC: ", count_sig_pos3, "\n")
```

    ## Significant (p<0.05), pos. PCC:  191

``` r
cat("Percent positive:               ", percent3, "%\n")
```

    ## Percent positive:                71.27 %

## 0.2 Peve

``` r
peve_miranda <- "../../E-Peve/output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve-miRanda-lncRNA-strict-parsed.txt"
peve_pcc_link <- "https://gannet.fish.washington.edu/kdurkin1/ravenbackups/deep-dive-expression/E-Peve/output/15-Peve-miRNA-lncRNA-PCC/PCC_miRNA_lncRNA.csv"
peve_pcc_loc <- "../../E-Peve/output/15-Peve-miRNA-lncRNA-PCC/PCC_miRNA_lncRNA.csv"
peve_miranda_pcc <- "../../E-Peve/output/15-Peve-miRNA-lncRNA-PCC/miranda_PCC_miRNA_lncRNA.csv"

# ---- Optionally PCC download file if missing ----
if (!file.exists(peve_pcc_loc)) {
  download.file(peve_pcc_link, peve_pcc_loc)
}


# ---- miRanda binding interactions ----
cat("\nNum. of miRNA-lncRNA binding interactions (miRanda):\n")
```

    ## 
    ## Num. of miRNA-lncRNA binding interactions (miRanda):

``` r
num_clusters <- sum(grepl("^>Cluster", readLines(peve_miranda)))
print(num_clusters)
```

    ## [1] 4116

``` r
cat("\nNum. of miRNA-lncRNA bindings >75% (miRanda):\n")
```

    ## 
    ## Num. of miRNA-lncRNA bindings >75% (miRanda):

``` r
miranda_df <- read.table(peve_miranda, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
miranda_df$V8 <- as.numeric(str_remove(miranda_df$V8, "%"))
num_75 <- sum(miranda_df$V8 >= 75, na.rm = TRUE)
print(num_75)
```

    ## [1] 2111

``` r
# ---- Significant PCC correlations ----
pcc_df <- read_csv(peve_pcc_loc, show_col_types = FALSE)
```

    ## New names:
    ## • `` -> `...1`

``` r
count_sig <- sum(pcc_df[[5]] < 0.05, na.rm = TRUE)
count_sig_pos <- sum(pcc_df[[5]] < 0.05 & pcc_df[[4]] > 0, na.rm = TRUE)
percent <- if (count_sig > 0) round((count_sig_pos / count_sig) * 100, 2) else 0

cat("\nSignificant PCC correlations + % positive\n")
```

    ## 
    ## Significant PCC correlations + % positive

``` r
cat("Significant PCC (p<0.05):       ", count_sig, "\n")
```

    ## Significant PCC (p<0.05):        18565

``` r
cat("Significant (p<0.05), pos. PCC: ", count_sig_pos, "\n")
```

    ## Significant (p<0.05), pos. PCC:  9835

``` r
cat("Percent positive:               ", percent, "%\n")
```

    ## Percent positive:                52.98 %

``` r
# ---- Overlap (binding AND sig. PCC) ----
miranda_pcc <- read_csv(peve_miranda_pcc, show_col_types = FALSE)
```

    ## New names:
    ## • `` -> `...1`

``` r
# Clean columns 9, 12, 13 (remove quotes/percent)
miranda_pcc <- miranda_pcc %>%
  mutate(
    across(c(9, 12, 13), ~ as.numeric(str_remove_all(as.character(.x), '["%]')))
  )

count_sig2 <- sum(miranda_pcc[[13]] < 0.05, na.rm = TRUE)
count_sig_pos2 <- sum(miranda_pcc[[13]] < 0.05 & miranda_pcc[[12]] > 0, na.rm = TRUE)
percent2 <- if (count_sig2 > 0) round((count_sig_pos2 / count_sig2) * 100, 2) else 0

cat("\nOverlap (binding AND sig. PCC) + % positive\n")
```

    ## 
    ## Overlap (binding AND sig. PCC) + % positive

``` r
cat("Significant PCC (p<0.05):       ", count_sig2, "\n")
```

    ## Significant PCC (p<0.05):        175

``` r
cat("Significant (p<0.05), pos. PCC: ", count_sig_pos2, "\n")
```

    ## Significant (p<0.05), pos. PCC:  106

``` r
cat("Percent positive:               ", percent2, "%\n")
```

    ## Percent positive:                60.57 %

``` r
# ---- Overlap >75% (binding >75% AND PCC) ----
count_sig3 <- sum(miranda_pcc[[13]] < 0.05 & miranda_pcc[[9]] > 75, na.rm = TRUE)
count_sig_pos3 <- sum(miranda_pcc[[13]] < 0.05 & miranda_pcc[[12]] > 0 & miranda_pcc[[9]] > 75, na.rm = TRUE)
percent3 <- if (count_sig3 > 0) round((count_sig_pos3 / count_sig3) * 100, 2) else 0

cat("\nOverlap >75% (binding >75% AND PCC) + % positive\n")
```

    ## 
    ## Overlap >75% (binding >75% AND PCC) + % positive

``` r
cat("Significant PCC (p<0.05):       ", count_sig3, "\n")
```

    ## Significant PCC (p<0.05):        86

``` r
cat("Significant (p<0.05), pos. PCC: ", count_sig_pos3, "\n")
```

    ## Significant (p<0.05), pos. PCC:  44

``` r
cat("Percent positive:               ", percent3, "%\n")
```

    ## Percent positive:                51.16 %

## 0.3 Ptuh

``` r
ptuh_miranda <- "../../F-Ptuh/output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh-miRanda-lncRNA-strict-parsed.txt"
ptuh_pcc_link <- "https://gannet.fish.washington.edu/kdurkin1/ravenbackups/deep-dive-expression/F-Ptuh/output/15-Ptuh-miRNA-lncRNA-PCC/PCC_miRNA_lncRNA.csv"
ptuh_pcc_loc <- "../../F-Ptuh/output/15-Ptuh-miRNA-lncRNA-PCC/PCC_miRNA_lncRNA.csv"
ptuh_miranda_pcc <- "../../F-Ptuh/output/15-Ptuh-miRNA-lncRNA-PCC/miranda_PCC_miRNA_lncRNA.csv"

# ---- Optionally PCC download file if missing ----
if (!file.exists(ptuh_pcc_loc)) {
  download.file(ptuh_pcc_link, ptuh_pcc_loc)
}

# ---- File line counts ----
cat("Num. of lines in miranda_PCC file:\n")
```

    ## Num. of lines in miranda_PCC file:

``` r
print(length(readLines(ptuh_miranda_pcc)))
```

    ## [1] 8342

``` r
# ---- miRanda binding interactions ----
cat("\nNum. of miRNA-lncRNA binding interactions (miRanda):\n")
```

    ## 
    ## Num. of miRNA-lncRNA binding interactions (miRanda):

``` r
num_clusters <- sum(grepl("^>Cluster", readLines(ptuh_miranda)))
print(num_clusters)
```

    ## [1] 8341

``` r
cat("\nNum. of miRNA-lncRNA bindings >75% (miRanda):\n")
```

    ## 
    ## Num. of miRNA-lncRNA bindings >75% (miRanda):

``` r
miranda_df <- read.table(ptuh_miranda, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
miranda_df$V8 <- as.numeric(str_remove(miranda_df$V8, "%"))
num_75 <- sum(miranda_df$V8 >= 75, na.rm = TRUE)
print(num_75)
```

    ## [1] 4246

``` r
# ---- Significant PCC correlations ----
pcc_df <- read_csv(ptuh_pcc_loc, show_col_types = FALSE)
```

    ## New names:
    ## • `` -> `...1`

``` r
count_sig <- sum(pcc_df[[5]] < 0.05, na.rm = TRUE)
count_sig_pos <- sum(pcc_df[[5]] < 0.05 & pcc_df[[4]] > 0, na.rm = TRUE)
percent <- if (count_sig > 0) round((count_sig_pos / count_sig) * 100, 2) else 0

cat("\nSignificant PCC correlations + % positive\n")
```

    ## 
    ## Significant PCC correlations + % positive

``` r
cat("Significant PCC (p<0.05):       ", count_sig, "\n")
```

    ## Significant PCC (p<0.05):        41931

``` r
cat("Significant (p<0.05), pos. PCC: ", count_sig_pos, "\n")
```

    ## Significant (p<0.05), pos. PCC:  28344

``` r
cat("Percent positive:               ", percent, "%\n")
```

    ## Percent positive:                67.6 %

``` r
# ---- Overlap (binding AND sig. PCC) ----
miranda_pcc <- read_csv(ptuh_miranda_pcc, show_col_types = FALSE)
```

    ## New names:
    ## • `` -> `...1`

``` r
# Clean columns 9, 12, 13 (remove quotes/percent)
miranda_pcc <- miranda_pcc %>%
  mutate(
    across(c(9, 12, 13), ~ as.numeric(str_remove_all(as.character(.x), '["%]')))
  )

count_sig2 <- sum(miranda_pcc[[13]] < 0.05, na.rm = TRUE)
count_sig_pos2 <- sum(miranda_pcc[[13]] < 0.05 & miranda_pcc[[12]] > 0, na.rm = TRUE)
percent2 <- if (count_sig2 > 0) round((count_sig_pos2 / count_sig2) * 100, 2) else 0

cat("\nOverlap (binding AND sig. PCC) + % positive\n")
```

    ## 
    ## Overlap (binding AND sig. PCC) + % positive

``` r
cat("Significant PCC (p<0.05):       ", count_sig2, "\n")
```

    ## Significant PCC (p<0.05):        564

``` r
cat("Significant (p<0.05), pos. PCC: ", count_sig_pos2, "\n")
```

    ## Significant (p<0.05), pos. PCC:  311

``` r
cat("Percent positive:               ", percent2, "%\n")
```

    ## Percent positive:                55.14 %

``` r
# ---- Overlap >75% (binding >75% AND PCC) ----
count_sig3 <- sum(miranda_pcc[[13]] < 0.05 & miranda_pcc[[9]] > 75, na.rm = TRUE)
count_sig_pos3 <- sum(miranda_pcc[[13]] < 0.05 & miranda_pcc[[12]] > 0 & miranda_pcc[[9]] > 75, na.rm = TRUE)
percent3 <- if (count_sig3 > 0) round((count_sig_pos3 / count_sig3) * 100, 2) else 0

cat("\nOverlap >75% (binding >75% AND PCC) + % positive\n")
```

    ## 
    ## Overlap >75% (binding >75% AND PCC) + % positive

``` r
cat("Significant PCC (p<0.05):       ", count_sig3, "\n")
```

    ## Significant PCC (p<0.05):        216

``` r
cat("Significant (p<0.05), pos. PCC: ", count_sig_pos3, "\n")
```

    ## Significant (p<0.05), pos. PCC:  130

``` r
cat("Percent positive:               ", percent3, "%\n")
```

    ## Percent positive:                60.19 %
