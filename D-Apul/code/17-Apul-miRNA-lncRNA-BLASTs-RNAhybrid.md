17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid
================
Kathleen Durkin
2024-12-11

- <a href="#1-prep-for-blasts" id="toc-1-prep-for-blasts">1 Prep for
  BLASTs</a>
  - <a href="#11-isolate-the-pre-mirna-and-mature-mirna-sequences"
    id="toc-11-isolate-the-pre-mirna-and-mature-mirna-sequences">1.1 Isolate
    the pre-mirna and mature mirna sequences</a>
  - <a href="#12-check-mirna-lengths" id="toc-12-check-mirna-lengths">1.2
    Check miRNA lengths</a>
  - <a href="#13-check-lncrnas" id="toc-13-check-lncrnas">1.3 check
    lncRNAs</a>
- <a href="#2-blasts" id="toc-2-blasts">2 BLASTs</a>
  - <a href="#21-make-databases" id="toc-21-make-databases">2.1 Make
    databases</a>
  - <a href="#22-run-blastn" id="toc-22-run-blastn">2.2 Run BLASTn</a>
- <a href="#3-examine-blast-tables" id="toc-3-examine-blast-tables">3
  Examine BLAST tables</a>
  - <a href="#31-lncrnas-as-mirna-precursors"
    id="toc-31-lncrnas-as-mirna-precursors">3.1 LncRNAs as miRNA
    precursors</a>
  - <a href="#32-lncrnas-as-mirna-sponges"
    id="toc-32-lncrnas-as-mirna-sponges">3.2 LncRNAs as miRNA sponges</a>
- <a href="#4-miranda" id="toc-4-miranda">4 miRanda</a>
  - <a href="#41-run-miranda" id="toc-41-run-miranda">4.1 Run miRanda</a>
- <a href="#5-summarize-results" id="toc-5-summarize-results">5 Summarize
  results</a>
- <a href="#6-rnahybrid" id="toc-6-rnahybrid">6 RNAhybrid</a>
  - <a href="#61-get-lncrna-gtf" id="toc-61-get-lncrna-gtf">6.1 Get lncRNA
    gtf</a>
  - <a href="#62-break-up-100bp-sequences"
    id="toc-62-break-up-100bp-sequences">6.2 Break up &gt;100bp
    sequences</a>
  - <a href="#63-get-fasta-of-broken-up-lncrna-gtf"
    id="toc-63-get-fasta-of-broken-up-lncrna-gtf">6.3 Get fasta of broken-up
    lncRNA gtf</a>
  - <a href="#64-run-rnahybrid" id="toc-64-run-rnahybrid">6.4 Run
    RNAhybrid</a>
  - <a href="#65-summarize-rnahybrid-results"
    id="toc-65-summarize-rnahybrid-results">6.5 Summarize RNAhybrid
    results</a>
- <a href="#7-summarize-final-results"
  id="toc-7-summarize-final-results">7 Summarize final results</a>
- <a href="#8-references" id="toc-8-references">8 References</a>

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
knitr::opts_chunk$set(
  echo = TRUE,         # Display code chunks
  eval = FALSE         # Evaluate code chunks
)
```

Two possible interactions between miRNA and lncRNA are:

1)  lncRNA acting as a precursor molecule for miRNA(s), so that the
    lncRNA contains one or many pre-miRNA sequences and will be broken
    down into pre-miRNAs molecules, which will then be processed into
    mature miRNAs.

2)  lncRNA acting as a “sponge” for miRNAs, so that an miRNA will bind
    to the lncRNA instead of being incorporated into an RISC complex to
    alter gene expression.

In situation 1 we would expect one or several **pre-miRNA sequences to
appear inside of a lncRNA**. This should be identifiable via BLASTn.

In situation 2 we would expect the **mature miRNA sequence to appear
inside a lncRNA**. Note that situation 2 is a bit more complicated,
because we can’t say for certain what sequence similarity is required
for binding. In cnidarians, miRNAs seem to act, like plants, through
complementarity of the full mature miRNA (this is in contrast to
e.g. mammals, where only binding of a short seed region is required)
(Moran et al. ([2014](#ref-moran_cnidarian_2014)), Admoni et al.
([2023](#ref-admoni_target_2023))). However, for lncRNA acting as
sponges, I don’t know whether to expect complementarity of the full
mature miRNA or only a section, and I don’t know what degree of
complementarity is required. **Work to identify lncRNA sponges could use
BLASTn, but will likely need to include additional methods like miranda
or RNAhybrid to identify potential binding.**

# 1 Prep for BLASTs

## 1.1 Isolate the pre-mirna and mature mirna sequences

``` bash
full_mirna_fasta="../output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out/mir.fasta"
premirna_fasta="../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_precursor.fasta"
mature_mirna_fasta="../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_mature.fasta"
star_mirna_fasta="../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_star.fasta"

# Pull out all sequences that DON'T contain "mature" or "star" in sequence name
# Note the pre-miRNAs have sequences for both strands
awk '
    # If the line starts with ">", check the header
    /^>/ {
        if ($0 ~ /mature/ || $0 ~ /star/) {
            print_seq = 0  # Skip sequences with "mature" or "star" in the header
        } else {
            print_seq = 1  # Mark sequences for printing
        }
    }
    # Print the header and the next two lines if marked for printing
    print_seq {
        print
        if (!/^>/) { getline; print }  # Capture second sequence line
    }
' "$full_mirna_fasta" > "$premirna_fasta"

# Pull out all sequences that contain "mature" in sequence name
grep -A 1 "mature" $full_mirna_fasta | grep -v "^--$" > $mature_mirna_fasta

# Pull out all sequences that contain "star" in sequence name
grep -A 1 "star" $full_mirna_fasta | grep -v "^--$" > $star_mirna_fasta
```

``` bash
premirna_fasta="../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_precursor.fasta"
mature_mirna_fasta="../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_mature.fasta"
star_mirna_fasta="../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_star.fasta"

# Check we have appropriate headers, same number of sequences in each
grep "^>" $premirna_fasta | head -2
echo ""
grep "^>" $mature_mirna_fasta | head -2
echo ""
grep "^>" $star_mirna_fasta | head -2
echo ""
grep "^>" $premirna_fasta | wc -l
echo ""
grep "^>" $mature_mirna_fasta | wc -l
echo ""
grep "^>" $star_mirna_fasta | wc -l
echo ""
```

    ## >Cluster_1826::ntLink_6:4847443-4847535(-)
    ## >Cluster_1832::ntLink_6:5157537-5157626(+)
    ## 
    ## >Cluster_1826.mature::ntLink_6:4847465-4847486(-)
    ## >Cluster_1832.mature::ntLink_6:5157559-5157579(+)
    ## 
    ## >Cluster_1826.star::ntLink_6:4847494-4847515(-)
    ## >Cluster_1832.star::ntLink_6:5157586-5157606(+)
    ## 
    ## 39
    ## 
    ## 39
    ## 
    ## 39

## 1.2 Check miRNA lengths

``` bash
# Extract sequence lengths for precursors
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_precursor.fasta > ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_precursor_lengths.txt

# Sequence lengths for matures
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_mature.fasta > ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_mature_lengths.txt
```

``` r
# Summary stats of precursor and mature lengths

precursor_lengths <- read.table("../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_precursor_lengths.txt", sep = " ", header = FALSE, col.names = c("seqID", "length"))
mature_lengths <- read.table("../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_mature_lengths.txt", sep = " ", header = FALSE, col.names = c("seqID", "length"))

cat("Average pre-miRNA length: ", mean(precursor_lengths$length))
```

    ## Average pre-miRNA length:  94.20513

``` r
cat("\n")
```

``` r
cat("Range of pre-miRNA lengths: ", range(precursor_lengths$length))
```

    ## Range of pre-miRNA lengths:  90 98

``` r
cat("\n")
```

``` r
cat("Average mature miRNA length: ", mean(mature_lengths$length))
```

    ## Average mature miRNA length:  22.25641

``` r
cat("\n")
```

``` r
cat("Range of mature miRNA lengths: ", range(mature_lengths$length))
```

    ## Range of mature miRNA lengths:  21 24

## 1.3 check lncRNAs

LncRNAs were identified from Apul RNA-seq data in
`deep-dive-expression/D-Apul/code/31-Apul-lncRNA` – see details there.
Fasta of Apul lncRNAs stored at
`deep-dive-expression/D-Apul/output/31-Apul-lncRNA/Apul_lncRNA.fasta`

``` bash
curl -L https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/D-Apul/output/31-Apul-lncRNA/Apul-lncRNA.fasta -o ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.fasta

echo "Number of lncRNAs:"
grep "^>" ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.fasta | wc -l
```

    ##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
    ##                                  Dload  Upload   Total   Spent    Left  Speed
    ##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  0     0    0     0    0     0      0      0 --:--:--  0:00:01 --:--:--     0  1 73.4M    1 1275k    0     0   598k      0  0:02:05  0:00:02  0:02:03  598k100 73.4M  100 73.4M    0     0  26.3M      0  0:00:02  0:00:02 --:--:-- 26.3M
    ## Number of lncRNAs:
    ## 31491

31491 total lncRNA in Apul.

``` bash
# Extract sequence lengths for precursors
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.fasta > ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_lengths.txt
```

``` r
# Summary stats of lncRNA lengths

lncRNA_lengths <- read.table("../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_lengths.txt", sep = " ", header = FALSE, col.names = c("seqID", "length"))

cat("Average lncRNA length: ", mean(lncRNA_lengths$length))
```

    ## Average lncRNA length:  2397.015

``` r
cat("\n")
```

``` r
cat("Range of lncRNA lengths: ", range(lncRNA_lengths$length))
```

    ## Range of lncRNA lengths:  201 146972

``` r
ggplot(lncRNA_lengths, aes(x = length)) +
  geom_histogram(binwidth = 500) +
  labs(title = "A. pulchra lncRNA sequence lengths",
       x = "Sequence Length [nucleotides]",
       y = "Frequency") +
  xlim(200, 150000) +
  ylim(0, 1000) +
  theme_minimal()
```

    ## Warning: Removed 7 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# 2 BLASTs

## 2.1 Make databases

Database of pre-miRNAs:

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_precursor.fasta \
-dbtype nucl \
-out ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/Apul-db/Apul_ShortStack_4.1.0_precursor
```

Database of mature miRNAs:

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_mature.fasta \
-dbtype nucl \
-out ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/Apul-db/Apul_ShortStack_4.1.0_mature
```

## 2.2 Run BLASTn

Generate a list of blast results. It seems plausible that a single
lncRNA, which would be hundreds or thousands of nucleotides long, could
interact with multiple miRNAs, so I will allow up to 10 hits (\~25% of
Apul miRNAs) for each lncRNA. I want to see the top hits no matter how
poor the match is, so I will not filter by e-value at this stage. I’ll
also include the “-word_size 4” option, which reduces the required
length of the initial match.

Full pre-miRNAs:

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.fasta \
-db ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/Apul-db/Apul_ShortStack_4.1.0_precursor \
-out ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/lncRNA_to_precursor_blastn.tab \
-num_threads 40 \
-word_size 4 \
-max_target_seqs 10 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/lncRNA_to_precursor_blastn.tab
```

    ## 295472 ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/lncRNA_to_precursor_blastn.tab

Note we have less than (10 \* $$# of lncRNAs$$) output alignments
because, while I did not set an evalue threshold, the default evalue
threshold of evalue=10 is still in place. That means extremely poor
matches were still excluded by default.

Mature miRNAs:

Note that I’m using the blastn-short option here because all of our
mature miRNAs are less than 30 nucleotides long (recommended by [BLAST
user
manual](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/))

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.fasta \
-db ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/Apul-db/Apul_ShortStack_4.1.0_mature \
-out ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/lncRNA_to_mature_blastn.tab \
-num_threads 40 \
-word_size 4 \
-max_target_seqs 10 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/lncRNA_to_mature_blastn.tab
```

    ## 301978 ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/lncRNA_to_mature_blastn.tab

# 3 Examine BLAST tables

Read into R and assign informative column labels

``` r
precursor_lncRNA_BLASTn <- read.table("../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/lncRNA_to_precursor_blastn.tab", sep="\t", header=FALSE)
mature_lncRNA_BLASTn <- read.table("../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/lncRNA_to_mature_blastn.tab", sep="\t", header=FALSE)

colnames(precursor_lncRNA_BLASTn) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
colnames(mature_lncRNA_BLASTn) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
```

## 3.1 LncRNAs as miRNA precursors

Are there any alignments of the full precursor miRNA to a lncRNA? Our
precursor sequences are 90-98 nucleotides long, so let’s look for any
alignments of at least 90 nucleotides with 0 mismatches.

``` r
precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  unique() %>%
  nrow()
```

    ## [1] 7

``` r
precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  select(qseqid) %>%
  unique() %>%
  nrow()
```

    ## [1] 7

``` r
precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  select(sseqid) %>%
  unique() %>%
  nrow()
```

    ## [1] 2

``` r
precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  unique() 
```

    ##              qseqid                                         sseqid pident
    ## 1 Apul_lncRNA_20880 Cluster_14768::ptg000023l:37965243-37965338(+)    100
    ## 2 Apul_lncRNA_20881 Cluster_14768::ptg000023l:37965243-37965338(+)    100
    ## 3 Apul_lncRNA_28983   Cluster_18728::ptg000035l:5346032-5346127(+)    100
    ## 4 Apul_lncRNA_28986   Cluster_18728::ptg000035l:5346032-5346127(+)    100
    ## 5 Apul_lncRNA_28987   Cluster_18728::ptg000035l:5346032-5346127(+)    100
    ## 6 Apul_lncRNA_28989   Cluster_18728::ptg000035l:5346032-5346127(+)    100
    ## 7 Apul_lncRNA_28990   Cluster_18728::ptg000035l:5346032-5346127(+)    100
    ##   length mismatch gapopen qstart  qend sstart send   evalue bitscore
    ## 1     96        0       0   1065  1160      1   96 4.46e-46      174
    ## 2     96        0       0   1070  1165      1   96 5.23e-46      174
    ## 3     96        0       0  10871 10966      1   96 1.17e-45      174
    ## 4     96        0       0  10398 10493      1   96 1.12e-45      174
    ## 5     96        0       0   9862  9957      1   96 1.07e-45      174
    ## 6     96        0       0   9807  9902      1   96 1.07e-45      174
    ## 7     96        0       0   9801  9896      1   96 1.07e-45      174

We have 7 alignments of a full pre-miRNA to a lncRNA with no mismatches.
7 lncRNA and 2 miRNA are represented.

Two things are notable. First, the lncRNA which contain full pre-miRNA
sequences also sit in the same genomic region as the pre-miRNA. This
isn’t particularly surprising, but it’s important to note that these
lncRNA precursors are not *in addition to* the pre-miRNA regions.

Second, there are instances where multiple lncRNA not only contain the
same pre-miRNA, but are also overlapping (in terms of genomic position).
This suggests to me that there are instances of lncRNA isoforms –
multiple variants of the same molecule. For the purposes of summary,
these will likely be treated as distinct lncRNA, since we don’t have a
good idea of what degree of difference is necessary to alter lncRNA
function.

Save these results

``` r
precursor_lncRNAs <- precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0)

write.table(precursor_lncRNAs, "../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/lncRNAs_as_miRNA_precursors.txt")
```

## 3.2 LncRNAs as miRNA sponges

I’m not sure whether to expect lncRNAs to bind miRNAs in the same way
cnidarian miRNA-mRNA binding occurs (nearly perfect complementarity of
mature sequence), or whether the mechanism could differ (e.g., requires
only a complementary seed region, as in vertebrate miRNA-mRNA binding).
that means I don’t know what alignment parameters to require for our
BLAST results.

For now let’s say the aligned region must be at least 8 nucleotides (the
expected length of an miRNA seed region), and let’s require a low evalue
of 1e-3, to generally restrict results to those with high
complementarity.

``` r
mature_lncRNA_BLASTn %>%
  filter(length >= 8) %>%
  filter(evalue <= 0.001)
```

    ##               qseqid                                                sseqid
    ## 1    Apul_lncRNA_698   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 2   Apul_lncRNA_1389   Cluster_10228.mature::ptg000017l:7471168-7471190(+)
    ## 3   Apul_lncRNA_1438   Cluster_18723.mature::ptg000035l:4808391-4808412(+)
    ## 4   Apul_lncRNA_1446   Cluster_18723.mature::ptg000035l:4808391-4808412(+)
    ## 5   Apul_lncRNA_1903   Cluster_18723.mature::ptg000035l:4808391-4808412(+)
    ## 6   Apul_lncRNA_2081  Cluster_2859.mature::ptg000001l:20063094-20063116(+)
    ## 7   Apul_lncRNA_2357  Cluster_5012.mature::ptg000008l:10754789-10754809(-)
    ## 8   Apul_lncRNA_2391      Cluster_4220.mature::ptg000007l:915927-915948(-)
    ## 9   Apul_lncRNA_2550   Cluster_18723.mature::ptg000035l:4808391-4808412(+)
    ## 10  Apul_lncRNA_4126       Cluster_19193.mature::ptg000039l:35786-35807(-)
    ## 11  Apul_lncRNA_4243    Cluster_2463.mature::ptg000001l:5548893-5548914(-)
    ## 12  Apul_lncRNA_5168      Cluster_5900.mature::ptg000009l:616601-616622(-)
    ## 13  Apul_lncRNA_5175      Cluster_5900.mature::ptg000009l:616601-616622(-)
    ## 14  Apul_lncRNA_5176      Cluster_5900.mature::ptg000009l:616601-616622(-)
    ## 15  Apul_lncRNA_5244 Cluster_14768.mature::ptg000023l:37965298-37965318(+)
    ## 16  Apul_lncRNA_5245 Cluster_14768.mature::ptg000023l:37965298-37965318(+)
    ## 17  Apul_lncRNA_5248    Cluster_2463.mature::ptg000001l:5548893-5548914(-)
    ## 18  Apul_lncRNA_5394   Cluster_15775.mature::ptg000025l:7472581-7472603(-)
    ## 19  Apul_lncRNA_5759 Cluster_14768.mature::ptg000023l:37965298-37965318(+)
    ## 20  Apul_lncRNA_6233   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 21  Apul_lncRNA_6586   Cluster_17791.mature::ptg000031l:6751957-6751979(-)
    ## 22  Apul_lncRNA_7460   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 23  Apul_lncRNA_8205 Cluster_15854.mature::ptg000025l:10668923-10668945(-)
    ## 24  Apul_lncRNA_8963  Cluster_2859.mature::ptg000001l:20063094-20063116(+)
    ## 25  Apul_lncRNA_9453   Cluster_17791.mature::ptg000031l:6751957-6751979(-)
    ## 26  Apul_lncRNA_9748   Cluster_15775.mature::ptg000025l:7472581-7472603(-)
    ## 27 Apul_lncRNA_10923    Cluster_3437.mature::ptg000004l:1859911-1859933(-)
    ## 28 Apul_lncRNA_11310   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 29 Apul_lncRNA_11475  Cluster_2859.mature::ptg000001l:20063094-20063116(+)
    ## 30 Apul_lncRNA_13360 Cluster_14768.mature::ptg000023l:37965298-37965318(+)
    ## 31 Apul_lncRNA_13628   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 32 Apul_lncRNA_14018      Cluster_4220.mature::ptg000007l:915927-915948(-)
    ## 33 Apul_lncRNA_14942    Cluster_2463.mature::ptg000001l:5548893-5548914(-)
    ## 34 Apul_lncRNA_15590   Cluster_10057.mature::ptg000016l:8599884-8599905(-)
    ## 35 Apul_lncRNA_16146      Cluster_5900.mature::ptg000009l:616601-616622(-)
    ## 36 Apul_lncRNA_16148      Cluster_5900.mature::ptg000009l:616601-616622(-)
    ## 37 Apul_lncRNA_16392  Cluster_2859.mature::ptg000001l:20063094-20063116(+)
    ## 38 Apul_lncRNA_16428   Cluster_15775.mature::ptg000025l:7472581-7472603(-)
    ## 39 Apul_lncRNA_16907  Cluster_3366.mature::ptg000002l:14046285-14046308(+)
    ## 40 Apul_lncRNA_17028      Cluster_1826.mature::ntLink_6:4847465-4847486(-)
    ## 41 Apul_lncRNA_17179      Cluster_5900.mature::ptg000009l:616601-616622(-)
    ## 42 Apul_lncRNA_17188      Cluster_5900.mature::ptg000009l:616601-616622(-)
    ## 43 Apul_lncRNA_17221      Cluster_5900.mature::ptg000009l:616601-616622(-)
    ## 44 Apul_lncRNA_17400      Cluster_5900.mature::ptg000009l:616601-616622(-)
    ## 45 Apul_lncRNA_17409      Cluster_5900.mature::ptg000009l:616601-616622(-)
    ## 46 Apul_lncRNA_17410      Cluster_5900.mature::ptg000009l:616601-616622(-)
    ## 47 Apul_lncRNA_17670   Cluster_15316.mature::ptg000024l:4086254-4086275(+)
    ## 48 Apul_lncRNA_19553    Cluster_1951.mature::ntLink_6:13351801-13351822(-)
    ## 49 Apul_lncRNA_20416 Cluster_14768.mature::ptg000023l:37965298-37965318(+)
    ## 50 Apul_lncRNA_20880 Cluster_14768.mature::ptg000023l:37965298-37965318(+)
    ## 51 Apul_lncRNA_20881 Cluster_14768.mature::ptg000023l:37965298-37965318(+)
    ## 52 Apul_lncRNA_21034    Cluster_3250.mature::ptg000002l:7337560-7337581(+)
    ## 53 Apul_lncRNA_21035    Cluster_3250.mature::ptg000002l:7337560-7337581(+)
    ## 54 Apul_lncRNA_21894   Cluster_15775.mature::ptg000025l:7472581-7472603(-)
    ## 55 Apul_lncRNA_22649    Cluster_3437.mature::ptg000004l:1859911-1859933(-)
    ## 56 Apul_lncRNA_22763      Cluster_4220.mature::ptg000007l:915927-915948(-)
    ## 57 Apul_lncRNA_24440   Cluster_18711.mature::ptg000035l:4339600-4339621(-)
    ## 58 Apul_lncRNA_24738   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 59 Apul_lncRNA_24929   Cluster_10051.mature::ptg000016l:7795530-7795551(+)
    ## 60 Apul_lncRNA_24935    Cluster_4254.mature::ptg000007l:3377335-3377356(+)
    ## 61 Apul_lncRNA_25032   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 62 Apul_lncRNA_25152  Cluster_5012.mature::ptg000008l:10754789-10754809(-)
    ## 63 Apul_lncRNA_26488   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 64 Apul_lncRNA_26624   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 65 Apul_lncRNA_26634   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 66 Apul_lncRNA_26678   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 67 Apul_lncRNA_26734 Cluster_15854.mature::ptg000025l:10668923-10668945(-)
    ## 68 Apul_lncRNA_26805    Cluster_3437.mature::ptg000004l:1859911-1859933(-)
    ## 69 Apul_lncRNA_27067   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 70 Apul_lncRNA_27124   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 71 Apul_lncRNA_27153   Cluster_10057.mature::ptg000016l:8599884-8599905(-)
    ## 72 Apul_lncRNA_27160    Cluster_5981.mature::ptg000009l:4940537-4940559(-)
    ## 73 Apul_lncRNA_27754 Cluster_15854.mature::ptg000025l:10668923-10668945(-)
    ## 74 Apul_lncRNA_28983   Cluster_18728.mature::ptg000035l:5346054-5346075(+)
    ## 75 Apul_lncRNA_28986   Cluster_18728.mature::ptg000035l:5346054-5346075(+)
    ## 76 Apul_lncRNA_28987   Cluster_18728.mature::ptg000035l:5346054-5346075(+)
    ## 77 Apul_lncRNA_28989   Cluster_18728.mature::ptg000035l:5346054-5346075(+)
    ## 78 Apul_lncRNA_28990   Cluster_18728.mature::ptg000035l:5346054-5346075(+)
    ## 79 Apul_lncRNA_29074   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ## 80 Apul_lncRNA_29390  Cluster_2859.mature::ptg000001l:20063094-20063116(+)
    ## 81 Apul_lncRNA_29391  Cluster_2859.mature::ptg000001l:20063094-20063116(+)
    ## 82 Apul_lncRNA_29392  Cluster_2859.mature::ptg000001l:20063094-20063116(+)
    ## 83 Apul_lncRNA_29393  Cluster_2859.mature::ptg000001l:20063094-20063116(+)
    ## 84 Apul_lncRNA_29394  Cluster_2859.mature::ptg000001l:20063094-20063116(+)
    ## 85 Apul_lncRNA_29395  Cluster_2859.mature::ptg000001l:20063094-20063116(+)
    ## 86 Apul_lncRNA_29894  Cluster_3367.mature::ptg000002l:14046591-14046614(+)
    ## 87 Apul_lncRNA_29894  Cluster_3366.mature::ptg000002l:14046285-14046308(+)
    ## 88 Apul_lncRNA_30021  Cluster_3367.mature::ptg000002l:14046591-14046614(+)
    ## 89 Apul_lncRNA_30021  Cluster_3366.mature::ptg000002l:14046285-14046308(+)
    ## 90 Apul_lncRNA_30388  Cluster_2859.mature::ptg000001l:20063094-20063116(+)
    ## 91 Apul_lncRNA_30521  Cluster_3366.mature::ptg000002l:14046285-14046308(+)
    ## 92 Apul_lncRNA_30636 Cluster_14768.mature::ptg000023l:37965298-37965318(+)
    ## 93 Apul_lncRNA_30761   Cluster_10419.mature::ptg000018l:2286829-2286850(+)
    ##     pident length mismatch gapopen qstart  qend sstart send   evalue bitscore
    ## 1   94.444     18        1       0    272   289      1   18 1.00e-03     29.2
    ## 2   90.909     22        1       1    853   873     22    1 1.00e-03     28.3
    ## 3   95.238     21        1       0   9861  9881     22    2 2.00e-04     34.6
    ## 4   90.476     21        2       0   3310  3330     22    2 1.00e-03     30.1
    ## 5   95.238     21        1       0  12463 12483      2   22 3.31e-04     34.6
    ## 6  100.000     16        0       0   2902  2917      6   21 1.00e-03     30.1
    ## 7   94.444     18        1       0     78    95      2   19 3.99e-04     29.2
    ## 8   91.304     23        1       1     48    70     22    1 7.90e-05     30.1
    ## 9   90.909     22        1       1    240   261      2   22 1.00e-03     28.3
    ## 10  94.737     19        1       0   2028  2046     21    3 1.00e-03     31.0
    ## 11  90.000     20        2       0   1087  1106     22    3 1.00e-03     28.3
    ## 12  90.476     21        2       0   1103  1123     22    2 5.89e-04     30.1
    ## 13  90.476     21        2       0   1764  1784     22    2 5.62e-04     30.1
    ## 14  90.476     21        2       0    870   890     22    2 2.89e-04     30.1
    ## 15 100.000     19        0       0   3056  3074     19    1 2.66e-05     35.6
    ## 16 100.000     19        0       0   3056  3074     19    1 2.66e-05     35.6
    ## 17  94.737     19        1       0    287   305     21    3 1.29e-04     31.0
    ## 18  86.364     22        3       0    250   271     23    2 1.00e-03     27.4
    ## 19  94.737     19        1       0    103   121     19    1 7.98e-05     31.0
    ## 20  94.444     18        1       0    188   205     18    1 8.98e-04     29.2
    ## 21 100.000     15        0       0    156   170      3   17 5.70e-04     28.3
    ## 22  94.444     18        1       0     33    50     18    1 7.54e-04     29.2
    ## 23 100.000     15        0       0    585   599     16    2 8.99e-04     28.3
    ## 24 100.000     15        0       0    903   917     15    1 1.00e-03     28.3
    ## 25 100.000     14        0       0    308   321     14    1 1.00e-03     26.5
    ## 26 100.000     14        0       0     19    32     17    4 1.00e-03     26.5
    ## 27 100.000     14        0       0     85    98      1   14 1.00e-03     26.5
    ## 28 100.000     18        0       0   1051  1068      1   18 2.99e-05     33.7
    ## 29 100.000     15        0       0    235   249     15    1 3.56e-04     28.3
    ## 30 100.000     19        0       0    824   842      1   19 3.17e-05     35.6
    ## 31  94.444     18        1       0    217   234     18    1 3.69e-04     29.2
    ## 32  94.444     18        1       0    115   132      3   20 1.00e-03     29.2
    ## 33 100.000     16        0       0   2650  2665      7   22 7.76e-04     30.1
    ## 34  89.474     19        2       0    214   232      4   22 1.00e-03     26.5
    ## 35  90.476     21        2       0    347   367      2   22 1.00e-03     30.1
    ## 36  90.476     21        2       0    818   838      2   22 1.00e-03     30.1
    ## 37 100.000     15        0       0    131   145      1   15 5.43e-04     28.3
    ## 38  86.364     22        3       0     19    40     23    2 8.77e-04     27.4
    ## 39  91.304     23        1       1    248   270      2   23 4.38e-04     30.1
    ## 40  94.118     17        1       0    123   139      6   22 1.00e-03     27.4
    ## 41  90.476     21        2       0      7    27      2   22 1.32e-04     30.1
    ## 42  90.476     21        2       0    114   134     22    2 4.77e-04     30.1
    ## 43  90.476     21        2       0    373   393      2   22 1.74e-04     30.1
    ## 44  90.476     21        2       0    131   151      2   22 4.58e-04     30.1
    ## 45  90.476     21        2       0    636   656     22    2 2.30e-04     30.1
    ## 46  90.476     21        2       0    963   983     22    2 3.15e-04     30.1
    ## 47 100.000     14        0       0     56    69     16    3 1.00e-03     26.5
    ## 48  90.476     21        2       0     58    78      2   22 9.68e-05     30.1
    ## 49  94.737     19        1       0    612   630      1   19 2.55e-04     31.0
    ## 50 100.000     21        0       0   1120  1140      1   21 2.48e-06     39.2
    ## 51 100.000     21        0       0   1125  1145      1   21 2.59e-06     39.2
    ## 52  95.000     20        1       0   2913  2932      1   20 2.59e-04     32.8
    ## 53  95.000     20        1       0   2641  2660      1   20 2.38e-04     32.8
    ## 54  94.444     18        1       0     95   112     21    4 4.28e-04     29.2
    ## 55  94.118     17        1       0    155   171     17    1 1.00e-03     27.4
    ## 56 100.000     15        0       0    265   279     16    2 8.05e-04     28.3
    ## 57  94.737     19        1       0     32    50      4   22 1.00e-03     31.0
    ## 58  94.444     18        1       0    418   435      1   18 7.66e-04     29.2
    ## 59  94.444     18        1       0    258   275      4   21 6.51e-04     29.2
    ## 60  90.000     20        2       0    794   813      1   20 1.00e-03     28.3
    ## 61  94.118     17        1       0     83    99     17    1 1.00e-03     27.4
    ## 62  94.444     18        1       0   1082  1099      4   21 1.00e-03     29.2
    ## 63  94.444     18        1       0    793   810     18    1 1.00e-03     29.2
    ## 64  94.444     18        1       0    609   626      1   18 1.00e-03     29.2
    ## 65  94.444     18        1       0    609   626      1   18 1.00e-03     29.2
    ## 66  94.444     18        1       0    572   589      1   18 1.00e-03     29.2
    ## 67  90.000     20        2       0    610   629     21    2 1.00e-03     28.3
    ## 68  94.118     17        1       0    124   140      1   17 1.00e-03     27.4
    ## 69  94.444     18        1       0    192   209     18    1 4.20e-04     29.2
    ## 70  94.444     18        1       0    629   646      1   18 1.00e-03     29.2
    ## 71 100.000     15        0       0    566   580      3   17 1.00e-03     28.3
    ## 72 100.000     14        0       0    271   284     21    8 1.00e-03     26.5
    ## 73 100.000     15        0       0    822   836     21    7 1.00e-03     28.3
    ## 74 100.000     22        0       0  10893 10914      1   22 1.67e-06     41.0
    ## 75 100.000     22        0       0  10420 10441      1   22 1.61e-06     41.0
    ## 76 100.000     22        0       0   9884  9905      1   22 1.54e-06     41.0
    ## 77 100.000     22        0       0   9829  9850      1   22 1.53e-06     41.0
    ## 78 100.000     22        0       0   9823  9844      1   22 1.53e-06     41.0
    ## 79 100.000     18        0       0  16964 16981     18    1 3.39e-04     33.7
    ## 80 100.000     15        0       0    144   158     15    1 6.60e-04     28.3
    ## 81 100.000     15        0       0    144   158     15    1 6.76e-04     28.3
    ## 82 100.000     15        0       0    144   158     15    1 8.12e-04     28.3
    ## 83 100.000     15        0       0    585   599     15    1 1.00e-03     28.3
    ## 84 100.000     15        0       0    144   158     15    1 6.73e-04     28.3
    ## 85 100.000     15        0       0    144   158     15    1 6.74e-04     28.3
    ## 86 100.000     16        0       0   1111  1126     21    6 5.00e-04     30.1
    ## 87 100.000     16        0       0   1111  1126     21    6 5.00e-04     30.1
    ## 88 100.000     16        0       0    130   145     21    6 1.42e-04     30.1
    ## 89 100.000     16        0       0    130   145     21    6 1.42e-04     30.1
    ## 90 100.000     15        0       0   1028  1042     15    1 1.00e-03     28.3
    ## 91 100.000     16        0       0    366   381     24    9 2.92e-04     30.1
    ## 92  94.737     19        1       0    114   132      1   19 1.06e-04     31.0
    ## 93  94.444     18        1       0    357   374      1   18 1.00e-03     29.2

93 putative lncRNA sponges with these parameters.

Ultimately though these results are insufficient to determine lncRNA
sponging. We need to evaluate miRNA-lncRNA binding.

# 4 miRanda

miRanda is a target prediction software, used to identify likely
miRNA-mRNA interactions.

Inputs:

- FASTA of A.pulchra lncRNAs

- FASTA of A.pulchra mature miRNAs

## 4.1 Run miRanda

``` bash

# score cutoff >100
# energy cutoff <-20
# strict binding

/home/shared/miRanda-3.3a/src/miranda \
../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_ShortStack_4.1.0_mature.fasta \
../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.fasta \
-sc 100 \
-en -20 \
-strict \
-out ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict_all.tab
```

# 5 Summarize results

Let’s look at the output

``` bash

echo "miranda run finished!"
echo "Counting number of interacting miRNA-lncRNA pairs"

zgrep -c "Performing Scan" ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict_all.tab

echo "Parsing output"
grep -A 1 "Scores for this hit:" ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict_all.tab | sort | grep '>' > ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict-parsed.txt

echo "counting number of putative interactions predicted (can include multiple interactions between single miRNA-lncRNA pair)"
wc -l ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict_all.tab
```

    ## miranda run finished!
    ## Counting number of interacting miRNA-lncRNA pairs
    ## 1228149
    ## Parsing output
    ## counting number of putative interactions predicted (can include multiple interactions between single miRNA-lncRNA pair)
    ## 9983050 ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict_all.tab

This is a lot of putative interactions! We can probably narrow it down
though. In vertebrates, miRNA-mRNA binding only requires complementarity
of an miRNA seed region of \~8 nucleotides. This requirement is built in
to miRanda target prediction. In cnidarians, however, miRNA-mRNA binding
is believed to require near-complete complementarity of the full mature
miRNA, similarly to plants ( Admoni et al.
([2023](#ref-admoni_target_2023)) , Admoni et al.
([2025](#ref-admoni_mirna-target_2025)) ). While I couldn’t find any
information on expected requirements for miRNA-lncRNA sponges, its
possible the binding will function similarly to miRNA-mRNA binding.
Let’s look at how many putative interactions are predicted for a binding
length of at least 21 nucleotides (the length of our smallest mature
miRNA).

``` bash
echo "number of putative interactions of at least 21 nucleotides"
awk -F'\t' '$7 >= 21' ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict-parsed.txt | wc -l
echo ""
echo "check some:"
awk -F'\t' '$7 >= 21' ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict-parsed.txt | head -5
```

    ## number of putative interactions of at least 21 nucleotides
    ## 3117
    ## 
    ## check some:
    ## >Cluster_10051.mature::ptg000016l:7795530-7795551(+) Apul_lncRNA_10490   154.00  -21.36  2 21    13001 13024 21  61.90%  85.71%
    ## >Cluster_10051.mature::ptg000016l:7795530-7795551(+) Apul_lncRNA_10685   163.00  -21.86  2 21    190 213 21  71.43%  80.95%
    ## >Cluster_10051.mature::ptg000016l:7795530-7795551(+) Apul_lncRNA_12228   151.00  -20.84  2 21    363 388 23  65.22%  69.57%
    ## >Cluster_10051.mature::ptg000016l:7795530-7795551(+) Apul_lncRNA_15174   162.00  -23.23  2 21    2539 2562   21  71.43%  85.71%
    ## >Cluster_10051.mature::ptg000016l:7795530-7795551(+) Apul_lncRNA_16240   150.00  -21.01  2 21    23682 23709 25  64.00%  72.00%

The header for this output is formatted as:

mirna Target Score Energy-Kcal/Mol Query-Aln(start-end)
Subjetct-Al(Start-End) Al-Len Subject-Identity Query-Identity

We can see from the percent identities (last 2 entries) that this number
includes alignments with multiple mismatches. Let’s filter again to
reduce the number of permissible mismatches. Let’s say we want no more
than 3 mismatches (a gap is counted as a mismatch). For an alignment of
21 nucleotides, this would be an percent identity of (21-3)/21 = 85.7%.
The miRNA is our “subject”, so we will filter by column 8.

``` bash
echo "number of putative interactions of at least 21 nucleotides, with at most 3 mismatches"
awk -F'\t' '$7 >= 21' ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict-parsed.txt | awk -F'\t' '$8 >= 85' | wc -l
echo ""
echo "check some:"
awk -F'\t' '$7 >= 21' ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict-parsed.txt | awk -F'\t' '$8 >= 85' | head -5
```

    ## number of putative interactions of at least 21 nucleotides, with at most 3 mismatches
    ## 29
    ## 
    ## check some:
    ## >Cluster_15316.mature::ptg000024l:4086254-4086275(+) Apul_lncRNA_9485    183.00  -26.73  2 21    266 289 21  85.71%  90.48%
    ## >Cluster_15775.mature::ptg000025l:7472581-7472603(-) Apul_lncRNA_11550   174.00  -22.10  2 22    3281 3303   21  85.71%  85.71%
    ## >Cluster_15775.mature::ptg000025l:7472581-7472603(-) Apul_lncRNA_19619   174.00  -23.17  2 22    382 404 21  85.71%  85.71%
    ## >Cluster_15775.mature::ptg000025l:7472581-7472603(-) Apul_lncRNA_23595   178.00  -21.90  2 22    3859 3881   21  85.71%  90.48%
    ## >Cluster_15775.mature::ptg000025l:7472581-7472603(-) Apul_lncRNA_23596   178.00  -21.90  2 22    3816 3838   21  85.71%  90.48%

This is a dramatically smaller number – only 29 interactions are at
least 21 nucleotides with \<=3 mismatches

# 6 RNAhybrid

**NOTE: the below code for RNAhybrid has NOT been rerun using updated
lncRNA files, since we have decided not to use the tool in final
downstream analyses (in favor of miRanda)**

RNAhybrid is another miRNA-mRNA target prediction tool, which bases its
predictions primarily on thermodynamic binding stability (unlike
miRanda, which considers sequence features expected of miRNA targets).
While the tool is normally used to predict miRNA-mRNA binding, it should
also work for miRNA-lncRNA binding

First we need to format our lncRNA and mature miRNA data. RNAhybrid
requires a query fasta file of mature miRNAs, and a target fasta file
(in this case, of lncRNAs). The problem is that RNAhybrid can only
handle fastas that contain sequences of 1000 nucleotides or fewer. Some
of our lncRNAs are thousands of nucleotides long, so we’ll need to
reformat this file.

I need to:

1.  Get a gff/gtf/bed file of our lncRNAs

2.  Use a bash script to modify the gff so that any sequences of \>1000
    nucleotides are broken up into multiple sub-sequences (and
    appropriately annotated as such)

3.  Convert this modified gff back into a fasta file.

## 6.1 Get lncRNA gtf

We have a *candidate* lncRNA gtf that then underwent some filtering and
was converted to our final Apul_lncRNA.fasta. Let’s filter the gtf to
retain only the lncRNAs that made it into our final Apul_lncRNA.fasta.

``` bash
lncRNAfasta=../output/31-Apul-lncRNA/Apul_lncRNA.fasta
lncRNAcoordinates=../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_coordinates.txt
candidategtf=../output/31-Apul-lncRNA/Apul_lncRNA_candidates.gtf
lncRNAgtf=../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_unformatted.gtf

# Step 1: Extract coordinates from FASTA headers
grep "^>" $lncRNAfasta | \
sed 's/>//' | \
awk -F'[:\\-]' '{print $3, $4+1, $5}' OFS="\t" \
> $lncRNAcoordinates

# Step 2: Keep only the candidate gtf entries whose coordinates
# exactly match those included in the lncRNAfasta coordinates
awk 'NR==FNR {ref[$1,$2,$3]; next} ($1,$4,$5) in ref' \
$lncRNAcoordinates \
$candidategtf \
> $lncRNAgtf
```

``` bash
lncRNAfasta=../output/31-Apul-lncRNA/Apul_lncRNA.fasta
lncRNAgtf=../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_unformatted.gtf

# Check
echo "some lncRNA fasta sequences: "
grep "^>" $lncRNAfasta | head -100 |tail -3
echo ""
echo "same index of filtered lncRNA gtf sequences: "
head -100 $lncRNAgtf | tail -3
echo ""
echo "number of lncRNA fasta sequences: "
grep "^>" $lncRNAfasta | wc -l
echo "number of filtered lncRNA gtf sequences: "
wc -l $lncRNAgtf
```

    ## some lncRNA fasta sequences: 
    ## grep: ../output/31-Apul-lncRNA/Apul_lncRNA.fasta: No such file or directory
    ## 
    ## same index of filtered lncRNA gtf sequences: 
    ## ntLink_6 StringTie   transcript  12654169    12654858    .   +   .   transcript_id "MSTRG.1611.1"; gene_id "MSTRG.1611"; xloc "XLOC_000892"; class_code "u"; tss_id "TSS1327";
    ## ntLink_6 StringTie   transcript  13210103    13212000    .   +   .   transcript_id "MSTRG.1688.1"; gene_id "MSTRG.1688"; xloc "XLOC_000920"; class_code "u"; tss_id "TSS1373";
    ## ntLink_6 StringTie   transcript  13210123    13212000    .   +   .   transcript_id "MSTRG.1688.2"; gene_id "MSTRG.1688"; xloc "XLOC_000920"; class_code "u"; tss_id "TSS1373";
    ## 
    ## number of lncRNA fasta sequences: 
    ## grep: ../output/31-Apul-lncRNA/Apul_lncRNA.fasta: No such file or directory
    ## 0
    ## number of filtered lncRNA gtf sequences: 
    ## 24183 ../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_unformatted.gtf

Looks like we’re good!

Before we proceed I also just want to fix the gtf formatting. Right now
it looks like, instead of being contained as single column 9, all the
extra info (transcript ID, gene ID, etc.) is in separate tab-delimited
columns. Let’s get it all correctly formatted inside of the 9th column.

``` bash
unformatted=../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_unformatted.gtf
formatted=../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.gtf
awk -F'\t' '{
    combined = $9
    for (i = 10; i <= 18; i++) {
        combined = combined $i
    }
    gsub(/ /, "", combined)  # Remove spaces from the combined column
    $9 = combined
    for (i = 10; i <= 18; i++) {
        $i = ""
    }
    $0 = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9
    print $0
}' OFS='\t' $unformatted > $formatted

# Check
head -3 $formatted | awk -F'\t' '{print $9}'
```

## 6.2 Break up \>100bp sequences

``` bash

# mRNA-only genome gff
# Count total sequences in lncRNA gtf
wc -l ../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.gtf

# Count the number of sequences that contain >1000 bp
awk '{if ($5 - $4 > 1000) count++} END {print count}' ../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.gtf

# Check how the sequence names are formatted
head -2 ../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.gtf
```

    ## 24183 ../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.gtf
    ## 9601
    ## ntLink_0 StringTie   transcript  84516   93551   .   +   .   transcript_id"MSTRG.10.1";gene_id"MSTRG.10";xloc"XLOC_000009";class_code"u";tss_id"TSS16";
    ## ntLink_0 StringTie   transcript  15629   19151   .   -   .   transcript_id"MSTRG.3.1";gene_id"MSTRG.3";xloc"XLOC_000010";class_code"u";tss_id"TSS17";

about 40% of our lncRNAs are too long, so we’ll need to break them up

I want to break up any sequence \>1000bp into 1000bp chunks, adding a
line to the gff for each chunk.

(I also want there to be overlap among the chunks, in case the break
between two chunks falls in the middle of an miRNA binding site. Let’s
say a 25bp overlap, since that is just over the maximum expected miRNA
length.)

for now though let’s not worry about the overlap.

The below code checks every sequence in the gtf and, for sequences over
1000 nucleotides long, breaks them up iteratively into 1000bp chunks.
When it breaks up a sequence, it also appends to the final column of the
line a “parent ID” showing the original lncRNA ID.

``` bash

awk -v chunk_size=1000 '
BEGIN {OFS="\t"}
{
    seq_length = $5 - $4
    parent_id = $1 ":" $4 "-" $5
    if (seq_length > chunk_size) {
        start = $4
        ogend = $5
        while (start < ogend) {
            end = start + chunk_size
            if (end > ogend) end = ogend
            $4 = start
            $5 = end
            temp_col9 = $9 "parent_id\"" parent_id "\""  # Preserve the existing content and append parent_id
            print $1, $2, $3, $4, $5, $6, $7, $8, temp_col9
            start = end
        }
    } else {
        $9 = $9 "parent_id\"" parent_id "\""  # Append parent_id to the existing content in $9
        print
    }
}' "../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.gtf" > "../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_MAX1000.gtf"
```

``` bash
MAX1000gtf=../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_MAX1000.gtf
# mRNA-only genome gff
# Count total sequences in genome gff
wc -l $MAX1000gtf

# Count the number of sequences that contain >1000 bp
awk '{if ($5 - $4 > 1000) count++} END {print count}' $MAX1000gtf

# Check how the sequence names are formatted
head -5 $MAX1000gtf
```

    ## 65494 ../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_MAX1000.gtf
    ## 
    ## ntLink_0 StringTie   transcript  84516   85516   .   +   .   transcript_id"MSTRG.10.1";gene_id"MSTRG.10";xloc"XLOC_000009";class_code"u";tss_id"TSS16";parent_id"ntLink_0:84516-93551"
    ## ntLink_0 StringTie   transcript  85516   86516   .   +   .   transcript_id"MSTRG.10.1";gene_id"MSTRG.10";xloc"XLOC_000009";class_code"u";tss_id"TSS16";parent_id"ntLink_0:84516-93551"
    ## ntLink_0 StringTie   transcript  86516   87516   .   +   .   transcript_id"MSTRG.10.1";gene_id"MSTRG.10";xloc"XLOC_000009";class_code"u";tss_id"TSS16";parent_id"ntLink_0:84516-93551"
    ## ntLink_0 StringTie   transcript  87516   88516   .   +   .   transcript_id"MSTRG.10.1";gene_id"MSTRG.10";xloc"XLOC_000009";class_code"u";tss_id"TSS16";parent_id"ntLink_0:84516-93551"
    ## ntLink_0 StringTie   transcript  88516   89516   .   +   .   transcript_id"MSTRG.10.1";gene_id"MSTRG.10";xloc"XLOC_000009";class_code"u";tss_id"TSS16";parent_id"ntLink_0:84516-93551"

Looks good!

## 6.3 Get fasta of broken-up lncRNA gtf

``` bash

# Use lncRNA gtf and genome fasta to extract lncRNA fastas

/home/shared/bedtools2/bin/bedtools getfasta \
-fi "../data/Apulchra-genome.fa" \
-bed "../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_MAX1000.gtf" \
-fo "../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_MAX1000.fa"
```

## 6.4 Run RNAhybrid

Now we can run RNAhybrid! I was getting a weird issue with all-zero
pvalues when I used RNAcalibrate-generated shape distribution parameters
in `16-Apul-RNAhybrid`, so I’ll just use the built-in 3utr_worm
parameter again.

I have RNAhybrid installed on a miniconda environment

    # Check path to the conda environment I'm using
    which conda

    # Install RNAhybrid if neccessary
    conda install -y -c genomedk rnahybrid

    # Check installation
    conda list rnahybrid

``` bash
#Start time: 12/12/2024 15:36
#

RNAhybrid \
-s 3utr_worm \
-e -20 \
-p 0.05 \
-c \
-t ../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_MAX1000.fa \
-q ../data/16-Apul-RNAhybrid/miRNA_mature-Apul-ShortStack_4.1.0-pulchra_genome.fasta \
> ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-RNAhybrid-lncRNA-compact_3utrworm.txt
```

## 6.5 Summarize RNAhybrid results

``` bash

# How many significant hybridizations predicted for each input?
wc -l ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-RNAhybrid-lncRNA-compact_3utrworm.txt
echo ""
head -3 ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-RNAhybrid-lncRNA-compact_3utrworm.txt
```

    ## 1566 ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-RNAhybrid-lncRNA-compact_3utrworm.txt
    ## 
    ## ntLink_6:15566281-15566670:389:Cluster_1826.mature::ntLink_6:4847465-4847486(-):22:-26.2:0.020079:150:   U              C   :    GGAAAGUGCUAUGA    :    UCUUUCACGAUACU    :UGUU              AGUA
    ## ntLink_8:10040388-10040616:228:Cluster_1826.mature::ntLink_6:4847465-4847486(-):22:-27.2:0.004032:123:G               A      : ACAAAGAGAGUGCUA       : UGUUUCUUUCACGAU       :                ACUAGUA
    ## ptg000001l:25192-26193:1000:Cluster_1826.mature::ntLink_6:4847465-4847486(-):22:-29.5:0.014814:926: G                A   :  GAAGAAAGUGCUAUGA    :  UUUCUUUCACGAUACU    :UG                AGUA

Now let’s read our RNAhybrid results into R for visualization. Note this
is going to be slightly more complicated than it sounds because the
RNAhybrid compact output is colon-delimited and our target- and
query-IDs contain intentional colons than could get confused with column
delimiters.

``` r
RNAhybrid_lncRNA <- read.table("../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-RNAhybrid-lncRNA-compact_3utrworm.txt", sep=":")

# Recombine Columns 1 and 2 (fix incorrect separation of target ID components)
RNAhybrid_lncRNA$V1 <- paste(RNAhybrid_lncRNA$V1, RNAhybrid_lncRNA$V2, sep = ":")
RNAhybrid_lncRNA$V2 <- NULL

# Do the same for Columns 4-7 (query ID components)
RNAhybrid_lncRNA$V4 <- paste(RNAhybrid_lncRNA$V4, RNAhybrid_lncRNA$V5, RNAhybrid_lncRNA$V6, RNAhybrid_lncRNA$V7 , sep = ":")
RNAhybrid_lncRNA$V4 <- gsub(":NA:", "::", RNAhybrid_lncRNA$V4)
RNAhybrid_lncRNA$V5 <- NULL
RNAhybrid_lncRNA$V6 <- NULL
RNAhybrid_lncRNA$V7 <- NULL

# Rename all columns for readability/accessibility 
colnames(RNAhybrid_lncRNA) <- c("target_name", "target_length", "query_name", "query_length",
                              "mfe", "pval", "position",
                              "noncomp_target_seq", "comp_target_seq", "comp_query_seq", "noncomp_query_seq")
```

Right now the “target” names are, for many lncRNAs, the broken-up
“chunks” of 1000bp. Let’s associate all of these chunks back to their
original “parent” lncRNAs.

``` r
# Read in current gtf
MAX1000gtf <- read.table("../data/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_MAX1000.gtf", sep="\t", head=FALSE)

# make new column for the full sequence ID (scaffold:beginning location-end location)
MAX1000gtf$target_name <- paste(MAX1000gtf$V1, paste(MAX1000gtf$V4-1, MAX1000gtf$V5, sep = "-"), sep = ":")


# separate out all of the extra info contained in column 9

# Define the ID types
id_types <- c("transcript_id", "gene_id", "xloc", "class_code", "cmp_ref_gene", "tss_id", "parent_id")

# Function to extract values for all ID types
extract_ids <- function(row, id_types) {
  # Split the row by ';'
  entries <- strsplit(row, ";")[[1]]
  
  # Initialize a named list with NA for all fields
  result <- setNames(rep(NA, length(id_types)), id_types)
  
  # Populate the list with actual values from the row
  for (entry in entries) {
    for (id_type in id_types) {
      if (startsWith(entry, id_type)) {
        result[[id_type]] <- sub(paste0("^", id_type), "", entry)
      }
    }
  }
  
  # Return the result as a named vector
  return(result)
}

# Apply the function to each row in column V9
parsed_data <- t(sapply(MAX1000gtf$V9, extract_ids, id_types = id_types))

# Convert the result into a data frame
parsed_df <- as.data.frame(parsed_data, stringsAsFactors = FALSE)
rownames(parsed_df) <- NULL  # Reset row names

# Combine the parsed data back into the original data frame
MAX1000gtf <- cbind(MAX1000gtf, parsed_df)

# Keep only the columns we may want to use
MAX1000gtf_reduced <- MAX1000gtf %>%
  select(target_name, parent_id, gene_id, cmp_ref_gene)
# remove duplicate rows (I believe stemming from isoforms)
MAX1000gtf_reduced <- MAX1000gtf_reduced[!duplicated(MAX1000gtf_reduced$target_name), ]

# Save these for later use
write.table(MAX1000gtf, "../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_MAX1000_expandedgtf_large.txt", sep="\t")
write.table(MAX1000gtf_reduced, "../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_MAX1000_expandedgtf.txt", sep="\t")
```

Now we can merge this gtf-based association table with the RNAhybrid
outputs

``` r
# merge and keep only columns of interest
RNAhybrid_lncRNA_annot <- left_join(RNAhybrid_lncRNA, MAX1000gtf_reduced, by = c("target_name" = "target_name")) %>%
  select(target_name, query_name, mfe, pval, parent_id, gene_id, cmp_ref_gene)

# Also, grab just the miRNA cluster name, for simplicity
RNAhybrid_lncRNA_annot$miRNA_cluster <- sub("\\..*", "", RNAhybrid_lncRNA_annot$query_name)

# move lncRNA parent ID and miRNA cluster name to first two columns
RNAhybrid_lncRNA_annot <- RNAhybrid_lncRNA_annot %>% select(parent_id, miRNA_cluster, everything())

# take a look
head(RNAhybrid_lncRNA_annot)
```

    ##                      parent_id miRNA_cluster                  target_name
    ## 1   ntLink_6:15566282-15566670  Cluster_1826   ntLink_6:15566281-15566670
    ## 2   ntLink_8:10040389-10040616  Cluster_1826   ntLink_8:10040388-10040616
    ## 3       ptg000001l:25193-27088  Cluster_1826       ptg000001l:25192-26193
    ## 4 ptg000001l:10926731-10928207  Cluster_1826 ptg000001l:10927730-10928207
    ## 5   ptg000002l:5112686-5113594  Cluster_1826   ptg000002l:5112685-5113594
    ## 6   ptg000007l:7556635-7557971  Cluster_1826   ptg000007l:7557634-7557971
    ##                                         query_name   mfe     pval     gene_id
    ## 1 Cluster_1826.mature::ntLink_6:4847465-4847486(-) -26.2 0.020079  MSTRG.1961
    ## 2 Cluster_1826.mature::ntLink_6:4847465-4847486(-) -27.2 0.004032  MSTRG.4113
    ## 3 Cluster_1826.mature::ntLink_6:4847465-4847486(-) -29.5 0.014814  MSTRG.7321
    ## 4 Cluster_1826.mature::ntLink_6:4847465-4847486(-) -29.1 0.005259  MSTRG.8409
    ## 5 Cluster_1826.mature::ntLink_6:4847465-4847486(-) -27.5 0.037381 MSTRG.10224
    ## 6 Cluster_1826.mature::ntLink_6:4847465-4847486(-) -25.8 0.019897 MSTRG.14688
    ##   cmp_ref_gene
    ## 1         <NA>
    ## 2         <NA>
    ## 3         <NA>
    ## 4         <NA>
    ## 5         <NA>
    ## 6         <NA>

``` r
# save
write.table(RNAhybrid_lncRNA_annot, "../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-RNAhybrid-lncRNA-compact_3utrworm-annot.txt", sep = "\t")
```

# 7 Summarize final results

LncRNA as miRNA precursors:

``` r
cat("Number of putative lncRNA precursors: ", length(filter(filter(precursor_lncRNA_BLASTn, length >= 90), mismatch == 0)$qseqid), "\n",
    "Number of miRNA whose precursors are lncRNA: ", length(unique(filter(filter(precursor_lncRNA_BLASTn, length >= 90), mismatch == 0)$sseqid)))
```

    ## Number of putative lncRNA precursors:  7 
    ##  Number of miRNA whose precursors are lncRNA:  2

Note: So there are 3 instances of a *unique* lncRNA containing a full
pre-miRNA sequence (and one of those instances occurs in 5 lncRNA
isoforms)

LncRNA as miRNA sponges:

``` bash

echo "miRanda Results:"
echo ""
echo "Number of putative interactions:"
wc -l  ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict_all.tab
echo ""
echo "Number of putative interactions of at least 21 nucleotides"
awk -F'\t' '$7 >= 21' ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict-parsed.txt | wc -l
echo ""
echo "Number of putative interactions of at least 21 nucleotides, with at most 3 mismatches"
awk -F'\t' '$7 >= 21' ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict-parsed.txt | awk -F'\t' '$8 >= 85' | wc -l
```

    ## miRanda Results:
    ## 
    ## Number of putative interactions:
    ## 9983050 ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict_all.tab
    ## 
    ## Number of putative interactions of at least 21 nucleotides
    ## 3117
    ## 
    ## Number of putative interactions of at least 21 nucleotides, with at most 3 mismatches
    ## 29

**NOTE: below RNAhybrid results used lncRNA files that are now
outdated**

``` r
cat("RNAhybrid Results: ", "\n", "\n",
    "(as a reminder)", "\n",
    "Number of A.pulchra lncRNAs: ", length(lncRNA_lengths$seqID), "\n",
    "Number of A.pulchra miRNAs: ", length(mature_lengths$seqID), "\n",
    "~~~~~~~~~~~~~~~~~~~~~~", "\n",
    "for p < 0.05: ", "\n",
    "Number of significant lncRNA-miRNA hybridizations: ", length(RNAhybrid_lncRNA_annot$parent_id), "\n",
    "Number of putative lncRNA sponges: ", length(unique(RNAhybrid_lncRNA_annot$parent_id)), "\n",
    "Number of miRNA putatively sequestered by lncRNA: ", length(unique(RNAhybrid_lncRNA_annot$miRNA_cluster)), "\n",
    "~~~~~~~~~~~~~~~~~~~~~~", "\n",
    "for p < 0.01: ", "\n",
    "Number of lncRNA-miRNA hybridizations: ", length(filter(RNAhybrid_lncRNA_annot, pval < 0.01)$parent_id), "\n",
    "Number of putative lncRNA sponges: ", length(unique(filter(RNAhybrid_lncRNA_annot, pval < 0.01)$parent_id)), "\n",
    "Number of miRNA putatively sequestered by lncRNA sponges: ", length(unique(filter(RNAhybrid_lncRNA_annot, pval < 0.01)$miRNA_cluster)))
```

    ## RNAhybrid Results:  
    ##  
    ##  (as a reminder) 
    ##  Number of A.pulchra lncRNAs:  31491 
    ##  Number of A.pulchra miRNAs:  39 
    ##  ~~~~~~~~~~~~~~~~~~~~~~ 
    ##  for p < 0.05:  
    ##  Number of significant lncRNA-miRNA hybridizations:  1566 
    ##  Number of putative lncRNA sponges:  1265 
    ##  Number of miRNA putatively sequestered by lncRNA:  39 
    ##  ~~~~~~~~~~~~~~~~~~~~~~ 
    ##  for p < 0.01:  
    ##  Number of lncRNA-miRNA hybridizations:  313 
    ##  Number of putative lncRNA sponges:  259 
    ##  Number of miRNA putatively sequestered by lncRNA sponges:  35

# 8 References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-admoni_target_2023" class="csl-entry">

Admoni, Yael, Arie Fridrich, Talya Razin, Miguel Salinas-Saavedra,
Michal Rabani, Uri Frank, and Yehu Moran. 2023. “Target Complementarity
in Cnidarians Supports a Common Origin for Animal and Plant <span
class="nocase">microRNAs</span>.” bioRxiv.
<https://doi.org/10.1101/2023.01.08.523153>.

</div>

<div id="ref-admoni_mirna-target_2025" class="csl-entry">

Admoni, Yael, Arie Fridrich, Paris K Weavers, Reuven Aharoni, Talya
Razin, Miguel Salinas-Saavedra, Michal Rabani, Uri Frank, and Yehu
Moran. 2025. “<span class="nocase">miRNA</span>-Target Complementarity
in Cnidarians Resembles Its Counterpart in Plants.” *EMBO Reports*,
January, 1–24. <https://doi.org/10.1038/s44319-024-00350-z>.

</div>

<div id="ref-moran_cnidarian_2014" class="csl-entry">

Moran, Yehu, David Fredman, Daniela Praher, Xin Li, Liang Wee, Fabian
Rentzsch, Phillip Zamore, Ulrich Technau, and Hervé Seitz. 2014.
“Cnidarian <span class="nocase">microRNAs</span> Frequently Regulate
Targets by Cleavage.” *Genome Research* 24 (March).
<https://doi.org/10.1101/gr.162503.113>.

</div>

</div>
