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
`deep-dive-expression/D-Apul/code/10-Apul-lncRNA` – see details there.
Fasta of Apul lncRNAs stored at
`deep-dive-expression/D-Apul/output/10-Apul-lncRNA/Apul_lncRNA.fasta`

``` bash
curl -L https://gannet.fish.washington.edu/acropora/E5-deep-dive-expression/output/01.6-lncRNA-pipline/Apul_lncRNA.fasta -o ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.fasta

echo "Number of lncRNAs:"
grep "^>" ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.fasta | wc -l
```

    ##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
    ##                                  Dload  Upload   Total   Spent    Left  Speed
    ##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0 12  392M   12 47.8M    0     0  85.0M      0  0:00:04 --:--:--  0:00:04 85.0M 40  392M   40  157M    0     0   101M      0  0:00:03  0:00:01  0:00:02  101M 68  392M   68  268M    0     0   104M      0  0:00:03  0:00:02  0:00:01  104M 90  392M   90  353M    0     0  99.2M      0  0:00:03  0:00:03 --:--:-- 99.2M100  392M  100  392M    0     0  97.6M      0  0:00:04  0:00:04 --:--:-- 97.6M
    ## Number of lncRNAs:
    ## 77200

``` bash
# Extract sequence lengths for precursors
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA.fasta > ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_lengths.txt
```

``` r
# Summary stats of lncRNA lengths

lncRNA_lengths <- read.table("../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_lengths.txt", sep = " ", header = FALSE, col.names = c("seqID", "length"))

cat("Average lncRNA length: ", mean(lncRNA_lengths$length))
```

    ## Average lncRNA length:  5201.144

``` r
cat("\n")
```

``` r
cat("Range of lncRNA lengths: ", range(lncRNA_lengths$length))
```

    ## Range of lncRNA lengths:  201 181838

``` r
ggplot(lncRNA_lengths, aes(x = length)) +
  geom_histogram(binwidth = 500) +
  labs(title = "A. pulchra lncRNA sequence lengths",
       x = "Sequence Length [nucleotides]",
       y = "Frequency") +
  xlim(200, 200000) +
  ylim(0, 1000) +
  theme_minimal()
```

    ## Warning: Removed 18 rows containing missing values or values outside the scale range
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

    ## 738169 ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/lncRNA_to_precursor_blastn.tab

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

    ## 751007 ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/blasts/lncRNA_to_mature_blastn.tab

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

    ## [1] 34

``` r
precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  select(qseqid) %>%
  unique() %>%
  nrow()
```

    ## [1] 31

``` r
precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  select(sseqid) %>%
  unique() %>%
  nrow()
```

    ## [1] 23

``` r
precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  unique() 
```

    ##                                      qseqid
    ## 1      transcript::ntLink_6:5155137-5159371
    ## 3      transcript::ntLink_6:5156340-5159371
    ## 4      transcript::ntLink_6:4834719-4876081
    ## 5    transcript::ntLink_6:13329954-13369257
    ## 6    transcript::ptg000001l:5544306-5580019
    ## 7  transcript::ptg000002l:14041039-14052409
    ## 8  transcript::ptg000002l:14041039-14052409
    ## 9      transcript::ptg000007l:908509-921806
    ## 10     transcript::ptg000009l:611774-618362
    ## 11     transcript::ptg000009l:611774-618362
    ## 12     transcript::ptg000009l:613601-618362
    ## 13     transcript::ptg000009l:613601-618362
    ## 14   transcript::ptg000009l:4929027-4948436
    ## 15   transcript::ptg000016l:7787959-7823075
    ## 16   transcript::ptg000016l:8593882-8603986
    ## 17 transcript::ptg000016l:11748348-11760789
    ## 18 transcript::ptg000016l:11748351-11760789
    ## 19 transcript::ptg000016l:11748671-11760568
    ## 20 transcript::ptg000016l:11748715-11760789
    ## 21   transcript::ptg000018l:2280996-2290234
    ## 22 transcript::ptg000023l:18698522-18723762
    ## 23 transcript::ptg000023l:37964178-37968910
    ## 24 transcript::ptg000023l:38035263-38040884
    ## 25   transcript::ptg000024l:4085668-4087287
    ## 26   transcript::ptg000025l:1976289-1990446
    ## 27 transcript::ptg000025l:10491944-10502595
    ## 29 transcript::ptg000025l:10668087-10677231
    ## 30   transcript::ptg000035l:5335161-5347816
    ## 31   transcript::ptg000035l:5335634-5347816
    ## 32   transcript::ptg000035l:5336170-5347816
    ## 33   transcript::ptg000035l:5336225-5347816
    ## 34   transcript::ptg000035l:5336231-5347816
    ## 35   transcript::ptg000035l:4337555-4340587
    ## 36       transcript::ptg000039l:33382-41881
    ##                                            sseqid pident length mismatch
    ## 1       Cluster_1832::ntLink_6:5157537-5157626(+)    100     90        0
    ## 3       Cluster_1832::ntLink_6:5157537-5157626(+)    100     90        0
    ## 4       Cluster_1826::ntLink_6:4847443-4847535(-)    100     93        0
    ## 5     Cluster_1951::ntLink_6:13351746-13351842(-)    100     97        0
    ## 6     Cluster_2463::ptg000001l:5548841-5548934(-)    100     94        0
    ## 7   Cluster_3367::ptg000002l:14046541-14046634(+)    100     94        0
    ## 8   Cluster_3366::ptg000002l:14046235-14046328(+)    100     94        0
    ## 9       Cluster_4220::ptg000007l:915905-916002(-)    100     98        0
    ## 10      Cluster_5900::ptg000009l:616579-616673(-)    100     95        0
    ## 11      Cluster_5899::ptg000009l:616285-616377(-)    100     93        0
    ## 12      Cluster_5900::ptg000009l:616579-616673(-)    100     95        0
    ## 13      Cluster_5899::ptg000009l:616285-616377(-)    100     93        0
    ## 14    Cluster_5981::ptg000009l:4940515-4940609(-)    100     95        0
    ## 15   Cluster_10051::ptg000016l:7795508-7795602(+)    100     95        0
    ## 16   Cluster_10057::ptg000016l:8599833-8599925(-)    100     93        0
    ## 17 Cluster_10093::ptg000016l:11751358-11751448(-)    100     91        0
    ## 18 Cluster_10093::ptg000016l:11751358-11751448(-)    100     91        0
    ## 19 Cluster_10093::ptg000016l:11751358-11751448(-)    100     91        0
    ## 20 Cluster_10093::ptg000016l:11751358-11751448(-)    100     91        0
    ## 21   Cluster_10419::ptg000018l:2286780-2286870(+)    100     91        0
    ## 22 Cluster_14402::ptg000023l:18708903-18708996(+)    100     94        0
    ## 23 Cluster_14768::ptg000023l:37965243-37965338(+)    100     96        0
    ## 24 Cluster_14768::ptg000023l:37965243-37965338(+)    100     96        0
    ## 25   Cluster_15316::ptg000024l:4086202-4086295(+)    100     94        0
    ## 26   Cluster_15671::ptg000025l:1981692-1981783(+)    100     92        0
    ## 27 Cluster_15851::ptg000025l:10501030-10501125(+)    100     96        0
    ## 29 Cluster_15854::ptg000025l:10668901-10668998(-)    100     98        0
    ## 30   Cluster_18728::ptg000035l:5346032-5346127(+)    100     96        0
    ## 31   Cluster_18728::ptg000035l:5346032-5346127(+)    100     96        0
    ## 32   Cluster_18728::ptg000035l:5346032-5346127(+)    100     96        0
    ## 33   Cluster_18728::ptg000035l:5346032-5346127(+)    100     96        0
    ## 34   Cluster_18728::ptg000035l:5346032-5346127(+)    100     96        0
    ## 35   Cluster_18711::ptg000035l:4339550-4339641(-)    100     92        0
    ## 36       Cluster_19193::ptg000039l:35764-35856(-)    100     93        0
    ##    gapopen qstart  qend sstart send   evalue bitscore
    ## 1        0   2400  2489      1   90 7.22e-43      163
    ## 3        0   1197  1286      1   90 5.16e-43      163
    ## 4        0  12724 12816     93    1 1.60e-43      168
    ## 5        0  21792 21888     97    1 1.02e-45      176
    ## 6        0   4535  4628     94    1 3.96e-44      170
    ## 7        0   5502  5595      1   94 1.28e-44      170
    ## 8        0   5196  5289      1   94 1.28e-44      170
    ## 9        0   7396  7493     98    1 1.00e-46      178
    ## 10       0   4805  4899     95    1 2.14e-45      172
    ## 11       0   4511  4603     93    1 2.61e-44      168
    ## 12       0   2978  3072     95    1 1.57e-45      172
    ## 13       0   2684  2776     93    1 1.91e-44      168
    ## 14       0  11488 11582     95    1 6.24e-45      172
    ## 15       0   7549  7643      1   95 1.11e-44      172
    ## 16       0   5951  6043     93    1 4.01e-44      168
    ## 17       0   3010  3100     91    1 5.93e-43      165
    ## 18       0   3007  3097     91    1 5.93e-43      165
    ## 19       0   2687  2777     91    1 5.67e-43      165
    ## 20       0   2643  2733     91    1 5.76e-43      165
    ## 21       0   5784  5874      1   91 4.46e-43      165
    ## 22       0  10381 10474      1   94 2.80e-44      170
    ## 23       0   1065  1160      1   96 4.46e-46      174
    ## 24       0   1070  1165      1   96 5.23e-46      174
    ## 25       0    534   627      1   94 1.87e-45      170
    ## 26       0   5403  5494      1   92 1.93e-43      167
    ## 27       0   9086  9181      1   96 9.93e-46      174
    ## 29       0    814   911     98    1 7.00e-47      178
    ## 30       0  10871 10966      1   96 1.17e-45      174
    ## 31       0  10398 10493      1   96 1.12e-45      174
    ## 32       0   9862  9957      1   96 1.07e-45      174
    ## 33       0   9807  9902      1   96 1.07e-45      174
    ## 34       0   9801  9896      1   96 1.07e-45      174
    ## 35       0   1995  2086     92    1 4.24e-44      167
    ## 36       0   2382  2474     93    1 3.37e-44      168

We have 34 alignments of a full pre-miRNA to a lncRNA with no
mismatches. 31 lncRNA and 23 miRNA are represented.

Two things are notable. First, the lncRNA which contain full pre-miRNA
sequences also sit in the same genomic region as the pre-miRNA. This
isn’t particularly surprising, but it’s important to note that these
lncRNA precursors are not *in addition to* the pre-miRNA regions.

Second, there are several instances where multiple lncRNA not only
contain the same pre-miRNA, but are also overlapping (in terms of
genomic position). For example, the Cluster_10093 precursor is present
in the following lncRNA: transcript::ptg000016l:11748348-11760789,
transcript::ptg000016l:11748351-11760789,
transcript::ptg000016l:11748671-11760568, and
transcript::ptg000016l:11748715-11760789. These four lncRNA sit on the
same scaffold and share the same end position. This suggests to me that
there are instances of lncRNA isoforms – multiple variants of the same
molecule. For the purposes of summary, these will likely be treated as
distinct lncRNA, since we don’t have a god idea of what degree of
difference is necessary to alter lncRNA function.

It is also possible that a lncRNA may be directly processed into a
mature miRNA, without first being processed into pre-miRNA. Let’s look
for those by searching for alignments of mature miRNAs to lncRNAs. Our
mature miRNAs range 21-24 nucleotides in length. Let’s look for
alignments of at least 21 nucleotides in length with 0 mismatches and 0
gaps.

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

    ##                                       qseqid
    ## 1       transcript::ntLink_6:5155137-5159371
    ## 2       transcript::ntLink_6:5155137-5159371
    ## 3       transcript::ntLink_6:5156340-5159371
    ## 4       transcript::ntLink_6:4834719-4876081
    ## 5       transcript::ntLink_6:4863260-4876081
    ## 6     transcript::ntLink_6:13329954-13369257
    ## 7     transcript::ntLink_6:16464141-16469269
    ## 8       transcript::ntLink_6:6543820-6544705
    ## 9         transcript::ntLink_8:486960-488023
    ## 10      transcript::ntLink_8:3186326-3196541
    ## 11      transcript::ntLink_8:3305138-3308767
    ## 12    transcript::ntLink_8:22685890-22707953
    ## 13    transcript::ntLink_8:22685890-22707953
    ## 14        transcript::ntLink_8:668902-669217
    ## 15      transcript::ntLink_8:4437374-4477709
    ## 16      transcript::ntLink_8:4437376-4456957
    ## 17      transcript::ntLink_8:4439933-4477709
    ## 18    transcript::ntLink_8:10028596-10038034
    ## 19    transcript::ntLink_8:10028632-10044386
    ## 20    transcript::ntLink_8:10028669-10038034
    ## 21    transcript::ntLink_8:10029178-10037663
    ## 22    transcript::ntLink_8:19660513-19664606
    ## 23    transcript::ntLink_8:26538658-26550794
    ## 24    transcript::ntLink_8:26538658-26550794
    ## 25    transcript::ntLink_8:26539821-26550794
    ## 26      transcript::ntLink_8:2019087-2019406
    ## 27      transcript::ntLink_8:5011621-5011845
    ## 28      transcript::ntLink_8:5882049-5882621
    ## 29      transcript::ntLink_8:7264347-7265303
    ## 30      transcript::ntLink_8:7710817-7712117
    ## 31    transcript::ntLink_8:35371419-35371867
    ## 32  transcript::ptg000001l:12164041-12187052
    ## 33  transcript::ptg000001l:14329678-14339113
    ## 34  transcript::ptg000001l:20300959-20305253
    ## 35    transcript::ptg000001l:4724869-4729321
    ## 36    transcript::ptg000001l:5544306-5580019
    ## 37    transcript::ptg000001l:5998212-6000651
    ## 38  transcript::ptg000001l:10644454-10649597
    ## 39  transcript::ptg000001l:16246677-16264746
    ## 40  transcript::ptg000001l:17411559-17412717
    ## 41  transcript::ptg000002l:14041039-14052409
    ## 42  transcript::ptg000002l:14041039-14052409
    ## 43    transcript::ptg000002l:4720943-4723122
    ## 44    transcript::ptg000002l:5121476-5123331
    ## 45    transcript::ptg000002l:5122370-5123331
    ## 46    transcript::ptg000002l:8695162-8699329
    ## 47    transcript::ptg000002l:8695162-8699329
    ## 48    transcript::ptg000002l:8709640-8710034
    ## 49    transcript::ptg000002l:8946097-8947707
    ## 50    transcript::ptg000002l:9549210-9564066
    ## 51    transcript::ptg000002l:9549247-9564066
    ## 52    transcript::ptg000002l:9549354-9564066
    ## 53    transcript::ptg000002l:9549733-9562715
    ## 54  transcript::ptg000002l:12978526-12982461
    ## 55    transcript::ptg000002l:1383381-1383717
    ## 56  transcript::ptg000002l:12980969-12981195
    ## 57      transcript::ptg000004l:567157-571808
    ## 58  transcript::ptg000004l:13570684-13574992
    ## 59    transcript::ptg000004l:6702800-6706683
    ## 60    transcript::ptg000004l:1098398-1099254
    ## 61  transcript::ptg000004l:13769536-13770031
    ## 62      transcript::ptg000007l:908509-921806
    ## 63    transcript::ptg000007l:2791406-2791715
    ## 64    transcript::ptg000007l:6664571-6674626
    ## 65    transcript::ptg000007l:7178407-7179102
    ## 66    transcript::ptg000007l:7299212-7307849
    ## 67    transcript::ptg000007l:9586839-9587490
    ## 68  transcript::ptg000008l:22187313-22188105
    ## 69  transcript::ptg000008l:24309893-24315279
    ## 70        transcript::ptg000008l:43441-43948
    ## 71    transcript::ptg000008l:1574209-1575066
    ## 72  transcript::ptg000008l:22810033-22817561
    ## 73  transcript::ptg000008l:22810033-22817561
    ## 74  transcript::ptg000008l:23154675-23160297
    ## 75  transcript::ptg000008l:23156767-23159540
    ## 76  transcript::ptg000008l:24476657-24479596
    ## 77  transcript::ptg000008l:11694626-11695619
    ## 78  transcript::ptg000008l:25572824-25573194
    ## 79  transcript::ptg000008l:29764952-29765232
    ## 80    transcript::ptg000009l:3067182-3070131
    ## 81    transcript::ptg000009l:3067366-3069923
    ## 82    transcript::ptg000009l:3067618-3070131
    ## 83      transcript::ptg000009l:611774-618362
    ## 84      transcript::ptg000009l:611774-618362
    ## 85      transcript::ptg000009l:613601-618362
    ## 86      transcript::ptg000009l:613601-618362
    ## 87    transcript::ptg000009l:4929027-4948436
    ## 88    transcript::ptg000009l:9757122-9763155
    ## 89    transcript::ptg000009l:6485786-6486105
    ## 90      transcript::ptg000010l:496924-544737
    ## 91      transcript::ptg000010l:511224-544737
    ## 92      transcript::ptg000010l:511755-544737
    ## 93      transcript::ptg000010l:511755-544737
    ## 94      transcript::ptg000010l:543528-544737
    ## 95      transcript::ptg000010l:644795-661260
    ## 96      transcript::ptg000010l:646771-661704
    ## 97      transcript::ptg000010l:646832-660585
    ## 98    transcript::ptg000011l:6306487-6306773
    ## 99    transcript::ptg000011l:1427911-1431394
    ## 100   transcript::ptg000011l:4169626-4173486
    ## 101   transcript::ptg000011l:5552041-5556134
    ## 102 transcript::ptg000011l:11748489-11751325
    ## 103 transcript::ptg000011l:13267609-13268257
    ## 104   transcript::ptg000012l:8690502-8694128
    ## 105     transcript::ptg000012l:209105-221118
    ## 106     transcript::ptg000012l:209105-230917
    ## 107     transcript::ptg000012l:209113-224136
    ## 108     transcript::ptg000012l:210636-230917
    ## 109   transcript::ptg000012l:8696749-8701622
    ## 110   transcript::ptg000012l:8696749-8719892
    ## 111   transcript::ptg000012l:8699097-8711423
    ## 112 transcript::ptg000012l:12816099-12821132
    ## 113   transcript::ptg000015l:6682048-6682420
    ## 114 transcript::ptg000015l:10500419-10505388
    ## 115 transcript::ptg000015l:11083780-11084101
    ## 116   transcript::ptg000015l:2489793-2490942
    ## 117   transcript::ptg000015l:8949731-8952903
    ## 118   transcript::ptg000015l:8951319-8952903
    ## 119 transcript::ptg000015l:10716991-10719083
    ## 120   transcript::ptg000015l:1837134-1837430
    ## 121   transcript::ptg000015l:5493755-5494425
    ## 122   transcript::ptg000015l:6877278-6878452
    ## 123   transcript::ptg000016l:7787959-7823075
    ## 124   transcript::ptg000016l:8593882-8603986
    ## 125 transcript::ptg000016l:11748348-11760789
    ## 126 transcript::ptg000016l:11748351-11760789
    ## 127 transcript::ptg000016l:11748671-11760568
    ## 128 transcript::ptg000016l:11748715-11760789
    ## 129   transcript::ptg000016l:7924625-7927489
    ## 130   transcript::ptg000017l:6610157-6614214
    ## 131   transcript::ptg000017l:6663205-6666987
    ## 132   transcript::ptg000017l:6663205-6666987
    ## 133   transcript::ptg000017l:7791158-7794702
    ## 134   transcript::ptg000017l:7791466-7794702
    ## 135   transcript::ptg000017l:8069280-8069537
    ## 136   transcript::ptg000018l:2280996-2290234
    ## 137   transcript::ptg000018l:1178196-1178577
    ## 138   transcript::ptg000018l:3750411-3754055
    ## 139   transcript::ptg000018l:3767996-3771672
    ## 140   transcript::ptg000018l:4041685-4045953
    ## 141   transcript::ptg000018l:4269165-4272768
    ## 142   transcript::ptg000018l:4472290-4475949
    ## 143   transcript::ptg000018l:4536497-4540765
    ## 144   transcript::ptg000018l:4763977-4767580
    ## 145   transcript::ptg000018l:9541756-9542015
    ## 146   transcript::ptg000018l:5706103-5706575
    ## 147   transcript::ptg000018l:7640158-7640363
    ## 148     transcript::ptg000019l:454985-456434
    ## 149   transcript::ptg000019l:4600832-4609460
    ## 150   transcript::ptg000019l:2123324-2123600
    ## 151   transcript::ptg000020l:1314938-1315340
    ## 152   transcript::ptg000020l:1768626-1770203
    ## 153   transcript::ptg000020l:4217558-4218084
    ## 154     transcript::ptg000020l:615954-621880
    ## 155     transcript::ptg000020l:617421-621014
    ## 156   transcript::ptg000020l:1313648-1318092
    ## 157   transcript::ptg000020l:1314814-1316327
    ## 158   transcript::ptg000020l:1337771-1342704
    ## 159   transcript::ptg000020l:1423894-1427670
    ## 160   transcript::ptg000020l:1746901-1747592
    ## 161   transcript::ptg000020l:1746922-1747966
    ## 162   transcript::ptg000020l:4216792-4221114
    ## 163   transcript::ptg000020l:4248313-4252564
    ## 164   transcript::ptg000020l:5844736-5860466
    ## 165   transcript::ptg000020l:5845880-5860466
    ## 166   transcript::ptg000020l:5850706-5857299
    ## 167   transcript::ptg000020l:6512247-6517004
    ## 168   transcript::ptg000020l:4098333-4098592
    ## 169 transcript::ptg000020l:12568237-12568621
    ## 170 transcript::ptg000021l:13502731-13514774
    ## 171 transcript::ptg000021l:13503260-13514774
    ## 172 transcript::ptg000021l:13505533-13524949
    ## 173 transcript::ptg000021l:13505815-13514774
    ## 174 transcript::ptg000021l:13822331-13873009
    ## 175 transcript::ptg000021l:13842888-13873009
    ## 176 transcript::ptg000021l:13843080-13866689
    ## 177 transcript::ptg000021l:13843252-13877113
    ## 178 transcript::ptg000021l:16927850-16951190
    ## 179 transcript::ptg000021l:19195511-19200021
    ## 180 transcript::ptg000021l:19207344-19212121
    ## 181 transcript::ptg000021l:12156895-12157167
    ## 182 transcript::ptg000021l:13553132-13555030
    ## 183   transcript::ptg000023l:2635562-2636412
    ## 184   transcript::ptg000023l:2690489-2693741
    ## 185 transcript::ptg000023l:18698522-18723762
    ## 186 transcript::ptg000023l:22895865-22909910
    ## 187 transcript::ptg000023l:22896062-22909910
    ## 188 transcript::ptg000023l:22897078-22909910
    ## 189 transcript::ptg000023l:24874248-24874516
    ## 190 transcript::ptg000023l:25121113-25135457
    ## 191 transcript::ptg000023l:25121316-25136081
    ## 192 transcript::ptg000023l:26751020-26756152
    ## 193 transcript::ptg000023l:37964178-37968910
    ## 194 transcript::ptg000023l:38035263-38040884
    ## 195   transcript::ptg000023l:3633190-3636532
    ## 196   transcript::ptg000023l:3633462-3636532
    ## 197 transcript::ptg000023l:14476955-14492814
    ## 198 transcript::ptg000023l:15900851-15906610
    ## 199 transcript::ptg000023l:29261981-29262762
    ## 200   transcript::ptg000023l:7881897-7882272
    ## 201 transcript::ptg000023l:22827886-22828248
    ## 202 transcript::ptg000023l:27548133-27548827
    ## 203   transcript::ptg000024l:3071153-3077218
    ## 204   transcript::ptg000024l:4085668-4087287
    ## 205   transcript::ptg000025l:1976289-1990446
    ## 206 transcript::ptg000025l:10491944-10502595
    ## 207 transcript::ptg000025l:10491944-10502595
    ## 208 transcript::ptg000025l:11114996-11126477
    ## 209 transcript::ptg000025l:11115026-11124835
    ## 210   transcript::ptg000025l:1794902-1798840
    ## 211   transcript::ptg000025l:3300700-3324958
    ## 212 transcript::ptg000025l:10668087-10677231
    ## 213 transcript::ptg000025l:10687372-10693526
    ## 214 transcript::ptg000025l:15677553-15683336
    ## 215 transcript::ptg000025l:17390422-17391830
    ## 216   transcript::ptg000025l:7437716-7438377
    ## 217 transcript::ptg000025l:15784341-15784905
    ## 218 transcript::ptg000025l:15805797-15806833
    ## 219 transcript::ptg000025l:18567820-18568066
    ## 220   transcript::ptg000026l:1848361-1849652
    ## 221   transcript::ptg000026l:1848748-1850032
    ## 222   transcript::ptg000026l:6257379-6258705
    ## 223 transcript::ptg000026l:10976950-10994993
    ## 224 transcript::ptg000026l:10977767-10998419
    ## 225 transcript::ptg000026l:13542914-13544592
    ## 226 transcript::ptg000027l:12369353-12370324
    ## 227 transcript::ptg000027l:14957936-14959033
    ## 228 transcript::ptg000027l:15254638-15255702
    ## 229 transcript::ptg000027l:15899910-15900948
    ## 230     transcript::ptg000028l:607498-608599
    ## 231       transcript::ptg000029c:30378-30679
    ## 232   transcript::ptg000030l:1796595-1796930
    ## 233   transcript::ptg000030l:2571981-2573130
    ## 234   transcript::ptg000030l:2983299-2984288
    ## 235   transcript::ptg000030l:3114544-3114851
    ## 236   transcript::ptg000031l:6673205-6674498
    ## 237     transcript::ptg000033l:812213-813203
    ## 238   transcript::ptg000035l:5335161-5347816
    ## 239   transcript::ptg000035l:5335634-5347816
    ## 240   transcript::ptg000035l:5336170-5347816
    ## 241   transcript::ptg000035l:5336225-5347816
    ## 242   transcript::ptg000035l:5336231-5347816
    ## 243   transcript::ptg000035l:3568395-3588450
    ## 244   transcript::ptg000035l:4179611-4183830
    ## 245   transcript::ptg000035l:4337555-4340587
    ## 246   transcript::ptg000035l:4735323-4735894
    ## 247   transcript::ptg000035l:4745735-4746320
    ## 248   transcript::ptg000035l:4748344-4749044
    ## 249   transcript::ptg000035l:4758339-4759387
    ## 250   transcript::ptg000035l:4761389-4761971
    ## 251   transcript::ptg000035l:4763998-4764581
    ## 252   transcript::ptg000036l:1647657-1649308
    ## 253   transcript::ptg000036l:1647657-1649308
    ## 254   transcript::ptg000036l:3107888-3108320
    ## 255   transcript::ptg000036l:3107888-3108320
    ## 256       transcript::ptg000039l:33382-41881
    ## 257     transcript::ptg000047l:944557-966962
    ## 258     transcript::ptg000047l:944557-967239
    ## 259     transcript::ptg000047l:944577-966962
    ## 260     transcript::ptg000047l:947376-966854
    ## 261   transcript::ptg000047l:8057376-8058484
    ## 262   transcript::ptg000047l:8253896-8258722
    ## 263     transcript::ptg000047l:164714-165595
    ## 264   transcript::ptg000047l:2314993-2323437
    ## 265   transcript::ptg000047l:2314993-2323437
    ## 266   transcript::ptg000047l:2314993-2323437
    ## 267   transcript::ptg000047l:2315100-2323437
    ## 268   transcript::ptg000047l:2315108-2323437
    ## 269   transcript::ptg000047l:2315768-2322229
    ## 270   transcript::ptg000047l:6228838-6229808
    ## 271   transcript::ptg000047l:1248153-1248450
    ## 272   transcript::ptg000047l:9202194-9203168
    ## 273         transcript::ptg000056l:8228-8735
    ## 274       transcript::ptg000057l:12147-12654
    ## 275       transcript::ptg000060c:17547-22922
    ## 276       transcript::ptg000065l:44680-45187
    ## 277         transcript::ptg000099l:1121-2503
    ## 278       transcript::ptg000099l:43788-44669
    ## 279       transcript::ptg000116l:32415-32922
    ## 280       transcript::ptg000121l:37294-37801
    ## 281       transcript::ptg000159l:17542-18049
    ##                                                    sseqid  pident length
    ## 1        Cluster_1832.mature::ntLink_6:5157559-5157579(+) 100.000     21
    ## 2        Cluster_1832.mature::ntLink_6:5157559-5157579(+) 100.000     21
    ## 3        Cluster_1832.mature::ntLink_6:5157559-5157579(+) 100.000     21
    ## 4        Cluster_1826.mature::ntLink_6:4847465-4847486(-) 100.000     22
    ## 5        Cluster_1826.mature::ntLink_6:4847465-4847486(-) 100.000     22
    ## 6      Cluster_1951.mature::ntLink_6:13351801-13351822(-) 100.000     22
    ## 7        Cluster_4220.mature::ptg000007l:915927-915948(-)  90.476     21
    ## 8     Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 9     Cluster_10228.mature::ptg000017l:7471168-7471190(+)  90.909     22
    ## 10    Cluster_18723.mature::ptg000035l:4808391-4808412(+)  95.238     21
    ## 11    Cluster_18723.mature::ptg000035l:4808391-4808412(+)  90.476     21
    ## 12    Cluster_15775.mature::ptg000025l:7472581-7472603(-)  90.909     22
    ## 13    Cluster_15775.mature::ptg000025l:7472581-7472603(-)  90.909     22
    ## 14   Cluster_5012.mature::ptg000008l:10754789-10754809(-)  94.737     19
    ## 15    Cluster_18723.mature::ptg000035l:4808391-4808412(+)  95.238     21
    ## 16    Cluster_18723.mature::ptg000035l:4808391-4808412(+)  95.238     21
    ## 17    Cluster_18723.mature::ptg000035l:4808391-4808412(+)  95.238     21
    ## 18     Cluster_3250.mature::ptg000002l:7337560-7337581(+)  90.909     22
    ## 19     Cluster_3250.mature::ptg000002l:7337560-7337581(+)  90.909     22
    ## 20     Cluster_3250.mature::ptg000002l:7337560-7337581(+)  90.909     22
    ## 21     Cluster_3250.mature::ptg000002l:7337560-7337581(+)  90.909     22
    ## 22   Cluster_2859.mature::ptg000001l:20063094-20063116(+) 100.000     16
    ## 23       Cluster_4220.mature::ptg000007l:915927-915948(-)  95.455     22
    ## 24       Cluster_4220.mature::ptg000007l:915927-915948(-)  95.455     22
    ## 25       Cluster_4220.mature::ptg000007l:915927-915948(-)  95.455     22
    ## 26   Cluster_5012.mature::ptg000008l:10754789-10754809(-)  94.444     18
    ## 27       Cluster_4220.mature::ptg000007l:915927-915948(-)  91.304     23
    ## 28       Cluster_5899.mature::ptg000009l:616336-616357(-)  94.444     18
    ## 29    Cluster_18723.mature::ptg000035l:4808391-4808412(+)  90.909     22
    ## 30    Cluster_18723.mature::ptg000035l:4808391-4808412(+)  90.909     22
    ## 31     Cluster_3437.mature::ptg000004l:1859911-1859933(-) 100.000     15
    ## 32       Cluster_4220.mature::ptg000007l:915927-915948(-)  95.000     20
    ## 33     Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.909     22
    ## 34     Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.476     21
    ## 35   Cluster_5012.mature::ptg000008l:10754789-10754809(-)  95.000     20
    ## 36     Cluster_2463.mature::ptg000001l:5548893-5548914(-) 100.000     22
    ## 37  Cluster_14402.mature::ptg000023l:18708925-18708946(+)  90.476     21
    ## 38        Cluster_19193.mature::ptg000039l:35786-35807(-)  94.737     19
    ## 39       Cluster_4220.mature::ptg000007l:915927-915948(-)  90.909     22
    ## 40     Cluster_2463.mature::ptg000001l:5548893-5548914(-)  90.000     20
    ## 41   Cluster_3367.mature::ptg000002l:14046591-14046614(+) 100.000     24
    ## 42   Cluster_3366.mature::ptg000002l:14046285-14046308(+) 100.000     24
    ## 43       Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 44       Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 45       Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 46  Cluster_14768.mature::ptg000023l:37965298-37965318(+) 100.000     19
    ## 47  Cluster_14768.mature::ptg000023l:37965298-37965318(+) 100.000     19
    ## 48     Cluster_2463.mature::ptg000001l:5548893-5548914(-)  94.737     19
    ## 49     Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.476     21
    ## 50   Cluster_5012.mature::ptg000008l:10754789-10754809(-) 100.000     17
    ## 51   Cluster_5012.mature::ptg000008l:10754789-10754809(-) 100.000     17
    ## 52   Cluster_5012.mature::ptg000008l:10754789-10754809(-) 100.000     17
    ## 53   Cluster_5012.mature::ptg000008l:10754789-10754809(-) 100.000     17
    ## 54  Cluster_14768.mature::ptg000023l:37965298-37965318(+)  94.737     19
    ## 55    Cluster_15775.mature::ptg000025l:7472581-7472603(-)  86.364     22
    ## 56  Cluster_14768.mature::ptg000023l:37965298-37965318(+)  94.737     19
    ## 57    Cluster_10228.mature::ptg000017l:7471168-7471190(+)  94.737     19
    ## 58       Cluster_1862.mature::ntLink_6:7263537-7263560(-) 100.000     16
    ## 59     Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.909     22
    ## 60    Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 61    Cluster_17791.mature::ptg000031l:6751957-6751979(-) 100.000     15
    ## 62       Cluster_4220.mature::ptg000007l:915927-915948(-) 100.000     22
    ## 63    Cluster_15775.mature::ptg000025l:7472581-7472603(-)  89.474     19
    ## 64  Cluster_14402.mature::ptg000023l:18708925-18708946(+) 100.000     17
    ## 65   Cluster_5012.mature::ptg000008l:10754789-10754809(-)  94.444     18
    ## 66  Cluster_14402.mature::ptg000023l:18708925-18708946(+) 100.000     17
    ## 67    Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 68    Cluster_16409.mature::ptg000026l:8745562-8745583(-)  95.238     21
    ## 69  Cluster_10093.mature::ptg000016l:11751407-11751428(-) 100.000     16
    ## 70     Cluster_4254.mature::ptg000007l:3377335-3377356(+) 100.000     15
    ## 71  Cluster_15854.mature::ptg000025l:10668923-10668945(-) 100.000     15
    ## 72   Cluster_3366.mature::ptg000002l:14046285-14046308(+)  91.667     24
    ## 73   Cluster_3367.mature::ptg000002l:14046591-14046614(+)  90.909     22
    ## 74   Cluster_3366.mature::ptg000002l:14046285-14046308(+)  87.500     24
    ## 75   Cluster_3366.mature::ptg000002l:14046285-14046308(+)  87.500     24
    ## 76       Cluster_4220.mature::ptg000007l:915927-915948(-) 100.000     16
    ## 77   Cluster_2859.mature::ptg000001l:20063094-20063116(+) 100.000     15
    ## 78    Cluster_17791.mature::ptg000031l:6751957-6751979(-) 100.000     14
    ## 79    Cluster_15775.mature::ptg000025l:7472581-7472603(-) 100.000     14
    ## 80     Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.476     21
    ## 81     Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.476     21
    ## 82     Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.476     21
    ## 83       Cluster_5900.mature::ptg000009l:616601-616622(-) 100.000     22
    ## 84       Cluster_5899.mature::ptg000009l:616336-616357(-) 100.000     22
    ## 85       Cluster_5900.mature::ptg000009l:616601-616622(-) 100.000     22
    ## 86       Cluster_5899.mature::ptg000009l:616336-616357(-) 100.000     22
    ## 87     Cluster_5981.mature::ptg000009l:4940537-4940559(-) 100.000     23
    ## 88  Cluster_14768.mature::ptg000023l:37965298-37965318(+)  94.737     19
    ## 89     Cluster_3437.mature::ptg000004l:1859911-1859933(-) 100.000     14
    ## 90    Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     18
    ## 91    Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     18
    ## 92    Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     18
    ## 93    Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     18
    ## 94    Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     18
    ## 95    Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     18
    ## 96    Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     18
    ## 97    Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     18
    ## 98   Cluster_2859.mature::ptg000001l:20063094-20063116(+) 100.000     15
    ## 99  Cluster_14402.mature::ptg000023l:18708925-18708946(+)  90.476     21
    ## 100    Cluster_1951.mature::ntLink_6:13351801-13351822(-)  95.000     20
    ## 101 Cluster_15851.mature::ptg000025l:10501052-10501073(+)  94.737     19
    ## 102 Cluster_14768.mature::ptg000023l:37965298-37965318(+)  94.737     19
    ## 103  Cluster_5012.mature::ptg000008l:10754789-10754809(-) 100.000     18
    ## 104 Cluster_14768.mature::ptg000023l:37965298-37965318(+)  94.737     19
    ## 105   Cluster_18723.mature::ptg000035l:4808391-4808412(+)  95.455     22
    ## 106   Cluster_18723.mature::ptg000035l:4808391-4808412(+)  95.455     22
    ## 107   Cluster_18723.mature::ptg000035l:4808391-4808412(+)  95.455     22
    ## 108   Cluster_18723.mature::ptg000035l:4808391-4808412(+)  95.455     22
    ## 109 Cluster_14768.mature::ptg000023l:37965298-37965318(+) 100.000     17
    ## 110 Cluster_14768.mature::ptg000023l:37965298-37965318(+) 100.000     17
    ## 111 Cluster_14768.mature::ptg000023l:37965298-37965318(+) 100.000     17
    ## 112 Cluster_14768.mature::ptg000023l:37965298-37965318(+)  94.737     19
    ## 113   Cluster_10228.mature::ptg000017l:7471168-7471190(+) 100.000     15
    ## 114 Cluster_14768.mature::ptg000023l:37965298-37965318(+) 100.000     19
    ## 115 Cluster_15854.mature::ptg000025l:10668923-10668945(-)  94.118     17
    ## 116  Cluster_5012.mature::ptg000008l:10754789-10754809(-) 100.000     20
    ## 117  Cluster_5012.mature::ptg000008l:10754789-10754809(-) 100.000     18
    ## 118  Cluster_5012.mature::ptg000008l:10754789-10754809(-) 100.000     18
    ## 119       Cluster_19193.mature::ptg000039l:35786-35807(-) 100.000     16
    ## 120   Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 121    Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.476     21
    ## 122      Cluster_4220.mature::ptg000007l:915927-915948(-)  94.444     18
    ## 123   Cluster_10051.mature::ptg000016l:7795530-7795551(+) 100.000     22
    ## 124   Cluster_10057.mature::ptg000016l:8599884-8599905(-) 100.000     22
    ## 125 Cluster_10093.mature::ptg000016l:11751407-11751428(-) 100.000     22
    ## 126 Cluster_10093.mature::ptg000016l:11751407-11751428(-) 100.000     22
    ## 127 Cluster_10093.mature::ptg000016l:11751407-11751428(-) 100.000     22
    ## 128 Cluster_10093.mature::ptg000016l:11751407-11751428(-) 100.000     22
    ## 129    Cluster_2463.mature::ptg000001l:5548893-5548914(-) 100.000     16
    ## 130  Cluster_5012.mature::ptg000008l:10754789-10754809(-) 100.000     16
    ## 131  Cluster_3366.mature::ptg000002l:14046285-14046308(+)  94.737     19
    ## 132  Cluster_3367.mature::ptg000002l:14046591-14046614(+) 100.000     16
    ## 133 Cluster_14402.mature::ptg000023l:18708925-18708946(+)  90.476     21
    ## 134 Cluster_14402.mature::ptg000023l:18708925-18708946(+)  90.476     21
    ## 135   Cluster_10057.mature::ptg000016l:8599884-8599905(-)  89.474     19
    ## 136   Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     22
    ## 137      Cluster_4220.mature::ptg000007l:915927-915948(-)  90.476     21
    ## 138      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 139      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 140      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 141      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 142      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 143      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 144      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 145    Cluster_4254.mature::ptg000007l:3377335-3377356(+) 100.000     14
    ## 146  Cluster_2859.mature::ptg000001l:20063094-20063116(+) 100.000     15
    ## 147   Cluster_15775.mature::ptg000025l:7472581-7472603(-)  86.364     22
    ## 148  Cluster_3366.mature::ptg000002l:14046285-14046308(+)  91.304     23
    ## 149   Cluster_17776.mature::ptg000031l:5461327-5461348(-) 100.000     17
    ## 150      Cluster_1826.mature::ntLink_6:4847465-4847486(-)  94.118     17
    ## 151      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 152      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 153      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 154   Cluster_10228.mature::ptg000017l:7471168-7471190(+)  90.476     21
    ## 155   Cluster_10228.mature::ptg000017l:7471168-7471190(+)  90.476     21
    ## 156      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 157      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 158      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 159      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 160      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 161      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 162      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 163      Cluster_5900.mature::ptg000009l:616601-616622(-)  90.476     21
    ## 164   Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     17
    ## 165   Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     17
    ## 166   Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     17
    ## 167 Cluster_14768.mature::ptg000023l:37965298-37965318(+)  94.737     19
    ## 168   Cluster_15316.mature::ptg000024l:4086254-4086275(+) 100.000     14
    ## 169   Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 170      Cluster_4220.mature::ptg000007l:915927-915948(-)  90.909     22
    ## 171      Cluster_4220.mature::ptg000007l:915927-915948(-)  90.909     22
    ## 172      Cluster_4220.mature::ptg000007l:915927-915948(-)  90.909     22
    ## 173      Cluster_4220.mature::ptg000007l:915927-915948(-)  90.909     22
    ## 174    Cluster_1951.mature::ntLink_6:13351801-13351822(-)  95.238     21
    ## 175    Cluster_1951.mature::ntLink_6:13351801-13351822(-)  95.238     21
    ## 176    Cluster_1951.mature::ntLink_6:13351801-13351822(-)  95.238     21
    ## 177    Cluster_1951.mature::ntLink_6:13351801-13351822(-)  95.238     21
    ## 178      Cluster_4220.mature::ptg000007l:915927-915948(-)  90.909     22
    ## 179    Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.476     21
    ## 180    Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.476     21
    ## 181    Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.476     21
    ## 182 Cluster_14768.mature::ptg000023l:37965298-37965318(+)  94.737     19
    ## 183 Cluster_14768.mature::ptg000023l:37965298-37965318(+)  94.737     19
    ## 184 Cluster_14768.mature::ptg000023l:37965298-37965318(+)  94.737     19
    ## 185 Cluster_14402.mature::ptg000023l:18708925-18708946(+) 100.000     22
    ## 186    Cluster_3250.mature::ptg000002l:7337560-7337581(+)  95.000     20
    ## 187    Cluster_3250.mature::ptg000002l:7337560-7337581(+)  95.000     20
    ## 188    Cluster_3250.mature::ptg000002l:7337560-7337581(+)  95.000     20
    ## 189   Cluster_10057.mature::ptg000016l:8599884-8599905(-)  94.118     17
    ## 190 Cluster_14768.mature::ptg000023l:37965298-37965318(+) 100.000     19
    ## 191 Cluster_14768.mature::ptg000023l:37965298-37965318(+) 100.000     19
    ## 192    Cluster_3437.mature::ptg000004l:1859911-1859933(-)  90.476     21
    ## 193 Cluster_14768.mature::ptg000023l:37965298-37965318(+) 100.000     21
    ## 194 Cluster_14768.mature::ptg000023l:37965298-37965318(+) 100.000     21
    ## 195    Cluster_3250.mature::ptg000002l:7337560-7337581(+)  95.000     20
    ## 196    Cluster_3250.mature::ptg000002l:7337560-7337581(+)  95.000     20
    ## 197 Cluster_14402.mature::ptg000023l:18708925-18708946(+)  95.000     20
    ## 198    Cluster_3437.mature::ptg000004l:1859911-1859933(-)  91.304     23
    ## 199      Cluster_1832.mature::ntLink_6:5157559-5157579(+) 100.000     15
    ## 200   Cluster_15775.mature::ptg000025l:7472581-7472603(-)  94.444     18
    ## 201    Cluster_3437.mature::ptg000004l:1859911-1859933(-)  94.118     17
    ## 202      Cluster_4220.mature::ptg000007l:915927-915948(-) 100.000     15
    ## 203    Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.476     21
    ## 204   Cluster_15316.mature::ptg000024l:4086254-4086275(+) 100.000     22
    ## 205   Cluster_15671.mature::ptg000025l:1981742-1981763(+) 100.000     22
    ## 206 Cluster_15851.mature::ptg000025l:10501052-10501073(+) 100.000     22
    ## 207 Cluster_15851.mature::ptg000025l:10501052-10501073(+) 100.000     22
    ## 208   Cluster_18723.mature::ptg000035l:4808391-4808412(+) 100.000     17
    ## 209   Cluster_18723.mature::ptg000035l:4808391-4808412(+) 100.000     17
    ## 210 Cluster_14768.mature::ptg000023l:37965298-37965318(+)  90.476     21
    ## 211      Cluster_4220.mature::ptg000007l:915927-915948(-)  90.909     22
    ## 212 Cluster_15854.mature::ptg000025l:10668923-10668945(-) 100.000     23
    ## 213   Cluster_18711.mature::ptg000035l:4339600-4339621(-)  94.737     19
    ## 214      Cluster_1826.mature::ntLink_6:4847465-4847486(-)  94.737     19
    ## 215    Cluster_1951.mature::ntLink_6:13351801-13351822(-) 100.000     17
    ## 216   Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 217   Cluster_10051.mature::ptg000016l:7795530-7795551(+)  94.444     18
    ## 218    Cluster_4254.mature::ptg000007l:3377335-3377356(+)  90.000     20
    ## 219   Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.118     17
    ## 220  Cluster_5012.mature::ptg000008l:10754789-10754809(-)  94.444     18
    ## 221  Cluster_5012.mature::ptg000008l:10754789-10754809(-)  94.444     18
    ## 222  Cluster_5012.mature::ptg000008l:10754789-10754809(-)  94.444     18
    ## 223      Cluster_4220.mature::ptg000007l:915927-915948(-)  90.909     22
    ## 224      Cluster_4220.mature::ptg000007l:915927-915948(-)  90.909     22
    ## 225  Cluster_5012.mature::ptg000008l:10754789-10754809(-)  94.737     19
    ## 226   Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 227   Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 228   Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 229   Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 230 Cluster_15854.mature::ptg000025l:10668923-10668945(-)  90.000     20
    ## 231    Cluster_3437.mature::ptg000004l:1859911-1859933(-)  94.118     17
    ## 232   Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 233   Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 234   Cluster_10057.mature::ptg000016l:8599884-8599905(-) 100.000     15
    ## 235    Cluster_5981.mature::ptg000009l:4940537-4940559(-) 100.000     14
    ## 236 Cluster_15854.mature::ptg000025l:10668923-10668945(-) 100.000     15
    ## 237   Cluster_17791.mature::ptg000031l:6751957-6751979(-)  90.000     20
    ## 238   Cluster_18728.mature::ptg000035l:5346054-5346075(+) 100.000     22
    ## 239   Cluster_18728.mature::ptg000035l:5346054-5346075(+) 100.000     22
    ## 240   Cluster_18728.mature::ptg000035l:5346054-5346075(+) 100.000     22
    ## 241   Cluster_18728.mature::ptg000035l:5346054-5346075(+) 100.000     22
    ## 242   Cluster_18728.mature::ptg000035l:5346054-5346075(+) 100.000     22
    ## 243   Cluster_10419.mature::ptg000018l:2286829-2286850(+) 100.000     18
    ## 244      Cluster_4220.mature::ptg000007l:915927-915948(-)  90.909     22
    ## 245   Cluster_18711.mature::ptg000035l:4339600-4339621(-) 100.000     22
    ## 246  Cluster_2859.mature::ptg000001l:20063094-20063116(+) 100.000     15
    ## 247  Cluster_2859.mature::ptg000001l:20063094-20063116(+) 100.000     15
    ## 248  Cluster_2859.mature::ptg000001l:20063094-20063116(+) 100.000     15
    ## 249  Cluster_2859.mature::ptg000001l:20063094-20063116(+) 100.000     15
    ## 250  Cluster_2859.mature::ptg000001l:20063094-20063116(+) 100.000     15
    ## 251  Cluster_2859.mature::ptg000001l:20063094-20063116(+) 100.000     15
    ## 252  Cluster_3367.mature::ptg000002l:14046591-14046614(+) 100.000     16
    ## 253  Cluster_3366.mature::ptg000002l:14046285-14046308(+) 100.000     16
    ## 254  Cluster_3367.mature::ptg000002l:14046591-14046614(+) 100.000     16
    ## 255  Cluster_3366.mature::ptg000002l:14046285-14046308(+) 100.000     16
    ## 256       Cluster_19193.mature::ptg000039l:35786-35807(-) 100.000     22
    ## 257    Cluster_5981.mature::ptg000009l:4940537-4940559(-) 100.000     17
    ## 258    Cluster_5981.mature::ptg000009l:4940537-4940559(-) 100.000     17
    ## 259    Cluster_5981.mature::ptg000009l:4940537-4940559(-) 100.000     17
    ## 260    Cluster_5981.mature::ptg000009l:4940537-4940559(-) 100.000     17
    ## 261  Cluster_2859.mature::ptg000001l:20063094-20063116(+) 100.000     15
    ## 262    Cluster_1951.mature::ntLink_6:13351801-13351822(-)  90.476     21
    ## 263    Cluster_4254.mature::ptg000007l:3377335-3377356(+) 100.000     15
    ## 264 Cluster_15928.mature::ptg000025l:15156437-15156458(+) 100.000     21
    ## 265 Cluster_15928.mature::ptg000025l:15156437-15156458(+) 100.000     21
    ## 266 Cluster_15928.mature::ptg000025l:15156437-15156458(+) 100.000     21
    ## 267 Cluster_15928.mature::ptg000025l:15156437-15156458(+) 100.000     21
    ## 268 Cluster_15928.mature::ptg000025l:15156437-15156458(+) 100.000     21
    ## 269 Cluster_15928.mature::ptg000025l:15156437-15156458(+) 100.000     21
    ## 270  Cluster_3366.mature::ptg000002l:14046285-14046308(+) 100.000     16
    ## 271 Cluster_14768.mature::ptg000023l:37965298-37965318(+)  94.737     19
    ## 272   Cluster_10419.mature::ptg000018l:2286829-2286850(+)  94.444     18
    ## 273    Cluster_4254.mature::ptg000007l:3377335-3377356(+) 100.000     15
    ## 274    Cluster_4254.mature::ptg000007l:3377335-3377356(+) 100.000     15
    ## 275  Cluster_5012.mature::ptg000008l:10754789-10754809(-) 100.000     17
    ## 276    Cluster_4254.mature::ptg000007l:3377335-3377356(+) 100.000     15
    ## 277    Cluster_4254.mature::ptg000007l:3377335-3377356(+) 100.000     15
    ## 278    Cluster_4254.mature::ptg000007l:3377335-3377356(+) 100.000     15
    ## 279    Cluster_4254.mature::ptg000007l:3377335-3377356(+) 100.000     15
    ## 280    Cluster_4254.mature::ptg000007l:3377335-3377356(+) 100.000     15
    ## 281    Cluster_4254.mature::ptg000007l:3377335-3377356(+) 100.000     15
    ##     mismatch gapopen qstart  qend sstart send   evalue bitscore
    ## 1          0       0   2422  2442      1   21 2.22e-06     39.2
    ## 2          0       0   2422  2442      1   21 2.22e-06     39.2
    ## 3          0       0   1219  1239      1   21 1.58e-06     39.2
    ## 4          0       0  36848 36869     22    1 3.96e-06     41.0
    ## 5          0       0   8307  8328     22    1 1.46e-06     41.0
    ## 6          0       0  21847 21868     22    1 3.76e-06     41.0
    ## 7          2       0   4291  4311     21    1 1.00e-03     30.1
    ## 8          1       0    272   289      1   18 1.00e-03     29.2
    ## 9          1       1    853   873     22    1 1.00e-03     28.3
    ## 10         1       0   9861  9881     22    2 2.00e-04     34.6
    ## 11         2       0   3310  3330     22    2 1.00e-03     30.1
    ## 12         2       0    696   717     23    2 1.00e-03     31.9
    ## 13         2       0    696   717     23    2 1.00e-03     31.9
    ## 14         1       0    112   130      3   21 1.13e-04     31.0
    ## 15         1       0  12465 12485      2   22 5.73e-04     34.6
    ## 16         1       0  12463 12483      2   22 3.31e-04     34.6
    ## 17         1       0   9906  9926      2   22 5.37e-04     34.6
    ## 18         2       0   6293  6314     22    1 6.46e-04     31.9
    ## 19         2       0   6257  6278     22    1 1.00e-03     31.9
    ## 20         2       0   6220  6241     22    1 6.41e-04     31.9
    ## 21         2       0   5711  5732     22    1 5.80e-04     31.9
    ## 22         0       0   2902  2917      6   21 1.00e-03     30.1
    ## 23         1       0   4468  4489      1   22 6.82e-05     36.5
    ## 24         1       0   4468  4489      1   22 6.82e-05     36.5
    ## 25         1       0   3305  3326      1   22 6.16e-05     36.5
    ## 26         1       0     78    95      2   19 3.99e-04     29.2
    ## 27         1       1     48    70     22    1 7.90e-05     30.1
    ## 28         1       0    116   133      3   20 6.61e-04     29.2
    ## 29         1       1     32    53      2   22 1.00e-03     28.3
    ## 30         1       1    240   261      2   22 1.00e-03     28.3
    ## 31         0       0    236   250      4   18 5.14e-04     28.3
    ## 32         1       0  12132 12151      2   21 1.00e-03     32.8
    ## 33         2       0   6545  6566     22    1 6.45e-04     31.9
    ## 34         2       0   1831  1851     22    2 1.00e-03     30.1
    ## 35         1       0    873   892      1   20 3.46e-04     32.8
    ## 36         0       0   4587  4608     22    1 3.42e-06     41.0
    ## 37         2       0    203   223     21    1 6.60e-04     30.1
    ## 38         1       0   2028  2046     21    3 1.00e-03     31.0
    ## 39         2       0  15249 15270      1   22 1.00e-03     31.9
    ## 40         2       0   1087  1106     22    3 1.00e-03     28.3
    ## 41         0       0   5552  5575      1   24 1.23e-07     44.6
    ## 42         0       0   5246  5269      1   24 1.23e-07     44.6
    ## 43         2       0   1103  1123     22    2 5.89e-04     30.1
    ## 44         2       0   1764  1784     22    2 5.62e-04     30.1
    ## 45         2       0    870   890     22    2 2.89e-04     30.1
    ## 46         0       0   3056  3074     19    1 2.66e-05     35.6
    ## 47         0       0   3056  3074     19    1 2.66e-05     35.6
    ## 48         1       0    287   305     21    3 1.29e-04     31.0
    ## 49         2       0   1082  1102      2   22 4.87e-04     30.1
    ## 50         0       0  11889 11905     20    4 8.77e-04     31.9
    ## 51         0       0  11852 11868     20    4 8.75e-04     31.9
    ## 52         0       0  11745 11761     20    4 8.68e-04     31.9
    ## 53         0       0  11366 11382     20    4 7.66e-04     31.9
    ## 54         1       0   2546  2564     19    1 1.00e-03     31.0
    ## 55         3       0    250   271     23    2 1.00e-03     27.4
    ## 56         1       0    103   121     19    1 7.98e-05     31.0
    ## 57         1       0   3199  3217      1   19 1.00e-03     31.0
    ## 58         0       0    351   366      1   16 1.00e-03     30.1
    ## 59         2       0   3383  3404     22    1 3.02e-04     31.9
    ## 60         1       0    188   205     18    1 8.98e-04     29.2
    ## 61         0       0    156   170      3   17 5.70e-04     28.3
    ## 62         0       0   7418  7439     22    1 1.51e-06     41.0
    ## 63         2       0    274   292     21    3 1.00e-03     26.5
    ## 64         0       0    644   660      5   21 6.88e-04     31.9
    ## 65         1       0    131   148      1   18 8.06e-04     29.2
    ## 66         0       0    644   660      5   21 5.91e-04     31.9
    ## 67         1       0     33    50     18    1 7.54e-04     29.2
    ## 68         1       0    282   302      1   21 2.16e-05     34.6
    ## 69         0       0   3360  3375      7   22 1.00e-03     30.1
    ## 70         0       0     68    82      8   22 5.84e-04     28.3
    ## 71         0       0    585   599     16    2 8.99e-04     28.3
    ## 72         2       0   2517  2540      1   24 4.23e-05     35.6
    ## 73         2       0   2517  2538      1   22 5.15e-04     31.9
    ## 74         3       0   2361  2384      1   24 1.00e-03     31.0
    ## 75         3       0    269   292      1   24 7.51e-04     31.0
    ## 76         0       0   1793  1808      1   16 7.96e-04     30.1
    ## 77         0       0    903   917     15    1 1.00e-03     28.3
    ## 78         0       0    308   321     14    1 1.00e-03     26.5
    ## 79         0       0     19    32     17    4 1.00e-03     26.5
    ## 80         2       0    911   931     22    2 7.99e-04     30.1
    ## 81         2       0    727   747     22    2 6.92e-04     30.1
    ## 82         2       0    475   495     22    2 6.80e-04     30.1
    ## 83         0       0   4827  4848     22    1 8.69e-07     41.0
    ## 84         0       0   4562  4583     22    1 8.69e-07     41.0
    ## 85         0       0   3000  3021     22    1 7.14e-07     41.0
    ## 86         0       0   2735  2756     22    1 7.14e-07     41.0
    ## 87         0       0  11510 11532     23    1 6.34e-07     42.8
    ## 88         1       0   2798  2816      1   19 1.00e-03     31.0
    ## 89         0       0     85    98      1   14 1.00e-03     26.5
    ## 90         0       0  47655 47672      1   18 6.80e-04     33.7
    ## 91         0       0  33355 33372      1   18 4.76e-04     33.7
    ## 92         0       0  32824 32841      1   18 5.58e-04     33.7
    ## 93         0       0  32824 32841      1   18 5.58e-04     33.7
    ## 94         0       0   1051  1068      1   18 2.99e-05     33.7
    ## 95         0       0  15524 15541      1   18 2.78e-04     33.7
    ## 96         0       0  13548 13565      1   18 2.52e-04     33.7
    ## 97         0       0  13487 13504      1   18 2.33e-04     33.7
    ## 98         0       0    235   249     15    1 3.56e-04     28.3
    ## 99         2       0    998  1018     22    2 1.00e-03     30.1
    ## 100        1       0   1322  1341     21    2 3.00e-04     32.8
    ## 101        1       0   3477  3495      2   20 1.00e-03     31.0
    ## 102        1       0   2230  2248     19    1 7.68e-04     31.0
    ## 103        0       0    153   170      2   19 1.76e-05     33.7
    ## 104        1       0    424   442      1   19 1.00e-03     31.0
    ## 105        0       1   5421  5441      1   22 8.22e-04     32.8
    ## 106        0       1   5421  5441      1   22 1.00e-03     32.8
    ## 107        0       1   5413  5433      1   22 8.87e-04     32.8
    ## 108        0       1   3890  3910      1   22 1.00e-03     32.8
    ## 109        0       0   3790  3806     21    5 3.79e-04     31.9
    ## 110        0       0   3790  3806     21    5 1.00e-03     31.9
    ## 111        0       0   1442  1458     21    5 8.44e-04     31.9
    ## 112        1       0   2678  2696      1   19 1.00e-03     31.0
    ## 113        0       0     57    71     20    6 4.25e-04     28.3
    ## 114        0       0    824   842      1   19 3.17e-05     35.6
    ## 115        1       0    110   126      5   21 1.00e-03     27.4
    ## 116        0       0    729   748      2   21 2.34e-06     37.4
    ## 117        0       0   2677  2694      2   19 7.05e-05     33.7
    ## 118        0       0   1089  1106      2   19 3.93e-05     33.7
    ## 119        0       0    599   614     21    6 5.65e-04     30.1
    ## 120        1       0    217   234     18    1 3.69e-04     29.2
    ## 121        2       0    601   621      2   22 2.22e-04     30.1
    ## 122        1       0    115   132      3   20 1.00e-03     29.2
    ## 123        0       0   7571  7592      1   22 3.36e-06     41.0
    ## 124        0       0   6002  6023     22    1 1.33e-06     41.0
    ## 125        0       0   3059  3080     22    1 1.64e-06     41.0
    ## 126        0       0   3056  3077     22    1 1.64e-06     41.0
    ## 127        0       0   2736  2757     22    1 1.57e-06     41.0
    ## 128        0       0   2692  2713     22    1 1.60e-06     41.0
    ## 129        0       0   2650  2665      7   22 7.76e-04     30.1
    ## 130        0       0   1638  1653      2   17 1.00e-03     30.1
    ## 131        1       0   1617  1635      6   24 1.00e-03     31.0
    ## 132        0       0   1617  1632      6   21 1.00e-03     30.1
    ## 133        2       0   2279  2299      1   21 1.00e-03     30.1
    ## 134        2       0   1971  1991      1   21 8.77e-04     30.1
    ## 135        2       0    214   232      4   22 1.00e-03     26.5
    ## 136        0       0   5833  5854      1   22 1.22e-06     41.0
    ## 137        2       0    138   158     21    1 1.25e-04     30.1
    ## 138        2       0    347   367      2   22 1.00e-03     30.1
    ## 139        2       0    818   838      2   22 1.00e-03     30.1
    ## 140        2       0   1391  1411      2   22 1.00e-03     30.1
    ## 141        2       0    660   680      2   22 1.00e-03     30.1
    ## 142        2       0    231   251      2   22 1.00e-03     30.1
    ## 143        2       0   1391  1411      2   22 1.00e-03     30.1
    ## 144        2       0    660   680      2   22 1.00e-03     30.1
    ## 145        0       0     22    35      5   18 1.00e-03     26.5
    ## 146        0       0    131   145      1   15 5.43e-04     28.3
    ## 147        3       0     19    40     23    2 8.77e-04     27.4
    ## 148        1       1    248   270      2   23 4.38e-04     30.1
    ## 149        0       0   6150  6166      2   18 5.90e-04     31.9
    ## 150        1       0    123   139      6   22 1.00e-03     27.4
    ## 151        2       0      7    27      2   22 1.32e-04     30.1
    ## 152        2       0    114   134     22    2 4.77e-04     30.1
    ## 153        2       0    373   393      2   22 1.74e-04     30.1
    ## 154        2       0   3392  3412     21    1 1.00e-03     30.1
    ## 155        2       0   1925  1945     21    1 1.00e-03     30.1
    ## 156        2       0   1123  1143      2   22 1.00e-03     30.1
    ## 157        2       0    131   151      2   22 4.58e-04     30.1
    ## 158        2       0    663   683      2   22 1.00e-03     30.1
    ## 159        2       0     46    66      2   22 1.00e-03     30.1
    ## 160        2       0    636   656     22    2 2.30e-04     30.1
    ## 161        2       0    963   983     22    2 3.15e-04     30.1
    ## 162        2       0   1139  1159      2   22 1.00e-03     30.1
    ## 163        2       0    141   161      2   22 1.00e-03     30.1
    ## 164        0       0   9218  9234      6   22 1.00e-03     31.9
    ## 165        0       0   8074  8090      6   22 8.61e-04     31.9
    ## 166        0       0   3248  3264      6   22 4.51e-04     31.9
    ## 167        1       0   3967  3985     19    1 1.00e-03     31.0
    ## 168        0       0     56    69     16    3 1.00e-03     26.5
    ## 169        1       0    349   366      1   18 4.39e-04     29.2
    ## 170        2       0   8594  8615      1   22 8.24e-04     31.9
    ## 171        2       0   8065  8086      1   22 7.88e-04     31.9
    ## 172        2       0   5792  5813      1   22 1.00e-03     31.9
    ## 173        2       0   5510  5531      1   22 6.13e-04     31.9
    ## 174        1       0  30903 30923     22    2 7.20e-04     34.6
    ## 175        1       0  10346 10366     22    2 5.10e-04     34.6
    ## 176        1       0  10154 10174     22    2 3.99e-04     34.6
    ## 177        1       0   9982 10002     22    2 4.81e-04     34.6
    ## 178        2       0    214   235      1   22 1.00e-03     31.9
    ## 179        2       0   3828  3848      2   22 1.00e-03     30.1
    ## 180        2       0   3829  3849      2   22 1.00e-03     30.1
    ## 181        2       0     58    78      2   22 9.68e-05     30.1
    ## 182        1       0    358   376     19    1 5.75e-04     31.0
    ## 183        1       0    612   630      1   19 2.55e-04     31.0
    ## 184        1       0    868   886     19    1 8.81e-04     31.0
    ## 185        0       0  10403 10424      1   22 2.88e-06     41.0
    ## 186        1       0   8344  8363      1   20 8.29e-04     32.8
    ## 187        1       0   8147  8166      1   20 8.17e-04     32.8
    ## 188        1       0   7131  7150      1   20 7.57e-04     32.8
    ## 189        1       0    112   128     17    1 1.00e-03     27.4
    ## 190        0       0  10626 10644      1   19 6.95e-05     35.6
    ## 191        0       0  10423 10441      1   19 7.15e-05     35.6
    ## 192        2       0   3588  3608      1   21 1.00e-03     30.1
    ## 193        0       0   1120  1140      1   21 2.48e-06     39.2
    ## 194        0       0   1125  1145      1   21 2.59e-06     39.2
    ## 195        1       0   2913  2932      1   20 2.59e-04     32.8
    ## 196        1       0   2641  2660      1   20 2.38e-04     32.8
    ## 197        1       0  12881 12900     20    1 1.00e-03     32.8
    ## 198        1       1   5447  5469      1   22 1.00e-03     30.1
    ## 199        0       0    572   586     15    1 1.00e-03     28.3
    ## 200        1       0     95   112     21    4 4.28e-04     29.2
    ## 201        1       0    155   171     17    1 1.00e-03     27.4
    ## 202        0       0    265   279     16    2 8.05e-04     28.3
    ## 203        2       0   1408  1428      2   22 1.00e-03     30.1
    ## 204        0       0    586   607      1   22 2.71e-07     41.0
    ## 205        0       0   5453  5474      1   22 1.61e-06     41.0
    ## 206        0       0   9108  9129      1   22 1.41e-06     41.0
    ## 207        0       0   9108  9129      1   22 1.41e-06     41.0
    ## 208        0       0   2733  2749     22    6 7.86e-04     31.9
    ## 209        0       0   2703  2719     22    6 6.71e-04     31.9
    ## 210        2       0   3441  3461     21    1 1.00e-03     30.1
    ## 211        2       0  15984 16005     22    1 1.00e-03     31.9
    ## 212        0       0    836   858     23    1 3.46e-07     42.8
    ## 213        1       0     32    50      4   22 1.00e-03     31.0
    ## 214        1       0   2110  2128     19    1 1.00e-03     31.0
    ## 215        0       0    279   295      3   19 1.22e-04     31.9
    ## 216        1       0    418   435      1   18 7.66e-04     29.2
    ## 217        1       0    258   275      4   21 6.51e-04     29.2
    ## 218        2       0    794   813      1   20 1.00e-03     28.3
    ## 219        1       0     83    99     17    1 1.00e-03     27.4
    ## 220        1       0   1082  1099      4   21 1.00e-03     29.2
    ## 221        1       0    695   712      4   21 1.00e-03     29.2
    ## 222        1       0    736   753      4   21 1.00e-03     29.2
    ## 223        2       0  13088 13109      1   22 1.00e-03     31.9
    ## 224        2       0  12271 12292      1   22 1.00e-03     31.9
    ## 225        1       0    252   270      2   20 5.08e-04     31.0
    ## 226        1       0    793   810     18    1 1.00e-03     29.2
    ## 227        1       0    609   626      1   18 1.00e-03     29.2
    ## 228        1       0    609   626      1   18 1.00e-03     29.2
    ## 229        1       0    572   589      1   18 1.00e-03     29.2
    ## 230        2       0    610   629     21    2 1.00e-03     28.3
    ## 231        1       0    124   140      1   17 1.00e-03     27.4
    ## 232        1       0    192   209     18    1 4.20e-04     29.2
    ## 233        1       0    629   646      1   18 1.00e-03     29.2
    ## 234        0       0    566   580      3   17 1.00e-03     28.3
    ## 235        0       0    271   284     21    8 1.00e-03     26.5
    ## 236        0       0    822   836     21    7 1.00e-03     28.3
    ## 237        2       0    849   868     20    1 1.00e-03     28.3
    ## 238        0       0  10893 10914      1   22 1.67e-06     41.0
    ## 239        0       0  10420 10441      1   22 1.61e-06     41.0
    ## 240        0       0   9884  9905      1   22 1.54e-06     41.0
    ## 241        0       0   9829  9850      1   22 1.53e-06     41.0
    ## 242        0       0   9823  9844      1   22 1.53e-06     41.0
    ## 243        0       0  16964 16981     18    1 3.39e-04     33.7
    ## 244        2       0   3798  3819     22    1 3.28e-04     31.9
    ## 245        0       0   2045  2066     22    1 4.54e-07     41.0
    ## 246        0       0    144   158     15    1 6.60e-04     28.3
    ## 247        0       0    144   158     15    1 6.76e-04     28.3
    ## 248        0       0    144   158     15    1 8.12e-04     28.3
    ## 249        0       0    585   599     15    1 1.00e-03     28.3
    ## 250        0       0    144   158     15    1 6.73e-04     28.3
    ## 251        0       0    144   158     15    1 6.74e-04     28.3
    ## 252        0       0   1111  1126     21    6 5.00e-04     30.1
    ## 253        0       0   1111  1126     21    6 5.00e-04     30.1
    ## 254        0       0    130   145     21    6 1.42e-04     30.1
    ## 255        0       0    130   145     21    6 1.42e-04     30.1
    ## 256        0       0   2404  2425     22    1 1.12e-06     41.0
    ## 257        0       0   9948  9964     22    6 1.00e-03     31.9
    ## 258        0       0   9948  9964     22    6 1.00e-03     31.9
    ## 259        0       0   9928  9944     22    6 1.00e-03     31.9
    ## 260        0       0   7129  7145     22    6 1.00e-03     31.9
    ## 261        0       0   1028  1042     15    1 1.00e-03     28.3
    ## 262        2       0   1495  1515      2   22 1.00e-03     30.1
    ## 263        0       0     68    82      8   22 1.00e-03     28.3
    ## 264        0       0   3611  3631      1   21 3.89e-06     39.2
    ## 265        0       0   3611  3631      1   21 3.89e-06     39.2
    ## 266        0       0   3611  3631      1   21 3.89e-06     39.2
    ## 267        0       0   3504  3524      1   21 3.84e-06     39.2
    ## 268        0       0   3496  3516      1   21 3.84e-06     39.2
    ## 269        0       0   2836  2856      1   21 2.98e-06     39.2
    ## 270        0       0    366   381     24    9 2.92e-04     30.1
    ## 271        1       0    114   132      1   19 1.06e-04     31.0
    ## 272        1       0    357   374      1   18 1.00e-03     29.2
    ## 273        0       0     68    82      8   22 5.84e-04     28.3
    ## 274        0       0     68    82      8   22 5.84e-04     28.3
    ## 275        0       0   2836  2852     20    4 3.67e-04     31.9
    ## 276        0       0     68    82      8   22 5.84e-04     28.3
    ## 277        0       0     68    82      8   22 1.00e-03     28.3
    ## 278        0       0     68    82      8   22 1.00e-03     28.3
    ## 279        0       0     68    82      8   22 5.84e-04     28.3
    ## 280        0       0     68    82      8   22 5.84e-04     28.3
    ## 281        0       0     68    82      8   22 5.84e-04     28.3

281 putative lncRNA sponges with these parameters.

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
    ## 3010800
    ## Parsing output
    ## counting number of putative interactions predicted (can include multiple interactions between single miRNA-lncRNA pair)
    ## 24962563 ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict_all.tab

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
    ## 19024
    ## 
    ## check some:
    ## >Cluster_10051.mature::ptg000016l:7795530-7795551(+) transcript::ntLink_6:10103156-10120487  155.00  -21.12  2 21    9484 9507   21  66.67%  76.19%
    ## >Cluster_10051.mature::ptg000016l:7795530-7795551(+) transcript::ntLink_6:10103183-10120487  155.00  -21.12  2 21    9457 9480   21  66.67%  76.19%
    ## >Cluster_10051.mature::ptg000016l:7795530-7795551(+) transcript::ntLink_6:10103183-10120487  155.00  -21.12  2 21    9457 9480   21  66.67%  76.19%
    ## >Cluster_10051.mature::ptg000016l:7795530-7795551(+) transcript::ntLink_6:10103183-10120487  155.00  -21.12  2 21    9457 9480   21  66.67%  76.19%
    ## >Cluster_10051.mature::ptg000016l:7795530-7795551(+) transcript::ntLink_6:10104556-10114046  155.00  -21.12  2 21    8084 8107   21  66.67%  76.19%

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
    ## 163
    ## 
    ## check some:
    ## >Cluster_10228.mature::ptg000017l:7471168-7471190(+) transcript::ptg000019l:1503591-1507652  188.00  -31.70  2 22    3554 3577   21  85.71%  95.24%
    ## >Cluster_10228.mature::ptg000017l:7471168-7471190(+) transcript::ptg000019l:1505407-1507652  188.00  -31.70  2 22    1738 1761   21  85.71%  95.24%
    ## >Cluster_10228.mature::ptg000017l:7471168-7471190(+) transcript::ptg000019l:1505407-1507652  188.00  -31.70  2 22    1738 1761   21  85.71%  95.24%
    ## >Cluster_10228.mature::ptg000017l:7471168-7471190(+) transcript::ptg000039l:790534-797321    188.00  -25.92  2 22    5787 5810   21  90.48%  90.48%
    ## >Cluster_10228.mature::ptg000017l:7471168-7471190(+) transcript::ptg000039l:790555-797321    188.00  -25.92  2 22    5766 5789   21  90.48%  90.48%

This is a dramatically smaller number – only 163 interactions are at
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
lncRNAfasta=../output/10-Apul-lncRNA/Apul_lncRNA.fasta
lncRNAcoordinates=../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul_lncRNA_coordinates.txt
candidategtf=../output/10-Apul-lncRNA/Apul_lncRNA_candidates.gtf
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
lncRNAfasta=../output/10-Apul-lncRNA/Apul_lncRNA.fasta
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
    ## >transcript::ntLink_6:12654168-12654858
    ## >transcript::ntLink_6:13210102-13212000
    ## >transcript::ntLink_6:13210122-13212000
    ## 
    ## same index of filtered lncRNA gtf sequences: 
    ## ntLink_6 StringTie   transcript  12654169    12654858    .   +   .   transcript_id "MSTRG.1611.1"; gene_id "MSTRG.1611"; xloc "XLOC_000892"; class_code "u"; tss_id "TSS1327";
    ## ntLink_6 StringTie   transcript  13210103    13212000    .   +   .   transcript_id "MSTRG.1688.1"; gene_id "MSTRG.1688"; xloc "XLOC_000920"; class_code "u"; tss_id "TSS1373";
    ## ntLink_6 StringTie   transcript  13210123    13212000    .   +   .   transcript_id "MSTRG.1688.2"; gene_id "MSTRG.1688"; xloc "XLOC_000920"; class_code "u"; tss_id "TSS1373";
    ## 
    ## number of lncRNA fasta sequences: 
    ## 24183
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

    ## Number of putative lncRNA precursors:  36 
    ##  Number of miRNA whose precursors are lncRNA:  23

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
    ## 24962563 ../output/17-Apul-miRNA-lncRNA-BLASTs-RNAhybrid/Apul-miRanda-lncRNA-strict_all.tab
    ## 
    ## Number of putative interactions of at least 21 nucleotides
    ## 19024
    ## 
    ## Number of putative interactions of at least 21 nucleotides, with at most 3 mismatches
    ## 163

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
    ##  Number of A.pulchra lncRNAs:  77200 
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
