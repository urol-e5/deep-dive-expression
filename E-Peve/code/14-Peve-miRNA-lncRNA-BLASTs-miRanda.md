14-Peve-miRNA-lncRNA-BLASTs-miRanda
================
Kathleen Durkin
2025-06-17

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
BLASTn, but will likely need to include additional methods like miRanda
to identify potential binding.**

# 1 Prep for BLASTs

## 1.1 Isolate the pre-mirna and mature mirna sequences

``` bash
full_mirna_fasta="../output/05-Peve-sRNA-ShortStack_4.1.0/ShortStack_out/mir.fasta"
premirna_fasta="../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_precursor.fasta"
mature_mirna_fasta="../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_mature.fasta"
star_mirna_fasta="../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_star.fasta"

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
premirna_fasta="../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_precursor.fasta"
mature_mirna_fasta="../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_mature.fasta"
star_mirna_fasta="../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_star.fasta"

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

    ## >Cluster_29::Porites_evermani_scaffold_1:1404250-1404342(-)
    ## >Cluster_589::Porites_evermani_scaffold_16:383386-383478(-)
    ## 
    ## >Cluster_29.mature::Porites_evermani_scaffold_1:1404272-1404293(-)
    ## >Cluster_589.mature::Porites_evermani_scaffold_16:383437-383458(-)
    ## 
    ## >Cluster_29.star::Porites_evermani_scaffold_1:1404301-1404322(-)
    ## >Cluster_589.star::Porites_evermani_scaffold_16:383406-383427(-)
    ## 
    ## 45
    ## 
    ## 45
    ## 
    ## 45

## 1.2 Check miRNA lengths

``` bash
# Extract sequence lengths for precursors
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_precursor.fasta > ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_precursor_lengths.txt

# Sequence lengths for matures
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_mature.fasta > ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_mature_lengths.txt
```

``` r
# Summary stats of precursor and mature lengths

precursor_lengths <- read.table("../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_precursor_lengths.txt", sep = " ", header = FALSE, col.names = c("seqID", "length"))
mature_lengths <- read.table("../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_mature_lengths.txt", sep = " ", header = FALSE, col.names = c("seqID", "length"))

cat("Average pre-miRNA length: ", mean(precursor_lengths$length))
```

    ## Average pre-miRNA length:  93.51111

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

    ## Average mature miRNA length:  21.91111

``` r
cat("\n")
```

``` r
cat("Range of mature miRNA lengths: ", range(mature_lengths$length))
```

    ## Range of mature miRNA lengths:  21 23

## 1.3 check lncRNAs

LncRNAs were identified from Peve RNA-seq data

Fasta of Peve lncRNAs stored at
`https://gannet.fish.washington.edu/acropora/E5-deep-dive-expression/output/01.6-lncRNA-pipline/Peve_lncRNA.fasta`

``` bash
curl -L https://gannet.fish.washington.edu/acropora/E5-deep-dive-expression/output/01.6-lncRNA-pipline/Peve_lncRNA.fasta -o ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_lncRNA.fasta

echo "Number of lncRNAs:"
grep "^>" ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_lncRNA.fasta | wc -l
```

    ##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
    ##                                  Dload  Upload   Total   Spent    Left  Speed
    ##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0 18  298M   18 55.4M    0     0  92.9M      0  0:00:03 --:--:--  0:00:03 92.8M 55  298M   55  165M    0     0   103M      0  0:00:02  0:00:01  0:00:01  103M 91  298M   91  274M    0     0   105M      0  0:00:02  0:00:02 --:--:--  105M100  298M  100  298M    0     0   106M      0  0:00:02  0:00:02 --:--:--  105M
    ## Number of lncRNAs:
    ## 49811

``` bash
# Extract sequence lengths for precursors
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_lncRNA.fasta > ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_lncRNA_lengths.txt
```

``` r
# Summary stats of lncRNA lengths

lncRNA_lengths <- read.table("../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_lncRNA_lengths.txt", sep = " ", header = FALSE, col.names = c("seqID", "length"))

cat("Average lncRNA length: ", mean(lncRNA_lengths$length))
```

    ## Average lncRNA length:  6134.399

``` r
cat("\n")
```

``` r
cat("Range of lncRNA lengths: ", range(lncRNA_lengths$length))
```

    ## Range of lncRNA lengths:  201 107091

``` r
ggplot(lncRNA_lengths, aes(x = length)) +
  geom_histogram(binwidth = 500) +
  labs(title = "A. pulchra lncRNA sequence lengths",
       x = "Sequence Length [nucleotides]",
       y = "Frequency") +
  xlim(200, 110000) +
  ylim(0, 1000) +
  theme_minimal()
```

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](14-Peve-miRNA-lncRNA-BLASTs-miRanda_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# 2 BLASTs

## 2.1 Make databases

Database of pre-miRNAs:

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_precursor.fasta \
-dbtype nucl \
-out ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/Peve-db/Peve_ShortStack_4.1.0_precursor
```

Database of mature miRNAs:

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_mature.fasta \
-dbtype nucl \
-out ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/Peve-db/Peve_ShortStack_4.1.0_mature
```

## 2.2 Run BLASTn

Generate a list of blast results. It seems plausible that a single
lncRNA, which would be hundreds or thousands of nucleotides long, could
interact with multiple miRNAs, so I will allow up to 10 hits (\~25% of
Peve miRNAs) for each lncRNA. I want to see the top hits no matter how
poor the match is, so I will not filter by e-value at this stage. I’ll
also include the “-word_size 4” option, which reduces the required
length of the initial match.

Full pre-miRNAs:

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_lncRNA.fasta \
-db ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/Peve-db/Peve_ShortStack_4.1.0_precursor \
-out ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_precursor_blastn.tab \
-num_threads 40 \
-word_size 4 \
-max_target_seqs 10 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_precursor_blastn.tab
```

    ## 479117 ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_precursor_blastn.tab

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
-query ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_lncRNA.fasta \
-db ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/Peve-db/Peve_ShortStack_4.1.0_mature \
-out ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_mature_blastn.tab \
-num_threads 40 \
-word_size 4 \
-max_target_seqs 10 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_mature_blastn.tab
```

    ## 476143 ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_mature_blastn.tab

# 3 Examine BLAST tables

Read into R and assign informative column labels

``` r
precursor_lncRNA_BLASTn <- read.table("../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_precursor_blastn.tab", sep="\t", header=FALSE)
mature_lncRNA_BLASTn <- read.table("../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_mature_blastn.tab", sep="\t", header=FALSE)

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

    ## [1] 37

``` r
precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  select(qseqid) %>%
  unique() %>%
  nrow()
```

    ## [1] 36

``` r
precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  select(sseqid) %>%
  unique() %>%
  nrow()
```

    ## [1] 21

``` r
precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  unique()
```

    ##                                                      qseqid
    ## 1   transcript::Porites_evermani_scaffold_1:1380412-1417060
    ## 2   transcript::Porites_evermani_scaffold_1:1380748-1417060
    ## 3  transcript::Porites_evermani_scaffold_1128:109236-127279
    ## 4  transcript::Porites_evermani_scaffold_1128:109400-127279
    ## 5     transcript::Porites_evermani_scaffold_1159:3706-16094
    ## 6     transcript::Porites_evermani_scaffold_1159:3706-16094
    ## 7    transcript::Porites_evermani_scaffold_1429:32475-53408
    ## 8    transcript::Porites_evermani_scaffold_1429:46011-48479
    ## 9     transcript::Porites_evermani_scaffold_148:92027-98639
    ## 10   transcript::Porites_evermani_scaffold_16:363080-391056
    ## 11   transcript::Porites_evermani_scaffold_1732:69845-79966
    ## 12   transcript::Porites_evermani_scaffold_1732:72129-79966
    ## 13   transcript::Porites_evermani_scaffold_26:377983-392680
    ## 14   transcript::Porites_evermani_scaffold_26:378057-392680
    ## 15    transcript::Porites_evermani_scaffold_316:75149-91320
    ## 16    transcript::Porites_evermani_scaffold_316:75149-91356
    ## 18    transcript::Porites_evermani_scaffold_316:75223-97019
    ## 19  transcript::Porites_evermani_scaffold_334:150792-156024
    ## 20  transcript::Porites_evermani_scaffold_334:151379-156024
    ## 21  transcript::Porites_evermani_scaffold_461:195334-219192
    ## 22  transcript::Porites_evermani_scaffold_461:195351-223079
    ## 23  transcript::Porites_evermani_scaffold_461:210393-223079
    ## 24   transcript::Porites_evermani_scaffold_47:472593-482329
    ## 25   transcript::Porites_evermani_scaffold_49:147530-165944
    ## 26  transcript::Porites_evermani_scaffold_613:139660-180006
    ## 28  transcript::Porites_evermani_scaffold_613:139663-179070
    ## 29  transcript::Porites_evermani_scaffold_613:139663-180006
    ## 30  transcript::Porites_evermani_scaffold_613:139962-180006
    ## 31    transcript::Porites_evermani_scaffold_866:21937-24095
    ## 32   transcript::Porites_evermani_scaffold_910:95262-102133
    ## 33   transcript::Porites_evermani_scaffold_910:95900-102133
    ## 34  transcript::Porites_evermani_scaffold_910:108432-120843
    ## 35  transcript::Porites_evermani_scaffold_910:129984-141109
    ## 36  transcript::Porites_evermani_scaffold_910:134861-141109
    ## 37  transcript::Porites_evermani_scaffold_942:121381-138291
    ## 38    transcript::Porites_evermani_scaffold_984:51531-58363
    ## 39    transcript::Porites_evermani_scaffold_984:51539-59053
    ##                                                           sseqid pident length
    ## 1     Cluster_29::Porites_evermani_scaffold_1:1404250-1404342(-)    100     93
    ## 2     Cluster_29::Porites_evermani_scaffold_1:1404250-1404342(-)    100     93
    ## 3  Cluster_9983::Porites_evermani_scaffold_1128:125526-125618(-)    100     93
    ## 4  Cluster_9983::Porites_evermani_scaffold_1128:125526-125618(-)    100     93
    ## 5     Cluster_10061::Porites_evermani_scaffold_1159:7406-7500(+)    100     95
    ## 6     Cluster_10060::Porites_evermani_scaffold_1159:6653-6743(+)    100     91
    ## 7   Cluster_10965::Porites_evermani_scaffold_1429:47235-47326(-)    100     92
    ## 8   Cluster_10965::Porites_evermani_scaffold_1429:47235-47326(-)    100     92
    ## 9     Cluster_2882::Porites_evermani_scaffold_148:94531-94624(-)    100     94
    ## 10    Cluster_589::Porites_evermani_scaffold_16:383386-383478(-)    100     93
    ## 11  Cluster_11997::Porites_evermani_scaffold_1732:76482-76574(+)    100     93
    ## 12  Cluster_11997::Porites_evermani_scaffold_1732:76482-76574(+)    100     93
    ## 13    Cluster_796::Porites_evermani_scaffold_26:382550-382645(-)    100     96
    ## 14    Cluster_796::Porites_evermani_scaffold_26:382550-382645(-)    100     96
    ## 15    Cluster_4629::Porites_evermani_scaffold_316:88415-88506(-)    100     92
    ## 16    Cluster_4629::Porites_evermani_scaffold_316:88415-88506(-)    100     92
    ## 18    Cluster_4629::Porites_evermani_scaffold_316:88415-88506(-)    100     92
    ## 19  Cluster_4735::Porites_evermani_scaffold_334:153554-153646(-)    100     93
    ## 20  Cluster_4735::Porites_evermani_scaffold_334:153554-153646(-)    100     93
    ## 21  Cluster_5882::Porites_evermani_scaffold_461:215454-215549(+)    100     96
    ## 22  Cluster_5882::Porites_evermani_scaffold_461:215454-215549(+)    100     96
    ## 23  Cluster_5882::Porites_evermani_scaffold_461:215454-215549(+)    100     96
    ## 24   Cluster_1140::Porites_evermani_scaffold_47:475972-476066(-)    100     95
    ## 25   Cluster_1167::Porites_evermani_scaffold_49:151587-151681(-)    100     95
    ## 26  Cluster_7053::Porites_evermani_scaffold_613:156453-156546(+)    100     94
    ## 28  Cluster_7053::Porites_evermani_scaffold_613:156453-156546(+)    100     94
    ## 29  Cluster_7053::Porites_evermani_scaffold_613:156453-156546(+)    100     94
    ## 30  Cluster_7053::Porites_evermani_scaffold_613:156453-156546(+)    100     94
    ## 31    Cluster_8634::Porites_evermani_scaffold_866:22803-22895(-)    100     93
    ## 32    Cluster_8884::Porites_evermani_scaffold_910:99233-99322(+)    100     90
    ## 33    Cluster_8884::Porites_evermani_scaffold_910:99233-99322(+)    100     90
    ## 34  Cluster_8887::Porites_evermani_scaffold_910:118720-118809(+)    100     90
    ## 35  Cluster_8888::Porites_evermani_scaffold_910:139331-139420(+)    100     90
    ## 36  Cluster_8888::Porites_evermani_scaffold_910:139331-139420(+)    100     90
    ## 37  Cluster_8988::Porites_evermani_scaffold_942:133648-133739(+)    100     92
    ## 38    Cluster_9149::Porites_evermani_scaffold_984:51832-51924(-)    100     93
    ## 39    Cluster_9149::Porites_evermani_scaffold_984:51832-51924(-)    100     93
    ##    mismatch gapopen qstart  qend sstart send   evalue bitscore
    ## 1         0       0  23838 23930     93    1 1.62e-43      168
    ## 2         0       0  23502 23594     93    1 1.60e-43      168
    ## 3         0       0  16290 16382     93    1 8.08e-44      168
    ## 4         0       0  16126 16218     93    1 8.01e-44      168
    ## 5         0       0   3700  3794      1   95 4.55e-45      172
    ## 6         0       0   2947  3037      1   91 6.75e-43      165
    ## 7         0       0  14760 14851     92    1 3.27e-43      167
    ## 8         0       0   1224  1315     92    1 3.94e-44      167
    ## 9         0       0   2504  2597     94    1 8.58e-45      170
    ## 10        0       0  20306 20398     93    1 1.24e-43      168
    ## 11        0       0   6637  6729      1   93 4.53e-44      168
    ## 12        0       0   4353  4445      1   93 3.55e-44      168
    ## 13        0       0   4567  4662     96    1 1.55e-45      174
    ## 14        0       0   4493  4588     96    1 1.54e-45      174
    ## 15        0       0  13266 13357     92    1 2.53e-43      167
    ## 16        0       0  13266 13357     92    1 2.53e-43      167
    ## 18        0       0  13192 13283     92    1 3.36e-43      167
    ## 19        0       0   2762  2854     93    1 2.37e-44      168
    ## 20        0       0   2175  2267     93    1 2.10e-44      168
    ## 21        0       0  20120 20215      1   96 2.48e-45      174
    ## 22        0       0  20103 20198      1   96 2.88e-45      174
    ## 23        0       0   5061  5156      1   96 1.34e-45      174
    ## 24        0       0   3379  3473     95    1 3.58e-45      172
    ## 25        0       0   4057  4151     95    1 6.77e-45      172
    ## 26        0       0  16793 16886      1   94 5.11e-44      170
    ## 28        0       0  16790 16883      1   94 4.99e-44      170
    ## 29        0       0  16790 16883      1   94 5.11e-44      170
    ## 30        0       0  16491 16584      1   94 5.07e-44      170
    ## 31        0       0    866   958     93    1 9.85e-45      168
    ## 32        0       0   3971  4060      1   90 1.32e-42      163
    ## 33        0       0   3333  3422      1   90 1.20e-42      163
    ## 34        0       0  10288 10377      1   90 2.36e-42      163
    ## 35        0       0   9347  9436      1   90 2.12e-42      163
    ## 36        0       0   4470  4559      1   90 1.20e-42      163
    ## 37        0       0  12267 12358      1   92 2.64e-43      167
    ## 38        0       0    301   393     93    1 3.10e-44      168
    ## 39        0       0    293   385     93    1 3.40e-44      168

We have 37 alignments of a full pre-miRNA to a lncRNA with no
mismatches. 36 lncRNA and 21 miRNA are represented.

Note that, as in A.pulchra, there are instances of a single pre-miRNA
matching to several overlapping lncRNA. For example, Cluster_29 is
contained within transcript::Porites_evermani_scaffold_1:1380412-1417060
and transcript::Porites_evermani_scaffold_1:1380748-1417060. This may
represent multiple isoforms of a single lncRNA gene.

Save these results

``` r
precursor_lncRNAs <- precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  unique()

write.table(precursor_lncRNAs, "../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/lncRNAs_as_miRNA_precursors.txt")
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

    ##                                                       qseqid
    ## 1      transcript::Porites_evermani_scaffold_1:920565-973715
    ## 2      transcript::Porites_evermani_scaffold_1:920565-973715
    ## 3      transcript::Porites_evermani_scaffold_1:943644-974204
    ## 4      transcript::Porites_evermani_scaffold_1:943644-974204
    ## 5    transcript::Porites_evermani_scaffold_1:1585708-1591493
    ## 6    transcript::Porites_evermani_scaffold_1:1380412-1417060
    ## 7    transcript::Porites_evermani_scaffold_1:1380748-1417060
    ## 8      transcript::Porites_evermani_scaffold_1:708835-709158
    ## 9   transcript::Porites_evermani_scaffold_1022:140030-141073
    ## 10  transcript::Porites_evermani_scaffold_1022:140176-141073
    ## 11   transcript::Porites_evermani_scaffold_106:125081-125454
    ## 12   transcript::Porites_evermani_scaffold_107:236779-275156
    ## 13   transcript::Porites_evermani_scaffold_107:236779-275156
    ## 14   transcript::Porites_evermani_scaffold_107:237948-275156
    ## 15   transcript::Porites_evermani_scaffold_107:237948-275156
    ## 16  transcript::Porites_evermani_scaffold_1079:109501-123509
    ## 17  transcript::Porites_evermani_scaffold_1079:109501-123509
    ## 18   transcript::Porites_evermani_scaffold_108:201235-201493
    ## 19     transcript::Porites_evermani_scaffold_111:73252-73642
    ## 20  transcript::Porites_evermani_scaffold_1128:109236-127279
    ## 21  transcript::Porites_evermani_scaffold_1128:109400-127279
    ## 22   transcript::Porites_evermani_scaffold_113:467406-480376
    ## 23   transcript::Porites_evermani_scaffold_113:467406-480376
    ## 24   transcript::Porites_evermani_scaffold_113:467406-480376
    ## 25   transcript::Porites_evermani_scaffold_113:467409-480376
    ## 26   transcript::Porites_evermani_scaffold_113:467409-480376
    ## 27   transcript::Porites_evermani_scaffold_113:467454-480376
    ## 28     transcript::Porites_evermani_scaffold_1159:3706-16094
    ## 29     transcript::Porites_evermani_scaffold_1159:3706-16094
    ## 30  transcript::Porites_evermani_scaffold_1162:119103-143883
    ## 31  transcript::Porites_evermani_scaffold_1162:119125-143883
    ## 32  transcript::Porites_evermani_scaffold_1162:132444-143883
    ## 33  transcript::Porites_evermani_scaffold_1162:129693-137434
    ## 34    transcript::Porites_evermani_scaffold_1165:11453-23290
    ## 35   transcript::Porites_evermani_scaffold_117:390567-391381
    ## 36    transcript::Porites_evermani_scaffold_1176:78697-83146
    ## 37  transcript::Porites_evermani_scaffold_1176:126058-131550
    ## 38    transcript::Porites_evermani_scaffold_1204:91416-92378
    ## 39    transcript::Porites_evermani_scaffold_1210:30896-42143
    ## 40    transcript::Porites_evermani_scaffold_1210:30896-42143
    ## 41    transcript::Porites_evermani_scaffold_1210:39391-42143
    ## 42    transcript::Porites_evermani_scaffold_1210:39391-42143
    ## 43  transcript::Porites_evermani_scaffold_1216:141298-141568
    ## 44    transcript::Porites_evermani_scaffold_1223:17191-17507
    ## 45    transcript::Porites_evermani_scaffold_1223:17191-17507
    ## 46     transcript::Porites_evermani_scaffold_1248:2544-18263
    ## 47    transcript::Porites_evermani_scaffold_125:91490-102059
    ## 48    transcript::Porites_evermani_scaffold_125:97286-102059
    ## 49    transcript::Porites_evermani_scaffold_125:97286-102059
    ## 50     transcript::Porites_evermani_scaffold_126:12826-13147
    ## 51    transcript::Porites_evermani_scaffold_1260:48211-49261
    ## 52  transcript::Porites_evermani_scaffold_1287:113596-114691
    ## 53  transcript::Porites_evermani_scaffold_1293:134814-135814
    ## 54    transcript::Porites_evermani_scaffold_1304:36120-39924
    ## 55    transcript::Porites_evermani_scaffold_1304:36120-39924
    ## 56    transcript::Porites_evermani_scaffold_1334:25998-34574
    ## 57  transcript::Porites_evermani_scaffold_1357:102815-104229
    ## 58  transcript::Porites_evermani_scaffold_1357:113383-115785
    ## 59   transcript::Porites_evermani_scaffold_136:421102-424046
    ## 60    transcript::Porites_evermani_scaffold_1410:20206-21088
    ## 61    transcript::Porites_evermani_scaffold_1410:20209-21088
    ## 62    transcript::Porites_evermani_scaffold_1429:32475-53408
    ## 63    transcript::Porites_evermani_scaffold_1429:46011-48479
    ## 64    transcript::Porites_evermani_scaffold_1442:76759-98681
    ## 65    transcript::Porites_evermani_scaffold_1442:77189-98681
    ## 66     transcript::Porites_evermani_scaffold_145:12752-17204
    ## 67     transcript::Porites_evermani_scaffold_145:12752-17204
    ## 68   transcript::Porites_evermani_scaffold_145:185297-190281
    ## 69      transcript::Porites_evermani_scaffold_1478:5635-5956
    ## 70     transcript::Porites_evermani_scaffold_148:76525-89252
    ## 71     transcript::Porites_evermani_scaffold_148:76525-89252
    ## 72     transcript::Porites_evermani_scaffold_148:76565-91947
    ## 73     transcript::Porites_evermani_scaffold_148:76565-91947
    ## 74     transcript::Porites_evermani_scaffold_148:76565-91947
    ## 75     transcript::Porites_evermani_scaffold_148:76565-91947
    ## 76     transcript::Porites_evermani_scaffold_148:92027-98639
    ## 77    transcript::Porites_evermani_scaffold_1501:27352-27790
    ## 78    transcript::Porites_evermani_scaffold_1501:28694-29132
    ## 79       transcript::Porites_evermani_scaffold_1513:140-3165
    ## 80    transcript::Porites_evermani_scaffold_1555:48165-50112
    ## 81    transcript::Porites_evermani_scaffold_1555:48165-50112
    ## 82    transcript::Porites_evermani_scaffold_1565:50666-51323
    ## 83    transcript::Porites_evermani_scaffold_1565:50666-51323
    ## 84       transcript::Porites_evermani_scaffold_1589:574-2877
    ## 85    transcript::Porites_evermani_scaffold_16:363080-391056
    ## 86    transcript::Porites_evermani_scaffold_1632:13812-17817
    ## 87    transcript::Porites_evermani_scaffold_1656:13542-15454
    ## 88   transcript::Porites_evermani_scaffold_1662:99010-103815
    ## 89    transcript::Porites_evermani_scaffold_1687:82104-85648
    ## 90    transcript::Porites_evermani_scaffold_1687:82104-85648
    ## 91    transcript::Porites_evermani_scaffold_1710:39465-42045
    ## 92   transcript::Porites_evermani_scaffold_1729:94708-105887
    ## 93    transcript::Porites_evermani_scaffold_1732:69845-79966
    ## 94    transcript::Porites_evermani_scaffold_1732:72129-79966
    ## 95    transcript::Porites_evermani_scaffold_1750:72607-87206
    ## 96    transcript::Porites_evermani_scaffold_1750:72607-87206
    ## 97     transcript::Porites_evermani_scaffold_1752:3816-18499
    ## 98     transcript::Porites_evermani_scaffold_1752:3816-18499
    ## 99     transcript::Porites_evermani_scaffold_1752:3828-26260
    ## 100    transcript::Porites_evermani_scaffold_1752:3828-26260
    ## 101    transcript::Porites_evermani_scaffold_1752:3828-52470
    ## 102    transcript::Porites_evermani_scaffold_1752:3828-52470
    ## 103   transcript::Porites_evermani_scaffold_1752:11122-18499
    ## 104   transcript::Porites_evermani_scaffold_1752:11122-18499
    ## 105   transcript::Porites_evermani_scaffold_1752:68629-75825
    ## 106   transcript::Porites_evermani_scaffold_1752:68629-75825
    ## 107  transcript::Porites_evermani_scaffold_1752:86211-100604
    ## 108  transcript::Porites_evermani_scaffold_1752:86211-100604
    ## 109  transcript::Porites_evermani_scaffold_1752:88696-103089
    ## 110  transcript::Porites_evermani_scaffold_1752:88696-103089
    ## 111   transcript::Porites_evermani_scaffold_1760:92340-97040
    ## 112   transcript::Porites_evermani_scaffold_1760:92340-97040
    ## 113   transcript::Porites_evermani_scaffold_1767:77997-78697
    ## 114  transcript::Porites_evermani_scaffold_1797:82322-103281
    ## 115   transcript::Porites_evermani_scaffold_1808:77443-84801
    ## 116   transcript::Porites_evermani_scaffold_1812:73190-84212
    ## 117   transcript::Porites_evermani_scaffold_1820:67344-68289
    ## 118   transcript::Porites_evermani_scaffold_1824:13748-14042
    ## 119  transcript::Porites_evermani_scaffold_184:326895-327198
    ## 120   transcript::Porites_evermani_scaffold_1860:47400-53630
    ## 121   transcript::Porites_evermani_scaffold_1860:47400-53630
    ## 122   transcript::Porites_evermani_scaffold_1890:38294-38620
    ## 123   transcript::Porites_evermani_scaffold_1999:33239-33569
    ## 124     transcript::Porites_evermani_scaffold_20:25200-25456
    ## 125     transcript::Porites_evermani_scaffold_20:25200-25456
    ## 126   transcript::Porites_evermani_scaffold_20:970835-972209
    ## 127   transcript::Porites_evermani_scaffold_2003:27353-28808
    ## 128  transcript::Porites_evermani_scaffold_206:175814-183681
    ## 129  transcript::Porites_evermani_scaffold_206:175814-183681
    ## 130   transcript::Porites_evermani_scaffold_2062:59875-90309
    ## 131   transcript::Porites_evermani_scaffold_21:839076-840328
    ## 132   transcript::Porites_evermani_scaffold_2113:30641-42751
    ## 133   transcript::Porites_evermani_scaffold_2113:30641-42751
    ## 134     transcript::Porites_evermani_scaffold_2119:6599-9631
    ## 135  transcript::Porites_evermani_scaffold_212:281783-282086
    ## 136   transcript::Porites_evermani_scaffold_2120:53722-55495
    ## 137    transcript::Porites_evermani_scaffold_2153:4231-29177
    ## 138    transcript::Porites_evermani_scaffold_2153:4231-29177
    ## 139    transcript::Porites_evermani_scaffold_2153:4265-27976
    ## 140    transcript::Porites_evermani_scaffold_2153:4265-27976
    ## 141   transcript::Porites_evermani_scaffold_2156:83498-86012
    ## 142   transcript::Porites_evermani_scaffold_2168:31980-37149
    ## 143   transcript::Porites_evermani_scaffold_2168:31980-37149
    ## 144   transcript::Porites_evermani_scaffold_2190:23736-31713
    ## 145   transcript::Porites_evermani_scaffold_2190:23736-31713
    ## 146   transcript::Porites_evermani_scaffold_22:587709-588435
    ## 147   transcript::Porites_evermani_scaffold_2239:54195-70401
    ## 148   transcript::Porites_evermani_scaffold_2239:63552-73656
    ## 149   transcript::Porites_evermani_scaffold_2270:70952-76996
    ## 150   transcript::Porites_evermani_scaffold_2270:73520-76996
    ## 151   transcript::Porites_evermani_scaffold_2283:33781-44642
    ## 152   transcript::Porites_evermani_scaffold_2317:39765-39992
    ## 153   transcript::Porites_evermani_scaffold_2317:39765-39992
    ## 154   transcript::Porites_evermani_scaffold_2327:62378-69825
    ## 155   transcript::Porites_evermani_scaffold_2335:62492-66862
    ## 156   transcript::Porites_evermani_scaffold_2335:62492-66862
    ## 157  transcript::Porites_evermani_scaffold_236:342625-343105
    ## 158   transcript::Porites_evermani_scaffold_2379:48679-52878
    ## 159   transcript::Porites_evermani_scaffold_2382:49533-52477
    ## 160    transcript::Porites_evermani_scaffold_2435:4070-15311
    ## 161    transcript::Porites_evermani_scaffold_2435:4070-15311
    ## 162     transcript::Porites_evermani_scaffold_2488:5213-6052
    ## 163  transcript::Porites_evermani_scaffold_252:209301-218430
    ## 164  transcript::Porites_evermani_scaffold_252:209301-218430
    ## 165  transcript::Porites_evermani_scaffold_252:213930-218430
    ## 166  transcript::Porites_evermani_scaffold_252:213930-218430
    ## 167   transcript::Porites_evermani_scaffold_2568:13571-18297
    ## 168   transcript::Porites_evermani_scaffold_2568:13571-18297
    ## 169   transcript::Porites_evermani_scaffold_2568:11405-20423
    ## 170   transcript::Porites_evermani_scaffold_2568:11405-20423
    ## 171   transcript::Porites_evermani_scaffold_2568:12242-18576
    ## 172   transcript::Porites_evermani_scaffold_2568:12242-18576
    ## 173  transcript::Porites_evermani_scaffold_257:235205-236309
    ## 174   transcript::Porites_evermani_scaffold_2591:16325-23286
    ## 175   transcript::Porites_evermani_scaffold_2593:25149-32107
    ## 176   transcript::Porites_evermani_scaffold_26:377983-392680
    ## 177   transcript::Porites_evermani_scaffold_26:378057-392680
    ## 178   transcript::Porites_evermani_scaffold_26:775167-778497
    ## 179   transcript::Porites_evermani_scaffold_26:775167-778497
    ## 180   transcript::Porites_evermani_scaffold_2602:14879-46627
    ## 181   transcript::Porites_evermani_scaffold_2602:14879-46627
    ## 182     transcript::Porites_evermani_scaffold_2631:1389-5334
    ## 183  transcript::Porites_evermani_scaffold_268:172494-188842
    ## 184  transcript::Porites_evermani_scaffold_272:316431-319943
    ## 185  transcript::Porites_evermani_scaffold_272:316431-319943
    ## 186   transcript::Porites_evermani_scaffold_2764:53473-62366
    ## 187   transcript::Porites_evermani_scaffold_2764:53473-62366
    ## 188  transcript::Porites_evermani_scaffold_277:293140-300568
    ## 189  transcript::Porites_evermani_scaffold_281:314069-317936
    ## 190  transcript::Porites_evermani_scaffold_281:314069-317936
    ## 191   transcript::Porites_evermani_scaffold_2896:21191-49288
    ## 192   transcript::Porites_evermani_scaffold_2896:21191-49288
    ## 193  transcript::Porites_evermani_scaffold_290:300981-301743
    ## 194        transcript::Porites_evermani_scaffold_2905:8-7590
    ## 195        transcript::Porites_evermani_scaffold_2905:8-7590
    ## 196   transcript::Porites_evermani_scaffold_3048:11325-16022
    ## 197   transcript::Porites_evermani_scaffold_3048:12156-20308
    ## 198   transcript::Porites_evermani_scaffold_3048:12156-20308
    ## 199   transcript::Porites_evermani_scaffold_3048:12734-29284
    ## 200   transcript::Porites_evermani_scaffold_3048:12953-29284
    ## 201   transcript::Porites_evermani_scaffold_3086:14248-17291
    ## 202    transcript::Porites_evermani_scaffold_311:17026-17374
    ## 203  transcript::Porites_evermani_scaffold_314:175295-204995
    ## 204  transcript::Porites_evermani_scaffold_314:175561-184949
    ## 205   transcript::Porites_evermani_scaffold_3152:38970-42202
    ## 206   transcript::Porites_evermani_scaffold_3152:38970-42202
    ## 207   transcript::Porites_evermani_scaffold_3152:38970-44378
    ## 208   transcript::Porites_evermani_scaffold_3152:38970-44378
    ## 209   transcript::Porites_evermani_scaffold_3152:38970-44378
    ## 210   transcript::Porites_evermani_scaffold_3152:38970-44378
    ## 211   transcript::Porites_evermani_scaffold_3152:39500-41967
    ## 212   transcript::Porites_evermani_scaffold_3152:39500-41967
    ## 213    transcript::Porites_evermani_scaffold_316:75149-91320
    ## 214    transcript::Porites_evermani_scaffold_316:75149-91356
    ## 215    transcript::Porites_evermani_scaffold_316:75149-91356
    ## 216    transcript::Porites_evermani_scaffold_316:75223-97019
    ## 217    transcript::Porites_evermani_scaffold_3160:3773-15980
    ## 218      transcript::Porites_evermani_scaffold_319:5554-5839
    ## 219   transcript::Porites_evermani_scaffold_3234:42922-44745
    ## 220  transcript::Porites_evermani_scaffold_327:132066-132634
    ## 221   transcript::Porites_evermani_scaffold_33:740925-750087
    ## 222   transcript::Porites_evermani_scaffold_33:740970-745664
    ## 223   transcript::Porites_evermani_scaffold_33:741157-744277
    ## 224     transcript::Porites_evermani_scaffold_3333:3563-5615
    ## 225  transcript::Porites_evermani_scaffold_334:150792-156024
    ## 226  transcript::Porites_evermani_scaffold_334:151379-156024
    ## 227   transcript::Porites_evermani_scaffold_3391:34834-45998
    ## 228   transcript::Porites_evermani_scaffold_3391:34834-45998
    ## 229   transcript::Porites_evermani_scaffold_3423:34033-34470
    ## 230   transcript::Porites_evermani_scaffold_3440:19793-20881
    ## 231    transcript::Porites_evermani_scaffold_345:49260-88511
    ## 232    transcript::Porites_evermani_scaffold_3489:8889-10033
    ## 233  transcript::Porites_evermani_scaffold_359:222469-223666
    ## 234  transcript::Porites_evermani_scaffold_359:222718-223410
    ## 235  transcript::Porites_evermani_scaffold_359:228982-229674
    ## 236  transcript::Porites_evermani_scaffold_359:245685-246649
    ## 237  transcript::Porites_evermani_scaffold_359:273157-274132
    ## 238      transcript::Porites_evermani_scaffold_364:7600-9612
    ## 239   transcript::Porites_evermani_scaffold_3647:25425-25821
    ## 240   transcript::Porites_evermani_scaffold_3647:25425-25821
    ## 241  transcript::Porites_evermani_scaffold_368:271584-272798
    ## 242    transcript::Porites_evermani_scaffold_373:30065-33959
    ## 243   transcript::Porites_evermani_scaffold_374:78782-107060
    ## 244   transcript::Porites_evermani_scaffold_374:89953-106967
    ## 245   transcript::Porites_evermani_scaffold_374:89953-106967
    ## 246   transcript::Porites_evermani_scaffold_374:91523-104913
    ## 247   transcript::Porites_evermani_scaffold_374:91523-104913
    ## 248     transcript::Porites_evermani_scaffold_3741:1103-4144
    ## 249   transcript::Porites_evermani_scaffold_38:547028-549268
    ## 250  transcript::Porites_evermani_scaffold_382:277890-283288
    ## 251   transcript::Porites_evermani_scaffold_3838:20314-36737
    ## 252   transcript::Porites_evermani_scaffold_3847:32813-33303
    ## 253   transcript::Porites_evermani_scaffold_3864:34838-35702
    ## 254   transcript::Porites_evermani_scaffold_39:167734-187936
    ## 255  transcript::Porites_evermani_scaffold_4:1045694-1050668
    ## 256  transcript::Porites_evermani_scaffold_4:1045694-1050668
    ## 257    transcript::Porites_evermani_scaffold_4:639186-642339
    ## 258     transcript::Porites_evermani_scaffold_4043:1544-3156
    ## 259    transcript::Porites_evermani_scaffold_4101:1896-12342
    ## 260    transcript::Porites_evermani_scaffold_4101:1896-12342
    ## 261  transcript::Porites_evermani_scaffold_411:137524-142355
    ## 262  transcript::Porites_evermani_scaffold_411:137524-142355
    ## 263  transcript::Porites_evermani_scaffold_412:146639-153893
    ## 264  transcript::Porites_evermani_scaffold_412:146639-153893
    ## 265    transcript::Porites_evermani_scaffold_4193:1635-30585
    ## 266   transcript::Porites_evermani_scaffold_4324:22984-23725
    ## 267    transcript::Porites_evermani_scaffold_433:61169-62558
    ## 268  transcript::Porites_evermani_scaffold_433:125352-128534
    ## 269  transcript::Porites_evermani_scaffold_433:126887-135958
    ## 270  transcript::Porites_evermani_scaffold_433:130234-135615
    ## 271  transcript::Porites_evermani_scaffold_433:153096-169835
    ## 272  transcript::Porites_evermani_scaffold_433:179394-181437
    ## 273  transcript::Porites_evermani_scaffold_433:225093-236913
    ## 274  transcript::Porites_evermani_scaffold_433:234332-236913
    ## 275   transcript::Porites_evermani_scaffold_4333:20267-20999
    ## 276  transcript::Porites_evermani_scaffold_434:105607-108878
    ## 277  transcript::Porites_evermani_scaffold_434:190592-190913
    ## 278   transcript::Porites_evermani_scaffold_4371:17765-21771
    ## 279   transcript::Porites_evermani_scaffold_4371:17765-21771
    ## 280   transcript::Porites_evermani_scaffold_4387:17082-22126
    ## 281        transcript::Porites_evermani_scaffold_4431:0-1044
    ## 282  transcript::Porites_evermani_scaffold_450:235626-241004
    ## 283  transcript::Porites_evermani_scaffold_452:101945-106670
    ## 284  transcript::Porites_evermani_scaffold_452:101950-146632
    ## 285        transcript::Porites_evermani_scaffold_4581:0-6394
    ## 286       transcript::Porites_evermani_scaffold_4581:14-7951
    ## 287       transcript::Porites_evermani_scaffold_4581:14-7951
    ## 288      transcript::Porites_evermani_scaffold_459:3715-4628
    ## 289  transcript::Porites_evermani_scaffold_461:195334-219192
    ## 290  transcript::Porites_evermani_scaffold_461:195351-223079
    ## 291  transcript::Porites_evermani_scaffold_461:210393-223079
    ## 292  transcript::Porites_evermani_scaffold_469:174357-176910
    ## 293   transcript::Porites_evermani_scaffold_47:472593-482329
    ## 294    transcript::Porites_evermani_scaffold_4714:1999-16191
    ## 295    transcript::Porites_evermani_scaffold_4714:1999-16191
    ## 296   transcript::Porites_evermani_scaffold_4809:14371-21317
    ## 297   transcript::Porites_evermani_scaffold_4809:16108-21317
    ## 298    transcript::Porites_evermani_scaffold_484:45687-52014
    ## 299    transcript::Porites_evermani_scaffold_484:45687-52014
    ## 300   transcript::Porites_evermani_scaffold_49:147530-165944
    ## 301   transcript::Porites_evermani_scaffold_49:531624-540967
    ## 302   transcript::Porites_evermani_scaffold_49:531624-540967
    ## 303     transcript::Porites_evermani_scaffold_5047:1877-8627
    ## 304    transcript::Porites_evermani_scaffold_509:33027-40843
    ## 305  transcript::Porites_evermani_scaffold_511:158791-178303
    ## 306  transcript::Porites_evermani_scaffold_511:158791-178303
    ## 307     transcript::Porites_evermani_scaffold_5122:393-16149
    ## 308     transcript::Porites_evermani_scaffold_5122:393-16149
    ## 309    transcript::Porites_evermani_scaffold_5122:7344-16149
    ## 310    transcript::Porites_evermani_scaffold_5122:7344-16149
    ## 311    transcript::Porites_evermani_scaffold_521:69491-72481
    ## 312    transcript::Porites_evermani_scaffold_521:69491-72481
    ## 313    transcript::Porites_evermani_scaffold_521:69992-72998
    ## 314   transcript::Porites_evermani_scaffold_5281:13425-15228
    ## 315  transcript::Porites_evermani_scaffold_532:146041-155315
    ## 316  transcript::Porites_evermani_scaffold_532:146154-155315
    ## 317  transcript::Porites_evermani_scaffold_532:100926-106341
    ## 318  transcript::Porites_evermani_scaffold_532:121147-123619
    ## 319  transcript::Porites_evermani_scaffold_538:195523-197740
    ## 320  transcript::Porites_evermani_scaffold_552:206218-215403
    ## 321  transcript::Porites_evermani_scaffold_552:206218-215403
    ## 322  transcript::Porites_evermani_scaffold_558:184879-227684
    ## 323  transcript::Porites_evermani_scaffold_558:184879-227684
    ## 324  transcript::Porites_evermani_scaffold_558:184879-227684
    ## 325  transcript::Porites_evermani_scaffold_558:184879-227684
    ## 326  transcript::Porites_evermani_scaffold_558:184898-227684
    ## 327  transcript::Porites_evermani_scaffold_558:184898-227684
    ## 328  transcript::Porites_evermani_scaffold_558:212607-227684
    ## 329  transcript::Porites_evermani_scaffold_558:212607-227684
    ## 330  transcript::Porites_evermani_scaffold_563:166346-175532
    ## 331  transcript::Porites_evermani_scaffold_563:166346-175532
    ## 332     transcript::Porites_evermani_scaffold_5647:5015-7597
    ## 333    transcript::Porites_evermani_scaffold_568:63403-93826
    ## 334    transcript::Porites_evermani_scaffold_568:63403-93826
    ## 335     transcript::Porites_evermani_scaffold_5697:252-10540
    ## 336      transcript::Porites_evermani_scaffold_5860:456-7132
    ## 337      transcript::Porites_evermani_scaffold_5860:456-7132
    ## 338     transcript::Porites_evermani_scaffold_5873:1061-1787
    ## 339     transcript::Porites_evermani_scaffold_6049:6069-8167
    ## 340     transcript::Porites_evermani_scaffold_6049:6143-8167
    ## 341    transcript::Porites_evermani_scaffold_609:22352-36231
    ## 342    transcript::Porites_evermani_scaffold_609:22358-36231
    ## 343  transcript::Porites_evermani_scaffold_613:139660-180006
    ## 344  transcript::Porites_evermani_scaffold_613:139660-180006
    ## 345  transcript::Porites_evermani_scaffold_613:139663-179070
    ## 346  transcript::Porites_evermani_scaffold_613:139663-180006
    ## 347  transcript::Porites_evermani_scaffold_613:139962-180006
    ## 348   transcript::Porites_evermani_scaffold_614:99032-101696
    ## 349   transcript::Porites_evermani_scaffold_614:99032-101696
    ## 350    transcript::Porites_evermani_scaffold_617:22652-25586
    ## 351  transcript::Porites_evermani_scaffold_632:130786-157955
    ## 352  transcript::Porites_evermani_scaffold_632:131313-140033
    ## 353  transcript::Porites_evermani_scaffold_632:173167-176852
    ## 354  transcript::Porites_evermani_scaffold_634:183915-188962
    ## 355     transcript::Porites_evermani_scaffold_6377:3294-3554
    ## 356       transcript::Porites_evermani_scaffold_640:522-4205
    ## 357    transcript::Porites_evermani_scaffold_644:44906-45221
    ## 358       transcript::Porites_evermani_scaffold_6458:17-2752
    ## 359       transcript::Porites_evermani_scaffold_6458:17-2752
    ## 360  transcript::Porites_evermani_scaffold_665:161804-162500
    ## 361       transcript::Porites_evermani_scaffold_6764:39-5101
    ## 362  transcript::Porites_evermani_scaffold_689:192514-200765
    ## 363     transcript::Porites_evermani_scaffold_7019:1762-3756
    ## 364   transcript::Porites_evermani_scaffold_71:580623-586006
    ## 365   transcript::Porites_evermani_scaffold_71:587785-592555
    ## 366    transcript::Porites_evermani_scaffold_710:46565-66234
    ## 367    transcript::Porites_evermani_scaffold_710:56999-58400
    ## 368        transcript::Porites_evermani_scaffold_7246:26-695
    ## 369    transcript::Porites_evermani_scaffold_729:41169-42951
    ## 370   transcript::Porites_evermani_scaffold_74:403052-408411
    ## 371    transcript::Porites_evermani_scaffold_763:21685-21993
    ## 372    transcript::Porites_evermani_scaffold_763:21685-21993
    ## 373    transcript::Porites_evermani_scaffold_763:87341-88332
    ## 374    transcript::Porites_evermani_scaffold_763:87341-88332
    ## 375  transcript::Porites_evermani_scaffold_768:164819-165869
    ## 376  transcript::Porites_evermani_scaffold_768:164876-165869
    ## 377   transcript::Porites_evermani_scaffold_79:199872-210892
    ## 378   transcript::Porites_evermani_scaffold_79:199927-210892
    ## 379   transcript::Porites_evermani_scaffold_79:282399-282762
    ## 380    transcript::Porites_evermani_scaffold_790:66155-76048
    ## 381    transcript::Porites_evermani_scaffold_790:66155-76048
    ## 382    transcript::Porites_evermani_scaffold_790:68425-76048
    ## 383    transcript::Porites_evermani_scaffold_790:68425-76048
    ## 384  transcript::Porites_evermani_scaffold_797:135670-135990
    ## 385   transcript::Porites_evermani_scaffold_81:526397-527406
    ## 386  transcript::Porites_evermani_scaffold_833:146401-165613
    ## 387  transcript::Porites_evermani_scaffold_833:146490-166167
    ## 388  transcript::Porites_evermani_scaffold_833:146490-166167
    ## 389  transcript::Porites_evermani_scaffold_833:146491-166167
    ## 390    transcript::Porites_evermani_scaffold_852:58686-80235
    ## 391    transcript::Porites_evermani_scaffold_866:21937-24095
    ## 392   transcript::Porites_evermani_scaffold_90:426819-436482
    ## 393   transcript::Porites_evermani_scaffold_90:426819-436482
    ## 394   transcript::Porites_evermani_scaffold_90:426819-439309
    ## 395   transcript::Porites_evermani_scaffold_90:426819-439309
    ## 396    transcript::Porites_evermani_scaffold_903:57641-77421
    ## 397    transcript::Porites_evermani_scaffold_903:57641-77421
    ## 398    transcript::Porites_evermani_scaffold_903:57641-77421
    ## 399    transcript::Porites_evermani_scaffold_903:57641-77421
    ## 400   transcript::Porites_evermani_scaffold_910:95262-102133
    ## 401   transcript::Porites_evermani_scaffold_910:95262-102133
    ## 402   transcript::Porites_evermani_scaffold_910:95262-102133
    ## 403   transcript::Porites_evermani_scaffold_910:95900-102133
    ## 404   transcript::Porites_evermani_scaffold_910:95900-102133
    ## 405   transcript::Porites_evermani_scaffold_910:95900-102133
    ## 406  transcript::Porites_evermani_scaffold_910:108432-120843
    ## 407  transcript::Porites_evermani_scaffold_910:108432-120843
    ## 408  transcript::Porites_evermani_scaffold_910:108432-120843
    ## 409  transcript::Porites_evermani_scaffold_910:129984-141109
    ## 410  transcript::Porites_evermani_scaffold_910:129984-141109
    ## 411  transcript::Porites_evermani_scaffold_910:129984-141109
    ## 412  transcript::Porites_evermani_scaffold_910:134861-141109
    ## 413  transcript::Porites_evermani_scaffold_910:134861-141109
    ## 414  transcript::Porites_evermani_scaffold_910:134861-141109
    ## 415  transcript::Porites_evermani_scaffold_926:111868-124526
    ## 416  transcript::Porites_evermani_scaffold_926:111868-132313
    ## 417  transcript::Porites_evermani_scaffold_935:170191-170800
    ## 418  transcript::Porites_evermani_scaffold_935:170191-170800
    ## 419  transcript::Porites_evermani_scaffold_942:121381-138291
    ## 420    transcript::Porites_evermani_scaffold_984:51531-58363
    ## 421    transcript::Porites_evermani_scaffold_984:51539-59053
    ##                                                                   sseqid
    ## 1    Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 2        Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 3    Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 4        Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 5    Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 6      Cluster_29.mature::Porites_evermani_scaffold_1:1404272-1404293(-)
    ## 7      Cluster_29.mature::Porites_evermani_scaffold_1:1404272-1404293(-)
    ## 8    Cluster_5882.mature::Porites_evermani_scaffold_461:215508-215529(+)
    ## 9      Cluster_15890.mature::Porites_evermani_scaffold_3893:3849-3870(-)
    ## 10     Cluster_15890.mature::Porites_evermani_scaffold_3893:3849-3870(-)
    ## 11     Cluster_29.mature::Porites_evermani_scaffold_1:1404272-1404293(-)
    ## 12   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 13       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 14   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 15       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 16   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 17       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 18   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 19   Cluster_4115.mature::Porites_evermani_scaffold_257:110361-110382(-)
    ## 20  Cluster_9983.mature::Porites_evermani_scaffold_1128:125548-125568(-)
    ## 21  Cluster_9983.mature::Porites_evermani_scaffold_1128:125548-125568(-)
    ## 22   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 23   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 24   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 25   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 26   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 27   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 28     Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 29     Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)
    ## 30   Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 31   Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 32   Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 33   Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 34   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 35   Cluster_7053.mature::Porites_evermani_scaffold_613:156505-156526(+)
    ## 36     Cluster_29.mature::Porites_evermani_scaffold_1:1404272-1404293(-)
    ## 37     Cluster_29.mature::Porites_evermani_scaffold_1:1404272-1404293(-)
    ## 38   Cluster_10934.mature::Porites_evermani_scaffold_1415:72071-72092(+)
    ## 39   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 40       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 41   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 42       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 43   Cluster_8988.mature::Porites_evermani_scaffold_942:133698-133719(+)
    ## 44   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 45       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 46   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 47   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 48   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 49   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 50   Cluster_4115.mature::Porites_evermani_scaffold_257:110361-110382(-)
    ## 51   Cluster_14500.mature::Porites_evermani_scaffold_2738:56885-56906(+)
    ## 52     Cluster_15890.mature::Porites_evermani_scaffold_3893:3849-3870(-)
    ## 53   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 54   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 55       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 56   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 57   Cluster_6914.mature::Porites_evermani_scaffold_594:158230-158250(+)
    ## 58   Cluster_6914.mature::Porites_evermani_scaffold_594:158230-158250(+)
    ## 59     Cluster_796.mature::Porites_evermani_scaffold_26:382572-382593(-)
    ## 60     Cluster_15890.mature::Porites_evermani_scaffold_3893:3849-3870(-)
    ## 61     Cluster_15890.mature::Porites_evermani_scaffold_3893:3849-3870(-)
    ## 62   Cluster_10965.mature::Porites_evermani_scaffold_1429:47285-47306(-)
    ## 63   Cluster_10965.mature::Porites_evermani_scaffold_1429:47285-47306(-)
    ## 64   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 65   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 66   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 67       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 68       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 69   Cluster_4115.mature::Porites_evermani_scaffold_257:110361-110382(-)
    ## 70   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 71       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 72   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 73       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 74   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 75       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 76     Cluster_2882.mature::Porites_evermani_scaffold_148:94583-94604(-)
    ## 77   Cluster_10934.mature::Porites_evermani_scaffold_1415:72071-72092(+)
    ## 78   Cluster_10934.mature::Porites_evermani_scaffold_1415:72071-72092(+)
    ## 79     Cluster_589.mature::Porites_evermani_scaffold_16:383437-383458(-)
    ## 80     Cluster_7658.mature::Porites_evermani_scaffold_730:82423-82444(-)
    ## 81     Cluster_7657.mature::Porites_evermani_scaffold_730:81385-81406(-)
    ## 82     Cluster_7658.mature::Porites_evermani_scaffold_730:82423-82444(-)
    ## 83     Cluster_7657.mature::Porites_evermani_scaffold_730:81385-81406(-)
    ## 84   Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 85     Cluster_589.mature::Porites_evermani_scaffold_16:383437-383458(-)
    ## 86   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 87       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 88   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 89     Cluster_7658.mature::Porites_evermani_scaffold_730:82423-82444(-)
    ## 90     Cluster_7657.mature::Porites_evermani_scaffold_730:81385-81406(-)
    ## 91   Cluster_2787.mature::Porites_evermani_scaffold_138:127966-127987(+)
    ## 92   Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 93   Cluster_11997.mature::Porites_evermani_scaffold_1732:76504-76525(+)
    ## 94   Cluster_11997.mature::Porites_evermani_scaffold_1732:76504-76525(+)
    ## 95   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 96       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 97   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 98       Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 99   Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 100      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 101  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 102      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 103  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 104      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 105  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 106      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 107  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 108      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 109  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 110      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 111  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 112      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 113    Cluster_29.mature::Porites_evermani_scaffold_1:1404272-1404293(-)
    ## 114  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 115  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 116  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 117    Cluster_29.mature::Porites_evermani_scaffold_1:1404272-1404293(-)
    ## 118      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 119  Cluster_16498.mature::Porites_evermani_scaffold_5010:12392-12413(+)
    ## 120  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 121      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 122  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 123  Cluster_6914.mature::Porites_evermani_scaffold_594:158230-158250(+)
    ## 124  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 125      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 126  Cluster_5882.mature::Porites_evermani_scaffold_461:215508-215529(+)
    ## 127  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 128  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 129      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 130  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 131  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 132  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 133      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 134    Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 135  Cluster_16498.mature::Porites_evermani_scaffold_5010:12392-12413(+)
    ## 136    Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)
    ## 137  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 138      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 139  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 140      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 141  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 142  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 143      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 144  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 145      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 146    Cluster_29.mature::Porites_evermani_scaffold_1:1404272-1404293(-)
    ## 147      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 148      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 149    Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 150    Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 151      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 152    Cluster_7658.mature::Porites_evermani_scaffold_730:82423-82444(-)
    ## 153    Cluster_7657.mature::Porites_evermani_scaffold_730:81385-81406(-)
    ## 154  Cluster_6914.mature::Porites_evermani_scaffold_594:158230-158250(+)
    ## 155  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 156      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 157  Cluster_10965.mature::Porites_evermani_scaffold_1429:47285-47306(-)
    ## 158      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 159    Cluster_796.mature::Porites_evermani_scaffold_26:382572-382593(-)
    ## 160  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 161      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 162  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 163  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 164      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 165  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 166      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 167    Cluster_7658.mature::Porites_evermani_scaffold_730:82423-82444(-)
    ## 168    Cluster_7657.mature::Porites_evermani_scaffold_730:81385-81406(-)
    ## 169    Cluster_7658.mature::Porites_evermani_scaffold_730:82423-82444(-)
    ## 170    Cluster_7657.mature::Porites_evermani_scaffold_730:81385-81406(-)
    ## 171    Cluster_7658.mature::Porites_evermani_scaffold_730:82423-82444(-)
    ## 172    Cluster_7657.mature::Porites_evermani_scaffold_730:81385-81406(-)
    ## 173    Cluster_796.mature::Porites_evermani_scaffold_26:382572-382593(-)
    ## 174  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 175  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 176    Cluster_796.mature::Porites_evermani_scaffold_26:382572-382593(-)
    ## 177    Cluster_796.mature::Porites_evermani_scaffold_26:382572-382593(-)
    ## 178  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 179      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 180  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 181      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 182  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 183  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 184      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 185  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 186      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 187  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 188  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 189      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 190  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 191  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 192      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 193    Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)
    ## 194      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 195  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 196  Cluster_11997.mature::Porites_evermani_scaffold_1732:76504-76525(+)
    ## 197  Cluster_11997.mature::Porites_evermani_scaffold_1732:76504-76525(+)
    ## 198  Cluster_11997.mature::Porites_evermani_scaffold_1732:76504-76525(+)
    ## 199  Cluster_11997.mature::Porites_evermani_scaffold_1732:76504-76525(+)
    ## 200  Cluster_11997.mature::Porites_evermani_scaffold_1732:76504-76525(+)
    ## 201  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 202  Cluster_16498.mature::Porites_evermani_scaffold_5010:12392-12413(+)
    ## 203  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 204  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 205      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 206  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 207      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 208  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 209      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 210  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 211      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 212  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 213    Cluster_4629.mature::Porites_evermani_scaffold_316:88465-88486(-)
    ## 214    Cluster_4629.mature::Porites_evermani_scaffold_316:88465-88486(-)
    ## 215    Cluster_4629.mature::Porites_evermani_scaffold_316:88465-88486(-)
    ## 216    Cluster_4629.mature::Porites_evermani_scaffold_316:88465-88486(-)
    ## 217  Cluster_2787.mature::Porites_evermani_scaffold_138:127966-127987(+)
    ## 218  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 219  Cluster_6914.mature::Porites_evermani_scaffold_594:158230-158250(+)
    ## 220  Cluster_16498.mature::Porites_evermani_scaffold_5010:12392-12413(+)
    ## 221  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 222  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 223  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 224  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 225  Cluster_4735.mature::Porites_evermani_scaffold_334:153605-153626(-)
    ## 226  Cluster_4735.mature::Porites_evermani_scaffold_334:153605-153626(-)
    ## 227  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 228      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 229    Cluster_589.mature::Porites_evermani_scaffold_16:383437-383458(-)
    ## 230  Cluster_6914.mature::Porites_evermani_scaffold_594:158230-158250(+)
    ## 231      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 232    Cluster_15890.mature::Porites_evermani_scaffold_3893:3849-3870(-)
    ## 233    Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 234    Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 235    Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 236    Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 237    Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 238  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 239  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 240      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 241    Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 242  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 243  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 244  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 245      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 246  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 247      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 248  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 249  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 250      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 251  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 252  Cluster_2787.mature::Porites_evermani_scaffold_138:127966-127987(+)
    ## 253  Cluster_14500.mature::Porites_evermani_scaffold_2738:56885-56906(+)
    ## 254  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 255      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 256  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 257  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 258    Cluster_589.mature::Porites_evermani_scaffold_16:383437-383458(-)
    ## 259  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 260      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 261      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 262  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 263  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 264      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 265  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 266  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 267  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 268  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 269  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 270  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 271  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 272  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 273  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 274  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 275  Cluster_6914.mature::Porites_evermani_scaffold_594:158230-158250(+)
    ## 276    Cluster_16738.mature::Porites_evermani_scaffold_6219:6549-6570(-)
    ## 277    Cluster_16738.mature::Porites_evermani_scaffold_6219:6549-6570(-)
    ## 278  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 279      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 280  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 281    Cluster_6255.mature::Porites_evermani_scaffold_502:58997-59018(-)
    ## 282    Cluster_589.mature::Porites_evermani_scaffold_16:383437-383458(-)
    ## 283  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 284  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 285  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 286  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 287  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 288  Cluster_5882.mature::Porites_evermani_scaffold_461:215508-215529(+)
    ## 289  Cluster_5882.mature::Porites_evermani_scaffold_461:215508-215529(+)
    ## 290  Cluster_5882.mature::Porites_evermani_scaffold_461:215508-215529(+)
    ## 291  Cluster_5882.mature::Porites_evermani_scaffold_461:215508-215529(+)
    ## 292  Cluster_4115.mature::Porites_evermani_scaffold_257:110361-110382(-)
    ## 293   Cluster_1140.mature::Porites_evermani_scaffold_47:475994-476015(-)
    ## 294  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 295      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 296  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 297  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 298  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 299      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 300   Cluster_1167.mature::Porites_evermani_scaffold_49:151640-151661(-)
    ## 301  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 302      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 303  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 304  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 305  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 306      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 307  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 308      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 309  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 310      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 311  Cluster_16498.mature::Porites_evermani_scaffold_5010:12392-12413(+)
    ## 312  Cluster_16498.mature::Porites_evermani_scaffold_5010:12392-12413(+)
    ## 313  Cluster_16498.mature::Porites_evermani_scaffold_5010:12392-12413(+)
    ## 314    Cluster_589.mature::Porites_evermani_scaffold_16:383437-383458(-)
    ## 315  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 316  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 317  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 318  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 319  Cluster_6914.mature::Porites_evermani_scaffold_594:158230-158250(+)
    ## 320  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 321      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 322  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 323      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 324  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 325      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 326  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 327      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 328  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 329      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 330  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 331      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 332    Cluster_9149.mature::Porites_evermani_scaffold_984:51883-51904(-)
    ## 333  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 334      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 335  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 336  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 337      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 338  Cluster_10965.mature::Porites_evermani_scaffold_1429:47285-47306(-)
    ## 339    Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 340    Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)
    ## 341  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 342  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 343  Cluster_7053.mature::Porites_evermani_scaffold_613:156505-156526(+)
    ## 344  Cluster_7053.mature::Porites_evermani_scaffold_613:156505-156526(+)
    ## 345  Cluster_7053.mature::Porites_evermani_scaffold_613:156505-156526(+)
    ## 346  Cluster_7053.mature::Porites_evermani_scaffold_613:156505-156526(+)
    ## 347  Cluster_7053.mature::Porites_evermani_scaffold_613:156505-156526(+)
    ## 348  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 349      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 350  Cluster_6914.mature::Porites_evermani_scaffold_594:158230-158250(+)
    ## 351  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 352  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 353  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 354    Cluster_589.mature::Porites_evermani_scaffold_16:383437-383458(-)
    ## 355  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 356  Cluster_6914.mature::Porites_evermani_scaffold_594:158230-158250(+)
    ## 357  Cluster_4115.mature::Porites_evermani_scaffold_257:110361-110382(-)
    ## 358      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 359  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 360  Cluster_10965.mature::Porites_evermani_scaffold_1429:47285-47306(-)
    ## 361  Cluster_15726.mature::Porites_evermani_scaffold_3707:36124-36145(-)
    ## 362  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 363    Cluster_29.mature::Porites_evermani_scaffold_1:1404272-1404293(-)
    ## 364  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 365  Cluster_7855.mature::Porites_evermani_scaffold_768:138025-138046(+)
    ## 366    Cluster_16738.mature::Porites_evermani_scaffold_6219:6549-6570(-)
    ## 367    Cluster_16738.mature::Porites_evermani_scaffold_6219:6549-6570(-)
    ## 368  Cluster_4115.mature::Porites_evermani_scaffold_257:110361-110382(-)
    ## 369  Cluster_11997.mature::Porites_evermani_scaffold_1732:76504-76525(+)
    ## 370    Cluster_589.mature::Porites_evermani_scaffold_16:383437-383458(-)
    ## 371    Cluster_7658.mature::Porites_evermani_scaffold_730:82423-82444(-)
    ## 372    Cluster_7657.mature::Porites_evermani_scaffold_730:81385-81406(-)
    ## 373    Cluster_7658.mature::Porites_evermani_scaffold_730:82423-82444(-)
    ## 374    Cluster_7657.mature::Porites_evermani_scaffold_730:81385-81406(-)
    ## 375    Cluster_15890.mature::Porites_evermani_scaffold_3893:3849-3870(-)
    ## 376    Cluster_15890.mature::Porites_evermani_scaffold_3893:3849-3870(-)
    ## 377    Cluster_16738.mature::Porites_evermani_scaffold_6219:6549-6570(-)
    ## 378    Cluster_16738.mature::Porites_evermani_scaffold_6219:6549-6570(-)
    ## 379    Cluster_29.mature::Porites_evermani_scaffold_1:1404272-1404293(-)
    ## 380  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 381      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 382  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 383      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 384    Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)
    ## 385   Cluster_1140.mature::Porites_evermani_scaffold_47:475994-476015(-)
    ## 386      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 387      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 388      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 389      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 390  Cluster_2787.mature::Porites_evermani_scaffold_138:127966-127987(+)
    ## 391    Cluster_8634.mature::Porites_evermani_scaffold_866:22854-22875(-)
    ## 392  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 393      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 394  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 395      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 396      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 397  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 398      Cluster_2854.mature::Porites_evermani_scaffold_145:5406-5427(+)
    ## 399  Cluster_13502.mature::Porites_evermani_scaffold_2259:40244-40265(+)
    ## 400  Cluster_8887.mature::Porites_evermani_scaffold_910:118742-118762(+)
    ## 401    Cluster_8884.mature::Porites_evermani_scaffold_910:99255-99275(+)
    ## 402  Cluster_8888.mature::Porites_evermani_scaffold_910:139353-139373(+)
    ## 403  Cluster_8887.mature::Porites_evermani_scaffold_910:118742-118762(+)
    ## 404    Cluster_8884.mature::Porites_evermani_scaffold_910:99255-99275(+)
    ## 405  Cluster_8888.mature::Porites_evermani_scaffold_910:139353-139373(+)
    ## 406  Cluster_8888.mature::Porites_evermani_scaffold_910:139353-139373(+)
    ## 407  Cluster_8887.mature::Porites_evermani_scaffold_910:118742-118762(+)
    ## 408    Cluster_8884.mature::Porites_evermani_scaffold_910:99255-99275(+)
    ## 409  Cluster_8888.mature::Porites_evermani_scaffold_910:139353-139373(+)
    ## 410  Cluster_8887.mature::Porites_evermani_scaffold_910:118742-118762(+)
    ## 411    Cluster_8884.mature::Porites_evermani_scaffold_910:99255-99275(+)
    ## 412  Cluster_8888.mature::Porites_evermani_scaffold_910:139353-139373(+)
    ## 413  Cluster_8887.mature::Porites_evermani_scaffold_910:118742-118762(+)
    ## 414    Cluster_8884.mature::Porites_evermani_scaffold_910:99255-99275(+)
    ## 415  Cluster_10965.mature::Porites_evermani_scaffold_1429:47285-47306(-)
    ## 416  Cluster_10965.mature::Porites_evermani_scaffold_1429:47285-47306(-)
    ## 417    Cluster_7658.mature::Porites_evermani_scaffold_730:82423-82444(-)
    ## 418    Cluster_7657.mature::Porites_evermani_scaffold_730:81385-81406(-)
    ## 419  Cluster_8988.mature::Porites_evermani_scaffold_942:133698-133719(+)
    ## 420    Cluster_9149.mature::Porites_evermani_scaffold_984:51883-51904(-)
    ## 421    Cluster_9149.mature::Porites_evermani_scaffold_984:51883-51904(-)
    ##      pident length mismatch gapopen qstart  qend sstart send   evalue bitscore
    ## 1   100.000     22        0       0  32822 32843     22    1 5.49e-06     41.0
    ## 2   100.000     22        0       0  32851 32872     22    1 5.49e-06     41.0
    ## 3   100.000     22        0       0   9774  9795      1   22 3.80e-06     41.0
    ## 4   100.000     22        0       0   9745  9766      1   22 3.80e-06     41.0
    ## 5    90.909     22        2       0   2953  2974      1   22 4.34e-04     31.9
    ## 6   100.000     22        0       0  23860 23881     22    1 3.78e-06     41.0
    ## 7   100.000     22        0       0  23524 23545     22    1 3.75e-06     41.0
    ## 8   100.000     14        0       0    130   143      9   22 1.00e-03     26.5
    ## 9   100.000     15        0       0    311   325      2   16 1.00e-03     28.3
    ## 10  100.000     15        0       0    165   179      2   16 1.00e-03     28.3
    ## 11  100.000     16        0       0    129   144      6   21 1.36e-04     30.1
    ## 12   95.455     22        1       0  37594 37615      1   22 1.68e-04     36.5
    ## 13   95.455     22        1       0  37592 37613     22    1 1.68e-04     36.5
    ## 14   95.455     22        1       0  36425 36446      1   22 1.63e-04     36.5
    ## 15   95.455     22        1       0  36423 36444     22    1 1.63e-04     36.5
    ## 16  100.000     22        0       0   5379  5400     22    1 1.74e-06     41.0
    ## 17  100.000     22        0       0   5381  5402      1   22 1.74e-06     41.0
    ## 18  100.000     14        0       0     81    94      1   14 1.00e-03     26.5
    ## 19  100.000     16        0       0     26    41     16    1 1.42e-04     30.1
    ## 20  100.000     21        0       0  16312 16332     21    1 7.82e-06     39.2
    ## 21  100.000     21        0       0  16148 16168     21    1 7.75e-06     39.2
    ## 22   95.455     22        1       0    381   402     22    1 6.85e-05     36.5
    ## 23   95.455     22        1       0    381   402     22    1 6.85e-05     36.5
    ## 24   95.455     22        1       0    381   402     22    1 6.85e-05     36.5
    ## 25   95.455     22        1       0    378   399     22    1 6.85e-05     36.5
    ## 26   95.455     22        1       0    378   399     22    1 6.85e-05     36.5
    ## 27   95.455     22        1       0    333   354     22    1 6.82e-05     36.5
    ## 28  100.000     22        0       0   3722  3743      1   22 1.54e-06     41.0
    ## 29  100.000     22        0       0   2996  3017      1   22 1.54e-06     41.0
    ## 30  100.000     19        0       0  11294 11312      4   22 1.31e-04     35.6
    ## 31  100.000     19        0       0  11272 11290      4   22 1.31e-04     35.6
    ## 32  100.000     19        0       0   4704  4722      4   22 7.06e-05     35.6
    ## 33  100.000     19        0       0    704   722      4   22 4.78e-05     35.6
    ## 34   95.455     22        1       0   6400  6421      1   22 6.25e-05     36.5
    ## 35  100.000     16        0       0     43    58     22    7 2.71e-04     30.1
    ## 36  100.000     16        0       0    782   797      4   19 1.00e-03     30.1
    ## 37  100.000     16        0       0   3146  3161      4   19 1.00e-03     30.1
    ## 38  100.000     15        0       0    682   696     19    5 1.00e-03     28.3
    ## 39   95.455     22        1       0  10047 10068     22    1 6.94e-05     36.5
    ## 40   95.455     22        1       0  10049 10070      1   22 6.94e-05     36.5
    ## 41   95.455     22        1       0   1552  1573     22    1 1.94e-05     36.5
    ## 42   95.455     22        1       0   1554  1575      1   22 1.94e-05     36.5
    ## 43   89.474     19        2       0    102   120      4   22 1.00e-03     26.5
    ## 44   90.909     22        2       0    177   198      1   22 3.63e-05     31.9
    ## 45   90.909     22        2       0    175   196     22    1 3.63e-05     31.9
    ## 46   90.909     22        2       0   1621  1642     22    1 1.00e-03     31.9
    ## 47   95.455     22        1       0   6711  6732      1   22 6.52e-05     36.5
    ## 48   95.455     22        1       0    915   936      1   22 2.94e-05     36.5
    ## 49   95.455     22        1       0    915   936      1   22 2.94e-05     36.5
    ## 50  100.000     14        0       0    258   271      2   15 1.00e-03     26.5
    ## 51   94.444     18        1       0    717   734     18    1 1.00e-03     29.2
    ## 52   90.000     20        2       0   1005  1024      2   21 1.00e-03     28.3
    ## 53  100.000     15        0       0    796   810      8   22 1.00e-03     28.3
    ## 54   90.909     22        2       0   1465  1486      1   22 3.27e-04     31.9
    ## 55   90.909     22        2       0   1463  1484     22    1 3.27e-04     31.9
    ## 56   90.909     22        2       0   1519  1540     22    1 6.45e-04     31.9
    ## 57   94.737     19        1       0    231   249      1   19 4.75e-04     31.0
    ## 58   94.737     19        1       0    231   249      1   19 7.18e-04     31.0
    ## 59  100.000     16        0       0    187   202      7   22 8.81e-04     30.1
    ## 60   94.444     18        1       0     48    65     21    4 1.00e-03     29.2
    ## 61   94.444     18        1       0     45    62     21    4 1.00e-03     29.2
    ## 62  100.000     22        0       0  14810 14831     22    1 2.60e-06     41.0
    ## 63  100.000     22        0       0   1274  1295     22    1 4.08e-07     41.0
    ## 64   95.238     21        1       0   3417  3437      1   21 4.04e-04     34.6
    ## 65   95.238     21        1       0   2987  3007      1   21 3.96e-04     34.6
    ## 66  100.000     22        0       0   1469  1490      1   22 7.38e-07     41.0
    ## 67  100.000     22        0       0   1440  1461      1   22 7.38e-07     41.0
    ## 68  100.000     16        0       0    546   561      7   22 1.00e-03     30.1
    ## 69  100.000     14        0       0     51    64     15    2 1.00e-03     26.5
    ## 70   90.909     22        2       0   2160  2181      1   22 8.19e-04     31.9
    ## 71   90.909     22        2       0   2158  2179     22    1 8.19e-04     31.9
    ## 72   90.909     22        2       0   2120  2141      1   22 1.00e-03     31.9
    ## 73   90.909     22        2       0   2118  2139     22    1 1.00e-03     31.9
    ## 74   90.909     22        2       0   2120  2141      1   22 1.00e-03     31.9
    ## 75   90.909     22        2       0   2118  2139     22    1 1.00e-03     31.9
    ## 76  100.000     22        0       0   2556  2577     22    1 9.59e-07     41.0
    ## 77  100.000     17        0       0    370   386      6   22 4.60e-05     31.9
    ## 78  100.000     17        0       0    370   386      6   22 4.60e-05     31.9
    ## 79   90.476     21        2       0   2902  2922     21    1 1.00e-03     30.1
    ## 80   94.737     19        1       0   1559  1577     20    2 5.82e-04     31.0
    ## 81   94.737     19        1       0   1559  1577     20    2 5.82e-04     31.0
    ## 82   94.444     18        1       0    329   346      3   20 8.49e-04     29.2
    ## 83   94.444     18        1       0    329   346      3   20 8.49e-04     29.2
    ## 84  100.000     19        0       0    967   985      4   22 1.62e-05     35.6
    ## 85  100.000     22        0       0  20357 20378     22    1 3.48e-06     41.0
    ## 86   90.476     21        2       0    730   750      1   21 1.00e-03     30.1
    ## 87   94.737     19        1       0    672   690      4   22 5.71e-04     31.0
    ## 88   90.909     22        2       0   4352  4373      1   22 3.61e-04     31.9
    ## 89   95.000     20        1       0    196   215      3   22 3.04e-04     32.8
    ## 90   95.000     20        1       0    196   215      3   22 3.04e-04     32.8
    ## 91  100.000     18        0       0   1767  1784      5   22 6.34e-05     33.7
    ## 92   90.909     22        2       0   8371  8392     22    1 8.41e-04     31.9
    ## 93  100.000     22        0       0   6659  6680      1   22 1.47e-06     41.0
    ## 94  100.000     22        0       0   4375  4396      1   22 1.14e-06     41.0
    ## 95  100.000     22        0       0   8439  8460      1   22 1.81e-06     41.0
    ## 96  100.000     22        0       0   8437  8458     22    1 1.81e-06     41.0
    ## 97  100.000     22        0       0  11424 11445     22    1 1.82e-06     41.0
    ## 98  100.000     22        0       0  11426 11447      1   22 1.82e-06     41.0
    ## 99  100.000     22        0       0  11412 11433     22    1 2.79e-06     41.0
    ## 100 100.000     22        0       0  11414 11435      1   22 2.79e-06     41.0
    ## 101 100.000     22        0       0  11412 11433     22    1 5.02e-06     41.0
    ## 102 100.000     22        0       0  11414 11435      1   22 5.02e-06     41.0
    ## 103 100.000     22        0       0   4118  4139     22    1 1.07e-06     41.0
    ## 104 100.000     22        0       0   4120  4141      1   22 1.07e-06     41.0
    ## 105  95.455     22        1       0   4789  4810     22    1 4.44e-05     36.5
    ## 106  95.455     22        1       0   4791  4812      1   22 4.44e-05     36.5
    ## 107 100.000     22        0       0   8227  8248     22    1 1.79e-06     41.0
    ## 108 100.000     22        0       0   8229  8250      1   22 1.79e-06     41.0
    ## 109 100.000     22        0       0   5742  5763     22    1 1.79e-06     41.0
    ## 110 100.000     22        0       0   5744  5765      1   22 1.79e-06     41.0
    ## 111  90.909     22        2       0   3182  3203      1   22 3.53e-04     31.9
    ## 112  90.909     22        2       0   3180  3201     22    1 3.53e-04     31.9
    ## 113 100.000     16        0       0    608   623     21    6 2.59e-04     30.1
    ## 114  95.455     22        1       0  12457 12478      1   22 1.11e-04     36.5
    ## 115  90.909     22        2       0   7213  7234     22    1 5.53e-04     31.9
    ## 116  95.455     22        1       0   9393  9414      1   22 6.80e-05     36.5
    ## 117 100.000     16        0       0    776   791     21    6 3.16e-04     30.1
    ## 118  94.118     17        1       0    186   202     21    5 1.00e-03     27.4
    ## 119  90.000     20        2       0     93   112     22    3 4.23e-04     28.3
    ## 120 100.000     22        0       0   1720  1741      1   22 9.03e-07     41.0
    ## 121 100.000     22        0       0   1691  1712      1   22 9.03e-07     41.0
    ## 122  90.000     20        2       0    117   136     22    3 4.13e-04     28.3
    ## 123  90.476     21        1       1    276   296     21    2 1.00e-03     26.5
    ## 124  94.118     17        1       0    158   174      5   21 1.00e-03     27.4
    ## 125  86.364     22        3       0    152   173     22    1 1.00e-03     27.4
    ## 126  94.737     19        1       0    630   648      4   22 4.61e-04     31.0
    ## 127 100.000     19        0       0    107   125      4   22 1.15e-05     35.6
    ## 128 100.000     22        0       0   4388  4409     22    1 1.14e-06     41.0
    ## 129 100.000     22        0       0   4390  4411      1   22 1.14e-06     41.0
    ## 130  95.455     22        1       0  15339 15360      1   22 1.61e-04     36.5
    ## 131 100.000     15        0       0    427   441     22    8 1.00e-03     28.3
    ## 132  90.909     22        2       0   1132  1153      1   22 7.79e-04     31.9
    ## 133  90.909     22        2       0   4401  4422     22    1 7.79e-04     31.9
    ## 134  94.737     19        1       0    393   411      2   20 1.00e-03     31.0
    ## 135  90.000     20        2       0    192   211      3   22 4.23e-04     28.3
    ## 136  95.000     20        1       0      1    20      3   22 1.71e-04     32.8
    ## 137 100.000     22        0       0  13147 13168     22    1 3.10e-06     41.0
    ## 138 100.000     22        0       0  12779 12800     22    1 3.10e-06     41.0
    ## 139 100.000     22        0       0  13113 13134     22    1 2.95e-06     41.0
    ## 140 100.000     22        0       0  12745 12766     22    1 2.95e-06     41.0
    ## 141  90.909     22        2       0   1941  1962     22    1 2.15e-04     31.9
    ## 142  90.909     22        2       0    154   175      1   22 3.88e-04     31.9
    ## 143  90.909     22        2       0    152   173     22    1 3.88e-04     31.9
    ## 144 100.000     19        0       0   3654  3672      4   22 4.92e-05     35.6
    ## 145  90.909     22        2       0   3266  3287      1   22 6.00e-04     31.9
    ## 146 100.000     15        0       0    132   146      7   21 1.00e-03     28.3
    ## 147 100.000     17        0       0  12521 12537      6   22 1.00e-03     31.9
    ## 148 100.000     17        0       0   3164  3180      6   22 7.60e-04     31.9
    ## 149  91.667     24        0       1   4164  4187     22    1 4.54e-04     32.8
    ## 150  91.667     24        0       1   1596  1619     22    1 2.98e-04     32.8
    ## 151 100.000     17        0       0   3552  3568      6   22 8.17e-04     31.9
    ## 152 100.000     14        0       0     23    36      5   18 1.00e-03     26.5
    ## 153 100.000     14        0       0     23    36      5   18 1.00e-03     26.5
    ## 154 100.000     17        0       0   6029  6045     18    2 5.60e-04     31.9
    ## 155 100.000     16        0       0   2409  2424      4   19 1.00e-03     30.1
    ## 156 100.000     16        0       0   2409  2424     17    2 1.00e-03     30.1
    ## 157 100.000     16        0       0    367   382      3   18 1.76e-04     30.1
    ## 158 100.000     17        0       0   2145  2161     22    6 3.61e-04     31.9
    ## 159 100.000     16        0       0    187   202      7   22 8.81e-04     30.1
    ## 160 100.000     22        0       0   6320  6341     22    1 1.63e-06     41.0
    ## 161 100.000     22        0       0   6349  6370     22    1 1.63e-06     41.0
    ## 162  94.444     18        1       0     96   113      1   18 1.00e-03     29.2
    ## 163 100.000     19        0       0   7701  7719      4   22 5.63e-05     35.6
    ## 164  90.909     22        2       0   7329  7350      1   22 6.86e-04     31.9
    ## 165 100.000     19        0       0   3072  3090      4   22 3.17e-05     35.6
    ## 166  90.909     22        2       0   2700  2721      1   22 3.87e-04     31.9
    ## 167 100.000     18        0       0   2894  2911      3   20 1.02e-04     33.7
    ## 168 100.000     18        0       0   2894  2911      3   20 1.02e-04     33.7
    ## 169 100.000     18        0       0    786   803      3   20 1.94e-04     33.7
    ## 170 100.000     18        0       0    786   803      3   20 1.94e-04     33.7
    ## 171 100.000     18        0       0   4223  4240      3   20 1.36e-04     33.7
    ## 172 100.000     18        0       0   4223  4240      3   20 1.36e-04     33.7
    ## 173 100.000     15        0       0    720   734     22    8 1.00e-03     28.3
    ## 174  95.455     22        1       0   2985  3006     22    1 4.29e-05     36.5
    ## 175  90.909     22        2       0   4475  4496      1   22 5.23e-04     31.9
    ## 176 100.000     22        0       0   4589  4610     22    1 1.83e-06     41.0
    ## 177 100.000     22        0       0   4515  4536     22    1 1.82e-06     41.0
    ## 178  95.455     22        1       0   2804  2825     22    1 2.35e-05     36.5
    ## 179  95.455     22        1       0   2806  2827      1   22 2.35e-05     36.5
    ## 180 100.000     22        0       0  16364 16385     22    1 3.28e-06     41.0
    ## 181 100.000     22        0       0  16393 16414     22    1 3.28e-06     41.0
    ## 182  95.000     20        1       0   3244  3263     22    3 3.39e-04     32.8
    ## 183  95.000     20        1       0   8815  8834      1   20 1.00e-03     32.8
    ## 184  95.455     22        1       0   2882  2903     22    1 2.47e-05     36.5
    ## 185  95.238     21        1       0   2884  2904      1   21 8.64e-05     34.6
    ## 186  95.238     21        1       0   2306  2326      2   22 1.92e-04     34.6
    ## 187  90.909     22        2       0   2303  2324     22    1 6.69e-04     31.9
    ## 188  90.909     22        2       0   1368  1389      1   22 5.58e-04     31.9
    ## 189  95.455     22        1       0    232   253     22    1 2.73e-05     36.5
    ## 190  95.238     21        1       0    234   254      1   21 9.51e-05     34.6
    ## 191  95.455     22        1       0    546   567      1   22 1.48e-04     36.5
    ## 192  95.455     22        1       0    517   538      1   22 1.48e-04     36.5
    ## 193  90.476     21        2       0     65    85      1   21 2.54e-04     30.1
    ## 194 100.000     21        0       0   4673  4693     22    2 3.84e-06     39.2
    ## 195  95.455     22        1       0   4675  4696      1   22 4.68e-05     36.5
    ## 196  90.909     22        2       0   3502  3523      1   22 3.53e-04     31.9
    ## 197  90.909     22        2       0   2671  2692      1   22 6.13e-04     31.9
    ## 198  90.909     22        2       0   2671  2692      1   22 6.13e-04     31.9
    ## 199  90.909     22        2       0   2093  2114      1   22 1.00e-03     31.9
    ## 200  90.909     22        2       0   1874  1895      1   22 1.00e-03     31.9
    ## 201  90.909     22        2       0   2570  2591      1   22 2.61e-04     31.9
    ## 202  90.000     20        2       0    237   256      3   22 4.42e-04     28.3
    ## 203  95.455     22        1       0   8065  8086     22    1 1.57e-04     36.5
    ## 204  95.455     22        1       0   7799  7820     22    1 5.79e-05     36.5
    ## 205  95.455     22        1       0   2066  2087      1   22 2.28e-05     36.5
    ## 206  95.238     21        1       0   2065  2085     21    1 7.95e-05     34.6
    ## 207  95.455     22        1       0   2066  2087      1   22 3.33e-05     36.5
    ## 208  95.238     21        1       0   2065  2085     21    1 1.16e-04     34.6
    ## 209  95.455     22        1       0   2066  2087      1   22 3.33e-05     36.5
    ## 210  95.238     21        1       0   2065  2085     21    1 1.16e-04     34.6
    ## 211  95.455     22        1       0   1536  1557      1   22 1.74e-05     36.5
    ## 212  95.238     21        1       0   1535  1555     21    1 6.06e-05     34.6
    ## 213 100.000     22        0       0  13316 13337     22    1 2.01e-06     41.0
    ## 214 100.000     22        0       0  13316 13337     22    1 2.01e-06     41.0
    ## 215 100.000     22        0       0  13316 13337     22    1 2.01e-06     41.0
    ## 216 100.000     22        0       0  13242 13263     22    1 2.71e-06     41.0
    ## 217  95.455     22        1       0  11279 11300     22    1 6.44e-05     36.5
    ## 218 100.000     15        0       0    130   144     22    8 3.97e-04     28.3
    ## 219 100.000     17        0       0    689   705     18    2 1.76e-04     31.9
    ## 220  90.000     20        2       0    405   424     22    3 7.31e-04     28.3
    ## 221  95.455     22        1       0   1503  1524      1   22 5.65e-05     36.5
    ## 222  95.455     22        1       0   1458  1479      1   22 2.89e-05     36.5
    ## 223  95.455     22        1       0   1271  1292      1   22 2.20e-05     36.5
    ## 224 100.000     19        0       0    704   722      4   22 1.44e-05     35.6
    ## 225 100.000     22        0       0   2813  2834     22    1 7.58e-07     41.0
    ## 226 100.000     22        0       0   2226  2247     22    1 6.73e-07     41.0
    ## 227  95.455     22        1       0   7579  7600      1   22 6.89e-05     36.5
    ## 228  95.455     22        1       0   7577  7598     22    1 6.89e-05     36.5
    ## 229  94.444     18        1       0    113   130      5   22 5.59e-04     29.2
    ## 230 100.000     16        0       0    908   923     17    2 3.64e-04     30.1
    ## 231  95.238     21        1       0   9549  9569      2   22 6.01e-04     34.6
    ## 232 100.000     17        0       0    637   653     18    2 1.10e-04     31.9
    ## 233 100.000     15        0       0    252   266     16    2 1.00e-03     28.3
    ## 234 100.000     15        0       0      3    17     16    2 8.95e-04     28.3
    ## 235 100.000     15        0       0      3    17     16    2 8.95e-04     28.3
    ## 236 100.000     15        0       0      3    17     16    2 1.00e-03     28.3
    ## 237 100.000     15        0       0      3    17     16    2 1.00e-03     28.3
    ## 238  90.909     22        2       0    410   431      1   22 1.72e-04     31.9
    ## 239  95.455     22        1       0    359   380     22    1 3.40e-06     36.5
    ## 240  95.455     22        1       0    361   382      1   22 3.40e-06     36.5
    ## 241 100.000     15        0       0    971   985      2   16 1.00e-03     28.3
    ## 242  90.909     22        2       0   3055  3076     22    1 3.34e-04     31.9
    ## 243 100.000     19        0       0  23638 23656      4   22 1.49e-04     35.6
    ## 244 100.000     19        0       0  12467 12485      4   22 8.99e-05     35.6
    ## 245  90.909     22        2       0  12094 12115      1   22 1.00e-03     31.9
    ## 246 100.000     19        0       0  10897 10915      4   22 7.07e-05     35.6
    ## 247  90.909     22        2       0  10524 10545      1   22 8.61e-04     31.9
    ## 248  95.455     22        1       0   1953  1974     22    1 2.14e-05     36.5
    ## 249  94.737     19        1       0   1485  1503     22    4 6.70e-04     31.0
    ## 250 100.000     17        0       0   2709  2725     22    6 4.05e-04     31.9
    ## 251  95.455     22        1       0  13530 13551     22    1 8.67e-05     36.5
    ## 252  95.455     22        1       0    361   382      1   22 4.24e-06     36.5
    ## 253 100.000     16        0       0    371   386     17    2 2.88e-04     30.1
    ## 254  90.909     22        2       0    319   340     22    1 1.00e-03     31.9
    ## 255  95.238     21        1       0    440   460     22    2 1.07e-04     34.6
    ## 256  90.909     22        2       0    442   463      1   22 3.73e-04     31.9
    ## 257  95.455     22        1       0    883   904     22    1 2.22e-05     36.5
    ## 258 100.000     16        0       0     99   114     18    3 5.42e-04     30.1
    ## 259 100.000     22        0       0   8004  8025      1   22 1.52e-06     41.0
    ## 260 100.000     22        0       0   8002  8023     22    1 1.52e-06     41.0
    ## 261  90.909     22        2       0   4145  4166     22    1 3.63e-04     31.9
    ## 262  90.476     21        2       0   4147  4167      1   21 1.00e-03     30.1
    ## 263  90.909     22        2       0   4681  4702      1   22 5.45e-04     31.9
    ## 264  90.909     22        2       0   4679  4700     22    1 5.45e-04     31.9
    ## 265  95.455     22        1       0  25233 25254     22    1 1.53e-04     36.5
    ## 266  94.737     19        1       0    656   674      4   22 2.75e-04     31.0
    ## 267 100.000     17        0       0     71    87      4   20 1.34e-04     31.9
    ## 268 100.000     19        0       0   2002  2020      4   22 2.24e-05     35.6
    ## 269 100.000     19        0       0    467   485      4   22 5.60e-05     35.6
    ## 270 100.000     19        0       0   1628  1646      4   22 3.32e-05     35.6
    ## 271 100.000     19        0       0    437   455      4   22 8.84e-05     35.6
    ## 272  94.737     19        1       0    701   719      4   22 6.10e-04     31.0
    ## 273 100.000     19        0       0  10469 10487      4   22 6.24e-05     35.6
    ## 274 100.000     19        0       0   1230  1248      4   22 1.82e-05     35.6
    ## 275 100.000     18        0       0    628   645      3   20 2.23e-05     33.7
    ## 276  95.000     20        1       0   3239  3258     20    1 2.81e-04     32.8
    ## 277  95.000     20        1       0    289   308     20    1 3.34e-05     32.8
    ## 278  90.909     22        2       0    858   879     22    1 3.44e-04     31.9
    ## 279  90.909     22        2       0    860   881      1   22 3.44e-04     31.9
    ## 280  90.909     22        2       0   1721  1742      1   22 3.79e-04     31.9
    ## 281  90.000     20        2       0    286   305      1   20 1.00e-03     28.3
    ## 282 100.000     16        0       0    364   379     18    3 1.00e-03     30.1
    ## 283  95.455     22        1       0   3888  3909     22    1 2.91e-05     36.5
    ## 284  95.455     22        1       0  13229 13250     22    1 1.96e-04     36.5
    ## 285  95.455     22        1       0   2737  2758      1   22 3.94e-05     36.5
    ## 286  95.455     22        1       0   2723  2744      1   22 4.90e-05     36.5
    ## 287  95.455     22        1       0   2723  2744      1   22 4.90e-05     36.5
    ## 288 100.000     15        0       0    813   827     15    1 1.00e-03     28.3
    ## 289 100.000     22        0       0  20174 20195      1   22 2.96e-06     41.0
    ## 290 100.000     22        0       0  20157 20178      1   22 3.45e-06     41.0
    ## 291 100.000     22        0       0   5115  5136      1   22 1.58e-06     41.0
    ## 292  94.737     19        1       0   2283  2301      2   20 7.64e-04     31.0
    ## 293 100.000     22        0       0   3401  3422     22    1 1.41e-06     41.0
    ## 294  95.455     22        1       0   7795  7816     22    1 7.49e-05     36.5
    ## 295  95.455     22        1       0   7797  7818      1   22 7.49e-05     36.5
    ## 296  95.455     22        1       0   3330  3351      1   22 4.28e-05     36.5
    ## 297  95.455     22        1       0   1593  1614      1   22 3.21e-05     36.5
    ## 298  90.909     22        2       0   5269  5290     22    1 4.75e-04     31.9
    ## 299  90.909     22        2       0   5271  5292      1   22 4.75e-04     31.9
    ## 300 100.000     22        0       0   4110  4131     22    1 2.29e-06     41.0
    ## 301 100.000     22        0       0   7065  7086     22    1 1.36e-06     41.0
    ## 302 100.000     22        0       0   7094  7115     22    1 1.36e-06     41.0
    ## 303  95.455     22        1       0   2334  2355      1   22 4.16e-05     36.5
    ## 304  90.909     22        2       0   1179  1200     22    1 5.87e-04     31.9
    ## 305 100.000     22        0       0  18278 18299     22    1 2.42e-06     41.0
    ## 306 100.000     22        0       0  18280 18301      1   22 2.42e-06     41.0
    ## 307 100.000     22        0       0  12575 12596      1   22 1.96e-06     41.0
    ## 308 100.000     21        0       0  12574 12594     21    1 6.83e-06     39.2
    ## 309 100.000     22        0       0   5624  5645      1   22 1.28e-06     41.0
    ## 310 100.000     21        0       0   5623  5643     21    1 4.46e-06     39.2
    ## 311  90.476     21        2       0   1812  1832     21    1 8.95e-04     30.1
    ## 312  90.476     21        2       0   1812  1832     21    1 8.95e-04     30.1
    ## 313  90.476     21        2       0   1311  1331     21    1 1.00e-03     30.1
    ## 314  90.476     21        2       0     84   104      2   22 6.07e-04     30.1
    ## 315 100.000     19        0       0    817   835      4   22 5.72e-05     35.6
    ## 316 100.000     19        0       0    704   722      4   22 5.65e-05     35.6
    ## 317 100.000     19        0       0    689   707      4   22 3.34e-05     35.6
    ## 318  94.737     19        1       0    704   722      4   22 7.39e-04     31.0
    ## 319 100.000     19        0       0   1173  1191      2   20 1.56e-05     35.6
    ## 320  95.455     22        1       0   6699  6720      1   22 5.67e-05     36.5
    ## 321  95.455     22        1       0   6697  6718     22    1 5.67e-05     36.5
    ## 322 100.000     22        0       0  41060 41081     22    1 4.42e-06     41.0
    ## 323 100.000     22        0       0  41062 41083      1   22 4.42e-06     41.0
    ## 324 100.000     22        0       0  41060 41081     22    1 4.42e-06     41.0
    ## 325 100.000     22        0       0  41062 41083      1   22 4.42e-06     41.0
    ## 326 100.000     22        0       0  41041 41062     22    1 4.42e-06     41.0
    ## 327 100.000     22        0       0  41043 41064      1   22 4.42e-06     41.0
    ## 328 100.000     22        0       0  13332 13353     22    1 1.87e-06     41.0
    ## 329 100.000     22        0       0  13334 13355      1   22 1.87e-06     41.0
    ## 330  95.455     22        1       0    809   830     22    1 5.67e-05     36.5
    ## 331  95.455     22        1       0    811   832      1   22 5.67e-05     36.5
    ## 332  90.909     22        2       0   2043  2064     22    1 2.21e-04     31.9
    ## 333 100.000     22        0       0   3245  3266      1   22 3.78e-06     41.0
    ## 334 100.000     22        0       0   3243  3264     22    1 3.78e-06     41.0
    ## 335  90.909     22        2       0   4507  4528     22    1 7.74e-04     31.9
    ## 336  95.455     22        1       0   3128  3149     22    1 4.12e-05     36.5
    ## 337  95.238     21        1       0   3130  3150      1   21 1.44e-04     34.6
    ## 338 100.000     15        0       0    467   481      6   20 1.00e-03     28.3
    ## 339 100.000     16        0       0    942   957      4   19 6.27e-04     30.1
    ## 340 100.000     16        0       0    868   883      4   19 6.05e-04     30.1
    ## 341  95.455     22        1       0   8214  8235      1   22 7.33e-05     36.5
    ## 342  95.455     22        1       0   8208  8229      1   22 7.33e-05     36.5
    ## 343 100.000     22        0       0  16845 16866      1   22 4.17e-06     41.0
    ## 344 100.000     22        0       0  16845 16866      1   22 4.17e-06     41.0
    ## 345 100.000     22        0       0  16842 16863      1   22 4.07e-06     41.0
    ## 346 100.000     22        0       0  16842 16863      1   22 4.17e-06     41.0
    ## 347 100.000     22        0       0  16543 16564      1   22 4.13e-06     41.0
    ## 348 100.000     22        0       0    742   763      1   22 4.41e-07     41.0
    ## 349 100.000     22        0       0    740   761     22    1 4.41e-07     41.0
    ## 350 100.000     20        0       0   1998  2017      2   21 5.92e-06     37.4
    ## 351  95.455     22        1       0   7902  7923      1   22 1.44e-04     36.5
    ## 352  95.455     22        1       0   7375  7396      1   22 5.38e-05     36.5
    ## 353 100.000     19        0       0   2016  2034      4   22 2.60e-05     35.6
    ## 354  94.737     19        1       0    394   412      4   22 1.00e-03     31.0
    ## 355  86.364     22        3       0    150   171     22    1 1.00e-03     27.4
    ## 356 100.000     18        0       0   3260  3277     19    2 9.06e-05     33.7
    ## 357 100.000     15        0       0     81    95      5   19 4.40e-04     28.3
    ## 358  95.455     22        1       0   1132  1153     22    1 1.93e-05     36.5
    ## 359  95.000     20        1       0   1134  1153      1   20 2.35e-04     32.8
    ## 360 100.000     16        0       0    412   427      3   18 2.58e-04     30.1
    ## 361 100.000     16        0       0    678   693      3   18 1.00e-03     30.1
    ## 362  95.455     22        1       0   2946  2967      1   22 5.09e-05     36.5
    ## 363 100.000     16        0       0     77    92      6   21 5.96e-04     30.1
    ## 364  95.455     22        1       0   1858  1879      1   22 3.32e-05     36.5
    ## 365  95.455     22        1       0   1270  1291      1   22 2.94e-05     36.5
    ## 366  95.000     20        1       0  11803 11822     20    1 1.00e-03     32.8
    ## 367  95.000     20        1       0   1369  1388     20    1 1.35e-04     32.8
    ## 368 100.000     16        0       0    128   143     16    1 2.48e-04     30.1
    ## 369 100.000     17        0       0    821   837      1   17 1.72e-04     31.9
    ## 370  90.909     22        2       0   2141  2162      1   22 4.02e-04     31.9
    ## 371  94.444     18        1       0    161   178      3   20 4.30e-04     29.2
    ## 372  94.444     18        1       0    161   178      3   20 4.30e-04     29.2
    ## 373  94.444     18        1       0    322   339      3   20 1.00e-03     29.2
    ## 374  94.444     18        1       0    322   339      3   20 1.00e-03     29.2
    ## 375 100.000     15        0       0    718   732     16    2 1.00e-03     28.3
    ## 376 100.000     15        0       0    661   675     16    2 1.00e-03     28.3
    ## 377  95.000     20        1       0   9003  9022     22    3 8.29e-04     32.8
    ## 378  95.000     20        1       0   8948  8967     22    3 8.25e-04     32.8
    ## 379 100.000     15        0       0     48    62     16    2 4.62e-04     28.3
    ## 380 100.000     22        0       0   8241  8262      1   22 1.44e-06     41.0
    ## 381 100.000     22        0       0   8239  8260     22    1 1.44e-06     41.0
    ## 382 100.000     22        0       0   5971  5992      1   22 1.11e-06     41.0
    ## 383 100.000     22        0       0   5969  5990     22    1 1.11e-06     41.0
    ## 384  94.737     19        0       1    261   279     19    2 1.00e-03     27.4
    ## 385  94.737     19        1       0    976   994     22    4 3.38e-04     31.0
    ## 386 100.000     17        0       0  15060 15076      6   22 1.00e-03     31.9
    ## 387 100.000     17        0       0  14971 14987      6   22 1.00e-03     31.9
    ## 388 100.000     17        0       0  14971 14987      6   22 1.00e-03     31.9
    ## 389 100.000     17        0       0  14970 14986      6   22 1.00e-03     31.9
    ## 390 100.000     18        0       0   1710  1727      5   22 3.97e-04     33.7
    ## 391 100.000     22        0       0    917   938     22    1 3.57e-07     41.0
    ## 392 100.000     22        0       0   3061  3082      1   22 1.40e-06     41.0
    ## 393 100.000     22        0       0   3032  3053      1   22 1.40e-06     41.0
    ## 394 100.000     22        0       0   3061  3082      1   22 1.55e-06     41.0
    ## 395 100.000     22        0       0   3032  3053      1   22 1.55e-06     41.0
    ## 396  95.455     22        1       0   8919  8940     22    1 1.04e-04     36.5
    ## 397  95.238     21        1       0   8921  8941      1   21 3.65e-04     34.6
    ## 398  95.455     22        1       0   8919  8940     22    1 1.04e-04     36.5
    ## 399  95.238     21        1       0   8921  8941      1   21 3.65e-04     34.6
    ## 400 100.000     21        0       0   4501  4521      1   21 3.48e-06     39.2
    ## 401 100.000     21        0       0   3993  4013      1   21 3.48e-06     39.2
    ## 402 100.000     20        0       0   4018  4037     21    2 1.21e-05     37.4
    ## 403 100.000     21        0       0   3863  3883      1   21 3.15e-06     39.2
    ## 404 100.000     21        0       0   3355  3375      1   21 3.15e-06     39.2
    ## 405 100.000     20        0       0   3380  3399     21    2 1.10e-05     37.4
    ## 406 100.000     21        0       0  10335 10355     21    1 5.38e-06     39.2
    ## 407 100.000     21        0       0  10310 10330      1   21 5.38e-06     39.2
    ## 408 100.000     20        0       0  10310 10329      1   20 1.88e-05     37.4
    ## 409 100.000     21        0       0   9369  9389      1   21 5.64e-06     39.2
    ## 410 100.000     20        0       0   9370  9389      2   21 1.97e-05     37.4
    ## 411 100.000     19        0       0   9370  9388      2   20 6.87e-05     35.6
    ## 412 100.000     21        0       0   4492  4512      1   21 3.16e-06     39.2
    ## 413 100.000     20        0       0   4493  4512      2   21 1.10e-05     37.4
    ## 414 100.000     19        0       0   4493  4511      2   20 3.85e-05     35.6
    ## 415 100.000     19        0       0   7329  7347      2   20 6.68e-05     35.6
    ## 416 100.000     19        0       0   7329  7347      2   20 1.08e-04     35.6
    ## 417  94.444     18        1       0     15    32      3   20 7.85e-04     29.2
    ## 418  94.444     18        1       0     15    32      3   20 7.85e-04     29.2
    ## 419 100.000     22        0       0  12317 12338      1   22 2.10e-06     41.0
    ## 420 100.000     22        0       0    352   373     22    1 9.91e-07     41.0
    ## 421 100.000     22        0       0    344   365     22    1 1.09e-06     41.0

421 putative lncRNA sponges with these parameters.

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
../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_ShortStack_4.1.0_mature.fasta \
../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve_lncRNA.fasta \
-sc 100 \
-en -20 \
-strict \
-out ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve-miRanda-lncRNA-strict_all.tab
```

# 5 Summarize results

Let’s look at the output

``` bash

echo "miranda run finished!"
echo "Counting number of interacting miRNA-lncRNA pairs"

zgrep -c "Performing Scan" ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve-miRanda-lncRNA-strict_all.tab

echo "Parsing output"
grep -A 1 "Scores for this hit:" ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve-miRanda-lncRNA-strict_all.tab | sort | grep '>' > ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve-miRanda-lncRNA-strict-parsed.txt

echo "counting number of putative interactions predicted (can include multiple interactions between single miRNA-lncRNA pair)"
wc -l ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve-miRanda-lncRNA-strict_all.tab
```

    ## miranda run finished!
    ## Counting number of interacting miRNA-lncRNA pairs
    ## 2241495
    ## Parsing output
    ## counting number of putative interactions predicted (can include multiple interactions between single miRNA-lncRNA pair)
    ## 18532367 ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve-miRanda-lncRNA-strict_all.tab

2690450 – This is a lot of putative interactions! We can probably narrow
it down though. In vertebrates, miRNA-mRNA binding only requires
complementarity of an miRNA seed region of \~8 nucleotides. This
requirement is built in to miRanda target prediction. In cnidarians,
however, miRNA-mRNA binding is believed to require near-complete
complementarity of the full mature miRNA, similarly to plants ( Admoni
et al. ([2023](#ref-admoni_target_2023)) , Admoni et al.
([2025](#ref-admoni_mirna-target_2025)) ). While I couldn’t find any
information on expected requirements for miRNA-lncRNA sponges, its
possible the binding will function similarly to miRNA-mRNA binding.
Let’s look at how many putative interactions are predicted for a binding
length of at least 21 nucleotides (the length of our smallest mature
miRNA).

``` bash
echo "number of putative interactions of at least 21 nucleotides"
awk -F'\t' '$7 >= 21' ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve-miRanda-lncRNA-strict-parsed.txt | wc -l
echo ""
echo "check some:"
awk -F'\t' '$7 >= 21' ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve-miRanda-lncRNA-strict-parsed.txt | head -5
```

    ## number of putative interactions of at least 21 nucleotides
    ## 9625
    ## 
    ## check some:
    ## >Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)   transcript::Porites_evermani_scaffold_1294:95643-100854 167.00  -21.39  2 21    3383 3407   22  72.73%  81.82%
    ## >Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)   transcript::Porites_evermani_scaffold_1294:95707-100854 167.00  -21.39  2 21    3319 3343   22  72.73%  81.82%
    ## >Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)   transcript::Porites_evermani_scaffold_1775:21921-37427  168.00  -20.04  2 21    14864 14886 21  85.71%  85.71%
    ## >Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)   transcript::Porites_evermani_scaffold_1847:9535-13170   146.00  -20.84  2 20    858 887 26  57.69%  69.23%
    ## >Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)   transcript::Porites_evermani_scaffold_2579:22-36631 163.00  -20.08  2 21    35697 35724 25  72.00%  72.00%

9625

The header for this output is formatted as:

mirna Target Score Energy-Kcal/Mol Query-Aln(start-end)
Subject-Al(Start-End) Al-Len Subject-Identity Query-Identity

We can see from the percent identities (last 2 entries) that this number
includes alignments with multiple mismatches. Let’s filter again to
reduce the number of permissible mismatches. Let’s say we want no more
than 3 mismatches (a gap is counted as a mismatch). For an alignment of
21 nucleotides, this would be an percent identity of (21-3)/21 = 85.7%.
The miRNA is our “subject”, so we will filter by column 8.

``` bash
echo "number of putative interactions of at least 21 nucleotides, with at most 3 mismatches"
awk -F'\t' '$7 >= 21' ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve-miRanda-lncRNA-strict-parsed.txt | awk -F'\t' '$8 >= 85' | wc -l
echo ""
echo "check some:"
awk -F'\t' '$7 >= 21' ../output/14-Peve-miRNA-lncRNA-BLASTs-miRanda/Peve-miRanda-lncRNA-strict-parsed.txt | awk -F'\t' '$8 >= 85' | head -5
```

    ## number of putative interactions of at least 21 nucleotides, with at most 3 mismatches
    ## 61
    ## 
    ## check some:
    ## >Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)   transcript::Porites_evermani_scaffold_1775:21921-37427  168.00  -20.04  2 21    14864 14886 21  85.71%  85.71%
    ## >Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)   transcript::Porites_evermani_scaffold_687:21178-58742   178.00  -22.78  2 21    31762 31785 21  85.71%  90.48%
    ## >Cluster_10060.mature::Porites_evermani_scaffold_1159:6702-6723(+)   transcript::Porites_evermani_scaffold_687:25412-58742   178.00  -22.78  2 21    27528 27551 21  85.71%  90.48%
    ## >Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)   transcript::Porites_evermani_scaffold_133:272228-275815 183.00  -26.00  2 21    1940 1964   22  86.36%  86.36%
    ## >Cluster_10061.mature::Porites_evermani_scaffold_1159:7428-7449(+)   transcript::Porites_evermani_scaffold_2270:70952-76996  187.00  -28.27  2 21    4164 4187   21  90.48%  90.48%

This is a dramatically smaller number – only 61 interactions were at
least 21 nucleotides with \<=3 mismatches

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
