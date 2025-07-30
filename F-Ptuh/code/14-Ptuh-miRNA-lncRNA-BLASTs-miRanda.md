14-Ptuh-miRNA-lncRNA-BLASTs-miRanda
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
full_mirna_fasta="../output/05-Ptuh-sRNA-ShortStack_4.1.0/ShortStack_out/mir.fasta"
premirna_fasta="../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_precursor.fasta"
mature_mirna_fasta="../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_mature.fasta"
star_mirna_fasta="../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_star.fasta"

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
premirna_fasta="../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_precursor.fasta"
mature_mirna_fasta="../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_mature.fasta"
star_mirna_fasta="../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_star.fasta"

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

    ## >Cluster_21::Pocillopora_meandrina_HIv1___Sc0000000:818027-818120(+)
    ## >Cluster_21.mature::Pocillopora_meandrina_HIv1___Sc0000000:818049-818070(+)
    ## 
    ## >Cluster_21.mature::Pocillopora_meandrina_HIv1___Sc0000000:818049-818070(+)
    ## >Cluster_36.mature::Pocillopora_meandrina_HIv1___Sc0000000:2872041-2872061(+)
    ## 
    ## >Cluster_21.star::Pocillopora_meandrina_HIv1___Sc0000000:818079-818100(+)
    ## >Cluster_36.star::Pocillopora_meandrina_HIv1___Sc0000000:2872070-2872090(+)
    ## 
    ## 111
    ## 
    ## 37
    ## 
    ## 37

## 1.2 Check miRNA lengths

``` bash
# Extract sequence lengths for precursors
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_precursor.fasta > ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_precursor_lengths.txt

# Sequence lengths for matures
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_mature.fasta > ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_mature_lengths.txt
```

``` r
# Summary stats of precursor and mature lengths

precursor_lengths <- read.table("../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_precursor_lengths.txt", sep = " ", header = FALSE, col.names = c("seqID", "length"))
mature_lengths <- read.table("../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_mature_lengths.txt", sep = " ", header = FALSE, col.names = c("seqID", "length"))

cat("Average pre-miRNA length: ", mean(precursor_lengths$length))
```

    ## Average pre-miRNA length:  46.25225

``` r
cat("\n")
```

``` r
cat("Range of pre-miRNA lengths: ", range(precursor_lengths$length))
```

    ## Range of pre-miRNA lengths:  21 101

``` r
cat("\n")
```

``` r
cat("Average mature miRNA length: ", mean(mature_lengths$length))
```

    ## Average mature miRNA length:  21.91892

``` r
cat("\n")
```

``` r
cat("Range of mature miRNA lengths: ", range(mature_lengths$length))
```

    ## Range of mature miRNA lengths:  21 23

## 1.3 check lncRNAs

LncRNAs were identified from Ptuh RNA-seq data – see details in
`F-Ptuh/code/17-Ptuh-lncRNA.Rmd`

Fasta of Ptuh lncRNAs stored at
`https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/F-Ptuh/output/17-Ptuh-lncRNA/Ptuh-lncRNA.fasta`

``` bash
curl -L https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/F-Ptuh/output/17-Ptuh-lncRNA/Ptuh-lncRNA.fasta -o ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_lncRNA.fasta

echo "Number of lncRNAs:"
grep "^>" ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_lncRNA.fasta | wc -l
```

    ##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
    ##                                  Dload  Upload   Total   Spent    Left  Speed
    ##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  0     0    0     0    0     0      0      0 --:--:--  0:00:01 --:--:--     0  0 49.0M    0     0    0     0      0      0 --:--:--  0:00:01 --:--:--     0100 49.0M  100 49.0M    0     0  23.6M      0  0:00:02  0:00:02 --:--:-- 23.6M
    ## Number of lncRNAs:
    ## 16153

16153 total lncRNA

``` bash
# Extract sequence lengths for precursors
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_lncRNA.fasta > ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_lncRNA_lengths.txt
```

``` r
# Summary stats of lncRNA lengths

lncRNA_lengths <- read.table("../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_lncRNA_lengths.txt", sep = " ", header = FALSE, col.names = c("seqID", "length"))

cat("Average lncRNA length: ", mean(lncRNA_lengths$length))
```

    ## Average lncRNA length:  3123.615

``` r
cat("\n")
```

``` r
cat("Range of lncRNA lengths: ", range(lncRNA_lengths$length))
```

    ## Range of lncRNA lengths:  201 227014

``` r
ggplot(lncRNA_lengths, aes(x = length)) +
  geom_histogram(binwidth = 500) +
  labs(title = "A. pulchra lncRNA sequence lengths",
       x = "Sequence Length [nucleotides]",
       y = "Frequency") +
  xlim(200, 250000) +
  ylim(0, 800) +
  theme_minimal()
```

    ## Warning: Removed 6 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](14-Ptuh-miRNA-lncRNA-BLASTs-miRanda_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# 2 BLASTs

## 2.1 Make databases

Database of pre-miRNAs:

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_precursor.fasta \
-dbtype nucl \
-out ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/Ptuh-db/Ptuh_ShortStack_4.1.0_precursor
```

Database of mature miRNAs:

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_mature.fasta \
-dbtype nucl \
-out ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/Ptuh-db/Ptuh_ShortStack_4.1.0_mature
```

## 2.2 Run BLASTn

Generate a list of blast results. It seems plausible that a single
lncRNA, which would be hundreds or thousands of nucleotides long, could
interact with multiple miRNAs, so I will allow up to 10 hits (\~25% of
Ptuh miRNAs) for each lncRNA. I want to see the top hits no matter how
poor the match is, so I will not filter by e-value at this stage. I’ll
also include the “-word_size 4” option, which reduces the required
length of the initial match.

Full pre-miRNAs:

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_lncRNA.fasta \
-db ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/Ptuh-db/Ptuh_ShortStack_4.1.0_precursor \
-out ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_precursor_blastn.tab \
-num_threads 40 \
-word_size 4 \
-max_target_seqs 10 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_precursor_blastn.tab
```

    ## 156432 ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_precursor_blastn.tab

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
-query ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_lncRNA.fasta \
-db ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/Ptuh-db/Ptuh_ShortStack_4.1.0_mature \
-out ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_mature_blastn.tab \
-num_threads 40 \
-word_size 4 \
-max_target_seqs 10 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_mature_blastn.tab
```

    ## 155901 ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_mature_blastn.tab

# 3 Examine BLAST tables

Read into R and assign informative column labels

``` r
precursor_lncRNA_BLASTn <- read.table("../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_precursor_blastn.tab", sep="\t", header=FALSE)
mature_lncRNA_BLASTn <- read.table("../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/blasts/lncRNA_to_mature_blastn.tab", sep="\t", header=FALSE)

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

    ## [1] 8

``` r
precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0) %>%
  select(qseqid) %>%
  unique() %>%
  nrow()
```

    ## [1] 8

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

    ##             qseqid
    ## 1 Ptuh_lncRNA_5896
    ## 2 Ptuh_lncRNA_5897
    ## 3 Ptuh_lncRNA_5898
    ## 4 Ptuh_lncRNA_7177
    ## 5 Ptuh_lncRNA_7178
    ## 6 Ptuh_lncRNA_7347
    ## 7 Ptuh_lncRNA_7348
    ## 8 Ptuh_lncRNA_7349
    ##                                                                    sseqid
    ## 1 Cluster_2859::Pocillopora_meandrina_HIv1___Sc0000008:5387738-5387830(-)
    ## 2 Cluster_2859::Pocillopora_meandrina_HIv1___Sc0000008:5387738-5387830(-)
    ## 3 Cluster_2859::Pocillopora_meandrina_HIv1___Sc0000008:5387738-5387830(-)
    ## 4 Cluster_3392::Pocillopora_meandrina_HIv1___Sc0000010:8012916-8013008(+)
    ## 5 Cluster_3392::Pocillopora_meandrina_HIv1___Sc0000010:8012916-8013008(+)
    ## 6 Cluster_3392::Pocillopora_meandrina_HIv1___Sc0000010:8012916-8013008(+)
    ## 7 Cluster_3392::Pocillopora_meandrina_HIv1___Sc0000010:8012916-8013008(+)
    ## 8 Cluster_3392::Pocillopora_meandrina_HIv1___Sc0000010:8012916-8013008(+)
    ##   pident length mismatch gapopen qstart qend sstart send   evalue bitscore
    ## 1    100     93        0       0   6176 6268     93    1 3.38e-44      168
    ## 2    100     93        0       0   5570 5662     93    1 3.13e-44      168
    ## 3    100     93        0       0   5169 5261     93    1 2.96e-44      168
    ## 4    100     93        0       0   4377 4469      1   93 2.33e-44      168
    ## 5    100     93        0       0   2315 2407      1   93 1.52e-44      168
    ## 6    100     93        0       0   4151 4243      1   93 2.13e-44      168
    ## 7    100     93        0       0   3877 3969      1   93 2.01e-44      168
    ## 8    100     93        0       0   3727 3819      1   93 2.02e-44      168

We have 8 alignments of a full pre-miRNA to a lncRNA with no mismatches.
8 lncRNA and 2 miRNA are represented.

Note that, as in A.pulchra and P.evermanni, there are instances of a
single pre-miRNA matching to several overlapping lncRNA, potentially
representing multiple isoforms of a single lncRNA gene.

Save these results

``` r
precursor_lncRNAs <- precursor_lncRNA_BLASTn %>% 
  filter(length >= 90) %>%
  filter(mismatch == 0)

write.table(precursor_lncRNAs, "../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/lncRNAs_as_miRNA_precursors.txt")
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
  filter(evalue <= 0.001) %>%
  unique()
```

    ##               qseqid
    ## 1    Ptuh_lncRNA_109
    ## 2    Ptuh_lncRNA_385
    ## 3    Ptuh_lncRNA_703
    ## 4    Ptuh_lncRNA_957
    ## 5    Ptuh_lncRNA_959
    ## 6    Ptuh_lncRNA_960
    ## 7    Ptuh_lncRNA_961
    ## 8    Ptuh_lncRNA_962
    ## 9    Ptuh_lncRNA_965
    ## 10   Ptuh_lncRNA_966
    ## 11   Ptuh_lncRNA_968
    ## 12   Ptuh_lncRNA_970
    ## 13   Ptuh_lncRNA_972
    ## 14   Ptuh_lncRNA_972
    ## 15   Ptuh_lncRNA_973
    ## 16   Ptuh_lncRNA_978
    ## 17   Ptuh_lncRNA_979
    ## 18  Ptuh_lncRNA_1120
    ## 19  Ptuh_lncRNA_1169
    ## 20  Ptuh_lncRNA_1169
    ## 21  Ptuh_lncRNA_3404
    ## 22  Ptuh_lncRNA_3514
    ## 23  Ptuh_lncRNA_4091
    ## 24  Ptuh_lncRNA_4407
    ## 25  Ptuh_lncRNA_4555
    ## 26  Ptuh_lncRNA_4950
    ## 27  Ptuh_lncRNA_5760
    ## 28  Ptuh_lncRNA_5896
    ## 29  Ptuh_lncRNA_5897
    ## 30  Ptuh_lncRNA_5898
    ## 31  Ptuh_lncRNA_6284
    ## 32  Ptuh_lncRNA_6285
    ## 33  Ptuh_lncRNA_6288
    ## 34  Ptuh_lncRNA_6297
    ## 35  Ptuh_lncRNA_6526
    ## 36  Ptuh_lncRNA_6621
    ## 37  Ptuh_lncRNA_6633
    ## 38  Ptuh_lncRNA_6643
    ## 39  Ptuh_lncRNA_6648
    ## 40  Ptuh_lncRNA_7177
    ## 41  Ptuh_lncRNA_7178
    ## 42  Ptuh_lncRNA_7179
    ## 43  Ptuh_lncRNA_7346
    ## 44  Ptuh_lncRNA_7347
    ## 45  Ptuh_lncRNA_7348
    ## 46  Ptuh_lncRNA_7349
    ## 47  Ptuh_lncRNA_8022
    ## 48  Ptuh_lncRNA_8416
    ## 49  Ptuh_lncRNA_8443
    ## 50  Ptuh_lncRNA_8464
    ## 51  Ptuh_lncRNA_8641
    ## 52  Ptuh_lncRNA_8884
    ## 53 Ptuh_lncRNA_10892
    ## 54 Ptuh_lncRNA_10961
    ## 55 Ptuh_lncRNA_10961
    ## 56 Ptuh_lncRNA_10962
    ## 57 Ptuh_lncRNA_10962
    ## 58 Ptuh_lncRNA_11583
    ## 59 Ptuh_lncRNA_12865
    ## 60 Ptuh_lncRNA_12866
    ## 61 Ptuh_lncRNA_13340
    ## 62 Ptuh_lncRNA_13982
    ## 63 Ptuh_lncRNA_14360
    ## 64 Ptuh_lncRNA_14969
    ## 65 Ptuh_lncRNA_15005
    ## 66 Ptuh_lncRNA_15006
    ## 67 Ptuh_lncRNA_16035
    ## 68 Ptuh_lncRNA_16036
    ## 69 Ptuh_lncRNA_16064
    ## 70 Ptuh_lncRNA_16065
    ## 71 Ptuh_lncRNA_16070
    ## 72 Ptuh_lncRNA_16079
    ##                                                                            sseqid
    ## 1  Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 2  Cluster_4435.mature::Pocillopora_meandrina_HIv1___Sc0000016:7549579-7549600(+)
    ## 3  Cluster_2839.mature::Pocillopora_meandrina_HIv1___Sc0000008:3619363-3619384(-)
    ## 4  Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 5  Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 6  Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 7  Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 8  Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 9  Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 10 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 11 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 12 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 13 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 14 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 15 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 16 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 17 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 18     Cluster_21.mature::Pocillopora_meandrina_HIv1___Sc0000000:818049-818070(+)
    ## 19 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 20     Cluster_21.mature::Pocillopora_meandrina_HIv1___Sc0000000:818049-818070(+)
    ## 21 Cluster_4118.mature::Pocillopora_meandrina_HIv1___Sc0000014:9365481-9365502(-)
    ## 22 Cluster_4118.mature::Pocillopora_meandrina_HIv1___Sc0000014:9365481-9365502(-)
    ## 23   Cluster_1793.mature::Pocillopora_meandrina_HIv1___Sc0000005:601626-601646(+)
    ## 24   Cluster_1793.mature::Pocillopora_meandrina_HIv1___Sc0000005:601626-601646(+)
    ## 25 Cluster_4754.mature::Pocillopora_meandrina_HIv1___Sc0000018:4564906-4564927(+)
    ## 26   Cluster_1793.mature::Pocillopora_meandrina_HIv1___Sc0000005:601626-601646(+)
    ## 27 Cluster_4437.mature::Pocillopora_meandrina_HIv1___Sc0000016:7550615-7550635(+)
    ## 28 Cluster_2859.mature::Pocillopora_meandrina_HIv1___Sc0000008:5387760-5387781(-)
    ## 29 Cluster_2859.mature::Pocillopora_meandrina_HIv1___Sc0000008:5387760-5387781(-)
    ## 30 Cluster_2859.mature::Pocillopora_meandrina_HIv1___Sc0000008:5387760-5387781(-)
    ## 31 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 32 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 33 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 34 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 35 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 36 Cluster_4439.mature::Pocillopora_meandrina_HIv1___Sc0000016:7555692-7555713(+)
    ## 37 Cluster_4439.mature::Pocillopora_meandrina_HIv1___Sc0000016:7555692-7555713(+)
    ## 38 Cluster_4439.mature::Pocillopora_meandrina_HIv1___Sc0000016:7555692-7555713(+)
    ## 39 Cluster_4439.mature::Pocillopora_meandrina_HIv1___Sc0000016:7555692-7555713(+)
    ## 40 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 41 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 42 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 43 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 44 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 45 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 46 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 47 Cluster_4118.mature::Pocillopora_meandrina_HIv1___Sc0000014:9365481-9365502(-)
    ## 48 Cluster_5612.mature::Pocillopora_meandrina_HIv1___Sc0000024:4808688-4808708(+)
    ## 49   Cluster_1793.mature::Pocillopora_meandrina_HIv1___Sc0000005:601626-601646(+)
    ## 50 Cluster_2839.mature::Pocillopora_meandrina_HIv1___Sc0000008:3619363-3619384(-)
    ## 51 Cluster_2839.mature::Pocillopora_meandrina_HIv1___Sc0000008:3619363-3619384(-)
    ## 52 Cluster_2859.mature::Pocillopora_meandrina_HIv1___Sc0000008:5387760-5387781(-)
    ## 53 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 54     Cluster_21.mature::Pocillopora_meandrina_HIv1___Sc0000000:818049-818070(+)
    ## 55 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 56 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 57     Cluster_21.mature::Pocillopora_meandrina_HIv1___Sc0000000:818049-818070(+)
    ## 58 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 59 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 60 Cluster_6382.mature::Pocillopora_meandrina_HIv1___Sc0000035:1989842-1989863(+)
    ## 61 Cluster_2859.mature::Pocillopora_meandrina_HIv1___Sc0000008:5387760-5387781(-)
    ## 62 Cluster_2839.mature::Pocillopora_meandrina_HIv1___Sc0000008:3619363-3619384(-)
    ## 63 Cluster_4118.mature::Pocillopora_meandrina_HIv1___Sc0000014:9365481-9365502(-)
    ## 64 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 65 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 66 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 67 Cluster_2839.mature::Pocillopora_meandrina_HIv1___Sc0000008:3619363-3619384(-)
    ## 68 Cluster_2839.mature::Pocillopora_meandrina_HIv1___Sc0000008:3619363-3619384(-)
    ## 69 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 70 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 71 Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)
    ## 72 Cluster_4118.mature::Pocillopora_meandrina_HIv1___Sc0000014:9365481-9365502(-)
    ##     pident length mismatch gapopen qstart  qend sstart send   evalue bitscore
    ## 1   95.000     20        1       0   1606  1625     20    1 4.94e-04     32.8
    ## 2  100.000     15        0       0    103   117      6   20 8.38e-04     28.3
    ## 3   90.476     21        1       1    101   121      3   22 8.90e-04     26.5
    ## 4   95.238     21        1       0  72355 72375      2   22 1.00e-03     34.6
    ## 5   95.238     21        1       0  72351 72371      2   22 1.00e-03     34.6
    ## 6   90.909     22        2       0  18795 18816      1   22 1.00e-03     31.9
    ## 7   90.909     22        2       0    611   632      1   22 1.96e-04     31.9
    ## 8   95.238     21        1       0  11097 11117      2   22 5.51e-04     34.6
    ## 9   95.238     21        1       0  29755 29775      2   22 4.98e-04     34.6
    ## 10  95.238     21        1       0  29755 29775      2   22 4.98e-04     34.6
    ## 11  95.238     21        1       0   5353  5373      2   22 4.09e-04     34.6
    ## 12  95.238     21        1       0   4774  4794      2   22 3.93e-04     34.6
    ## 13  95.238     21        1       0  13596 13616      2   22 3.98e-04     34.6
    ## 14  90.909     22        2       0  12464 12485      1   22 1.00e-03     31.9
    ## 15  95.238     21        1       0    106   126      2   22 2.26e-04     34.6
    ## 16  90.000     20        2       0   1109  1128      3   22 1.00e-03     28.3
    ## 17  90.000     20        2       0   1096  1115      3   22 1.00e-03     28.3
    ## 18 100.000     16        0       0     71    86      5   20 1.38e-04     30.1
    ## 19  86.364     22        3       0    268   289     22    1 1.00e-03     27.4
    ## 20  89.474     19        2       0    300   318     19    1 1.00e-03     26.5
    ## 21  90.000     20        2       0     46    65     20    1 1.00e-03     28.3
    ## 22 100.000     18        0       0  14602 14619     20    3 3.61e-04     33.7
    ## 23 100.000     15        0       0    709   723      1   15 1.00e-03     28.3
    ## 24 100.000     15        0       0    442   456      1   15 1.00e-03     28.3
    ## 25  90.476     21        2       0   3469  3489      1   21 1.00e-03     30.1
    ## 26  90.000     20        2       0    985  1004      2   21 1.00e-03     28.3
    ## 27 100.000     15        0       0   1071  1085      1   15 1.00e-03     28.3
    ## 28 100.000     22        0       0   6198  6219     22    1 9.65e-07     41.0
    ## 29 100.000     22        0       0   5592  5613     22    1 8.92e-07     41.0
    ## 30 100.000     22        0       0   5191  5212     22    1 8.44e-07     41.0
    ## 31  90.909     22        2       0  14178 14199      1   22 1.00e-03     31.9
    ## 32  95.455     22        1       0  23095 23116      1   22 1.05e-04     36.5
    ## 33  90.909     22        2       0  13691 13712      1   22 8.78e-04     31.9
    ## 34  90.476     21        2       0    328   348      2   22 3.14e-04     30.1
    ## 35  95.455     22        1       0    661   682      1   22 6.54e-06     36.5
    ## 36  90.476     21        2       0   5282  5302      2   22 1.00e-03     30.1
    ## 37  90.476     21        2       0     34    54      2   22 1.00e-03     30.1
    ## 38  90.476     21        2       0   2598  2618      2   22 7.39e-04     30.1
    ## 39  90.476     21        2       0    106   126      2   22 2.03e-04     30.1
    ## 40 100.000     22        0       0   4399  4420      1   22 7.61e-07     41.0
    ## 41 100.000     22        0       0   2337  2358      1   22 4.78e-07     41.0
    ## 42  95.238     21        1       0   1583  1603     22    2 4.89e-05     34.6
    ## 43  95.238     21        1       0   2425  2445     22    2 7.18e-05     34.6
    ## 44 100.000     22        0       0   4173  4194      1   22 6.95e-07     41.0
    ## 45 100.000     22        0       0   3899  3920      1   22 6.57e-07     41.0
    ## 46 100.000     22        0       0   3749  3770      1   22 6.37e-07     41.0
    ## 47  90.000     20        2       0     29    48      3   22 1.00e-03     28.3
    ## 48 100.000     14        0       0     72    85     14    1 1.00e-03     26.5
    ## 49  94.737     19        1       0   3069  3087      1   19 1.00e-03     31.0
    ## 50  94.737     19        1       0    281   299      1   19 8.05e-04     31.0
    ## 51  94.737     19        1       0   1469  1487     19    1 5.82e-04     31.0
    ## 52 100.000     15        0       0    435   449      3   17 1.00e-03     28.3
    ## 53  89.474     19        2       0    153   171     22    4 1.00e-03     26.5
    ## 54 100.000     19        0       0    127   145     19    1 2.96e-06     35.6
    ## 55  90.909     22        2       0     95   116     22    1 3.61e-05     31.9
    ## 56  90.909     22        2       0   1352  1373     22    1 1.58e-04     31.9
    ## 57 100.000     16        0       0   1384  1399     19    4 5.51e-04     30.1
    ## 58 100.000     17        0       0   4642  4658      5   21 1.00e-03     31.9
    ## 59  95.238     21        1       0     66    86     22    2 4.71e-05     34.6
    ## 60  95.238     21        1       0     66    86     22    2 4.71e-05     34.6
    ## 61 100.000     16        0       0    122   137      5   20 2.31e-04     30.1
    ## 62  94.118     17        1       0    243   259     20    4 1.00e-03     27.4
    ## 63 100.000     15        0       0     83    97      6   20 2.80e-04     28.3
    ## 64  95.455     22        1       0  15552 15573      1   22 8.32e-05     36.5
    ## 65  95.455     22        1       0  50963 50984      1   22 1.86e-04     36.5
    ## 66  95.455     22        1       0  36945 36966      1   22 1.63e-04     36.5
    ## 67  94.737     19        1       0    353   371      1   19 5.32e-04     31.0
    ## 68  94.737     19        1       0    321   339      1   19 5.23e-04     31.0
    ## 69  90.909     22        2       0  25202 25223      1   22 1.00e-03     31.9
    ## 70  90.909     22        2       0  23694 23715      1   22 1.00e-03     31.9
    ## 71  95.455     22        1       0    361   382      1   22 3.73e-06     36.5
    ## 72  90.000     20        2       0     29    48      3   22 1.00e-03     28.3

72 putative lncRNA sponges with these parameters.

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
../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_ShortStack_4.1.0_mature.fasta \
../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh_lncRNA.fasta \
-sc 100 \
-en -20 \
-strict \
-out ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh-miRanda-lncRNA-strict_all.tab
```

# 5 Summarize results

Let’s look at the output

``` bash

echo "miranda run finished!"
echo "Counting number of interacting miRNA-lncRNA pairs"

zgrep -c "Performing Scan" ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh-miRanda-lncRNA-strict_all.tab

echo "Parsing output"
grep -A 1 "Scores for this hit:" ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh-miRanda-lncRNA-strict_all.tab | sort | grep '>' > ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh-miRanda-lncRNA-strict-parsed.txt

echo "counting number of putative interactions predicted (can include multiple interactions between single miRNA-lncRNA pair)"
wc -l ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh-miRanda-lncRNA-strict_all.tab
```

    ## miranda run finished!
    ## Counting number of interacting miRNA-lncRNA pairs
    ## 597661
    ## Parsing output
    ## counting number of putative interactions predicted (can include multiple interactions between single miRNA-lncRNA pair)
    ## 4887694 ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh-miRanda-lncRNA-strict_all.tab

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
awk -F'\t' '$7 >= 21' ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh-miRanda-lncRNA-strict-parsed.txt | wc -l
echo ""
echo "check some:"
awk -F'\t' '$7 >= 21' ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh-miRanda-lncRNA-strict-parsed.txt | head -5
```

    ## number of putative interactions of at least 21 nucleotides
    ## 1831
    ## 
    ## check some:
    ## >Cluster_1015.mature::Pocillopora_meandrina_HIv1___Sc0000002:10845744-10845765(-)    Ptuh_lncRNA_008 159.00  -22.47  2 21    2138 2162   22  63.64%  81.82%
    ## >Cluster_1015.mature::Pocillopora_meandrina_HIv1___Sc0000002:10845744-10845765(-)    Ptuh_lncRNA_008 159.00  -22.47  2 21    290 314 22  63.64%  81.82%
    ## >Cluster_1015.mature::Pocillopora_meandrina_HIv1___Sc0000002:10845744-10845765(-)    Ptuh_lncRNA_009 159.00  -22.47  2 21    1902 1926   22  63.64%  81.82%
    ## >Cluster_1015.mature::Pocillopora_meandrina_HIv1___Sc0000002:10845744-10845765(-)    Ptuh_lncRNA_009 159.00  -22.47  2 21    54 78   22  63.64%  81.82%
    ## >Cluster_1015.mature::Pocillopora_meandrina_HIv1___Sc0000002:10845744-10845765(-)    Ptuh_lncRNA_010 154.00  -21.59  2 20    3803 3827   21  61.90%  80.95%

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
awk -F'\t' '$7 >= 21' ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh-miRanda-lncRNA-strict-parsed.txt | awk -F'\t' '$8 >= 85' | wc -l
echo ""
echo "check some:"
awk -F'\t' '$7 >= 21' ../output/14-Ptuh-miRNA-lncRNA-BLASTs-miRanda/Ptuh-miRanda-lncRNA-strict-parsed.txt | awk -F'\t' '$8 >= 85' | head -5
```

    ## number of putative interactions of at least 21 nucleotides, with at most 3 mismatches
    ## 9
    ## 
    ## check some:
    ## >Cluster_2859.mature::Pocillopora_meandrina_HIv1___Sc0000008:5387760-5387781(-)  Ptuh_lncRNA_5192    179.00  -22.10  2 21    438 461 21  85.71%  85.71%
    ## >Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)  Ptuh_lncRNA_14188   174.00  -20.99  2 21    18351 18374 21  85.71%  85.71%
    ## >Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)  Ptuh_lncRNA_14214   174.00  -20.99  2 21    43029 43052 21  85.71%  85.71%
    ## >Cluster_3392.mature::Pocillopora_meandrina_HIv1___Sc0000010:8012938-8012959(+)  Ptuh_lncRNA_14216   174.00  -20.99  2 21    27367 27390 21  85.71%  85.71%
    ## >Cluster_4118.mature::Pocillopora_meandrina_HIv1___Sc0000014:9365481-9365502(-)  Ptuh_lncRNA_13362   174.00  -22.17  2 21    1037 1060   21  85.71%  85.71%

This is a dramatically smaller number – only 9 interactions are at least
115 nucleotides with \<=3 mismatches

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
