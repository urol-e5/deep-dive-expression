23-Apul-miRNA-methylation-machinery
================
Kathleen Durkin
2025-01-27

- <a href="#01-which-genes-are-involved-in-dna-methylation"
  id="toc-01-which-genes-are-involved-in-dna-methylation">0.1 Which genes
  are involved in DNA methylation?</a>

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

One way miRNAs and DNA methylation can interact to influence expression
is through miRNAs influencing the transcription/translation of protein
machinery involved in DNA methylation. In humans, this has been
identified as a potential mechanism of carcinogenisis! (Reviewed in
Karimzadeh et al. (2021))

## 0.1 Which genes are involved in DNA methylation?

“DNA methylation machinery is composed of several enzymes to accurately
regulate the gene expression together. The main enzyme-encoding genes
comprise a small family of DNA cytosine-5 methyltransferase, which are
involved in established (DNA methyltransferase 3 alpha (DNMT3A) and DNA
methyltransferase 3 beta (DNMT3B)) and maintenance (DNMT1) of
methylation patterns. The other genes are methyl-CpG-binding proteins
(MBPs; methyl-CpG-binding domain (MBD)1–4, methyl-CpG-binding protein 2
(MECP2), UHRF1/2, ZBTB33, ZBTB4, and ZBTB38) that assist in expression
repression and 10–11 translocation (TET) enzymes (TET1, TET2, and TET3);
seemingly, are implicated in DNA demethylation $$10, 11$$. Some
experiments demonstrated that, based on the importance of oncogenes
and/or tumor suppressor genes alterations in the development of each
type of cancer, DNA methylation regulator genes may enhance or suppress
cancer development in each specific type of cancer $$12,13,14,15,16$$.”
Karimzadeh et al. (2021)

So we’re looking for miRNAs that target genes encoding: - DNA
methyltransferases (DNMTs) - methyl-CpG-binding proteins - 10–11
translocation (TET) enzymes

First let’s see what a quick search of our annotated A.pulchra genome
turns up. I want to find any annotations related to DNA
methyltransferase

``` bash

awk -F'\t' '$5 ~ /methyltransferase/ {print $2, $5}' OFS='\t' \
../output/02-Apul-reference-annotation/Apulcra-genome-mRNA-IDmapping-2024_12_12.tab \
| grep "DNA"
```

    ## "ntLink_8:23126215-23133100" "DNA methyltransferase 1-associated protein 1 (DNMAP1) (DNMT1-associated protein 1)"
    ## "ptg000001l:8684194-8690468" "Alkylated DNA repair protein alkB homolog 8 (Probable alpha-ketoglutarate-dependent dioxygenase ABH8) (S-adenosyl-L-methionine-dependent tRNA methyltransferase ABH8) (tRNA (carboxymethyluridine(34)-5-O)-methyltransferase ABH8) (EC 2.1.1.-, EC 2.1.1.229)"
    ## "ptg000001l:8684194-8690468" "Alkylated DNA repair protein alkB homolog 8 (Probable alpha-ketoglutarate-dependent dioxygenase ABH8) (S-adenosyl-L-methionine-dependent tRNA methyltransferase ABH8) (tRNA (carboxymethyluridine(34)-5-O)-methyltransferase ABH8) (EC 2.1.1.-, EC 2.1.1.229)"
    ## "ptg000023l:4341384-4351489" "N(6)-adenine-specific methyltransferase METTL4 (Methyltransferase-like protein 4) (N(6)-adenine-specific DNA methyltransferase METTL4) (EC 2.1.1.72) (snRNA (2'-O-methyladenosine-N(6)-)-methyltransferase METTL4) (EC 2.1.1.-)"
    ## "ptg000023l:30502755-30503307"   "Methylated-DNA--protein-cysteine methyltransferase (EC 2.1.1.63) (6-O-methylguanine-DNA methyltransferase) (MGMT) (O-6-methylguanine-DNA-alkyltransferase)"
    ## "ptg000047l:7043865-7069644" "DNA (cytosine-5)-methyltransferase 1 (Dnmt1) (EC 2.1.1.37) (DNA methyltransferase GgaI) (DNA MTase GgaI) (M.GgaI) (MCMT)"
    ## "ptg000047l:7043865-7069644" "DNA (cytosine-5)-methyltransferase 1 (Dnmt1) (EC 2.1.1.37) (DNA methyltransferase GgaI) (DNA MTase GgaI) (M.GgaI) (MCMT)"
    ## "ptg000047l:7043865-7069644" "DNA (cytosine-5)-methyltransferase 1 (Dnmt1) (EC 2.1.1.37) (DNA methyltransferase GgaI) (DNA MTase GgaI) (M.GgaI) (MCMT)"

There’s definitely DNMT machinery genes anontated in our A.pulchra
genome (as expected)! Now let’s create a data frame that associated each
annotated mRNA with it’s gene ID (i.e. FUN######) and do a more broad
search for different machinery involved in methylation

``` r
# Read in data
mRNA_annot <- read.table("../output/02-Apul-reference-annotation/Apulcra-genome-mRNA-IDmapping-2024_12_12.tab", header=TRUE, sep='\t') %>% select(-X)
mRNA_geneIDs <- read.table("../output/15-Apul-annotate-UTRs/Apul-mRNA-FUNids.txt", header=FALSE, sep='\t')

# Join
mRNA_annot_geneID <- left_join(mRNA_annot, mRNA_geneIDs, by="V1")
mRNA_annot_geneID$V4 <- gsub("Parent=", "", mRNA_annot_geneID$V4)

# To double-check that this worked correctly, check for NAs in the column containing gene IDs
# Every mRNA has a gene ID, so this should return FALSE
any(is.na(mRNA_annot_geneID$V3.y))
```

    ## [1] FALSE

Now let’s select genes that encode different protein machinery relevant
to DNA methylation

``` r
methyl_machinery <- subset(mRNA_annot_geneID, 
                           (grepl("methyltransferase", mRNA_annot_geneID$Protein.names) 
                            & grepl("DNA", mRNA_annot_geneID$Protein.names))
                           | grepl("DNMT", mRNA_annot_geneID$Protein.names)
                           | grepl("Dnmt", mRNA_annot_geneID$Protein.names)
                           | grepl("MBP", mRNA_annot_geneID$Protein.names)
                           | grepl("MBD", mRNA_annot_geneID$Protein.names)
                           | grepl("MECP", mRNA_annot_geneID$Protein.names)
                           | grepl("UHRF", mRNA_annot_geneID$Protein.names)
                           | grepl("ZBTB", mRNA_annot_geneID$Protein.names)
                           | grepl("TET", mRNA_annot_geneID$Protein.names)
                           )

nrow(methyl_machinery)
```

    ## [1] 23

We find 23 mRNAs in our A.pulchra genome that encode known DNA
methylation machinery!

Note that I stipulated the term “methyltransferase” must be accompanied
by “DNA” also present in the Protein.names column. I imposed this
restriction because, when I searched for “methyltransferase” alone,
methyltransferases related to histone modification and modification with
nucleotides other than cytosine came up.

Now let’s see whether any of these methylation machinery genes are
putative miRNA targets, using the output from our target-prediction
software, miRanda.

``` r
miranda_apul <- read.table("../output/18-Apul-interactions-functional-annotation/miRanda_miRNA_mRNA.txt", header=TRUE, sep='\t')

miranda_methyl_machinery <- left_join(methyl_machinery, miranda_apul, by=c("V4" = "mRNA_FUNid"))

length(na.omit(miranda_apul$miRNA_cluster))
```

    ## [1] 6109

``` r
length(unique(na.omit(miranda_apul$miRNA_cluster)))
```

    ## [1] 39

``` r
length(na.omit(miranda_methyl_machinery$miRNA_cluster))
```

    ## [1] 0

``` r
length(unique(na.omit(miranda_methyl_machinery$miRNA_cluster)))
```

    ## [1] 0

Miranda originally identified 6109 putative targets of 39 of our miRNAs,
but none of these targets are included in our list of DNA methylation
machinery mRNA. :(

So this means that our **A.pulchra miRNAs do not putatively bind to any
methylation machinery mRNA** and/or I’m missing some genes in my search.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-karimzadeh_regulation_2021" class="csl-entry">

Karimzadeh, Mohammad Reza, Peyman Pourdavoud, Naeim Ehtesham, Mohaddese
Qadbeigi, Masood Movahedi Asl, Behrang Alani, Meysam Mosallaei, and
Bahram Pakzad. 2021. “Regulation of DNA Methylation Machinery by
Epi-<span class="nocase">miRNAs</span> in Human Cancer: Emerging New
Targets in Cancer Therapy.” *Cancer Gene Therapy* 28 (3): 157–74.
<https://doi.org/10.1038/s41417-020-00210-7>.

</div>

</div>
