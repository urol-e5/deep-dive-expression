26-Apul-siRNA-mRNA-interactions
================
Kathleen Durkin
2024-02-27

- <a href="#1-obtain-sirna-fasta" id="toc-1-obtain-sirna-fasta">1 Obtain
  siRNA fasta</a>
  - <a href="#11-extract-sirna-regions-and-convert-to-bed"
    id="toc-11-extract-sirna-regions-and-convert-to-bed">1.1 Extract siRNA
    regions and convert to .bed</a>
  - <a href="#12-convert-sample-srnaseq-bam-files-to-bed"
    id="toc-12-convert-sample-srnaseq-bam-files-to-bed">1.2 Convert sample
    sRNAseq bam files to .bed</a>
  - <a href="#13-use-bedtools-multiintersect"
    id="toc-13-use-bedtools-multiintersect">1.3 Use bedtools
    <code>multiintersect</code></a>
  - <a href="#14-annotate-with-cluster-information-and-filter"
    id="toc-14-annotate-with-cluster-information-and-filter">1.4 Annotate
    with cluster information and filter</a>
- <a href="#2-prep-for-blasts" id="toc-2-prep-for-blasts">2 Prep for
  BLASTs</a>
  - <a href="#21-check-sirna-lengths" id="toc-21-check-sirna-lengths">2.1
    Check siRNA lengths</a>
- <a href="#3-blasts" id="toc-3-blasts">3 BLASTs</a>
  - <a href="#31-make-databases" id="toc-31-make-databases">3.1 Make
    databases</a>
  - <a href="#32-run-blastn" id="toc-32-run-blastn">3.2 Run BLASTn</a>
- <a href="#4-examine-blast-tables" id="toc-4-examine-blast-tables">4
  Examine BLAST tables</a>
- <a href="#5-transposable-elements" id="toc-5-transposable-elements">5
  Transposable Elements</a>

The vast majority to work with siRNAs appears to focus on their
application in therapeutic treatments, in the field of human health.
siRNAs are designed to bind to specific locations are part of these
therapeutics. Because of this, while some tools have been developed to
predict siRNA targets (Grešová, Alexiou, and Giassa
([2022](#ref-gresova_small_2022))), these tools focus on identifying
“off-targets”, or potential interactions *other* than the siRNAs primary
target. One of these tools, siFi (Lück et al.
([2019](#ref-luck_sirna-finder_2019))) could be useful to us, but is
unfortunately only available as a windows application, which we can’t
install on our servers.I may circle back to siFi later by running it
locally, but for now we’ll just try BLAST. siRNAs require full- or
near-full complementarity to bind to their targets, which should make
searching for targets through BLAST relatively straightforward.

# 1 Obtain siRNA fasta

ShortStack outputs a gff of genomic regions with accumulating sRNAs
(`Results.gff3`) and annotates those regions as “unkown sRNA”, “miRNA”,
or “siRNA”. However, these regions are generally several hundred
basepairs long, since they are “accumulation regions”, not the \~21nt
sRNA molecules themselves. While ShortStack also identifies and outputs
fastas for the miRNAs, it doesn’t do that for siRNAs. This means we need
to ID the precise siRNA sequences ourselves.

Approach:

1.  Extract regions of accumulating siRNA from ShortStack output

2.  Convert gff of siRNA regions to .bed

3.  Convert sample sRNAseq bam files to .bed

4.  Use bedtools `multiintersect` to find sRNA sequences that appear
    within the ShortStack siRNA regions and in our sample sRNAseq reads.

## 1.1 Extract siRNA regions and convert to .bed

``` bash
awk '$3 ~ /siRNA/' ../output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out/Results.gff3 > ../output/26-Apul-siRNA-mRNA-interactions/Apul-ShortStack-siRNA_regions.gff3

awk '$3 ~ /siRNA/ {print $1 "\t" $4-1 "\t" $5}' ../output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out/Results.gff3 > ../output/26-Apul-siRNA-mRNA-interactions/Apul-ShortStack-siRNA_regions.bed
```

How many siRNA regions are there? This should be our total number of
siRNA at the end of this process

``` bash
wc -l ../output/26-Apul-siRNA-mRNA-interactions/Apul-ShortStack-siRNA_regions.bed
```

    ## 142 ../output/26-Apul-siRNA-mRNA-interactions/Apul-ShortStack-siRNA_regions.bed

## 1.2 Convert sample sRNAseq bam files to .bed

``` bash
reads_dir="../output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out"
output_dir="../output/26-Apul-siRNA-mRNA-interactions"
/home/shared/bedtools2/bin/bamToBed \
-i ${reads_dir}/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam \
> ${output_dir}/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bed

/home/shared/bedtools2/bin/bamToBed \
-i ${reads_dir}/sRNA-ACR-145-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam \
> ${output_dir}/sRNA-ACR-145-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bed

/home/shared/bedtools2/bin/bamToBed \
-i ${reads_dir}/sRNA-ACR-150-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam \
> ${output_dir}/sRNA-ACR-150-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bed

/home/shared/bedtools2/bin/bamToBed \
-i ${reads_dir}/sRNA-ACR-173-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam \
> ${output_dir}/sRNA-ACR-173-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bed

/home/shared/bedtools2/bin/bamToBed \
-i ${reads_dir}/sRNA-ACR-178-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam \
> ${output_dir}/sRNA-ACR-178-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bed
```

## 1.3 Use bedtools `multiintersect`

``` bash
output_dir="../output/26-Apul-siRNA-mRNA-interactions"

/home/shared/bedtools2/bin/bedtools multiinter -i \
${output_dir}/Apul-ShortStack-siRNA_regions.bed \
${output_dir}/sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bed \
${output_dir}/sRNA-ACR-145-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bed \
${output_dir}/sRNA-ACR-150-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bed \
${output_dir}/sRNA-ACR-173-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bed \
${output_dir}/sRNA-ACR-178-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bed \
-header \
> ${output_dir}/Apul-shared-siRNA_clusters.bed
```

How many outputs?

``` bash
wc -l ../output/26-Apul-siRNA-mRNA-interactions/Apul-shared-siRNA_clusters.bed
```

    ## 2899994 ../output/26-Apul-siRNA-mRNA-interactions/Apul-shared-siRNA_clusters.bed

Ignoring the header line, there are 2899993 intersecting regions.

First, we’ll only keep regions identified as siRNA regions by ShortStack
(present in the first input file) and present in at least two of the
samples

``` bash
awk '$4 >= 3 && $6 > 0' ../output/26-Apul-siRNA-mRNA-interactions/Apul-shared-siRNA_clusters.bed > ../output/26-Apul-siRNA-mRNA-interactions/Apul-shared-siRNA_clusters-filtered.bed

wc -l ../output/26-Apul-siRNA-mRNA-interactions/Apul-shared-siRNA_clusters-filtered.bed
```

    ## 490 ../output/26-Apul-siRNA-mRNA-interactions/Apul-shared-siRNA_clusters-filtered.bed

There are more than 142 shared clusters here, so most of these are not
par of the true siRNA. Let’s use intersectBed to annotate each shared
segment with the ShortStack siRNA accumulation region it falls within.

``` bash
/home/shared/bedtools2/bin/bedtools intersect -wa -wb -a ../output/26-Apul-siRNA-mRNA-interactions/Apul-shared-siRNA_clusters-filtered.bed -b ../output/26-Apul-siRNA-mRNA-interactions/Apul-ShortStack-siRNA_regions.bed > ../output/26-Apul-siRNA-mRNA-interactions/Apul-shared-siRNA_clusters-annotated.bed
```

## 1.4 Annotate with cluster information and filter

Now read in to R and annotate with additional information on each siRNA
region (provided in the ShortStack gff3)

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

# Read in data
intersectBed <- read.table("../output/26-Apul-siRNA-mRNA-interactions/Apul-shared-siRNA_clusters-annotated.bed", sep="\t", header=FALSE)
siRNA_gff3 <- read.table("../output/26-Apul-siRNA-mRNA-interactions/Apul-ShortStack-siRNA_regions.gff3", sep="\t", header=FALSE)

# Create a column that contains full genomic location info
intersectBed$region <- paste(intersectBed$V12, "_", intersectBed$V13, "_", intersectBed$V14)
# For siRNA_gff, need to first convert genomic loction columns to bed format
siRNA_gff3$V4 <-siRNA_gff3$V4 - 1
siRNA_gff3$region <- paste(siRNA_gff3$V1, "_", siRNA_gff3$V4, "_", siRNA_gff3$V5)
# Remove unnecessary columns and rename columns
colnames(siRNA_gff3) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "region")

# Annotate intersectBed results with the gff info
intersectBed <- left_join(intersectBed, siRNA_gff3, by="region")

# Extract cluster and dicer call info from "attributes" column
intersectBed <- intersectBed %>% separate(attributes, into = c("ID", "DicerCall", "MIRNA"), sep = ";") %>%
  mutate(ID = gsub("ID=Cluster_", "Cluster_", ID), 
         DicerCall = gsub("DicerCall=", "", DicerCall)) %>%
  select(-MIRNA)

# Calculate lengths of intersecting regions
intersectBed$length <- intersectBed$V3 - intersectBed$V2
```

``` r
# Keep only intersections with lengths of at least the lowest DicerCall value (21)
intersectBed_filt <- intersectBed[intersectBed$length >= 21,]

# Remove intersections with lengths >=30. These may be menaingful (e.g. a precursor), but no siRNA should be longer than ~25nt.
intersectBed_filt <- intersectBed_filt[intersectBed_filt$length <= 30,]

nrow(intersectBed)
```

    ## [1] 490

``` r
nrow(intersectBed_filt)
```

    ## [1] 155

``` r
length(unique(intersectBed$ID))
```

    ## [1] 130

``` r
length(unique(intersectBed_filt$ID))
```

    ## [1] 116

We’ve filtered 490 intersecting regions to 155, and from coverage of 130
ShortStack siRNA “Clusters” to 116 (out of the 142 output by ShortStack)

In cases where there are multiple long intersections within a single
ShortStack siRNA region, I can’t think of a robust way to decide which
intersecting region is more likely to constitute a genuine siRNA.
Instead, I’ll just label both as different siRNAs of the same cluster
for now.

``` r
# Initialize updated_ID col
intersectBed_filt$updated_ID <- intersectBed_filt$ID

# Append _a, _b, etc. to the IDs that are repeated
for(i in 1:nrow(intersectBed_filt)) {
  # Find all rows with the same ID
  matching_rows <- which(intersectBed_filt$ID == intersectBed_filt$ID[i])
  
  # If the ID is repeated, append _a, _b, etc.
  if(length(matching_rows) > 1) {
    # Find the position of the current row within the matching rows
    position <- which(matching_rows == i)
    
    # Append the corresponding letter (_a, _b, etc.) to the ID
    intersectBed_filt$updated_ID[i] <- paste0(intersectBed_filt$ID[i], "_", letters[position])
  }
}

length(unique(intersectBed_filt$updated_ID))
```

    ## [1] 155

When accounting for multiple unique putative siRNA loci in a single
siRNA region, there are 155 siRNA loci (in 116 siRNA-accumulating
regions)

Ok, we now have a set of putative siRNA loci. Let’s save them as a gff
and extract fastas for them

``` r
intersectBed_filt$attributes <- paste0("siRNALociID=", intersectBed_filt$updated_ID, ";DicerCall=", intersectBed_filt$DicerCall, ";Length=", intersectBed_filt$length)
intersectBed_filt_gff <- intersectBed_filt %>% select(V1, source, type, V2, V3, score, strand, phase, attributes)

#convert genomic location to gff3 format
intersectBed_filt_gff$V2 <- intersectBed_filt_gff$V2 + 1

write.table(intersectBed_filt_gff, "../output/26-Apul-siRNA-mRNA-interactions/Apul-siRNA_loci.gff3", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

``` bash
/home/shared/bedtools2/bin/bedtools getfasta -fi ../data/Apulchra-genome.fa -bed ../output/26-Apul-siRNA-mRNA-interactions/Apul-siRNA_loci.gff3 -fo ../output/26-Apul-siRNA-mRNA-interactions/Apul-siRNA_loci.fa
```

We have a fasta file of putative siRNA loci! Now we can move on to
searching for siRNA targets

# 2 Prep for BLASTs

## 2.1 Check siRNA lengths

``` bash
awk '/^>/ {if (seqlen){print seqlen}; printf $0" " ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' ../output/26-Apul-siRNA-mRNA-interactions/Apul-siRNA_loci.fa > ../output/26-Apul-siRNA-mRNA-interactions/Apul-siRNA_loci-lengths.txt
```

``` r
# Summary stats of precursor and mature lengths

lengths <- read.table("../output/26-Apul-siRNA-mRNA-interactions/Apul-siRNA_loci-lengths.txt", sep = " ", header = FALSE, col.names = c("seqID", "length"))

cat("Average siRNA length: ", mean(lengths$length))
```

    ## Average siRNA length:  23.45161

``` r
cat("\n")
```

``` r
cat("Range of siRNA lengths: ", range(lengths$length))
```

    ## Range of siRNA lengths:  21 30

# 3 BLASTs

## 3.1 Make databases

Database of siRNAs:

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../output/26-Apul-siRNA-mRNA-interactions/Apul-siRNA_loci.fa \
-dbtype nucl \
-out ../output/26-Apul-siRNA-mRNA-interactions/blasts/Apul-db/Apul-siRNA_loci
```

## 3.2 Run BLASTn

Generate a list of blast results. It’s possible that a single siRNA
could bind to multiple places in the genome so I will allow up to 20
hits. Since siRNAs bind with high complementarity, I set minimum percent
identity to 85%. I’ll also include the “-word_size 4” option, which
reduces the required length of the initial match

``` bash
/home/shared/ncbi-blast-2.11.0+/bin/blastn \
-task blastn \
-query ../data/Apulchra-genome.fa \
-db ../output/26-Apul-siRNA-mRNA-interactions/blasts/Apul-db/Apul-siRNA_loci \
-out ../output/26-Apul-siRNA-mRNA-interactions/blasts/siRNA_to_genome_blastn.tab \
-num_threads 40 \
-word_size 4 \
-perc_identity 85 \
-max_target_seqs 20 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/26-Apul-siRNA-mRNA-interactions/blasts/siRNA_to_genome_blastn.tab
```

    ## 2201 ../output/26-Apul-siRNA-mRNA-interactions/blasts/siRNA_to_genome_blastn.tab

# 4 Examine BLAST tables

Read into R and assign informative column labels

``` r
siRNA_genome_BLASTn <- read.table("../output/26-Apul-siRNA-mRNA-interactions/blasts/siRNA_to_genome_blastn.tab", sep="\t", header=FALSE)

colnames(siRNA_genome_BLASTn) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
```

``` r
# Remove any matches of <21bp (the length of our shortest siRNA)
siRNA_genome_BLASTn <- siRNA_genome_BLASTn[siRNA_genome_BLASTn$length >= 21,]

# Remove instances of an siRNA aligning to itsself
siRNA_genome_BLASTn$genome_loc <- paste0(siRNA_genome_BLASTn$qseqid, ":", siRNA_genome_BLASTn$qstart-1, "-", siRNA_genome_BLASTn$qend)
siRNA_genome_BLASTn <- siRNA_genome_BLASTn[siRNA_genome_BLASTn$sseqid != siRNA_genome_BLASTn$genome_loc,]

nrow(siRNA_genome_BLASTn)
```

    ## [1] 829

``` r
length(unique(siRNA_genome_BLASTn$sseqid))
```

    ## [1] 89

``` r
length(unique(siRNA_genome_BLASTn$genome_loc))
```

    ## [1] 724

Of the 155 siRNA loci, 89 align to a 724 different locations in the
genome with very high complementarity. Some genome locations match
multiple siRNA

# 5 Transposable Elements

TE annotation for A.pulchra stored in project OSF:
<https://osf.io/y8963/files/osfstorage>

``` bash
awk '$10 ~ /transposon/' ../data/apul.hifiasm.s55_pa.p_ctg.fa.k32.w100.z1000.ntLink.5rounds.fa.out | wc -l
```

    ## 0

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-gresova_small_2022" class="csl-entry">

Grešová, Katarína, Panagiotis Alexiou, and Ilektra-Chara Giassa. 2022.
“Small RNA Targets: Advances in Prediction Tools and High-Throughput
Profiling.” *Biology* 11 (12): 1798.
<https://doi.org/10.3390/biology11121798>.

</div>

<div id="ref-luck_sirna-finder_2019" class="csl-entry">

Lück, Stefanie, Tino Kreszies, Marc Strickert, Patrick Schweizer, Markus
Kuhlmann, and Dimitar Douchkov. 2019. “<span
class="nocase">siRNA</span>-Finder (Si-Fi) Software for RNAi-Target
Design and Off-Target Prediction.” *Frontiers in Plant Science* 10
(August). <https://doi.org/10.3389/fpls.2019.01023>.

</div>

</div>
