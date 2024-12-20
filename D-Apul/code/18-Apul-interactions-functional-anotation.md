18-Apul-interactions-functional-annotation
================
Kathleen Durkin
2024-12-11

- <a href="#1-mirna-mrna-interactions"
  id="toc-1-mirna-mrna-interactions">1 miRNA-mRNA interactions</a>
  - <a href="#11-rnahybrid-output" id="toc-11-rnahybrid-output">1.1
    RNAhybrid output</a>
  - <a href="#12-jills-interaction-plot-stuff"
    id="toc-12-jills-interaction-plot-stuff">1.2 Jill’s interaction plot
    stuff</a>
- <a
  href="#2-how-do-miranda-and-rnahybrid-compare-for-mirna-mrna-target-prediction"
  id="toc-2-how-do-miranda-and-rnahybrid-compare-for-mirna-mrna-target-prediction">2
  How do miRanda and RNAhybrid compare for miRNA-mRNA target
  prediction?</a>

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
```

Perform functional annotation for miRNA-mRNA, lncRNA-mRNA, lncRNA-miRNA,
etc. putative interactions.

The reference genome was annotated using Uniprot/Swissprot in
`deep-dive-expression/D-Apul/code/02-Apul-reference-annotation`.

Load in reference annotations mapping

``` r
genome_IDmapping <- read.table("../output/02-Apul-reference-annotation/Apulcra-genome-mRNA-IDmapping-2024_12_12.tab", header=TRUE)
```

# 1 miRNA-mRNA interactions

## 1.1 RNAhybrid output

``` r
# Read in RNAhybrid results for miRNAs that bind in the 3'UTR region of an mRNA
RNAhybrid_3UTR <- read.table("../output/16-Apul-RNAhybrid/Apul-RNAhybrid-mRNA-compact_3utr_worm-formatted.txt", sep="\t", header=TRUE)

# Filter to only retain highly likely hybridizations
RNAhybrid_3UTR_p0.01 <- RNAhybrid_3UTR %>%
  filter(pval < 0.01)

# join with ID mapping to annotate
RNAhybrid_3UTR_p0.01_FA <- left_join(RNAhybrid_3UTR_p0.01, genome_IDmapping, by = c("target_name" = "V3"))
```

## 1.2 Jill’s interaction plot stuff

``` r
# read in data
miRanda_pcor <-  read.csv("../output/09-Apul-mRNA-miRNA-interactions/PCC_miRNA_mRNA.csv")

# read in FUN id - mRNA association table
mRNA_FUN_table <- read.table("../output/15-Apul-annotate-UTRs/Apul-mRNA-FUNids.txt")
# Remove "Parent=" prefix of the FUN IDs
mRNA_FUN_table$V4 <- gsub("Parent=", "", mRNA_FUN_table$V4)

# read in functional annotation mapping table and ensure unique rows
annot_tab <- read.table("../output/02-Apul-reference-annotation/Apulcra-genome-mRNA-IDmapping-2024_12_12.tab", header=TRUE) %>% 
  distinct(V1, .keep_all = TRUE)
```

``` r
# Use FUN ids to annotate with associated mRNA coordinates
miRanda_pcor_annot <- left_join(miRanda_pcor, mRNA_FUN_table, by = c("mRNA" = "V4")) %>%
  select(-V2, -V3, -V5)

# Use these mRNA coordinates to map to functional annotations
miRanda_pcor_annot <- left_join(miRanda_pcor_annot, annot_tab, by = c("V1" = "V1"))
```

# 2 How do miRanda and RNAhybrid compare for miRNA-mRNA target prediction?

``` r
# Read in the data
miRanda_miRNA_mRNA <- read.table("../data/09-Apul-mRNA-miRNA-interactions/miranda_strict_all_1kb_parsed_apul_updated.txt")
colnames(miRanda_miRNA_mRNA) <- c("mirna", "target",  "score", "energy", "miRNA_start", "miRNA_end", "target_start", "target_end", "aln_length", "subject_identity", "Qquery_identity")

# Separate the miRNA cluster names and locations
miRanda_miRNA_mRNA <- miRanda_miRNA_mRNA %>% separate(mirna, into = c("miRNA_cluster", "miRNA_location"), sep = "::")
miRanda_miRNA_mRNA$miRNA_cluster <- gsub("^>", "", miRanda_miRNA_mRNA$miRNA_cluster)  # Remove leading >
miRanda_miRNA_mRNA$miRNA_cluster <- gsub("\\.mature$", "", miRanda_miRNA_mRNA$miRNA_cluster)  # Remove trailing .mature

# Separate the mRNA FUN ids and locations
miRanda_miRNA_mRNA <- miRanda_miRNA_mRNA %>% separate(target, into = c("mRNA_FUNid", "mRNA_location"), sep = "::")

# Check
head(miRanda_miRNA_mRNA)
```

    ##   miRNA_cluster                miRNA_location mRNA_FUNid
    ## 1 Cluster_10051 ptg000016l:7795530-7795551(+) FUN_000580
    ## 2 Cluster_10051 ptg000016l:7795530-7795551(+) FUN_000640
    ## 3 Cluster_10051 ptg000016l:7795530-7795551(+) FUN_001049
    ## 4 Cluster_10051 ptg000016l:7795530-7795551(+) FUN_001564
    ## 5 Cluster_10051 ptg000016l:7795530-7795551(+) FUN_002479
    ## 6 Cluster_10051 ptg000016l:7795530-7795551(+) FUN_002556
    ##                mRNA_location score energy miRNA_start miRNA_end target_start
    ## 1   ntLink_6:4580655-4581655   162 -20.38           2        21           11
    ## 2   ntLink_6:5278186-5279186   156 -20.16           2        17          132
    ## 3   ntLink_6:7885174-7886174   169 -23.61           2        20          195
    ## 4 ntLink_6:13794539-13795539   163 -21.78           2        16          237
    ## 5   ntLink_8:2255969-2256969   166 -24.91           2        21          413
    ## 6   ntLink_8:3610299-3611299   159 -22.63           2        21          537
    ##   target_end aln_length subject_identity Qquery_identity
    ## 1         31         19           73.68%          84.21%
    ## 2        153         15           73.33%          86.67%
    ## 3        218         20           80.00%          90.00%
    ## 4        258         14           85.71%          92.86%
    ## 5        433         19           73.68%          89.47%
    ## 6        559         20           65.00%          85.00%

``` r
write.table(miRanda_miRNA_mRNA, "../output/18-Apul-interactions-functional-annotation/miRanda_miRNA_mRNA.txt", col.names = TRUE, row.names = FALSE,sep = "\t")
```

``` r
# Read in the data
RNAhybrid_miRNA_mRNA <- read.table("../output/16-Apul-RNAhybrid/Apul-RNAhybrid-3UTR-compact_3utr_worm-formatted.txt", sep="\t", header=TRUE)

# Separate the miRNA cluster names and locations
RNAhybrid_miRNA_mRNA <- RNAhybrid_miRNA_mRNA %>% separate(query_name, into = c("miRNA_cluster", "miRNA_location"), sep = "::")
RNAhybrid_miRNA_mRNA$miRNA_cluster <- gsub("\\.mature$", "", RNAhybrid_miRNA_mRNA$miRNA_cluster)  # Remove trailing .mature

write.table(RNAhybrid_miRNA_mRNA, "../output/18-Apul-interactions-functional-annotation/RNAhybrid_miRNA_mRNA.txt", col.names = TRUE, row.names = FALSE,sep = "\t")
```

Both miRanda and RNAhybrid were run using 3’UTR regions as the target.
Since there may be some formattign differences in identifying the 3UTr
regions, I don’t think a left_join will work, because the target regions
won’t match exactly.

Instead, maybe I can use bedtools intersect? Intersect will identify
regions where both the miRanda *and* RNAhybrid outputs show a
target-miRNA hybridization

Actually first, let’s try to look at the two in IGV

Make IGV inputs (super basic gff)

``` bash
# Format miRanda for IGV as a gff
#!/bin/bash

# Input and output file paths
INPUT_FILE="../output/18-Apul-interactions-functional-annotation/miRanda_miRNA_mRNA.txt"  # Replace with the path to your input file
OUTPUT_FILE="../output/18-Apul-interactions-functional-annotation/miranda_strict_all_1kb_parsed_apul_updated.gff"

# Write GFF3 header to the output file
echo "##gff-version 3" > "$OUTPUT_FILE"

# Process the input file, skipping the header line
tail -n +2 "$INPUT_FILE" | while IFS=$'\t' read -r miRNA_cluster    miRNA_location  mRNA_FUNid  mRNA_location   score   energy  miRNA_start miRNA_end   target_start    target_end  aln_length  subject_identity    Qquery_identity
do
  # Extract locus name and coordinates from target_name
  locus=$(echo "$mRNA_location" | cut -d'"' -f2 | cut -d':' -f1)
  start_coord=$(echo "$mRNA_location" | cut -d':' -f2 | cut -d'-' -f1)
  start_gff=$((start_coord + target_start))
  end_gff=$((start_gff + aln_length))

  # Extract strandedness from query_name
  strand=$(echo "$miRNA_location" | grep -o '(-\|+)' | tr -d '()')

  # Write the GFF3 line
  echo -e "$locus\tRNAhybrid\tmiRNA_binding\t$start_gff\t$end_gff\t.\t$strand\t.\tID=$miRNA_cluster;energy=$energy;score=$score" >> "$OUTPUT_FILE"
done
```

``` bash
/home/shared/bedtools2/bin/bedtools intersect \
-a ../output/18-Apul-interactions-functional-annotation/miranda_strict_all_1kb_parsed_apul_updated.gff \
-b ../output/16-Apul-RNAhybrid/Apul-RNAhybrid-3UTR-compact_3utr_worm.gff \
-wa


/home/shared/bedtools2/bin/bedtools intersect \
-a ../output/16-Apul-RNAhybrid/Apul-RNAhybrid-3UTR-compact_3utr_worm.gff \
-b ../output/18-Apul-interactions-functional-annotation/miranda_strict_all_1kb_parsed_apul_updated.gff \
-wa
```

    ## ptg000019l   RNAhybrid   miRNA_binding   3905988 3906006 .   +   .   ID="Cluster_10051";energy=-20.49;score=157
    ## ptg000039l   RNAhybrid   miRNA_binding   1038025 1038045 .   +   .   ID="Cluster_10051";energy=-27.39;score=167
    ## ntLink_8 RNAhybrid   miRNA_binding   17362798    17362810    .   -   .   ID="Cluster_10057";energy=-25.09;score=161
    ## ntLink_8 RNAhybrid   miRNA_binding   17362798    17362810    .   -   .   ID="Cluster_10057";energy=-25.09;score=161
    ## ntLink_8 RNAhybrid   miRNA_binding   31354508    31354525    .   -   .   ID="Cluster_10057";energy=-23.81;score=168
    ## ptg000002l   RNAhybrid   miRNA_binding   3411163 3411182 .   -   .   ID="Cluster_10057";energy=-23.78;score=166
    ## ptg000008l   RNAhybrid   miRNA_binding   6343721 6343733 .   -   .   ID="Cluster_10057";energy=-26.72;score=161
    ## ptg000011l   RNAhybrid   miRNA_binding   11819542    11819559    .   -   .   ID="Cluster_10057";energy=-23.81;score=168
    ## ptg000015l   RNAhybrid   miRNA_binding   2404071 2404083 .   -   .   ID="Cluster_10057";energy=-25.78;score=165
    ## ptg000015l   RNAhybrid   miRNA_binding   2404071 2404083 .   -   .   ID="Cluster_10057";energy=-25.78;score=165
    ## ptg000020l   RNAhybrid   miRNA_binding   16091584    16091603    .   -   .   ID="Cluster_10057";energy=-21.48;score=156
    ## ptg000020l   RNAhybrid   miRNA_binding   16493595    16493614    .   -   .   ID="Cluster_10057";energy=-24.36;score=152
    ## ptg000023l   RNAhybrid   miRNA_binding   24634237    24634258    .   -   .   ID="Cluster_10057";energy=-21.92;score=163
    ## ptg000025l   RNAhybrid   miRNA_binding   1210277 1210296 .   -   .   ID="Cluster_10057";energy=-21.68;score=152
    ## ptg000026l   RNAhybrid   miRNA_binding   556244  556256  .   -   .   ID="Cluster_10057";energy=-26.71;score=165
    ## ptg000026l   RNAhybrid   miRNA_binding   556244  556256  .   -   .   ID="Cluster_10057";energy=-26.71;score=165
    ## ptg000026l   RNAhybrid   miRNA_binding   896120  896132  .   -   .   ID="Cluster_10057";energy=-26.71;score=165
    ## ptg000026l   RNAhybrid   miRNA_binding   896120  896132  .   -   .   ID="Cluster_10057";energy=-26.71;score=165
    ## ptg000036l   RNAhybrid   miRNA_binding   859237  859254  .   -   .   ID="Cluster_10057";energy=-24.85;score=163
    ## ptg000047l   RNAhybrid   miRNA_binding   2819153 2819168 .   -   .   ID="Cluster_10093";energy=-26.08;score=168
    ## ptg000047l   RNAhybrid   miRNA_binding   2819153 2819168 .   -   .   ID="Cluster_10093";energy=-26.08;score=168
    ## ptg000047l   RNAhybrid   miRNA_binding   2819153 2819168 .   -   .   ID="Cluster_10093";energy=-26.08;score=168
    ## ptg000047l   RNAhybrid   miRNA_binding   2819153 2819168 .   -   .   ID="Cluster_10093";energy=-26.08;score=168
    ## ptg000021l   RNAhybrid   miRNA_binding   13553488    13553506    .   +   .   ID="Cluster_14768";energy=-39.74;score=191
    ## ntLink_6 RNAhybrid   miRNA_binding   15902887    15902912    .   +   .   ID="Cluster_15316";energy=-24.16;score=155
    ## ptg000034l   RNAhybrid   miRNA_binding   183625  183643  .   +   .   ID="Cluster_15316";energy=-27.08;score=179
    ## ptg000064l   RNAhybrid   miRNA_binding   99588   99606   .   +   .   ID="Cluster_15316";energy=-27.08;score=179
    ## ptg000009l   RNAhybrid   miRNA_binding   10224585    10224598    .   -   .   ID="Cluster_15340";energy=-22.98;score=170
    ## ptg000025l   RNAhybrid   miRNA_binding   5797264 5797281 .   -   .   ID="Cluster_15340";energy=-30.09;score=182
    ## ptg000025l   RNAhybrid   miRNA_binding   5797264 5797281 .   -   .   ID="Cluster_15340";energy=-30.09;score=182
    ## ptg000025l   RNAhybrid   miRNA_binding   5797264 5797281 .   -   .   ID="Cluster_15340";energy=-30.09;score=182
    ## ptg000025l   RNAhybrid   miRNA_binding   5797264 5797281 .   -   .   ID="Cluster_15340";energy=-30.09;score=182
    ## ptg000012l   RNAhybrid   miRNA_binding   17288350    17288364    .   +   .   ID="Cluster_15851";energy=-23.52;score=175
    ## ptg000008l   RNAhybrid   miRNA_binding   8603023 8603039 .   -   .   ID="Cluster_15854";energy=-24.28;score=181
    ## ptg000008l   RNAhybrid   miRNA_binding   8603023 8603039 .   -   .   ID="Cluster_15854";energy=-24.28;score=181
    ## ptg000008l   RNAhybrid   miRNA_binding   8603023 8603039 .   -   .   ID="Cluster_15854";energy=-24.28;score=181
    ## ptg000008l   RNAhybrid   miRNA_binding   8603023 8603039 .   -   .   ID="Cluster_15854";energy=-24.28;score=181
    ## ptg000023l   RNAhybrid   miRNA_binding   7799752 7799770 .   -   .   ID="Cluster_16409";energy=-24.64;score=177
    ## ntLink_6 RNAhybrid   miRNA_binding   4819493 4819515 .   -   .   ID="Cluster_17791";energy=-20.26;score=174
    ## ntLink_6 RNAhybrid   miRNA_binding   4819493 4819515 .   -   .   ID="Cluster_17791";energy=-20.26;score=174
    ## ntLink_6 RNAhybrid   miRNA_binding   4819493 4819515 .   -   .   ID="Cluster_17791";energy=-20.26;score=174
    ## ntLink_6 RNAhybrid   miRNA_binding   4819493 4819515 .   -   .   ID="Cluster_17791";energy=-20.26;score=174
    ## ptg000021l   RNAhybrid   miRNA_binding   13589804    13589824    .   -   .   ID="Cluster_1826";energy=-25.71;score=179
    ## ntLink_8 RNAhybrid   miRNA_binding   3046118 3046138 .   -   .   ID="Cluster_1862";energy=-20.09;score=177
    ## ptg000039l   RNAhybrid   miRNA_binding   483717  483733  .   -   .   ID="Cluster_19193";energy=-26.35;score=169
    ## ntLink_8 RNAhybrid   miRNA_binding   7813582 7813602 .   -   .   ID="Cluster_2463";energy=-25.95;score=179
    ## ptg000001l   RNAhybrid   miRNA_binding   10400679    10400699    .   -   .   ID="Cluster_2463";energy=-25.95;score=178
    ## ptg000004l   RNAhybrid   miRNA_binding   11697979    11697993    .   -   .   ID="Cluster_2463";energy=-24.38;score=171
    ## ptg000007l   RNAhybrid   miRNA_binding   9978668 9978688 .   -   .   ID="Cluster_2463";energy=-25.95;score=178
    ## ptg000018l   RNAhybrid   miRNA_binding   13812471    13812487    .   -   .   ID="Cluster_2463";energy=-24.03;score=169
    ## ptg000026l   RNAhybrid   miRNA_binding   7994104 7994124 .   -   .   ID="Cluster_2463";energy=-25.95;score=179
    ## ntLink_8 RNAhybrid   miRNA_binding   20869785    20869806    .   +   .   ID="Cluster_3366";energy=-22.64;score=169
    ## ntLink_8 RNAhybrid   miRNA_binding   20869785    20869806    .   +   .   ID="Cluster_3366";energy=-22.64;score=169
    ## ntLink_8 RNAhybrid   miRNA_binding   30064515    30064532    .   +   .   ID="Cluster_3366";energy=-25.67;score=170
    ## ntLink_8 RNAhybrid   miRNA_binding   30064515    30064532    .   +   .   ID="Cluster_3366";energy=-25.67;score=170
    ## ptg000020l   RNAhybrid   miRNA_binding   5862528 5862544 .   +   .   ID="Cluster_3366";energy=-23.94;score=169
    ## ptg000020l   RNAhybrid   miRNA_binding   5862528 5862544 .   +   .   ID="Cluster_3366";energy=-23.94;score=169
    ## ptg000024l   RNAhybrid   miRNA_binding   9175741 9175760 .   +   .   ID="Cluster_3366";energy=-21.37;score=174
    ## ptg000024l   RNAhybrid   miRNA_binding   9175741 9175760 .   +   .   ID="Cluster_3366";energy=-21.37;score=174
    ## ptg000024l   RNAhybrid   miRNA_binding   9175741 9175760 .   +   .   ID="Cluster_3366";energy=-21.37;score=174
    ## ptg000024l   RNAhybrid   miRNA_binding   9175741 9175760 .   +   .   ID="Cluster_3366";energy=-21.37;score=174
    ## ptg000036l   RNAhybrid   miRNA_binding   4523015 4523036 .   +   .   ID="Cluster_3366";energy=-24.55;score=176
    ## ntLink_8 RNAhybrid   miRNA_binding   20869785    20869806    .   +   .   ID="Cluster_3367";energy=-22.64;score=169
    ## ntLink_8 RNAhybrid   miRNA_binding   20869785    20869806    .   +   .   ID="Cluster_3367";energy=-22.64;score=169
    ## ntLink_8 RNAhybrid   miRNA_binding   30064515    30064532    .   +   .   ID="Cluster_3367";energy=-23.68;score=170
    ## ntLink_8 RNAhybrid   miRNA_binding   30064515    30064532    .   +   .   ID="Cluster_3367";energy=-23.68;score=170
    ## ptg000020l   RNAhybrid   miRNA_binding   5862528 5862544 .   +   .   ID="Cluster_3367";energy=-25.17;score=169
    ## ptg000020l   RNAhybrid   miRNA_binding   5862528 5862544 .   +   .   ID="Cluster_3367";energy=-25.17;score=169
    ## ptg000024l   RNAhybrid   miRNA_binding   9175741 9175760 .   +   .   ID="Cluster_3367";energy=-21.37;score=174
    ## ptg000024l   RNAhybrid   miRNA_binding   9175741 9175760 .   +   .   ID="Cluster_3367";energy=-21.37;score=174
    ## ptg000024l   RNAhybrid   miRNA_binding   9175741 9175760 .   +   .   ID="Cluster_3367";energy=-21.37;score=174
    ## ptg000024l   RNAhybrid   miRNA_binding   9175741 9175760 .   +   .   ID="Cluster_3367";energy=-21.37;score=174
    ## ptg000036l   RNAhybrid   miRNA_binding   4523015 4523036 .   +   .   ID="Cluster_3367";energy=-22.6;score=176
    ## ntLink_8 RNAhybrid   miRNA_binding   16931011    16931030    .   -   .   ID="Cluster_4220";energy=-30.98;score=188
    ## ntLink_8 RNAhybrid   miRNA_binding   16953818    16953837    .   -   .   ID="Cluster_4220";energy=-30.98;score=188
    ## ntLink_8 RNAhybrid   miRNA_binding   16978257    16978276    .   -   .   ID="Cluster_4220";energy=-30.98;score=188
    ## ptg000018l   RNAhybrid   miRNA_binding   1167302 1167318 .   -   .   ID="Cluster_4220";energy=-29.24;score=181
    ## ptg000035l   RNAhybrid   miRNA_binding   9011316 9011333 .   -   .   ID="Cluster_4220";energy=-30.63;score=182
    ## ntLink_6 RNAhybrid   miRNA_binding   5772709 5772728 .   +   .   ID="Cluster_4254";energy=-22.88;score=176
    ## ptg000030l   RNAhybrid   miRNA_binding   2128306 2128325 .   +   .   ID="Cluster_4254";energy=-24.43;score=178
    ## ptg000030l   RNAhybrid   miRNA_binding   2128306 2128325 .   +   .   ID="Cluster_4254";energy=-24.43;score=178
    ## ptg000030l   RNAhybrid   miRNA_binding   2128306 2128325 .   +   .   ID="Cluster_4254";energy=-24.43;score=178
    ## ptg000030l   RNAhybrid   miRNA_binding   2128306 2128325 .   +   .   ID="Cluster_4254";energy=-24.43;score=178
    ## ntLink_6 RNAhybrid   miRNA_binding   6374246 6374262 .   -   .   ID="Cluster_5012";energy=-26.78;score=185
    ## ptg000021l   RNAhybrid   miRNA_binding   16212974    16212995    .   -   .   ID="Cluster_5981";energy=-22.19;score=146
    ## ptg000021l   RNAhybrid   miRNA_binding   13589810    13589832    .   -   .   ID=Cluster_1826.mature::ntLink_6:4847465-4847486(-);MFE=-28.7;Pval=0.022695
    ## ntLink_8 RNAhybrid   miRNA_binding   7813588 7813610 .   -   .   ID=Cluster_2463.mature::ptg000001l:5548893-5548914(-);MFE=-27.8;Pval=0.039252
    ## ptg000001l   RNAhybrid   miRNA_binding   10400686    10400708    .   -   .   ID=Cluster_2463.mature::ptg000001l:5548893-5548914(-);MFE=-27.8;Pval=0.039252
    ## ptg000004l   RNAhybrid   miRNA_binding   11697985    11698007    .   -   .   ID=Cluster_2463.mature::ptg000001l:5548893-5548914(-);MFE=-28.3;Pval=0.030148
    ## ptg000007l   RNAhybrid   miRNA_binding   9978675 9978697 .   -   .   ID=Cluster_2463.mature::ptg000001l:5548893-5548914(-);MFE=-27.8;Pval=0.039252
    ## ptg000018l   RNAhybrid   miRNA_binding   13812477    13812499    .   -   .   ID=Cluster_2463.mature::ptg000001l:5548893-5548914(-);MFE=-27.8;Pval=0.039252
    ## ptg000026l   RNAhybrid   miRNA_binding   7994110 7994132 .   -   .   ID=Cluster_2463.mature::ptg000001l:5548893-5548914(-);MFE=-27.8;Pval=0.039252
    ## ptg000036l   RNAhybrid   miRNA_binding   4522999 4523022 .   +   .   ID=Cluster_2859.mature::ptg000001l:20063094-20063116(+);MFE=-26.4;Pval=0.003604
    ## ptg000036l   RNAhybrid   miRNA_binding   4522999 4523022 .   +   .   ID=Cluster_2859.mature::ptg000001l:20063094-20063116(+);MFE=-26.4;Pval=0.003604
    ## ntLink_8 RNAhybrid   miRNA_binding   20869798    20869822    .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-21.5;Pval=0.024451
    ## ntLink_8 RNAhybrid   miRNA_binding   20869798    20869822    .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-21.5;Pval=0.024451
    ## ntLink_8 RNAhybrid   miRNA_binding   30064524    30064548    .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-26.9;Pval=0.007283
    ## ntLink_8 RNAhybrid   miRNA_binding   30064524    30064548    .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-26.9;Pval=0.007283
    ## ptg000020l   RNAhybrid   miRNA_binding   5862537 5862561 .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-28.8;Pval=0.032792
    ## ptg000020l   RNAhybrid   miRNA_binding   5862537 5862561 .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-28.8;Pval=0.032792
    ## ptg000024l   RNAhybrid   miRNA_binding   9175751 9175775 .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-26.2;Pval=0.043413
    ## ptg000024l   RNAhybrid   miRNA_binding   9175751 9175775 .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-26.2;Pval=0.043413
    ## ptg000024l   RNAhybrid   miRNA_binding   9175751 9175775 .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-26.2;Pval=0.043413
    ## ptg000024l   RNAhybrid   miRNA_binding   9175751 9175775 .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-26.2;Pval=0.043413
    ## ptg000024l   RNAhybrid   miRNA_binding   9175751 9175775 .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-26.2;Pval=0.043413
    ## ptg000024l   RNAhybrid   miRNA_binding   9175751 9175775 .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-26.2;Pval=0.043413
    ## ptg000024l   RNAhybrid   miRNA_binding   9175751 9175775 .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-26.2;Pval=0.043413
    ## ptg000024l   RNAhybrid   miRNA_binding   9175751 9175775 .   +   .   ID=Cluster_3366.mature::ptg000002l:14046285-14046308(+);MFE=-26.2;Pval=0.043413
    ## ntLink_8 RNAhybrid   miRNA_binding   20869798    20869822    .   +   .   ID=Cluster_3367.mature::ptg000002l:14046591-14046614(+);MFE=-21.5;Pval=0.029415
    ## ntLink_8 RNAhybrid   miRNA_binding   20869798    20869822    .   +   .   ID=Cluster_3367.mature::ptg000002l:14046591-14046614(+);MFE=-21.5;Pval=0.029415
    ## ntLink_8 RNAhybrid   miRNA_binding   30064524    30064548    .   +   .   ID=Cluster_3367.mature::ptg000002l:14046591-14046614(+);MFE=-26.9;Pval=0.008933
    ## ntLink_8 RNAhybrid   miRNA_binding   30064524    30064548    .   +   .   ID=Cluster_3367.mature::ptg000002l:14046591-14046614(+);MFE=-26.9;Pval=0.008933
    ## ptg000020l   RNAhybrid   miRNA_binding   5862537 5862561 .   +   .   ID=Cluster_3367.mature::ptg000002l:14046591-14046614(+);MFE=-28.8;Pval=0.03925
    ## ptg000020l   RNAhybrid   miRNA_binding   5862537 5862561 .   +   .   ID=Cluster_3367.mature::ptg000002l:14046591-14046614(+);MFE=-28.8;Pval=0.03925
    ## ntLink_8 RNAhybrid   miRNA_binding   16931011    16931033    .   -   .   ID=Cluster_4220.mature::ptg000007l:915927-915948(-);MFE=-34.1;Pval=0.0027
    ## ntLink_8 RNAhybrid   miRNA_binding   16953818    16953840    .   -   .   ID=Cluster_4220.mature::ptg000007l:915927-915948(-);MFE=-34.1;Pval=0.0027
    ## ntLink_8 RNAhybrid   miRNA_binding   16978257    16978279    .   -   .   ID=Cluster_4220.mature::ptg000007l:915927-915948(-);MFE=-34.1;Pval=0.0027
    ## ptg000018l   RNAhybrid   miRNA_binding   1167306 1167328 .   -   .   ID=Cluster_4220.mature::ptg000007l:915927-915948(-);MFE=-30.4;Pval=0.017938
    ## ptg000035l   RNAhybrid   miRNA_binding   9011319 9011341 .   -   .   ID=Cluster_4220.mature::ptg000007l:915927-915948(-);MFE=-30.2;Pval=0.01986
    ## ntLink_6 RNAhybrid   miRNA_binding   5772715 5772737 .   +   .   ID=Cluster_4254.mature::ptg000007l:3377335-3377356(+);MFE=-26.2;Pval=0.025451
    ## ptg000030l   RNAhybrid   miRNA_binding   2128313 2128335 .   +   .   ID=Cluster_4254.mature::ptg000007l:3377335-3377356(+);MFE=-26;Pval=0.034916
    ## ptg000030l   RNAhybrid   miRNA_binding   2128313 2128335 .   +   .   ID=Cluster_4254.mature::ptg000007l:3377335-3377356(+);MFE=-26;Pval=0.034916
    ## ptg000030l   RNAhybrid   miRNA_binding   2128313 2128335 .   +   .   ID=Cluster_4254.mature::ptg000007l:3377335-3377356(+);MFE=-26;Pval=0.034916
    ## ptg000030l   RNAhybrid   miRNA_binding   2128313 2128335 .   +   .   ID=Cluster_4254.mature::ptg000007l:3377335-3377356(+);MFE=-26;Pval=0.034916
    ## ntLink_6 RNAhybrid   miRNA_binding   6374249 6374270 .   -   .   ID=Cluster_5012.mature::ptg000008l:10754789-10754809(-);MFE=-31;Pval=0.007748
    ## ptg000019l   RNAhybrid   miRNA_binding   3905999 3906021 .   +   .   ID=Cluster_10051.mature::ptg000016l:7795530-7795551(+);MFE=-25.5;Pval=0.013755
    ## ptg000039l   RNAhybrid   miRNA_binding   1038032 1038054 .   +   .   ID=Cluster_10051.mature::ptg000016l:7795530-7795551(+);MFE=-29;Pval=0.020136
    ## ntLink_8 RNAhybrid   miRNA_binding   17362806    17362828    .   -   .   ID=Cluster_10057.mature::ptg000016l:8599884-8599905(-);MFE=-28.6;Pval=0.036082
    ## ptg000008l   RNAhybrid   miRNA_binding   6343729 6343751 .   -   .   ID=Cluster_10057.mature::ptg000016l:8599884-8599905(-);MFE=-28.2;Pval=0.044287
    ## ptg000015l   RNAhybrid   miRNA_binding   2404079 2404101 .   -   .   ID=Cluster_10057.mature::ptg000016l:8599884-8599905(-);MFE=-28.2;Pval=0.044287
    ## ptg000026l   RNAhybrid   miRNA_binding   556252  556274  .   -   .   ID=Cluster_10057.mature::ptg000016l:8599884-8599905(-);MFE=-28.6;Pval=0.036082
    ## ptg000026l   RNAhybrid   miRNA_binding   896128  896150  .   -   .   ID=Cluster_10057.mature::ptg000016l:8599884-8599905(-);MFE=-28.6;Pval=0.036082
    ## ptg000047l   RNAhybrid   miRNA_binding   2819158 2819180 .   -   .   ID=Cluster_10093.mature::ptg000016l:11751407-11751428(-);MFE=-30.4;Pval=0.012598
    ## ptg000047l   RNAhybrid   miRNA_binding   2819158 2819180 .   -   .   ID=Cluster_10093.mature::ptg000016l:11751407-11751428(-);MFE=-30.4;Pval=0.012598
    ## ptg000047l   RNAhybrid   miRNA_binding   2819158 2819180 .   -   .   ID=Cluster_10093.mature::ptg000016l:11751407-11751428(-);MFE=-30.4;Pval=0.012598
    ## ptg000047l   RNAhybrid   miRNA_binding   2819158 2819180 .   -   .   ID=Cluster_10093.mature::ptg000016l:11751407-11751428(-);MFE=-30.4;Pval=0.012598
    ## ptg000021l   RNAhybrid   miRNA_binding   13553487    13553508    .   +   .   ID=Cluster_14768.mature::ptg000023l:37965298-37965318(+);MFE=-41.5;Pval=0.000033
    ## ptg000021l   RNAhybrid   miRNA_binding   16212956    16212977    .   +   .   ID=Cluster_14768.mature::ptg000023l:37965298-37965318(+);MFE=-27.9;Pval=0.009385
    ## ntLink_6 RNAhybrid   miRNA_binding   15902888    15902910    .   +   .   ID=Cluster_15316.mature::ptg000024l:4086254-4086275(+);MFE=-28.7;Pval=0.034534
    ## ptg000034l   RNAhybrid   miRNA_binding   183630  183652  .   +   .   ID=Cluster_15316.mature::ptg000024l:4086254-4086275(+);MFE=-30.8;Pval=0.011661
    ## ptg000064l   RNAhybrid   miRNA_binding   99593   99615   .   +   .   ID=Cluster_15316.mature::ptg000024l:4086254-4086275(+);MFE=-30.8;Pval=0.011661
    ## ntLink_8 RNAhybrid   miRNA_binding   17362805    17362827    .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-29.5;Pval=0.031681
    ## ntLink_8 RNAhybrid   miRNA_binding   31354517    31354539    .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-28.9;Pval=0.04276
    ## ptg000002l   RNAhybrid   miRNA_binding   3411170 3411192 .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-29;Pval=0.04068
    ## ptg000009l   RNAhybrid   miRNA_binding   10224592    10224614    .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-28.9;Pval=0.04276
    ## ptg000011l   RNAhybrid   miRNA_binding   11819551    11819573    .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-30.6;Pval=0.018216
    ## ptg000015l   RNAhybrid   miRNA_binding   2404078 2404100 .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-30.2;Pval=0.022286
    ## ptg000020l   RNAhybrid   miRNA_binding   16091592    16091614    .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-30.3;Pval=0.021191
    ## ptg000020l   RNAhybrid   miRNA_binding   16493605    16493627    .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-32.2;Pval=0.008105
    ## ptg000023l   RNAhybrid   miRNA_binding   24634245    24634267    .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-31.9;Pval=0.009436
    ## ptg000025l   RNAhybrid   miRNA_binding   1210284 1210306 .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-31.3;Pval=0.012788
    ## ptg000025l   RNAhybrid   miRNA_binding   5797267 5797289 .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-31.8;Pval=0.008009
    ## ptg000025l   RNAhybrid   miRNA_binding   5797267 5797289 .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-31.8;Pval=0.008009
    ## ptg000025l   RNAhybrid   miRNA_binding   5797267 5797289 .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-31.8;Pval=0.008009
    ## ptg000025l   RNAhybrid   miRNA_binding   5797267 5797289 .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-31.8;Pval=0.008009
    ## ptg000026l   RNAhybrid   miRNA_binding   556251  556273  .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-31.2;Pval=0.013451
    ## ptg000026l   RNAhybrid   miRNA_binding   896127  896149  .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-31.2;Pval=0.013451
    ## ptg000036l   RNAhybrid   miRNA_binding   859245  859267  .   -   .   ID=Cluster_15340.mature::ptg000024l:5256476-5256497(-);MFE=-31;Pval=0.014883
    ## ptg000012l   RNAhybrid   miRNA_binding   17288356    17288378    .   +   .   ID=Cluster_15851.mature::ptg000025l:10501052-10501073(+);MFE=-29;Pval=0.022003
    ## ptg000008l   RNAhybrid   miRNA_binding   8603028 8603051 .   -   .   ID=Cluster_15854.mature::ptg000025l:10668923-10668945(-);MFE=-27.3;Pval=0.021392
    ## ptg000008l   RNAhybrid   miRNA_binding   8603028 8603051 .   -   .   ID=Cluster_15854.mature::ptg000025l:10668923-10668945(-);MFE=-27.3;Pval=0.021392
    ## ptg000008l   RNAhybrid   miRNA_binding   8603028 8603051 .   -   .   ID=Cluster_15854.mature::ptg000025l:10668923-10668945(-);MFE=-27.3;Pval=0.021392
    ## ptg000008l   RNAhybrid   miRNA_binding   8603028 8603051 .   -   .   ID=Cluster_15854.mature::ptg000025l:10668923-10668945(-);MFE=-27.3;Pval=0.021392
    ## ptg000023l   RNAhybrid   miRNA_binding   7799756 7799778 .   -   .   ID=Cluster_16409.mature::ptg000026l:8745562-8745583(-);MFE=-27.7;Pval=0.04698
    ## ntLink_8 RNAhybrid   miRNA_binding   3046116 3046138 .   -   .   ID=Cluster_17776.mature::ptg000031l:5461327-5461348(-);MFE=-29.4;Pval=0.001937
    ## ntLink_6 RNAhybrid   miRNA_binding   4819483 4819506 .   -   .   ID=Cluster_17791.mature::ptg000031l:6751957-6751979(-);MFE=-21.9;Pval=0.022713
    ## ntLink_6 RNAhybrid   miRNA_binding   4819483 4819506 .   -   .   ID=Cluster_17791.mature::ptg000031l:6751957-6751979(-);MFE=-21.9;Pval=0.022713
    ## ntLink_6 RNAhybrid   miRNA_binding   4819483 4819506 .   -   .   ID=Cluster_17791.mature::ptg000031l:6751957-6751979(-);MFE=-21.9;Pval=0.022713
    ## ntLink_6 RNAhybrid   miRNA_binding   4819483 4819506 .   -   .   ID=Cluster_17791.mature::ptg000031l:6751957-6751979(-);MFE=-21.9;Pval=0.022713
    ## ptg000039l   RNAhybrid   miRNA_binding   483724  483746  .   -   .   ID=Cluster_19193.mature::ptg000039l:35786-35807(-);MFE=-28.3;Pval=0.04365

``` r
miRanda_RNAhybrid_miRNA_mRNA <- left_join(miRanda_miRNA_mRNA, RNAhybrid_miRNA_mRNA, by = c("miRNA_cluster" = "target_name"))
```
