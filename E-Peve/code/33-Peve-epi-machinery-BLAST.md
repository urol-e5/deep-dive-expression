33-Peve-epi-machinery-BLAST
================
Kathleen Durkin
2026-04-28

- [0.1 BLASTp](#01-blastp)
- [0.2 Reduce](#02-reduce)

Using Steven’s code from doing this in `timeseries`:
<https://github.com/urol-e5/timeseries_molecular/blob/dada661a07a45bf80b42953c4f8f8b7dad6824db/E-Peve/code/06-Peve-epimods-blast.qmd>

NOTE: Instead of using Steven’s e-value threshold of 1e-05, I’ll use the
more stringent threshold used in Ashey et al. 2025 in Jill’s check for
the existence of ncRNA machinery, for consistency: **1e-40**

Will use the protien database already generated in
32-Peve-ncRNA-machinery-BLAST

Epimachinery fasta provided by Hollie, obtain from timeseries repo if
necessary:

<https://raw.githubusercontent.com/urol-e5/timeseries_molecular/dada661a07a45bf80b42953c4f8f8b7dad6824db/D-Apul/data/Machinery.fasta>

``` bash
cd ../../data
curl -o Machinery.fasta https://raw.githubusercontent.com/urol-e5/timeseries_molecular/dada661a07a45bf80b42953c4f8f8b7dad6824db/D-Apul/data/Machinery.fasta
```

``` bash
head ../../data/Machinery.fasta
```

## 0.1 BLASTp

``` bash
fasta="../../data/Machinery.fasta"

/home/shared/ncbi-blast-2.15.0+/bin/blastp \
-query $fasta \
-db ../output/32-Peve-ncRNA-machinery-BLAST/Peve-proteins \
-out ../output/33-Peve-epi-machinery-BLAST/Mach-blastp-Peve_out.tab \
-evalue 1E-40 \
-num_threads 48 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/33-Peve-epi-machinery-BLAST/Mach-blastp-Peve_out.tab
```

    398 ../output/33-Peve-epi-machinery-BLAST/Mach-blastp-Peve_out.tab

This blast search using a 1e-40 evlaue threshold produces 398 hits,
compared to the 480 hits produced in [Steven’s `timeseries`
code](https://github.com/urol-e5/timeseries_molecular/blob/dada661a07a45bf80b42953c4f8f8b7dad6824db/E-Peve/code/06-Peve-epimods-blast.qmd)
using 1e-05

``` bash
head ../output/33-Peve-epi-machinery-BLAST/Mach-blastp-Peve_out.tab
```

    Dnmt1-201   Peve_00007510   50.452  1548    655 24  33  1538    65  1542    0.0 1437
    Dnmt3a-201  Peve_00023930   47.329  674 314 10  264 907 52  714 0.0 618
    Dnmt3b-201  Peve_00023930   45.104  674 303 11  228 860 68  715 0.0 583
    Dnmt3l-201  Peve_00033130   60.899  2491    910 22  456 2926    5   2451    0.0 3091
    USP9Y-201   Peve_00033130   59.044  2532    949 24  32  2521    13  2498    0.0 3009
    USP10-201   Peve_00029730   44.134  537 267 14  276 798 275 792 1.63e-129   405
    USP11-201   Peve_00004975   45.695  906 414 14  93  933 33  925 0.0 705
    USP12-201   Peve_00031882   76.630  368 71  2   17  369 1   368 0.0 591
    USP14-201   Peve_00036006   62.288  472 163 5   5   470 1   463 0.0 602
    USP15-201   Peve_00004975   58.245  661 244 11  286 934 283 923 0.0 734

## 0.2 Reduce

Let’s also make sure we have an easy-to-use csv associating each protein
with any matches.

``` r
raw <- read.table("../output/33-Peve-epi-machinery-BLAST/Mach-blastp-Peve_out.tab")

formatted <- raw  %>% dplyr::select(V1, V2)
formatted$target <- raw$V2
formatted$gene_name <- raw$V1
formatted <- formatted %>% dplyr::select(-V1, -V2)

# Remove trailing ##s after the target gene names (artifact of Ensembl)
formatted$gene_name <- sub("-2[0-9]{2}$", "", formatted$gene_name)

write.csv(formatted, "../output/33-Peve-epi-machinery-BLAST/Mach-blastp-Peve_reduced.csv", row.names = FALSE)
```
