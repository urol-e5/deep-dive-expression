33-Ptuh-epi-machinery-BLAST
================
Kathleen Durkin
2026-04-28

- [0.1 BLASTp](#01-blastp)

Using Steven’s code from doing this in `timeseries`:
<https://github.com/urol-e5/timeseries_molecular/blob/dada661a07a45bf80b42953c4f8f8b7dad6824db/F-Ptua/code/03-Ptua-epimods-blast.qmd>

NOTE: Instead of using Steven’s e-value threshold of 1e-05, I’ll use the
more stringent threshold used in Ashey et al. 2025 in Jill’s check for
the existence of ncRNA machinery, for consistency: **1e-40**

Will use the protein database already generated in
32-Ptuh-ncRNA-machinery-BLAST

Epimachinery fasta provided by Hollie, obtain from timeseries repo if
needed:

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
-db ../output/32-Ptuh-ncRNA-machinery-BLAST/Ptuh-proteins \
-out ../output/33-Ptuh-epi-machinery-BLAST/Mach-blastp-Ptuh_out.tab \
-evalue 1E-40 \
-num_threads 48 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/33-Ptuh-epi-machinery-BLAST/Mach-blastp-Ptuh_out.tab
```

    403 ../output/33-Ptuh-epi-machinery-BLAST/Mach-blastp-Ptuh_out.tab

This blast search using a 1e-40 evlaue threshold produces 403 hits,
compared to the 480 hits produced in [Steven’s `timeseries`
code](https://github.com/urol-e5/timeseries_molecular/blob/dada661a07a45bf80b42953c4f8f8b7dad6824db/F-Ptua/code/03-Ptua-epimods-blast.qmd)
using 1e-05

``` bash
head ../output/33-Ptuh-epi-machinery-BLAST/Mach-blastp-Ptuh_out.tab
```

    Dnmt1-201   Pocillopora_meandrina_HIv1___RNAseq.g8948.t1    57.918  1345    517 16  294 1598    191 1526    0.0 1543
    Dnmt3a-201  Pocillopora_meandrina_HIv1___RNAseq.g2960.t2    48.233  651 294 7   286 905 78  716 0.0 621
    Dnmt3b-201  Pocillopora_meandrina_HIv1___RNAseq.g2960.t2    43.776  715 322 12  181 857 44  716 0.0 596
    Dnmt3l-201  Pocillopora_meandrina_HIv1___RNAseq.g16036.t1   58.795  1859    717 17  453 2289    1   1832    0.0 2202
    USP9Y-201   Pocillopora_meandrina_HIv1___RNAseq.g16036.t1   58.297  1832    714 17  30  1838    27  1831    0.0 2143
    USP10-201   Pocillopora_meandrina_HIv1___RNAseq.g12276.t1   47.193  481 210 15  338 798 366 822 1.53e-127   401
    USP11-201   Pocillopora_meandrina_HIv1___RNAseq.g19447.t1   45.005  931 442 19  62  935 2   919 0.0 733
    USP12-201   Pocillopora_meandrina_HIv1___RNAseq.g18898.t1   74.255  369 78  3   16  369 1   367 0.0 563
    USP14-201   Pocillopora_meandrina_HIv1___RNAseq.g22603.t1   62.805  492 168 4   1   487 1   482 0.0 646
    USP15-201   Pocillopora_meandrina_HIv1___RNAseq.g19447.t1   53.340  943 387 15  5   934 13  915 0.0 954
