35-Apul-epi-machinery-BLAST
================
Kathleen Durkin
2026-04-28

- [0.1 BLASTp](#01-blastp)

Using Steven’s code from doing this in `timeseries`:
<https://github.com/urol-e5/timeseries_molecular/blob/main/D-Apul/code/25-Apul-epimods-blast.qmd>

NOTE: Instead of using Steven’s e-value threshold of 1e-05, I’ll use the
more stringent threshold used in Ashey et al. 2025 in Jill’s check for
the existence of ncRNA machinery, for consistency: **1e-40**

Will use the protien database already generated in
34-Apul-ncRNA-machinery-BLAST

Epimachinery fasta provided by Hollie, obtain from timeseries repo:

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
-db ../output/34-Apul-ncRNA-machinery-BLAST/Apul-proteins \
-out ../output/35-Apul-epi-machinery-BLAST/Mach-blastp-Apul_out.tab \
-evalue 1E-40 \
-num_threads 48 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/35-Apul-epi-machinery-BLAST/Mach-blastp-Apul_out.tab
```

    397 ../output/35-Apul-epi-machinery-BLAST/Mach-blastp-Apul_out.tab

This blast search using a 1e-40 evlaue threshold produces 397 hits,
compared to the 446 hits produced in [Steven’s `timeseries`
code](https://github.com/urol-e5/timeseries_molecular/blob/main/D-Apul/code/25-Apul-epimods-blast.qmd)
using 1e-05

``` bash
head ../output/35-Apul-epi-machinery-BLAST/Mach-blastp-Apul_out.tab
```

    Dnmt1-201   FUN_043321-T1   51.653  1603    653 29  48  1605    4   1529    0.0 1528
    Dnmt3a-201  FUN_021938-T1   46.097  679 317 9   267 907 51  718 0.0 614
    Dnmt3b-201  FUN_021938-T1   43.529  680 320 12  221 860 64  719 0.0 570
    Dnmt3l-201  FUN_041562-T1   60.977  2519    890 26  469 2954    35  2493    0.0 3088
    USP9Y-201   FUN_041562-T1   60.859  2468    890 26  67  2502    67  2490    0.0 3036
    USP10-201   FUN_034781-T1   44.492  463 242 8   341 798 343 795 1.20e-117   374
    USP11-201   FUN_005808-T1   47.303  890 401 15  93  931 35  907 0.0 751
    USP12-201   FUN_024199-T1   76.152  369 73  2   16  369 1   369 0.0 578
    USP14-201   FUN_014539-T1   61.122  499 173 4   1   487 1   490 0.0 635
    USP15-201   FUN_005808-T1   53.961  934 377 13  8   935 22  908 0.0 973
