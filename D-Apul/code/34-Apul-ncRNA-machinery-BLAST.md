34-Apul-ncRNA-machinery-BLAST
================
Kathleen Durkin
2026-04-14

- [1 Proteins](#1-proteins)
- [2 Make BLAST db](#2-make-blast-db)
- [3 BLASTp](#3-blastp)
- [4 Parse and annotate](#4-parse-and-annotate)

Using Steven’s code from doing this in `timeseries`:
<https://github.com/urol-e5/timeseries_molecular/blob/main/D-Apul/code/34-Apul-ncRNA-machinery-BLAST.qmd>

# 1 Proteins

Will be using A. pulchra protein sequences here:

<https://raw.githubusercontent.com/urol-e5/deep-dive-expression/main/D-Apul/data/Apulchra-genome.pep.faa>

``` bash
head ../data/Apulchra-genome.pep.faa
```

# 2 Make BLAST db

``` bash

/home/shared/ncbi-blast-2.15.0+/bin/makeblastdb \
-in ../data/Apulchra-genome.pep.faa \
-dbtype prot \
-out ../output/34-Apul-ncRNA-machinery-BLAST/Apul-proteins
```

``` bash
head ../../data/ncRNA_machinery.fasta
```

# 3 BLASTp

``` bash
fasta="../../data/ncRNA_machinery.fasta"

/home/shared/ncbi-blast-2.15.0+/bin/blastp \
-query $fasta \
-db ../output/34-Apul-ncRNA-machinery-BLAST/Apul-proteins \
-out ../output/34-Apul-ncRNA-machinery-BLAST/ncRNAmach-blastp-Apul_out.tab \
-evalue 1E-05 \
-num_threads 48 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/34-Apul-ncRNA-machinery-BLAST/ncRNAmach-blastp-Apul_out.tab
```

    865 ../output/34-Apul-ncRNA-machinery-BLAST/ncRNAmach-blastp-Apul_out.tab

``` bash
head ../output/34-Apul-ncRNA-machinery-BLAST/ncRNAmach-blastp-Apul_out.tab
```

    tr|A0A0C2N1X6|A0A0C2N1X6_THEKT  FUN_045726-T1   22.805  820 544 20  10  760 5   804 6.50e-50    188
    tr|A0A2B4SX82|A0A2B4SX82_STYPI  FUN_045726-T1   77.861  804 178 0   1   804 1   804 0.0 1334
    tr|A0A6P8J329|A0A6P8J329_ACTTE  FUN_045726-T1   64.144  806 281 8   1   800 1   804 0.0 1103
    tr|A0A6S7G2D7|A0A6S7G2D7_PARCT  FUN_045726-T1   50.435  805 392 4   1   798 1   805 0.0 886
    tr|A0A7M5V3C0|A0A7M5V3C0_9CNID  FUN_045726-T1   50.675  815 390 5   1   814 1   804 0.0 873
    tr|A0A913XIF0|A0A913XIF0_EXADI  FUN_045726-T1   63.433  804 289 5   1   799 1   804 0.0 1094
    tr|A0A9W9ZWX6|A0A9W9ZWX6_9CNID  FUN_045726-T1   81.491  805 148 1   1   805 1   804 0.0 1382
    tr|A0AAD9QN77|A0AAD9QN77_ACRCE  FUN_045726-T1   92.485  825 6   2   1   789 1   805 0.0 1563
    tr|A0AAU9VQE5|A0AAU9VQE5_9CNID  FUN_045726-T1   78.856  804 170 0   1   804 1   804 0.0 1348
    tr|A0ABM4DB04|A0ABM4DB04_HYDVU  FUN_045726-T1   51.870  802 373 7   14  810 11  804 0.0 868

# 4 Parse and annotate

This .tab file associates A. pulchra proteins with ncRNA machinery
proteins found in the reference db. However, the names are just taken
straight from the fasta headers, which ends up being mostly just the
Uniprot accession number. Let’s parse this into a more interpretable db.

parse the reference ncRNA machinery fasta headers

``` bash
cd ../../data/
  
grep "^>" ncRNA_machinery.fasta | \
sed 's/^>//' | \
awk -F'|' '{
  # Get accession
  acc = $2
  
  # Everything after the third pipe
  rest = $3
  
  # Strip entry_name (first token before space)
  sub(/^[^ ]+ /, "", rest)
  
  # Parse fields using OS=, OX=, GN=, CGN=, PE=
  protein = rest; sub(/ CGN=.*/, "", protein); sub(/ OS=.*/, "", protein)
  
  species = rest; sub(/.*OS=/, "", species); sub(/ OX=.*/, "", species)
  
  if (rest ~ /GN=/) {
    gene = rest; sub(/.*GN=/, "", gene); sub(/ PE=.*/, "", gene)
  } else {
    gene = "NA"
  }
  
  if (rest ~ /CGN=/) {
    cgn = rest; sub(/.*CGN=/, "", cgn); sub(/ OS=.*/, "", cgn)
  } else {
    cgn = "NA"
  }
  
  # Quote protein and species in case of commas
  gsub(/"/, "\"\"", protein)
  gsub(/"/, "\"\"", species)
  
  print acc "," "\"" protein "\"" "," gene "," cgn "," "\"" species "\""
}' | \
sed '1i accession,protein_name,gene_name,std_gene_name,species' \
> ncRNA_machinery_reference_table.csv
```

Now parse blastp output rows to reduce to accession numbers

``` bash
cd ../output/34-Apul-ncRNA-machinery-BLAST

awk -F'\t' 'BEGIN{OFS="\t"} {split($1,a,"|"); $1=a[2]; print}' ncRNAmach-blastp-Apul_out.tab > ncRNAmach-blastp-Apul_parsed.tab
```

Now read in and match with reference db, so that blast results are
associated with protein/gene names (not just accession numbers)

(Note that I can’t rely on just the gene names, because most of the
Uniprot entries used non-standard gene names)

``` r
library(dplyr)

# Read reference table
ref <- read.csv("../../data/ncRNA_machinery_reference_table.csv", stringsAsFactors = FALSE)

# Read parsed BLAST output
blast <- read.delim("../output/34-Apul-ncRNA-machinery-BLAST/ncRNAmach-blastp-Apul_parsed.tab", header = FALSE,
                    col.names = c("accession", "target", "pident", "length",
                                  "mismatch", "gapopen", "qstart", "qend",
                                  "sstart", "send", "evalue", "bitscore"))

# Join
blast_annotated <- blast %>%
  left_join(ref, by = "accession") %>%
  dplyr::select(target, protein_name, gene_name, std_gene_name, species, accession, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)

write.csv(blast_annotated, "../output/34-Apul-ncRNA-machinery-BLAST/ncRNAmach-blastp-Apul_annotated.csv", row.names = FALSE)
```

Finally, let’s reduce this to unique associations. In the ncNRA
machinery db I compiled, I sometimes included multiple sequences of a
given protein sequenced in different species. This was partly due to
convenience, but also because (a) I wanted “backups” in case one of the
sequences was mis-annotated or had other problems (these are all
unreviewed entries), and (b) in case a protein is highly divergent, and
only a subset of its entries matched our species.

Now that I’ve done the blast, though, I only really need to know the
$$apulchra mRNA$$ - $$standardized gene name$$ association. Let’s make
this reduced db.

``` r
blast_reduc <- blast_annotated %>% dplyr::select(target, std_gene_name) %>% distinct()

write.csv(blast_reduc, "../output/34-Apul-ncRNA-machinery-BLAST/ncRNAmach-blastp-Apul_annotated_reduced.csv", row.names = FALSE)
```
