32-Peve-ncRNA-machinery-BLAST
================
Kathleen Durkin
2026-04-15

- [1 Proteins](#1-proteins)
- [2 Make BLAST db](#2-make-blast-db)
- [3 BLASTp](#3-blastp)
- [4 Parse and annotate](#4-parse-and-annotate)

NOTE: Instead of using Steven’s e-value threshold of 1e-05, I’ll use the
more stringent threshold used in Ashey et al. 2025 in Jill’s check for
the existence of ncRNA machinery, for consistency: **1e-40**

# 1 Proteins

Will be using A. pulchra protein sequences here:

<https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.annot.pep.fa>

Dowload if necessary

``` bash
cd ../data
curl -o Porites_evermanni_v1.annot.pep.fa https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.annot.pep.fa
```

``` bash
head ../data/Porites_evermanni_v1.annot.pep.fa
```

# 2 Make BLAST db

``` bash

/home/shared/ncbi-blast-2.15.0+/bin/makeblastdb \
-in ../data/Porites_evermanni_v1.annot.pep.fa \
-dbtype prot \
-out ../output/32-Peve-ncRNA-machinery-BLAST/Peve-proteins
```

``` bash
head ../../data/ncRNA_machinery.fasta
```

# 3 BLASTp

``` bash
fasta="../../data/ncRNA_machinery.fasta"

/home/shared/ncbi-blast-2.15.0+/bin/blastp \
-query $fasta \
-db ../output/32-Peve-ncRNA-machinery-BLAST/Peve-proteins \
-out ../output/32-Peve-ncRNA-machinery-BLAST/ncRNAmach-blastp-Peve_out.tab \
-evalue 1E-40 \
-num_threads 48 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/32-Peve-ncRNA-machinery-BLAST/ncRNAmach-blastp-Peve_out.tab
```

    772 ../output/32-Peve-ncRNA-machinery-BLAST/ncRNAmach-blastp-Peve_out.tab

``` bash
head ../output/32-Peve-ncRNA-machinery-BLAST/ncRNAmach-blastp-Peve_out.tab
```

    tr|A0A224ASV7|A0A224ASV7_ACRTE  Peve_00019588   92.350  732 52  2   1   732 1   728 0.0 1276
    tr|Q967E1|Q967E1_DENKL  Peve_00019588   82.273  220 39  0   5   224 6   225 2.11e-124   388
    tr|A0A3M6U8K5|A0A3M6U8K5_POCDA  Peve_00006460   82.440  1082    168 6   1   1078    1   1064    0.0 1604
    tr|A0A6P8H2M2|A0A6P8H2M2_ACTTE  Peve_00006460   58.420  1152    379 22  1   1140    1   1064    0.0 1095
    tr|A0A7M5UZN6|A0A7M5UZN6_9CNID  Peve_00006460   52.361  953 348 19  1   898 1   902 0.0 850
    tr|A0A7M5V808|A0A7M5V808_9CNID  Peve_00006460   52.466  953 347 19  1   898 1   902 0.0 855
    tr|A0A913YAN3|A0A913YAN3_EXADI  Peve_00006460   63.417  1118    338 18  1   1101    1   1064    0.0 1197
    tr|A0AAD9UTI2|A0AAD9UTI2_ACRCE  Peve_00006460   80.863  1066    183 9   1   1053    1   1058    0.0 1558
    tr|A0AAU9WW55|A0AAU9WW55_9CNID  Peve_00006460   78.295  1138    167 9   1   1132    1   1064    0.0 1551
    tr|T2M958|T2M958_HYDVU  Peve_00006460   53.611  914 360 18  1   862 1   902 0.0 843

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
cd ../output/32-Peve-ncRNA-machinery-BLAST

awk -F'\t' 'BEGIN{OFS="\t"} {split($1,a,"|"); $1=a[2]; print}' ncRNAmach-blastp-Peve_out.tab > ncRNAmach-blastp-Peve_parsed.tab
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
blast <- read.delim("../output/32-Peve-ncRNA-machinery-BLAST/ncRNAmach-blastp-Peve_parsed.tab", header = FALSE,
                    col.names = c("accession", "target", "pident", "length",
                                  "mismatch", "gapopen", "qstart", "qend",
                                  "sstart", "send", "evalue", "bitscore"))

# Join
blast_annotated <- blast %>%
  left_join(ref, by = "accession") %>%
  dplyr::select(target, protein_name, gene_name, std_gene_name, species, accession, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)

write.csv(blast_annotated, "../output/32-Peve-ncRNA-machinery-BLAST/ncRNAmach-blastp-Peve_annotated.csv", row.names = FALSE)
```

Finally, let’s reduce this to unique associations. In the ncNRA
machinery db I compiled, I sometimes included multiple sequences of a
given protein sequenced in different species. This was partly due to
convenience, but also because (a) I wanted “backups” in case one of the
sequences was mis-annotated or had other problems (these are all
unreviewed entries), and (b) in case a protein is highly divergent, and
only a subset of its entries matched our species.

Now that I’ve done the blast, though, I only really need to know the
$$p.evermanni mRNA$$ - $$standardized gene name$$ association. Let’s
make this reduced db.

``` r
blast_reduc <- blast_annotated %>% dplyr::select(target, std_gene_name) %>% distinct()

write.csv(blast_reduc, "../output/32-Peve-ncRNA-machinery-BLAST/ncRNAmach-blastp-Peve_annotated_reduced.csv", row.names = FALSE)
```
