32-Ptuh-ncRNA-machinery-BLAST
================
Kathleen Durkin
2026-04-15

- [1 Proteins](#1-proteins)
- [2 Make BLAST db](#2-make-blast-db)
- [3 BLASTp](#3-blastp)
- [4 Parse and annotate](#4-parse-and-annotate)

# 1 Proteins

Will be using P. meandrina protein sequences here:
<https://gannet.fish.washington.edu/seashell/snaps/Pocillopora_meandrina_HIv1.genes.pep.faa>

Download if necessary

``` bash
cd ../data
curl -o Pocillopora_meandrina_HIv1.genes.pep.faa https://gannet.fish.washington.edu/seashell/snaps/Pocillopora_meandrina_HIv1.genes.pep.faa
```

``` bash
head ../data/Pocillopora_meandrina_HIv1.genes.pep.faa
```

# 2 Make BLAST db

``` bash

/home/shared/ncbi-blast-2.15.0+/bin/makeblastdb \
-in ../data/Pocillopora_meandrina_HIv1.genes.pep.faa \
-dbtype prot \
-out ../output/32-Ptuh-ncRNA-machinery-BLAST/Ptuh-proteins
```

``` bash
head ../../data/ncRNA_machinery.fasta
```

# 3 BLASTp

``` bash
fasta="../../data/ncRNA_machinery.fasta"

/home/shared/ncbi-blast-2.15.0+/bin/blastp \
-query $fasta \
-db ../output/32-Ptuh-ncRNA-machinery-BLAST/Ptuh-proteins \
-out ../output/32-Ptuh-ncRNA-machinery-BLAST/ncRNAmach-blastp-Ptuh_out.tab \
-evalue 1E-05 \
-num_threads 48 \
-max_target_seqs 1 \
-max_hsps 1 \
-outfmt 6
```

``` bash
wc -l ../output/32-Ptuh-ncRNA-machinery-BLAST/ncRNAmach-blastp-Ptuh_out.tab
```

    865 ../output/32-Ptuh-ncRNA-machinery-BLAST/ncRNAmach-blastp-Ptuh_out.tab

``` bash
head ../output/32-Ptuh-ncRNA-machinery-BLAST/ncRNAmach-blastp-Ptuh_out.tab
```

    tr|A0A0C2N1X6|A0A0C2N1X6_THEKT  Pocillopora_meandrina_HIv1___RNAseq.g6540.t1    20.798  827 556 23  10  762 5   806 4.15e-43    167
    tr|A0A2B4SX82|A0A2B4SX82_STYPI  Pocillopora_meandrina_HIv1___RNAseq.g6540.t1    95.161  806 39  0   1   806 1   806 0.0 1602
    tr|A0A6P8J329|A0A6P8J329_ACTTE  Pocillopora_meandrina_HIv1___RNAseq.g6540.t1    65.675  807 271 6   1   802 1   806 0.0 1139
    tr|A0A6S7G2D7|A0A6S7G2D7_PARCT  Pocillopora_meandrina_HIv1___RNAseq.g6540.t1    52.847  808 370 5   1   799 1   806 0.0 915
    tr|A0A7M5V3C0|A0A7M5V3C0_9CNID  Pocillopora_meandrina_HIv1___RNAseq.g6540.t1    50.428  817 393 5   1   816 1   806 0.0 876
    tr|A0A913XIF0|A0A913XIF0_EXADI  Pocillopora_meandrina_HIv1___RNAseq.g6540.t1    64.764  806 279 4   1   801 1   806 0.0 1129
    tr|A0A9W9ZWX6|A0A9W9ZWX6_9CNID  Pocillopora_meandrina_HIv1___RNAseq.g6540.t1    83.643  807 131 1   1   807 1   806 0.0 1429
    tr|A0AAD9QN77|A0AAD9QN77_ACRCE  Pocillopora_meandrina_HIv1___RNAseq.g6540.t1    73.908  824 159 2   1   788 1   804 0.0 1268
    tr|A0AAU9VQE5|A0AAU9VQE5_9CNID  Pocillopora_meandrina_HIv1___RNAseq.g6540.t1    100.000 806 0   0   1   806 1   806 0.0 1671
    tr|A0ABM4DB04|A0ABM4DB04_HYDVU  Pocillopora_meandrina_HIv1___RNAseq.g6540.t1    50.860  814 390 6   1   812 1   806 0.0 865

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
cd ../output/32-Ptuh-ncRNA-machinery-BLAST

awk -F'\t' 'BEGIN{OFS="\t"} {split($1,a,"|"); $1=a[2]; print}' ncRNAmach-blastp-Ptuh_out.tab > ncRNAmach-blastp-Ptuh_parsed.tab
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
blast <- read.delim("../output/32-Ptuh-ncRNA-machinery-BLAST/ncRNAmach-blastp-Ptuh_parsed.tab", header = FALSE,
                    col.names = c("accession", "target", "pident", "length",
                                  "mismatch", "gapopen", "qstart", "qend",
                                  "sstart", "send", "evalue", "bitscore"))

# Join
blast_annotated <- blast %>%
  left_join(ref, by = "accession") %>%
  dplyr::select(target, protein_name, gene_name, std_gene_name, species, accession, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore)

write.csv(blast_annotated, "../output/32-Ptuh-ncRNA-machinery-BLAST/ncRNAmach-blastp-Ptuh_annotated.csv", row.names = FALSE)
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

write.csv(blast_reduc, "../output/32-Ptuh-ncRNA-machinery-BLAST/ncRNAmach-blastp-Ptuh_annotated_reduced.csv", row.names = FALSE)
```
