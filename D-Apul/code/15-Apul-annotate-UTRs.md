11-Apul-annotate-UTRs
================
Kathleen Durkin
2024-11-20

For miRNA target prediction, we expect most binding to occur in the
3’UTR regions. However, our reference genome `Apulcra-genome.gff` is not
annotated with UTRs. We need to annotate those manually

First let’s take a look at what our reference `Apulcra-genome.gff` file
looks like

``` bash
grep -v '^#' ../data/Apulcra-genome.gff | cut -s -f 3 | sort | uniq -c | sort -rn
echo ""
head -10 ../data/Apulcra-genome.gff
```

    ##  209537 exon
    ##  201613 CDS
    ##   44371 gene
    ##   36447 mRNA
    ##    7924 tRNA
    ## 
    ## ##gff-version 3
    ## ntLink_0 funannotate gene    1105    7056    .   +   .   ID=FUN_000001;
    ## ntLink_0 funannotate mRNA    1105    7056    .   +   .   ID=FUN_000001-T1;Parent=FUN_000001;product=hypothetical protein;
    ## ntLink_0 funannotate exon    1105    1188    .   +   .   ID=FUN_000001-T1.exon1;Parent=FUN_000001-T1;
    ## ntLink_0 funannotate exon    1861    1941    .   +   .   ID=FUN_000001-T1.exon2;Parent=FUN_000001-T1;
    ## ntLink_0 funannotate exon    2762    2839    .   +   .   ID=FUN_000001-T1.exon3;Parent=FUN_000001-T1;
    ## ntLink_0 funannotate exon    5044    7056    .   +   .   ID=FUN_000001-T1.exon4;Parent=FUN_000001-T1;
    ## ntLink_0 funannotate CDS 1105    1188    .   +   0   ID=FUN_000001-T1.cds;Parent=FUN_000001-T1;
    ## ntLink_0 funannotate CDS 1861    1941    .   +   0   ID=FUN_000001-T1.cds;Parent=FUN_000001-T1;
    ## ntLink_0 funannotate CDS 2762    2839    .   +   0   ID=FUN_000001-T1.cds;Parent=FUN_000001-T1;

According to [NCBI](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/),
the standard gff3 format should annotate regions like this:

    gene1            ================================    ID=gene1
    mRNA1            ================================    ID=mRNA1;Parent=gene1
    exon1            ====                                Parent=mRNA1
    five_prime_UTR   ==                                  Parent=mRNA1
    CDS1               ==                                Parent=mRNA1
    CDS2                     ==========                  Parent=mRNA1
    CDS3                                  ====          Parent=mRNA1
    three_prime_UTR                            ======    Parent=mRNA1

Note how a region annotated as `mRNA` should contain both the 5’UTR and
3’UTR, as well as all CDS regions.

Unfortunately, looking at our gff, it seems that `Apulcra-genome.gff` is
formatted more like this:

    gene1            ================================    ID=gene1
    mRNA1            ================================    ID=mRNA1;Parent=gene1
    exon1            ======                              Parent=mRNA1
    CDS1             ======                              Parent=mRNA1
    exon2                          ==================    Parent=mRNA1
    CDS2                           ==================    Parent=mRNA1

Where an annotated mRNA contains only the CDS regions and internal UTR
regions. That means there is no straightforward way to identify the
actual 5’ and 3’ UTR regions from the gff.

Instead, we can just extract a best-guess region from immediately before
and immediately after the annotated mRNA. Jill Ashey validated during
our `deep-dive` work that a 1000bp region captures the 3’UTR well (see
post and code
[here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-06-15-e5-deepdive-miRNA-TargetPrediction.md)).
This is also convenient because RNAhybrid accepts sequences of maximum
1000bp in length.

So we want to 1) isolate mRNA regions in `Apulcra-genome.gff` 2)
determine sense of the annotated region (to inform which end is 5’ and
which is 3’) 3) make a new line in gff that covers the 1000bp
immediately preceding the mRNA, annotated appropriately as 5’/3’ UTR 4)
make new line that covers the 1000bp immediately following the mRNA,
annotated appropriately as 5’/3’ UTR 5) ensure none of the newly
annotated UTR regions overlap with an existing mRNA region – remove
overlapping region from UTR annotation if necessary

Note: while the majority of miRNA functional binding occurs in the
3’UTR, it is possible for functional binding to occur in other regions,
including the 5’UTR and CDS. I want to run target prediction for all
regions, at least at first, so I want to annotate both the 3’ and 5’ UTR
regions of our genome.

We’ll use a modified version of Jill’s code for this

``` bash

# extract the mRNAs
grep $'\tmRNA\t' ../data/Apulcra-genome.gff > ../output/15-Apul-annotate-UTRs/Apulcra-genome-mRNA_only.gff

# Let's also isolate the CDS while we're at it
grep $'\tCDS\t' ../data/Apulcra-genome.gff > ../output/15-Apul-annotate-UTRs/Apulcra-genome-CDS_only.gff


# check
wc -l ../output/15-Apul-annotate-UTRs/Apulcra-genome-mRNA_only.gff
echo ""
head -5 ../output/15-Apul-annotate-UTRs/Apulcra-genome-mRNA_only.gff
```

    ## 36447 ../output/15-Apul-annotate-UTRs/Apulcra-genome-mRNA_only.gff
    ## 
    ## ntLink_0 funannotate mRNA    1105    7056    .   +   .   ID=FUN_000001-T1;Parent=FUN_000001;product=hypothetical protein;
    ## ntLink_0 funannotate mRNA    10215   15286   .   +   .   ID=FUN_000002-T1;Parent=FUN_000002;product=hypothetical protein;
    ## ntLink_0 funannotate mRNA    32057   33275   .   +   .   ID=FUN_000003-T1;Parent=FUN_000003;product=hypothetical protein;
    ## ntLink_0 funannotate mRNA    34824   42794   .   +   .   ID=FUN_000004-T1;Parent=FUN_000004;product=hypothetical protein;
    ## ntLink_0 funannotate mRNA    45953   51024   .   +   .   ID=FUN_000005-T1;Parent=FUN_000005;product=hypothetical protein;

``` bash
cd ../output/15-Apul-annotate-UTRs

# Extract scaffold lengths
cat is ../../data/Apulchra-genome.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Apul.Chromosome_lengths.txt

# Extract scaffold names
awk -F" " '{print $1}' Apul.Chromosome_lengths.txt > Apul.Chromosome_names.txt

# Check
wc -l Apul.Chromosome_lengths.txt 
echo ""
head -3 Apul.Chromosome_lengths.txt 
echo ""
head -3 Apul.Chromosome_names.txt 
```

    ## cat: is: No such file or directory
    ## 174 Apul.Chromosome_lengths.txt
    ## 
    ## ntLink_7 182921
    ## ntLink_8 37517065
    ## ntLink_0 96579
    ## 
    ## ntLink_7
    ## ntLink_8
    ## ntLink_0

The following code will sort the mRNA gff, extract 1kb down the 3’ end
of mRNA, subtract portions of the 1kb flank (representing the 3’UTR)
from any overlapping mRNA, and make fasta file of the 3’UTRs. It will do
the same for 5’UTRs.

``` bash
cd ../output/15-Apul-annotate-UTRs
export PATH="/home/shared/bedtools2/bin:$PATH"

echo "Sorting gffs by chromosome" $(date)

sortBed -faidx Apul.Chromosome_names.txt -i Apulcra-genome-mRNA_only.gff > Apul_GFFannotation.mRNA_sorted.gff

echo "Sorting complete!" $(date)

echo "Extracting 1kb 3' UTRs" $(date)

bedtools flank -i Apul_GFFannotation.mRNA_sorted.gff -g Apul.Chromosome_lengths.txt -l 0 -r 1000 -s | \
awk '{gsub("mRNA","3prime_UTR",$3); print $0 }' | \
awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | \
tr ' ' '\t' > Apul.GFFannotation.3UTR_1kb.gff

echo "Subtract portions of UTRs that overlap nearby genes" $(date)

bedtools subtract -a Apul.GFFannotation.3UTR_1kb.gff -b Apul_GFFannotation.mRNA_sorted.gff > Apul.GFFannotation.3UTR_1kb_corrected.gff 
echo "3' UTRs identified!" $(date)

echo "Extracting 3' UTR sequences" $(date)

bedtools getfasta -fi ../../data/Apulchra-genome.fa -bed Apul.GFFannotation.3UTR_1kb_corrected.gff -fo Apul_3UTR_1kb.fasta -fullHeader


echo "Extracting 1kb 5' UTRs" $(date)

bedtools flank -i Apul_GFFannotation.mRNA_sorted.gff -g Apul.Chromosome_lengths.txt -l 1000 -r 0 -s | \
awk '{gsub("mRNA","five_prime_UTR",$3); print $0 }' | \
awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | \
tr ' ' '\t' > Apul.GFFannotation.5UTR_1kb.gff

echo "Subtract portions of 5' UTRs that overlap nearby genes" $(date)

bedtools subtract -a Apul.GFFannotation.5UTR_1kb.gff -b Apul_GFFannotation.mRNA_sorted.gff > Apul.GFFannotation.5UTR_1kb_corrected.gff

echo "5' UTRs identified!" $(date)

echo "Extracting 5' UTR sequences" $(date)

bedtools getfasta -fi ../../data/Apulchra-genome.fa -bed Apul.GFFannotation.5UTR_1kb_corrected.gff -fo Apul_5UTR_1kb.fasta -fullHeader

echo "Sequence extraction complete!" $(date)
```

    ## Sorting gffs by chromosome Wed Nov 20 17:37:39 PST 2024
    ## Sorting complete! Wed Nov 20 17:37:39 PST 2024
    ## Extracting 1kb 3' UTRs Wed Nov 20 17:37:39 PST 2024
    ## Subtract portions of UTRs that overlap nearby genes Wed Nov 20 17:37:39 PST 2024
    ## 3' UTRs identified! Wed Nov 20 17:37:40 PST 2024
    ## Extracting 3' UTR sequences Wed Nov 20 17:37:40 PST 2024
    ## Extracting 1kb 5' UTRs Wed Nov 20 17:37:40 PST 2024
    ## Subtract portions of 5' UTRs that overlap nearby genes Wed Nov 20 17:37:40 PST 2024
    ## 5' UTRs identified! Wed Nov 20 17:37:40 PST 2024
    ## Extracting 5' UTR sequences Wed Nov 20 17:37:40 PST 2024
    ## Sequence extraction complete! Wed Nov 20 17:37:41 PST 2024

Check We expect, for an mRNA in + sense (forward strand), the 3’UTR
region to be the 1000bp immediately following the mRNA, and the 5’UTR to
be the 1000bp immediately before the mRNA (For a - sense mRNA, we would
expect the opposite, with 3’ preceding the mRNA and 5’ following it.)

``` bash
cd ../output/15-Apul-annotate-UTRs

head -1 Apul_GFFannotation.mRNA_sorted.gff
head -1 Apul.GFFannotation.3UTR_1kb_corrected.gff
head -1 Apul.GFFannotation.5UTR_1kb_corrected.gff
echo ""
head -2 Apul_GFFannotation.mRNA_sorted.gff | tail -1
head -2 Apul.GFFannotation.3UTR_1kb_corrected.gff | tail -1
head -2 Apul.GFFannotation.5UTR_1kb_corrected.gff | tail -1
echo ""
head -3 Apul_GFFannotation.mRNA_sorted.gff | tail -1
head -3 Apul.GFFannotation.3UTR_1kb_corrected.gff | tail -1
head -3 Apul.GFFannotation.5UTR_1kb_corrected.gff | tail -1
```

    ## ntLink_7 funannotate mRNA    79  4679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical protein;
    ## ntLink_7 funannotate 3prime_UTR  4680    5679    .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical
    ## ntLink_7 funannotate five_prime_UTR  1   78  .   +   .   ID=FUN_002303-T1;Parent=FUN_002303;product=hypothetical
    ## 
    ## ntLink_7 funannotate mRNA    12385   16904   .   -   .   ID=FUN_002304-T1;Parent=FUN_002304;product=hypothetical protein;
    ## ntLink_7 funannotate 3prime_UTR  11385   12384   .   -   .   ID=FUN_002304-T1;Parent=FUN_002304;product=hypothetical
    ## ntLink_7 funannotate five_prime_UTR  16905   17904   .   -   .   ID=FUN_002304-T1;Parent=FUN_002304;product=hypothetical
    ## 
    ## ntLink_7 funannotate mRNA    18480   24187   .   +   .   ID=FUN_002305-T1;Parent=FUN_002305;product=hypothetical protein;
    ## ntLink_7 funannotate 3prime_UTR  24188   25187   .   +   .   ID=FUN_002305-T1;Parent=FUN_002305;product=hypothetical
    ## ntLink_7 funannotate five_prime_UTR  17480   18479   .   +   .   ID=FUN_002305-T1;Parent=FUN_002305;product=hypothetical

We’re good!