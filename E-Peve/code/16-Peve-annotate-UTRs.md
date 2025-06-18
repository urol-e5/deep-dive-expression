16-Peve-annotate-UTRs
================
Kathleen Durkin
2025-06-18

For miRNA target prediction, we expect most binding to occur in the
3’UTR regions. However, our reference genome is not annotated with UTRs.
We need to annotate those manually. I will also generate annotation
files for 5’UTR regions and CDS regions

First let’s take a look at what our reference
`Porites_evermanni_v1.annot.gff` file looks like

Download files if necessary

``` bash
curl -L https://gannet.fish.washington.edu/seashell/snaps/Porites_evermanni_v1.fa -o ../data/Porites_evermanni_v1.fa
```

    ##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
    ##                                  Dload  Upload   Total   Spent    Left  Speed
    ##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0 12  585M   12 76.1M    0     0  82.4M      0  0:00:07 --:--:--  0:00:07 82.4M 28  585M   28  165M    0     0  94.7M      0  0:00:06  0:00:01  0:00:05 94.7M 46  585M   46  274M    0     0   100M      0  0:00:05  0:00:02  0:00:03  100M 65  585M   65  381M    0     0   101M      0  0:00:05  0:00:03  0:00:02  101M 83  585M   83  489M    0     0  99.3M      0  0:00:05  0:00:04  0:00:01 99.3M 97  585M   97  573M    0     0  99.7M      0  0:00:05  0:00:05 --:--:--  103M100  585M  100  585M    0     0  99.9M      0  0:00:05  0:00:05 --:--:--  102M

``` bash
grep -v '^#' ../data/Porites_evermanni_v1.annot.gff | cut -s -f 3 | sort | uniq -c | sort -rn
echo ""
head -10 ../data/Porites_evermanni_v1.annot.gff
```

    ##  231320 CDS
    ##   40389 mRNA
    ##   15098 UTR
    ## 
    ## Porites_evermani_scaffold_1  Gmove   mRNA    3107    4488    543 -   .   ID=Peve_00000001;Name=Peve_00000001;start=0;stop=1;cds_size=543
    ## Porites_evermani_scaffold_1  Gmove   CDS 3107    3444    .   -   .   Parent=Peve_00000001
    ## Porites_evermani_scaffold_1  Gmove   CDS 4284    4488    .   -   .   Parent=Peve_00000001
    ## Porites_evermani_scaffold_1  Gmove   mRNA    424479  429034  2439.63 -   .   ID=Peve_00000002;Name=Peve_00000002;start=1;stop=1;cds_size=2019
    ## Porites_evermani_scaffold_1  Gmove   CDS 424479  425361  .   -   .   Parent=Peve_00000002
    ## Porites_evermani_scaffold_1  Gmove   CDS 426181  426735  .   -   .   Parent=Peve_00000002
    ## Porites_evermani_scaffold_1  Gmove   CDS 427013  427140  .   -   .   Parent=Peve_00000002
    ## Porites_evermani_scaffold_1  Gmove   CDS 427665  427724  .   -   .   Parent=Peve_00000002
    ## Porites_evermani_scaffold_1  Gmove   CDS 428642  429034  .   -   .   Parent=Peve_00000002
    ## Porites_evermani_scaffold_1  Gmove   mRNA    429394  438909  1570.66 +   .   ID=Peve_00000003;Name=Peve_00000003;start=1;stop=1;cds_size=1458

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

Unfortunately, looking at our gff, it seems that
`Porites_evermanni_v1.annot.gff` is formatted more like this:

    gene1            ================================    ID=gene1
    mRNA1            ================================    ID=mRNA1;Parent=gene1
    CDS1             ======                              Parent=mRNA1
    CDS2                           ==================    Parent=mRNA1

Where an annotated mRNA contains only the CDS regions and internal UTR
regions. That means there is no straightforward way to identify the
actual 5’ and 3’ UTR regions from the gff.

Instead, we can just extract a best-guess region from immediately before
and immediately after the annotated mRNA. Jill Ashey validated during
our `deep-dive` work that a 1000bp region captures the 3’UTR well (see
post and code
[here](https://github.com/JillAshey/JillAshey_Putnam_Lab_Notebook/blob/master/_posts/2024-06-15-e5-deepdive-miRNA-TargetPrediction.md)).

So we want to 1) isolate mRNA regions in
`Porites_evermanni_v1.annot.gff` 2) determine sense of the annotated
region (to inform which end is 5’ and which is 3’) 3) make new lines in
gff that covers the 1000bp immediately preceding/following the mRNA,
annotated appropriately as 5’/3’ UTR 4) ensure none of the newly
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
grep $'\tmRNA\t' ../data/Porites_evermanni_v1.annot.gff > ../output/16-Peve-annotate-UTRs/Peve-genome-mRNA_only.gff

# Let's also isolate the CDS while we're at it
grep $'\tCDS\t' ../data/Porites_evermanni_v1.annot.gff > ../output/16-Peve-annotate-UTRs/Peve-genome-CDS_only.gff


# check
wc -l ../output/16-Peve-annotate-UTRs/Peve-genome-mRNA_only.gff
echo ""
head -5 ../output/16-Peve-annotate-UTRs/Peve-genome-mRNA_only.gff
```

    ## 40389 ../output/16-Peve-annotate-UTRs/Peve-genome-mRNA_only.gff
    ## 
    ## Porites_evermani_scaffold_1  Gmove   mRNA    3107    4488    543 -   .   ID=Peve_00000001;Name=Peve_00000001;start=0;stop=1;cds_size=543
    ## Porites_evermani_scaffold_1  Gmove   mRNA    424479  429034  2439.63 -   .   ID=Peve_00000002;Name=Peve_00000002;start=1;stop=1;cds_size=2019
    ## Porites_evermani_scaffold_1  Gmove   mRNA    429394  438909  1570.66 +   .   ID=Peve_00000003;Name=Peve_00000003;start=1;stop=1;cds_size=1458
    ## Porites_evermani_scaffold_1  Gmove   mRNA    440255  447192  6037.5  -   .   ID=Peve_00000004;Name=Peve_00000004;start=1;stop=1;cds_size=1725
    ## Porites_evermani_scaffold_1  Gmove   mRNA    448045  460031  780.955 +   .   ID=Peve_00000005;Name=Peve_00000005;start=0;stop=0;cds_size=1173

``` bash
cd ../output/16-Peve-annotate-UTRs

# Extract scaffold lengths
cat is ../../data/Porites_evermanni_v1.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Peve.Chromosome_lengths.txt

# Extract scaffold names
awk -F" " '{print $1}' Peve.Chromosome_lengths.txt > Peve.Chromosome_names.txt

# Check
wc -l Peve.Chromosome_lengths.txt 
echo ""
head -3 Peve.Chromosome_lengths.txt 
echo ""
head -3 Peve.Chromosome_names.txt 
```

    ## cat: is: No such file or directory
    ## 8186 Peve.Chromosome_lengths.txt
    ## 
    ## Porites_evermani_scaffold_1  1802771
    ## Porites_evermani_scaffold_2  1721744
    ## Porites_evermani_scaffold_3  1532543
    ## 
    ## Porites_evermani_scaffold_1
    ## Porites_evermani_scaffold_2
    ## Porites_evermani_scaffold_3

The following code will sort the mRNA gff, extract 1kb down the 3’ end
of mRNA, subtract portions of the 1kb flank (representing the 3’UTR)
from any overlapping mRNA, and make fasta file of the 3’UTRs. It will do
the same for 5’UTRs.

``` bash
cd ../output/16-Peve-annotate-UTRs
#export PATH="/home/shared/bedtools2/bin:$PATH"
export PATH="/srlab/programs/bedtools:$PATH"

echo "Sorting gffs by chromosome" $(date)

bedtools sort -faidx Peve.Chromosome_names.txt -i Peve-genome-mRNA_only.gff > Peve_GFFannotation.mRNA_sorted.gff

echo "Sorting complete!" $(date)

echo "Extracting 1kb 3' UTRs" $(date)

bedtools flank -i Peve_GFFannotation.mRNA_sorted.gff -g Peve.Chromosome_lengths.txt -l 0 -r 1000 -s | \
awk '{gsub("mRNA","3prime_UTR",$3); print $0 }' | \
awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | \
tr ' ' '\t' > Peve.GFFannotation.3UTR_1kb.gff

echo "Subtract portions of UTRs that overlap nearby genes" $(date)

bedtools subtract -a Peve.GFFannotation.3UTR_1kb.gff -b Peve_GFFannotation.mRNA_sorted.gff > Peve.GFFannotation.3UTR_1kb_corrected.gff 
echo "3' UTRs identified!" $(date)

echo "Extracting 3' UTR sequences" $(date)

bedtools getfasta -fi ../../data/Porites_evermanni_v1.fa -bed Peve.GFFannotation.3UTR_1kb_corrected.gff -fo Peve_3UTR_1kb.fasta -fullHeader


echo "Extracting 1kb 5' UTRs" $(date)

bedtools flank -i Peve_GFFannotation.mRNA_sorted.gff -g Peve.Chromosome_lengths.txt -l 1000 -r 0 -s | \
awk '{gsub("mRNA","five_prime_UTR",$3); print $0 }' | \
awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | \
tr ' ' '\t' > Peve.GFFannotation.5UTR_1kb.gff

echo "Subtract portions of 5' UTRs that overlap nearby genes" $(date)

bedtools subtract -a Peve.GFFannotation.5UTR_1kb.gff -b Peve_GFFannotation.mRNA_sorted.gff > Peve.GFFannotation.5UTR_1kb_corrected.gff

echo "5' UTRs identified!" $(date)

echo "Extracting 5' UTR sequences" $(date)

bedtools getfasta -fi ../../data/Porites_evermanni_v1.fa -bed Peve.GFFannotation.5UTR_1kb_corrected.gff -fo Peve_5UTR_1kb.fasta -fullHeader

echo "Sequence extraction complete!" $(date)
```

    ## Sorting gffs by chromosome Wed Jun 18 03:15:10 PM PDT 2025
    ## Sorting complete! Wed Jun 18 03:15:10 PM PDT 2025
    ## Extracting 1kb 3' UTRs Wed Jun 18 03:15:10 PM PDT 2025
    ## Subtract portions of UTRs that overlap nearby genes Wed Jun 18 03:15:10 PM PDT 2025
    ## 3' UTRs identified! Wed Jun 18 03:15:10 PM PDT 2025
    ## Extracting 3' UTR sequences Wed Jun 18 03:15:10 PM PDT 2025
    ## Warning: the index file is older than the FASTA file.
    ## Extracting 1kb 5' UTRs Wed Jun 18 03:15:11 PM PDT 2025
    ## Subtract portions of 5' UTRs that overlap nearby genes Wed Jun 18 03:15:11 PM PDT 2025
    ## 5' UTRs identified! Wed Jun 18 03:15:11 PM PDT 2025
    ## Extracting 5' UTR sequences Wed Jun 18 03:15:11 PM PDT 2025
    ## Warning: the index file is older than the FASTA file.
    ## Sequence extraction complete! Wed Jun 18 03:15:12 PM PDT 2025

Check We expect, for an mRNA in + sense (forward strand), the 3’UTR
region to be the 1000bp immediately following the mRNA, and the 5’UTR to
be the 1000bp immediately before the mRNA (For a - sense mRNA, we would
expect the opposite, with 3’ preceding the mRNA and 5’ following it.)

``` bash
cd ../output/16-Peve-annotate-UTRs

head -1 Peve_GFFannotation.mRNA_sorted.gff
head -1 Peve.GFFannotation.3UTR_1kb_corrected.gff
head -1 Peve.GFFannotation.5UTR_1kb_corrected.gff
echo ""
head -2 Peve_GFFannotation.mRNA_sorted.gff | tail -1
head -2 Peve.GFFannotation.3UTR_1kb_corrected.gff | tail -1
head -2 Peve.GFFannotation.5UTR_1kb_corrected.gff | tail -1
echo ""
head -3 Peve_GFFannotation.mRNA_sorted.gff | tail -1
head -3 Peve.GFFannotation.3UTR_1kb_corrected.gff | tail -1
head -3 Peve.GFFannotation.5UTR_1kb_corrected.gff | tail -1
```

    ## Porites_evermani_scaffold_1  Gmove   mRNA    3107    4488    543 -   .   ID=Peve_00000001;Name=Peve_00000001;start=0;stop=1;cds_size=543
    ## Porites_evermani_scaffold_1  Gmove   3prime_UTR  2107    3106    543 -   .   ID=Peve_00000001;Name=Peve_00000001;start=0;stop=1;cds_size=543
    ## Porites_evermani_scaffold_1  Gmove   five_prime_UTR  4489    5488    543 -   .   ID=Peve_00000001;Name=Peve_00000001;start=0;stop=1;cds_size=543
    ## 
    ## Porites_evermani_scaffold_1  Gmove   mRNA    8492    11026   747 +   .   ID=Peve_00000039;Name=Peve_00000039;start=1;stop=1;cds_size=747
    ## Porites_evermani_scaffold_1  Gmove   3prime_UTR  11027   12026   747 +   .   ID=Peve_00000039;Name=Peve_00000039;start=1;stop=1;cds_size=747
    ## Porites_evermani_scaffold_1  Gmove   five_prime_UTR  7492    8491    747 +   .   ID=Peve_00000039;Name=Peve_00000039;start=1;stop=1;cds_size=747
    ## 
    ## Porites_evermani_scaffold_1  Gmove   mRNA    26134   26404   261 +   .   ID=Peve_00000112;Name=Peve_00000112;start=0;stop=0;cds_size=261
    ## Porites_evermani_scaffold_1  Gmove   3prime_UTR  26405   27404   261 +   .   ID=Peve_00000112;Name=Peve_00000112;start=0;stop=0;cds_size=261
    ## Porites_evermani_scaffold_1  Gmove   five_prime_UTR  25134   26133   261 +   .   ID=Peve_00000112;Name=Peve_00000112;start=0;stop=0;cds_size=261

We’re good!

Finally, we may want to use the simpler mRNA-associated gene ids to
denote 3UTR/5UTR regions, rather than the genomic location. Let’s create
helper files for making those associations

``` bash

FUNid_helper() {
    local input_gff="$1"
    local output_file="$2"

    # Check if the input file exists
    if [[ ! -f "$input_gff" ]]; then
        echo "Error: File $input_gff not found!"
        return 1
    fi

    # Process the GFF file
    awk 'BEGIN {OFS="\t"} {
        # Combine location
        location = $1 ":" $4-1 "-" $5;

        # Extract type
        type = $3;

        # Extract and format the attributes (last column)
        attributes = $9;
        gsub(";", "\t", attributes); # Replace semicolons with tabs

        # Print the desired output
        print location, type, attributes;
    }' "$input_gff" > "$output_file"

    echo "Processed $input_gff -> $output_file"
}

FUNid_helper "../output/16-Peve-annotate-UTRs/Peve-genome-mRNA_only.gff" "../output/16-Peve-annotate-UTRs/Peve-mRNA-geneIDs.txt"
FUNid_helper "../output/16-Peve-annotate-UTRs/Peve.GFFannotation.3UTR_1kb_corrected.gff" "../output/16-Peve-annotate-UTRs/Peve-3UTR-geneIDs.txt"
FUNid_helper "../output/16-Peve-annotate-UTRs/Peve.GFFannotation.5UTR_1kb_corrected.gff" "../output/16-Peve-annotate-UTRs/Peve-5UTR-geneIDs.txt"
```

    ## Processed ../output/16-Peve-annotate-UTRs/Peve-genome-mRNA_only.gff -> ../output/16-Peve-annotate-UTRs/Peve-mRNA-geneIDs.txt
    ## Processed ../output/16-Peve-annotate-UTRs/Peve.GFFannotation.3UTR_1kb_corrected.gff -> ../output/16-Peve-annotate-UTRs/Peve-3UTR-geneIDs.txt
    ## Processed ../output/16-Peve-annotate-UTRs/Peve.GFFannotation.5UTR_1kb_corrected.gff -> ../output/16-Peve-annotate-UTRs/Peve-5UTR-geneIDs.txt
