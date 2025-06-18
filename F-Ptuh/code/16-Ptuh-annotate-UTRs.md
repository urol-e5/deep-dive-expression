16-Ptuh-annotate-UTRs
================
Kathleen Durkin
2025-06-18

For miRNA target prediction, we expect most binding to occur in the
3’UTR regions. However, our reference genome is not annotated with UTRs.
We need to annotate those manually. I will also generate annotation
files for 5’UTR regions and CDS regions

First let’s take a look at what our reference
`Pocillopora_meandrina_HIv1.genes-validated.gff3.gff3` file looks like

Download files if necessary

``` bash
curl -L https://owl.fish.washington.edu/halfshell/genomic-databank/Pocillopora_meandrina_HIv1.assembly.fasta -o ../data/Pocillopora_meandrina_HIv1.assembly.fasta
```

    ##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
    ##                                  Dload  Upload   Total   Spent    Left  Speed
    ##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  5  359M    5 19.7M    0     0  20.8M      0  0:00:17 --:--:--  0:00:17 20.7M 12  359M   12 43.7M    0     0  22.1M      0  0:00:16  0:00:01  0:00:15 22.1M 38  359M   38  138M    0     0  44.2M      0  0:00:08  0:00:03  0:00:05 44.2M 61  359M   61  221M    0     0  56.0M      0  0:00:06  0:00:03  0:00:03 56.0M 85  359M   85  306M    0     0  59.8M      0  0:00:05  0:00:05 --:--:-- 59.8M100  359M  100  359M    0     0  62.0M      0  0:00:05  0:00:05 --:--:-- 70.1M

``` bash
grep -v '^#' ../data/Pocillopora_meandrina_HIv1.genes-validated.gff3 | cut -s -f 3 | sort | uniq -c | sort -rn
echo ""
head -10 ../data/Pocillopora_meandrina_HIv1.genes-validated.gff3
```

    ##  208535 exon
    ##  208535 CDS
    ##   31840 mRNA
    ##   31840 gene
    ## 
    ## ##gff-version 3
    ## ##sequence-region   Pocillopora_meandrina_HIv1___xfSc0000716 887 39392
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    gene    887 6811    .   -   .   ID=gene-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    mRNA    887 6811    .   -   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    CDS 887 973 .   -   0   ID=cds-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1;Parent=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    exon    887 973 .   -   0   ID=exon-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1-1;Parent=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    CDS 1828    1882    .   -   1   ID=cds-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1;Parent=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    exon    1828    1882    .   -   1   ID=exon-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1-2;Parent=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    CDS 2308    2371    .   -   2   ID=cds-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1;Parent=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    exon    2308    2371    .   -   2   ID=exon-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1-3;Parent=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1

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
`Pocillopora_meandrina_HIv1.genes-validated.gff3` is formatted more like
this:

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
`Pocillopora_meandrina_HIv1.genes-validated.gff3` 2) determine sense of
the annotated region (to inform which end is 5’ and which is 3’) 3) make
new lines in gff that covers the 1000bp immediately preceding/following
the mRNA, annotated appropriately as 5’/3’ UTR 4) ensure none of the
newly annotated UTR regions overlap with an existing mRNA region –
remove overlapping region from UTR annotation if necessary

Note: while the majority of miRNA functional binding occurs in the
3’UTR, it is possible for functional binding to occur in other regions,
including the 5’UTR and CDS. I want to run target prediction for all
regions, at least at first, so I want to annotate both the 3’ and 5’ UTR
regions of our genome.

We’ll use a modified version of Jill’s code for this

``` bash

# extract the mRNAs
grep $'\tmRNA\t' ../data/Pocillopora_meandrina_HIv1.genes-validated.gff3 > ../output/16-Ptuh-annotate-UTRs/Ptuh-genome-mRNA_only.gff

# Let's also isolate the CDS while we're at it
grep $'\tCDS\t' ../data/Pocillopora_meandrina_HIv1.genes-validated.gff3 > ../output/16-Ptuh-annotate-UTRs/Ptuh-genome-CDS_only.gff


# check
wc -l ../output/16-Ptuh-annotate-UTRs/Ptuh-genome-mRNA_only.gff
echo ""
head -5 ../output/16-Ptuh-annotate-UTRs/Ptuh-genome-mRNA_only.gff
```

    ## 31840 ../output/16-Ptuh-annotate-UTRs/Ptuh-genome-mRNA_only.gff
    ## 
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    mRNA    887 6811    .   -   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g29951.t1
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    mRNA    7305    12239   .   +   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29952.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g29952.t1
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    mRNA    21622   23112   .   -   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29953.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g29953.t1
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    mRNA    25313   26016   .   -   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29954.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g29954.t1
    ## Pocillopora_meandrina_HIv1___xfSc0000716 AUGUSTUS    mRNA    34771   39392   .   +   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g29955.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g29955.t1

``` bash
cd ../output/16-Ptuh-annotate-UTRs

# Extract scaffold lengths
cat is ../../data/Pocillopora_meandrina_HIv1.assembly.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > Ptuh.Chromosome_lengths.txt

# Extract scaffold names
awk -F" " '{print $1}' Ptuh.Chromosome_lengths.txt > Ptuh.Chromosome_names.txt

# Check
wc -l Ptuh.Chromosome_lengths.txt 
echo ""
head -3 Ptuh.Chromosome_lengths.txt 
echo ""
head -3 Ptuh.Chromosome_names.txt 
```

    ## cat: is: No such file or directory
    ## 212 Ptuh.Chromosome_lengths.txt
    ## 
    ## Pocillopora_meandrina_HIv1___Sc0000000   21651136
    ## Pocillopora_meandrina_HIv1___Sc0000001   21507942
    ## Pocillopora_meandrina_HIv1___Sc0000002   19741974
    ## 
    ## Pocillopora_meandrina_HIv1___Sc0000000
    ## Pocillopora_meandrina_HIv1___Sc0000001
    ## Pocillopora_meandrina_HIv1___Sc0000002

The following code will sort the mRNA gff, extract 1kb down the 3’ end
of mRNA, subtract portions of the 1kb flank (representing the 3’UTR)
from any overlapping mRNA, and make fasta file of the 3’UTRs. It will do
the same for 5’UTRs.

``` bash
cd ../output/16-Ptuh-annotate-UTRs
#export PATH="/home/shared/bedtools2/bin:$PATH"
export PATH="/srlab/programs/bedtools:$PATH"

echo "Sorting gffs by chromosome" $(date)

bedtools sort -faidx Ptuh.Chromosome_names.txt -i Ptuh-genome-mRNA_only.gff > Ptuh_GFFannotation.mRNA_sorted.gff

echo "Sorting complete!" $(date)

echo "Extracting 1kb 3' UTRs" $(date)

bedtools flank -i Ptuh_GFFannotation.mRNA_sorted.gff -g Ptuh.Chromosome_lengths.txt -l 0 -r 1000 -s | \
awk '{gsub("mRNA","3prime_UTR",$3); print $0 }' | \
awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | \
tr ' ' '\t' > Ptuh.GFFannotation.3UTR_1kb.gff

echo "Subtract portions of UTRs that overlap nearby genes" $(date)

bedtools subtract -a Ptuh.GFFannotation.3UTR_1kb.gff -b Ptuh_GFFannotation.mRNA_sorted.gff > Ptuh.GFFannotation.3UTR_1kb_corrected.gff 
echo "3' UTRs identified!" $(date)

echo "Extracting 3' UTR sequences" $(date)

bedtools getfasta -fi ../../data/Pocillopora_meandrina_HIv1.assembly.fasta -bed Ptuh.GFFannotation.3UTR_1kb_corrected.gff -fo Ptuh_3UTR_1kb.fasta -fullHeader


echo "Extracting 1kb 5' UTRs" $(date)

bedtools flank -i Ptuh_GFFannotation.mRNA_sorted.gff -g Ptuh.Chromosome_lengths.txt -l 1000 -r 0 -s | \
awk '{gsub("mRNA","five_prime_UTR",$3); print $0 }' | \
awk '{if($5-$4 > 3)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | \
tr ' ' '\t' > Ptuh.GFFannotation.5UTR_1kb.gff

echo "Subtract portions of 5' UTRs that overlap nearby genes" $(date)

bedtools subtract -a Ptuh.GFFannotation.5UTR_1kb.gff -b Ptuh_GFFannotation.mRNA_sorted.gff > Ptuh.GFFannotation.5UTR_1kb_corrected.gff

echo "5' UTRs identified!" $(date)

echo "Extracting 5' UTR sequences" $(date)

bedtools getfasta -fi ../../data/Pocillopora_meandrina_HIv1.assembly.fasta -bed Ptuh.GFFannotation.5UTR_1kb_corrected.gff -fo Ptuh_5UTR_1kb.fasta -fullHeader

echo "Sequence extraction complete!" $(date)
```

    ## Sorting gffs by chromosome Wed Jun 18 03:37:47 PM PDT 2025
    ## Sorting complete! Wed Jun 18 03:37:47 PM PDT 2025
    ## Extracting 1kb 3' UTRs Wed Jun 18 03:37:47 PM PDT 2025
    ## Subtract portions of UTRs that overlap nearby genes Wed Jun 18 03:37:47 PM PDT 2025
    ## 3' UTRs identified! Wed Jun 18 03:37:47 PM PDT 2025
    ## Extracting 3' UTR sequences Wed Jun 18 03:37:47 PM PDT 2025
    ## Warning: the index file is older than the FASTA file.
    ## Extracting 1kb 5' UTRs Wed Jun 18 03:37:48 PM PDT 2025
    ## Subtract portions of 5' UTRs that overlap nearby genes Wed Jun 18 03:37:48 PM PDT 2025
    ## 5' UTRs identified! Wed Jun 18 03:37:48 PM PDT 2025
    ## Extracting 5' UTR sequences Wed Jun 18 03:37:48 PM PDT 2025
    ## Warning: the index file is older than the FASTA file.
    ## Sequence extraction complete! Wed Jun 18 03:37:49 PM PDT 2025

Check. We expect, for an mRNA in + sense (forward strand), the 3’UTR
region to be the 1000bp immediately following the mRNA, and the 5’UTR to
be the 1000bp immediately before the mRNA (For a - sense mRNA, we would
expect the opposite, with 3’ preceding the mRNA and 5’ following it.)

``` bash
cd ../output/16-Ptuh-annotate-UTRs

head -1 Ptuh_GFFannotation.mRNA_sorted.gff
head -1 Ptuh.GFFannotation.3UTR_1kb_corrected.gff
head -1 Ptuh.GFFannotation.5UTR_1kb_corrected.gff
echo ""
head -2 Ptuh_GFFannotation.mRNA_sorted.gff | tail -1
head -2 Ptuh.GFFannotation.3UTR_1kb_corrected.gff | tail -1
head -2 Ptuh.GFFannotation.5UTR_1kb_corrected.gff | tail -1
echo ""
head -3 Ptuh_GFFannotation.mRNA_sorted.gff | tail -1
head -3 Ptuh.GFFannotation.3UTR_1kb_corrected.gff | tail -1
head -3 Ptuh.GFFannotation.5UTR_1kb_corrected.gff | tail -1
```

    ## Pocillopora_meandrina_HIv1___Sc0000000   AUGUSTUS    mRNA    10771   23652   .   +   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g20902.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g20902.t1
    ## Pocillopora_meandrina_HIv1___Sc0000000   AUGUSTUS    3prime_UTR  23653   24652   .   +   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g20902.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g20902.t1
    ## Pocillopora_meandrina_HIv1___Sc0000000   AUGUSTUS    five_prime_UTR  9771    10770   .   +   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g20902.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g20902.t1
    ## 
    ## Pocillopora_meandrina_HIv1___Sc0000000   AUGUSTUS    mRNA    25129   33764   .   -   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g20903.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g20903.t1
    ## Pocillopora_meandrina_HIv1___Sc0000000   AUGUSTUS    3prime_UTR  24129   25128   .   -   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g20903.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g20903.t1
    ## Pocillopora_meandrina_HIv1___Sc0000000   AUGUSTUS    five_prime_UTR  33765   33907   .   -   .   ID=mrna-Pocillopora_meandrina_HIv1___RNAseq.g20903.t1;Parent=gene-Pocillopora_meandrina_HIv1___RNAseq.g20903.t1
    ## 
    ## Pocillopora_meandrina_HIv1___Sc0000000   AUGUSTUS    mRNA    33908   39296   .   +   .   ID=mrna-Pocillopora_meandrina_HIv1___TS.g25664.t1;Parent=gene-Pocillopora_meandrina_HIv1___TS.g25664.t1
    ## Pocillopora_meandrina_HIv1___Sc0000000   AUGUSTUS    3prime_UTR  39297   40296   .   +   .   ID=mrna-Pocillopora_meandrina_HIv1___TS.g25664.t1;Parent=gene-Pocillopora_meandrina_HIv1___TS.g25664.t1
    ## Pocillopora_meandrina_HIv1___Sc0000000   AUGUSTUS    five_prime_UTR  33765   33907   .   +   .   ID=mrna-Pocillopora_meandrina_HIv1___TS.g25664.t1;Parent=gene-Pocillopora_meandrina_HIv1___TS.g25664.t1

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

FUNid_helper "../output/16-Ptuh-annotate-UTRs/Ptuh-genome-mRNA_only.gff" "../output/16-Ptuh-annotate-UTRs/Ptuh-mRNA-geneIDs.txt"
FUNid_helper "../output/16-Ptuh-annotate-UTRs/Ptuh.GFFannotation.3UTR_1kb_corrected.gff" "../output/16-Ptuh-annotate-UTRs/Ptuh-3UTR-geneIDs.txt"
FUNid_helper "../output/16-Ptuh-annotate-UTRs/Ptuh.GFFannotation.5UTR_1kb_corrected.gff" "../output/16-Ptuh-annotate-UTRs/Ptuh-5UTR-geneIDs.txt"
```

    ## Processed ../output/16-Ptuh-annotate-UTRs/Ptuh-genome-mRNA_only.gff -> ../output/16-Ptuh-annotate-UTRs/Ptuh-mRNA-geneIDs.txt
    ## Processed ../output/16-Ptuh-annotate-UTRs/Ptuh.GFFannotation.3UTR_1kb_corrected.gff -> ../output/16-Ptuh-annotate-UTRs/Ptuh-3UTR-geneIDs.txt
    ## Processed ../output/16-Ptuh-annotate-UTRs/Ptuh.GFFannotation.5UTR_1kb_corrected.gff -> ../output/16-Ptuh-annotate-UTRs/Ptuh-5UTR-geneIDs.txt
