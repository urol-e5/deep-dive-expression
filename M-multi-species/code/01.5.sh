#!/bin/bash

# Long non-coding RNA discovery pipeline with checkpoints
# Author: Steven Roberts (revised with checkpoints)

# Set Variables
DATA_DIR="../data/01.5-lncRNA"
OUTPUT_DIR="../output/01.5-lncRNA"
THREADS="42"
FASTQ_SOURCE="https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_meandrina/trimmed/"
FASTQ_SUFFIX="fastq.gz"
GENOME_SOURCE="https://owl.fish.washington.edu/halfshell/genomic-databank/Pocillopora_meandrina_HIv1.assembly.fasta"
GTF_SOURCE="https://raw.githubusercontent.com/urol-e5/timeseries_molecular/d5f546705e3df40558eeaa5c18b122c79d2f4453/F-Ptua/data/Pocillopora_meandrina_HIv1.genes-validated.gtf"
GFF_SOURCE="https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/F-Ptuh/data/Pocillopora_meandrina_HIv1.genes-validated.gff3"
GFFPATTERN='class_code "u"|class_code "x"|class_code "o"|class_code "i"'

# Tool Paths
HISAT2_DIR="/home/shared/hisat2-2.2.1/"
SAMTOOLS_DIR="/home/shared/samtools-1.12/"
STRINGTIE_DIR="/home/shared/stringtie-2.2.1.Linux_x86_64"
GFFCOMPARE_DIR="/home/shared/gffcompare-0.12.6.Linux_x86_64"
BEDTOOLS_DIR="/home/shared/bedtools2/bin"
CPC2_DIR="/home/shared/CPC2_standalone-1.0.1"
CONDA_PATH="/opt/anaconda/anaconda3/bin/conda"

GENOME_FASTA="$DATA_DIR/genome.fasta"
GENOME_GTF="$DATA_DIR/genome.gtf"
GENOME_GFF="$DATA_DIR/genome.gff"
FASTQ_DIR="$DATA_DIR/fastq"
GENOME_INDEX="$OUTPUT_DIR/genome.index"

# Create directories checkpoint
if [ ! -d "$DATA_DIR" ]; then
    mkdir -p "$DATA_DIR"
fi

if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

if [ ! -d "$FASTQ_DIR" ]; then
    mkdir -p "$FASTQ_DIR"
fi

echo "Directories created or already exist."

# Download Genome and Reads checkpoint
if [ ! -f "$GENOME_FASTA" ]; then
    echo "Downloading FASTQ files..."
    wget -nv -r --no-directories --no-parent -P "$FASTQ_DIR" -A "*${FASTQ_SUFFIX}" "$FASTQ_SOURCE"
    echo "Downloading genome and feature files..."
    curl -o "$GENOME_FASTA" "$GENOME_SOURCE"
    curl -o "$GENOME_GTF" "$GTF_SOURCE"
    curl -o "$GENOME_GFF" "$GFF_SOURCE"
else
    echo "Genome and reads already downloaded. Skipping download."
fi

# Validate downloads checkpoint
if [ ! -f "$OUTPUT_DIR/validate.done" ]; then
    output_fasta=$(head -1 "$GENOME_FASTA")
    output_gff=$(head -2 "$GENOME_GFF")
    output_gtf=$(head -1 "$GENOME_GTF")
    if [[ "$output_fasta" == *html* || "$output_gff" == *html* || "$output_gtf" == *html* ]]; then
        echo "FAIL - Downloaded an HTML instead of a genome feature file!"
        exit 1
    else
        echo "Download validation passed:"
        echo "$output_fasta"
        echo "$output_gff"
        echo "$output_gtf"
        touch "$OUTPUT_DIR/validate.done"
    fi
else
    echo "Download validation already complete. Skipping."
fi

# HISAT Indexing checkpoint
if [ ! -f "$OUTPUT_DIR/hisat_index.done" ]; then
    echo "Extracting exons and splice sites, and building HISAT2 index..."
    "${HISAT2_DIR}hisat2_extract_exons.py" "$GENOME_GTF" > "$OUTPUT_DIR/exon.txt"
    "${HISAT2_DIR}hisat2_extract_splice_sites.py" "$GENOME_GTF" > "$OUTPUT_DIR/splice_sites.txt"
    "${HISAT2_DIR}hisat2-build" -p "$THREADS" "$GENOME_FASTA" "$GENOME_INDEX" --exon "$OUTPUT_DIR/exon.txt" --ss "$OUTPUT_DIR/splice_sites.txt" 2> "$OUTPUT_DIR/hisat2-build_stats.txt"
    touch "$OUTPUT_DIR/hisat_index.done"
else
    echo "HISAT indexing already complete. Skipping."
fi

# Align Reads checkpoint
if [ ! -f "$OUTPUT_DIR/align_reads.done" ]; then
    echo "Aligning reads with HISAT2..."
    for r2 in "$FASTQ_DIR"/*_R2_*."$FASTQ_SUFFIX"; do
        base=$(basename "$r2")
        sample="${base%%_R2_*}"
        r1="${r2/_R2_/_R1_}"
        output="$OUTPUT_DIR/${sample}.sam"
        "${HISAT2_DIR}hisat2" -x "$GENOME_INDEX" -p "$THREADS" -1 "$r1" -2 "$r2" -S "$output"
    done
    touch "$OUTPUT_DIR/align_reads.done"
else
    echo "Reads alignment already complete. Skipping."
fi

# Convert SAM to BAM checkpoint
if [ ! -f "$OUTPUT_DIR/sam_to_bam.done" ]; then
    echo "Converting SAM files to sorted BAM files..."
    for samfile in "$OUTPUT_DIR"/*.sam; do
      bamfile="${samfile%.sam}.bam"
      sorted_bamfile="${samfile%.sam}.sorted.bam"
      "${SAMTOOLS_DIR}samtools" view -bS -@ "$THREADS" "$samfile" > "$bamfile"
      "${SAMTOOLS_DIR}samtools" sort -@ "$THREADS" "$bamfile" -o "$sorted_bamfile"
      "${SAMTOOLS_DIR}samtools" index -@ "$THREADS" "$sorted_bamfile"
    done
    touch "$OUTPUT_DIR/sam_to_bam.done"
else
    echo "SAM to BAM conversion already complete. Skipping."
fi

# StringTie Assembly checkpoint
if [ ! -f "$OUTPUT_DIR/stringtie_assembly.done" ]; then
    echo "Running StringTie assembly..."
    find "$OUTPUT_DIR" -name "*sorted.bam" \
    | xargs -n 1 basename -s .sorted.bam \
    | xargs -I{} "${STRINGTIE_DIR}/stringtie" \
        -p "$THREADS" \
        -G "$GENOME_GFF" \
        -o "$OUTPUT_DIR/{}.gtf" \
        "$OUTPUT_DIR/{}.sorted.bam"
    touch "$OUTPUT_DIR/stringtie_assembly.done"
else
    echo "StringTie assembly already complete. Skipping."
fi

# Merge GTFs checkpoint
if [ ! -f "$OUTPUT_DIR/merge_gtf.done" ]; then
    echo "Merging GTF files with StringTie..."
    "${STRINGTIE_DIR}/stringtie" --merge -G "$GENOME_GFF" -o "$OUTPUT_DIR/stringtie_merged.gtf" "$OUTPUT_DIR/"*.gtf
    touch "$OUTPUT_DIR/merge_gtf.done"
else
    echo "GTF merging already complete. Skipping."
fi

# GFFCompare checkpoint
if [ ! -f "$OUTPUT_DIR/gffcompare.done" ]; then
    echo "Running GFFCompare..."
    "${GFFCOMPARE_DIR}gffcompare" -r "$GENOME_GFF" -o "$OUTPUT_DIR/gffcompare_merged" "$OUTPUT_DIR/stringtie_merged.gtf"
    touch "$OUTPUT_DIR/gffcompare.done"
else
    echo "GFFCompare already complete. Skipping."
fi

# Filter GTF checkpoint
if [ ! -f "$OUTPUT_DIR/filter_gtf.done" ]; then
    echo "Filtering GTF for lncRNA candidates..."
    awk '$3 == "transcript" && $1 !~ /^#/' "$OUTPUT_DIR/gffcompare_merged.annotated.gtf" | \
    grep -E "$GFFPATTERN" | \
    awk '($5 - $4 > 199) || ($4 - $5 > 199)' > "$OUTPUT_DIR/lncRNA_candidates.gtf"
    touch "$OUTPUT_DIR/filter_gtf.done"
else
    echo "GTF filtering already complete. Skipping."
fi

# Extract sequences using Bedtools checkpoint
if [ ! -f "$OUTPUT_DIR/bedtools_done.done" ]; then
    echo "Extracting sequences with Bedtools..."
    "${BEDTOOLS_DIR}bedtools" getfasta -fi "$GENOME_FASTA" -bed "$OUTPUT_DIR/lncRNA_candidates.gtf" -fo "$OUTPUT_DIR/lncRNA_candidates.fasta" -name -split
    touch "$OUTPUT_DIR/bedtools_done.done"
else
    echo "Sequence extraction already complete. Skipping."
fi

# Run CPC2 checkpoint
if [ ! -f "$OUTPUT_DIR/cpc2.done" ]; then
    echo "Running CPC2..."
    "${CPC2_DIR}CPC2.py" -i "$OUTPUT_DIR/lncRNA_candidates.fasta" -o "$OUTPUT_DIR/CPC2"
    touch "$OUTPUT_DIR/cpc2.done"
else
    echo "CPC2 analysis already complete. Skipping."
fi

# Filter noncoding transcripts checkpoint
if [ ! -f "$OUTPUT_DIR/filter_noncoding.done" ]; then
    echo "Filtering noncoding transcripts..."
    awk '$8 == "noncoding" {print $1}' "$OUTPUT_DIR/CPC2.txt" > "$OUTPUT_DIR/noncoding_transcripts_ids.txt"
    touch "$OUTPUT_DIR/filter_noncoding.done"
else
    echo "Filtering noncoding transcripts already complete. Skipping."
fi

# Subset fasta checkpoint
if [ ! -f "$OUTPUT_DIR/subset_fasta.done" ]; then
    echo "Subsetting fasta for noncoding transcripts..."
    "${SAMTOOLS_DIR}samtools" faidx "$OUTPUT_DIR/lncRNA_candidates.fasta" -r "$OUTPUT_DIR/noncoding_transcripts_ids.txt" > "$OUTPUT_DIR/lncRNA.fasta"
    touch "$OUTPUT_DIR/subset_fasta.done"
else
    echo "Fasta subsetting already complete. Skipping."
fi

# Convert to GTF checkpoint (assumes lncRNA.bed is available)
if [ ! -f "$OUTPUT_DIR/convert_gtf.done" ]; then
    echo "Converting BED to GTF..."
    awk 'BEGIN{OFS="\t"; count=1} {printf "%s\t.\tlncRNA\t%d\t%d\t.\t+\t.\tgene_id \"lncRNA_%03d\";\n", $1, $2, $3, count++;}' "$OUTPUT_DIR/lncRNA.bed" > "$OUTPUT_DIR/lncRNA.gtf"
    touch "$OUTPUT_DIR/convert_gtf.done"
else
    echo "GTF conversion already complete. Skipping."
fi

# Summary checkpoint
if [ ! -f "$OUTPUT_DIR/summary.done" ]; then
    echo "Generating summary..."
    awk 'BEGIN { total_entries=0; min_length=1e9; max_length=0; sum_length=0; } { total_entries++; start=$4; end=$5; gene_length=end-start+1; if (gene_length<min_length) min_length=gene_length; if (gene_length>max_length) max_length=gene_length; sum_length+=gene_length; } END { print "Total entries:", total_entries; print "Min length:", min_length; print "Max length:", max_length; print "Avg length:", sum_length/total_entries; }' "$OUTPUT_DIR/lncRNA.gtf"
    touch "$OUTPUT_DIR/summary.done"
else
    echo "Summary already complete. Skipping."
fi