12-Peve-WGBS-bismark.Rmd
================
Zoe Dellaert
2025-04-10

- [0.1 This is the downstream methylation analysis of the WGBS data for
  *Porites
  evermanni*](#01-this-is-the-downstream-methylation-analysis-of-the-wgbs-data-for-porites-evermanni)
- [0.2 Important file locations:](#02-important-file-locations)
  - [0.2.1 Note: Most of this code is based on the E5 Time Series
    Molecular code by Steven Roberts
    here](#021-note-most-of-this-code-is-based-on-the-e5-time-series-molecular-code-by-steven-roberts-here)
- [0.3 Generate Bismark Bisulfite
  Genome](#03-generate-bismark-bisulfite-genome)
  - [0.3.1 output:](#031-output)
  - [0.3.2 Compress and generate md5](#032-compress-and-generate-md5)
  - [0.3.3 Output file location: Bismark
    Genome](#033-output-file-location-bismark-genome)
- [0.4 Test parameters](#04-test-parameters)
  - [0.4.1 Results from parameter
    tests:](#041-results-from-parameter-tests)
- [0.5 Align to genome](#05-align-to-genome)
  - [0.5.1 Output file location: All Bismark output
    files](#051-output-file-location-all-bismark-output-files)
- [0.6 Post-alignment code is based once again on Steven’s
  code](#06-post-alignment-code-is-based-once-again-on-stevens-code)
  - [0.6.1 Deduplication, Sorting, and methylation extraction &
    calling](#061-deduplication-sorting-and-methylation-extraction--calling)
  - [0.6.2 View output](#062-view-output)
  - [0.6.3 Make summary reports](#063-make-summary-reports)
  - [0.6.4 Sorting cov files and filtering for coverage and gene
    intersection](#064-sorting-cov-files-and-filtering-for-coverage-and-gene-intersection)
  - [0.6.5 Output file location: All Bismark output
    files](#065-output-file-location-all-bismark-output-files)

## 0.1 This is the downstream methylation analysis of the WGBS data for *Porites evermanni*

Reads were trimmed and QC’d in [this
code](https://github.com/urol-e5/deep-dive-expression/blob/main/E-Peve/code/01.00-E-Peve-WGBS-trimming-cutadapt-FastQC-MultiQC.md)

## 0.2 Important file locations:

1.  [Trimmed WGBS
    Reads](https://gannet.fish.washington.edu/gitrepos/urol-e5/deep-dive-expression/E-Peve/output/01.00-E-Peve-WGBS-trimming-cutadapt-FastQC-MultiQC/)
2.  [Bismark
    Genome](https://gannet.fish.washington.edu/gitrepos/urol-e5/deep-dive-expression/E-Peve/data/Bisulfite_Genome/)
3.  [All Bismark output files (BAMs, .cov files, .bedgraph
    files)](https://gannet.fish.washington.edu/gitrepos/urol-e5/deep-dive-expression/E-Peve/output/12-Peve-WGBS/bismark_cutadapt/)
4.  [Bismark and Qualimap
    MultiQC](https://gannet.fish.washington.edu/gitrepos/urol-e5/deep-dive-expression/E-Peve/output/12-Peve-WGBS/bismark_cutadapt/multiqc_report.html)
    and [base bismark
    report](https://gannet.fish.washington.edu/gitrepos/urol-e5/deep-dive-expression/E-Peve/output/12-Peve-WGBS/bismark_cutadapt/bismark_summary_report.html)

### 0.2.1 Note: Most of this code is based on the [E5 Time Series Molecular](https://github.com/urol-e5/timeseries_molecular) code by Steven Roberts [here](https://github.com/urol-e5/timeseries_molecular/blob/main/D-Apul/code/15.5-Apul-bismark.qmd)

## 0.3 Generate Bismark Bisulfite Genome

``` bash
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=200GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load bowtie2/2.5.2

cd ../data

bismark_genome_preparation --verbose --parallel 10 ./
```

### 0.3.1 output:

``` bash
Using 10 threads for the top and bottom strand indexing processes each, so using 20 cores in total
Writing bisulfite genomes out into a single MFA (multi FastA) file

Bisulfite Genome Indexer version v0.23.1 (last modified: 27 Jan 2021)

Step I - Prepare genome folders - completed



Step II - Genome bisulfite conversions - completed


Bismark Genome Preparation - Step III: Launching the Bowtie 2 indexer
Preparing indexing of CT converted genome in /scratch3/workspace/zdellaert_uri_edu-deep_dive/deep-dive-expression/E-Peve/data/Bisulfite_Genome/CT_conversion/
Building a SMALL index
Preparing indexing of GA converted genome in /scratch3/workspace/zdellaert_uri_edu-deep_dive/deep-dive-expression/E-Peve/data/Bisulfite_Genome/GA_conversion/
Building a SMALL index
Renaming BS_GA.3.bt2.tmp to BS_GA.3.bt2
Renaming BS_GA.4.bt2.tmp to BS_GA.4.bt2
Renaming BS_GA.1.bt2.tmp to BS_GA.1.bt2
Renaming BS_GA.2.bt2.tmp to BS_GA.2.bt2
Renaming BS_GA.rev.1.bt2.tmp to BS_GA.rev.1.bt2
Renaming BS_GA.rev.2.bt2.tmp to BS_GA.rev.2.bt2
Renaming BS_CT.3.bt2.tmp to BS_CT.3.bt2
Renaming BS_CT.4.bt2.tmp to BS_CT.4.bt2
Renaming BS_CT.1.bt2.tmp to BS_CT.1.bt2
Renaming BS_CT.2.bt2.tmp to BS_CT.2.bt2
Renaming BS_CT.rev.1.bt2.tmp to BS_CT.rev.1.bt2
Renaming BS_CT.rev.2.bt2.tmp to BS_CT.rev.2.bt2
```

### 0.3.2 Compress and generate md5

``` bash
cd ../data
tar -czvf Bisulfite_Genome.tar.gz Bisulfite_Genome
md5sum Bisulfite_Genome.tar.gz | tee Bisulfite_Genome.tar.gz.md5
```

``` bash
885c84c233ca313f23b8670903a18db6  Bisulfite_Genome.tar.gz
```

### 0.3.3 Output file location: [Bismark Genome](https://gannet.fish.washington.edu/gitrepos/urol-e5/deep-dive-expression/E-Peve/data/Bisulfite_Genome/)

## 0.4 Test parameters

``` bash
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=30 #split one task over multiple CPU
#SBATCH --array=0-4 #for 5 samples
#SBATCH --mem=100GB
#SBATCH -t 00:30:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load bowtie2/2.5.2

# Set directories and files
reads_dir="../output/01.00-E-Peve-WGBS-trimming-cutadapt-FastQC-MultiQC/"
genome_folder="../data/"

output_dir="../output/12-Peve-WGBS/bismark_paramtest_cutadapt"
checkpoint_file="${output_dir}/completed_samples.log"

# make output directory
mkdir -p ${output_dir}

# Create the checkpoint file if it doesn't exist
touch ${checkpoint_file}

# Get the list of sample files and corresponding sample names
files=(${reads_dir}*_R1_001.fastq.gz)
file="${files[$SLURM_ARRAY_TASK_ID]}"
sample_name=$(basename "$file" "_R1_001.fastq.gz")

# Check if the sample has already been processed
if grep -q "^${sample_name}$" ${checkpoint_file}; then
    echo "Sample ${sample_name} already processed. Skipping..."
    exit 0
fi

# Define log files for stdout and stderr
stdout_log="${output_dir}/${sample_name}_stdout.log"
stderr_log="${output_dir}/${sample_name}_stderr.log"

# Define the array of score_min parameters to test
score_min_params=(
    "L,0,-0.4"
    "L,0,-0.6"
    "L,0,-0.8"
    "L,0,-1.0"
    "L,-1,-0.6"
)

# Loop through each score_min parameter
for score_min in "${score_min_params[@]}"; do
    echo "Running Bismark for sample ${sample_name} with score_min ${score_min}"
    
    # Create a subdirectory for this parameter
    param_output_dir="${output_dir}/${sample_name}_score_${score_min//,/}"
    mkdir -p ${param_output_dir}

    # Run Bismark alignment
    bismark \
        -genome ${genome_folder} \
        -p 8 \
        -u 10000 \
        -score_min ${score_min} \
        --non_directional \
        -1 ${reads_dir}${sample_name}_R1_001.fastq.gz \
        -2 ${reads_dir}${sample_name}_R2_001.fastq.gz \
        -o ${param_output_dir} \
        --basename ${sample_name}_${score_min//,/} \
        2> "${param_output_dir}/${sample_name}-bismark_summary.txt"

    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Sample ${sample_name} with score_min ${score_min} processed successfully."
    else
        echo "Sample ${sample_name} with score_min ${score_min} failed. Check ${stderr_log} for details."
    fi
done

# Mark the sample as completed in the checkpoint file
if [ $? -eq 0 ]; then
    echo ${sample_name} >> ${checkpoint_file}
    echo "All tests for sample ${sample_name} completed."
else
    echo "Sample ${sample_name} encountered errors. Check logs for details."
fi

# Define directories
summary_file="${output_dir}/parameter_comparison_summary.csv"

# Initialize summary file
echo "Sample,Score_Min,Alignment_Rate" > ${summary_file}

# Loop through parameter output directories
for dir in ${output_dir}/*_score_*; do
    if [ -d "$dir" ]; then
        # Extract sample name and score_min parameter from directory name
        sample_name=$(basename "$dir" | cut -d '_' -f1-3)
        score_min=$(basename "$dir" | grep -o "score_.*" | sed 's/score_//; s/_/,/g')

        # Locate the summary file
        summary_file_path="${dir}/${sample_name}_${score_min}_PE_report.txt"

        # Extract metrics
        mapping=$(grep "Mapping efficiency:" ${summary_file_path} | awk '{print "mapping efficiency ", $3}')
        
        # Append to the summary file
        echo "${sample_name},${score_min},${mapping}" >> ${summary_file}
    fi
done
```

### 0.4.1 Results from parameter tests:

| Sample          | Score_Min | Alignment_Rate |
|:----------------|:----------|:---------------|
| trimmed_421_S8  | L0-0.4    | 29.20%         |
| trimmed_421_S8  | L0-0.6    | 36.00%         |
| trimmed_421_S8  | L0-0.8    | 41.00%         |
| trimmed_421_S8  | L0-1.0    | 47.30%         |
| trimmed_421_S8  | L-1-0.6   | 36.30%         |
| trimmed_471_S6  | L0-0.4    | 35.90%         |
| trimmed_471_S6  | L0-0.6    | 42.80%         |
| trimmed_471_S6  | L0-0.8    | 48.40%         |
| trimmed_471_S6  | L0-1.0    | 54.70%         |
| trimmed_471_S6  | L-1-0.6   | 43.00%         |
| trimmed_487_S9  | L0-0.4    | 34.60%         |
| trimmed_487_S9  | L0-0.6    | 41.40%         |
| trimmed_487_S9  | L0-0.8    | 47.50%         |
| trimmed_487_S9  | L0-1.0    | 54.20%         |
| trimmed_487_S9  | L-1-0.6   | 41.80%         |
| trimmed_489_S10 | L0-0.4    | 32.50%         |
| trimmed_489_S10 | L0-0.6    | 40.30%         |
| trimmed_489_S10 | L0-0.8    | 47.00%         |
| trimmed_489_S10 | L0-1.0    | 54.20%         |
| trimmed_489_S10 | L-1-0.6   | 40.70%         |
| trimmed_491_S7  | L0-0.4    | 42.00%         |
| trimmed_491_S7  | L0-0.6    | 50.50%         |
| trimmed_491_S7  | L0-0.8    | 56.60%         |
| trimmed_491_S7  | L0-1.0    | 62.10%         |
| trimmed_491_S7  | L-1-0.6   | 50.80%         |

I ran the parameter testing before fixing the file sample names. For
reference:

- 421 = POR-82
- 471 = POR-76
- 487 = POR-71
- 489 = POR-79
- 491 = POR-73

## 0.5 Align to genome

``` bash
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=48 #split one task over multiple CPU
#SBATCH --array=0-4 #for 5 samples
#SBATCH --mem=400GB
#SBATCH -t 48:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80 #email you when job stops and/or fails or is nearing its time limit
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load bowtie2/2.5.2

# Set directories and files
reads_dir="../output/01.00-E-Peve-WGBS-trimming-cutadapt-FastQC-MultiQC/"
genome_folder="../data/"

output_dir="../output/12-Peve-WGBS/bismark_cutadapt"
checkpoint_file="${output_dir}/completed_samples.log"

# make output directory
mkdir -p ${output_dir}

# Create the checkpoint file if it doesn't exist
touch ${checkpoint_file}

# Get the list of sample files and corresponding sample names
files=(${reads_dir}*_R1_001.fastq.gz)
file="${files[$SLURM_ARRAY_TASK_ID]}"
sample_name=$(basename "$file" "_R1_001.fastq.gz")

# Check if the sample has already been processed
if grep -q "^${sample_name}$" ${checkpoint_file}; then
    echo "Sample ${sample_name} already processed. Skipping..."
    exit 0
fi

# Define log files for stdout and stderr
stdout_log="${output_dir}/${sample_name}_stdout.log"
stderr_log="${output_dir}/${sample_name}_stderr.log"

    # Run Bismark alignment
    bismark \
        -genome ${genome_folder} \
        -p 48 \
        -score_min L,0,-1.0 \
        --non_directional \
        -1 ${reads_dir}${sample_name}_R1_001.fastq.gz \
        -2 ${reads_dir}${sample_name}_R2_001.fastq.gz \
        -o ${output_dir} \
        --basename ${sample_name} \
        2> "${output_dir}/${sample_name}-bismark_summary.txt"

# Check if the command was successful
if [ $? -eq 0 ]; then
    # Append the sample name to the checkpoint file
    echo ${sample_name} >> ${checkpoint_file}
    echo "Sample ${sample_name} processed successfully."
else
    echo "Sample ${sample_name} failed. Check ${stderr_log} for details."
fi

# Define directories
summary_file="${output_dir}/alignment_summary.csv"

# Initialize summary file
echo "Sample,Score_Min,Alignment_Rate" > ${summary_file}

# Loop through parameter output directories
for file in ${output_dir}/*_report.txt; do
    # Extract sample name and score_min parameter from directory name
    sample_name=$(basename "$file" | cut -d'_' -f1-3)
    score_min="L0-1.0"

    # Locate the summary file
    summary_file_path="${output_dir}/${sample_name}_PE_report.txt"

    # Extract metrics
    mapping=$(grep "Mapping efficiency:" ${summary_file_path} | awk '{gsub("%", "", $3); print $3}')

    # Append to the summary file
    echo "${sample_name},${score_min},${mapping}" >> ${summary_file}
done
```

### 0.5.1 Output file location: [All Bismark output files](https://gannet.fish.washington.edu/gitrepos/urol-e5/deep-dive-expression/E-Peve/output/12-Peve-WGBS/bismark_cutadapt/)

## 0.6 Post-alignment code is based once again on [Steven’s code](https://github.com/urol-e5/timeseries_molecular/blob/main/D-Apul/code/15.5-Apul-bismark.qmd)

### 0.6.1 Deduplication, Sorting, and methylation extraction & calling

``` bash
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=24
#SBATCH --mem=250GB
#SBATCH -t 24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file

# load modules needed
module load uri/main
module load Bismark/0.23.1-foss-2021b
module load parallel/20240822

# set directories
bismark_dir="../output/12-Peve-WGBS/bismark_cutadapt/"
genome_folder="../data/"

### deduplicate bams
find ${bismark_dir}*.bam | \
xargs -n 1 basename | \
sed 's/^trimmed_//' | sed 's/_pe.bam$//' | \
parallel -j 8 deduplicate_bismark \
--bam \
--paired \
--output_dir ${bismark_dir} \
${bismark_dir}trimmed_{}_pe.bam

### methylation extraction

find ${bismark_dir}*deduplicated.bam | xargs -n 1 -I{} \
bismark_methylation_extractor --bedGraph --counts --comprehensive --merge_non_CpG \
--multicore 24 --buffer_size 75% --output ${bismark_dir} "{}"

### methylation call

find ${bismark_dir}*deduplicated.bismark.cov.gz | \
xargs -n 1 basename | \
sed 's/^trimmed_//' | sed 's/_pe.deduplicated.bismark.cov.gz$//' | \
parallel -j 24 coverage2cytosine \
--genome_folder ${genome_folder} \
-o ${bismark_dir}{} \
--merge_CpG \
--zero_based \
${bismark_dir}trimmed_{}_pe.deduplicated.bismark.cov.gz

### sort bams

# change modules

module purge
module load samtools/1.19.2

find ${bismark_dir}*deduplicated.bam | \
xargs -n 1 basename | \
sed 's/^trimmed_//' | sed 's/_pe.deduplicated.bam$//' | \
xargs -I{} samtools \
sort --threads 24 \
${bismark_dir}trimmed_{}_pe.deduplicated.bam \
-o ${bismark_dir}{}.sorted.bam
```

This took just over 6 hours with max memory used per node as 249.99GiB.

### 0.6.2 View output

``` bash
head ${bismark_dir}*evidence.cov
```

### 0.6.3 Make summary reports

``` bash
#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --mem=250GB
#SBATCH -t 08:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails
#SBATCH --error=scripts/outs_errs/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=scripts/outs_errs/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file

module load uri/main
module load Bismark/0.23.1-foss-2021b
module load all/MultiQC/1.12-foss-2021b
module load qualimap/2.2.1

cd ../output/12-Peve-WGBS/bismark_cutadapt/

bam2nuc --genome_folder ../../../data/ *_pe.deduplicated.bam

mkdir -p qualimap/bamqc

for bamFile in *sorted.bam; do
    prefix=$(basename $bamFile .bam)

    qualimap \
    --java-mem-size=29491M \
    bamqc \
     \
    -bam ${bamFile}  \
     \
    -p non-strand-specific \
    --collect-overlap-pairs \
    -outdir qualimap/bamqc/${prefix} \
    -nt 6
done

bismark2report
bismark2summary *pe.bam

multiqc .
```

#### 0.6.3.1 Output file location: [Bismark and Qualimap MultiQC](https://gannet.fish.washington.edu/gitrepos/urol-e5/deep-dive-expression/E-Peve/output/12-Peve-WGBS/bismark_cutadapt/multiqc_report.html) and [base bismark report](https://gannet.fish.washington.edu/gitrepos/urol-e5/deep-dive-expression/E-Peve/output/12-Peve-WGBS/bismark_cutadapt/bismark_summary_report.html)

### 0.6.4 Sorting cov files and filtering for coverage and gene intersection

``` bash
# salloc -p cpu -c 8 --mem 16G

# load modules needed
module load bedtools2/2.31.1

# set directories and files
root_dir="/scratch3/workspace/zdellaert_uri_edu-deep_dive/deep-dive-expression/E-Peve/"
bismark_dir="${root_dir}/output/12-Peve-WGBS/bismark_cutadapt/"
genome_folder="${root_dir}/data/"
gtf_name="Porites_evermanni_validated"

# Sort .cov files
cd ${bismark_dir}
for file in *merged_CpG_evidence.cov
do
  sample=$(basename "${file}" .CpG_report.merged_CpG_evidence.cov)
  bedtools sort -i "${file}" \
  > "${sample}"_sorted.cov
done

# Create bedgraphs for 5X coverage 
for file in *_sorted.cov
do
  sample=$(basename "${file}" _sorted.cov)
  cat "${file}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4}}' \
  > "${sample}"_5x_sorted.bedgraph
done

# Create .tab files for 5x coverage
for file in *_sorted.cov
do
  sample=$(basename "${file}" _sorted.cov)
  cat "${file}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4, $5, $6}}' \
  > "${sample}"_5x_sorted.tab
done

# Intersection of all samples
multiIntersectBed -i *_5x_sorted.tab > CpG.all.samps.5x_sorted.bed

## change number after == to your number of samples
cat CpG.all.samps.5x_sorted.bed | awk '$4 ==5' > CpG.filt.all.samps.5x_sorted.bed 

# Intersection of all samples with gene bodies:
awk '{if ($3 == "transcript") {print}}' "${genome_folder}/${gtf_name}.gtf"  > "${genome_folder}/${gtf_name}_transcripts.gtf"

for i in *5x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b "${genome_folder}/${gtf_name}_transcripts.gtf" \
  > ${i}_gene
done

# Keep only loci intersecting with genes found in all samples
for i in *_5x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.all.samps.5x_sorted.bed \
  > ${i}_CpG_5x_enrichment.bed
done

# Global methylation levels 
for file in *_sorted.cov; do
  sample=$(basename "$file" _sorted.cov)
  awk '{methylated+=$5; unmethylated+=$6} END {print "'$sample'", methylated/(methylated+unmethylated)}' "$file"
done > global_methylation_levels.txt
```

### 0.6.5 Output file location: [All Bismark output files](https://gannet.fish.washington.edu/gitrepos/urol-e5/deep-dive-expression/E-Peve/output/12-Peve-WGBS/bismark_cutadapt/)
