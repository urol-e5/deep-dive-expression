12-Ptuh-WGBS-bismark.Rmd
================
Zoe Dellaert
2025-04-10

- [0.1 Generate Bismark Bisulfite
  Genome](#01-generate-bismark-bisulfite-genome)
  - [0.1.1 output:](#011-output)
  - [0.1.2 Compress and generate md5](#012-compress-and-generate-md5)
- [0.2 Test parameters](#02-test-parameters)
  - [0.2.1 Results from parameter
    tests:](#021-results-from-parameter-tests)
- [0.3 Align to genome](#03-align-to-genome)
- [0.4 Post-alignment code is based once again on Steven’s
  code](#04-post-alignment-code-is-based-once-again-on-stevens-code)
  - [0.4.1 Deduplication, Sorting, and methylation extraction &
    calling](#041-deduplication-sorting-and-methylation-extraction--calling)
  - [0.4.2 View output](#042-view-output)
  - [0.4.3 Make summary reports](#043-make-summary-reports)

## 0.1 Generate Bismark Bisulfite Genome

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

### 0.1.1 output:

``` bash
Using 10 threads for the top and bottom strand indexing processes each, so using 20 cores in total
Writing bisulfite genomes out into a single MFA (multi FastA) file

Bisulfite Genome Indexer version v0.23.1 (last modified: 27 Jan 2021)

Step I - Prepare genome folders - completed



Step II - Genome bisulfite conversions - completed


Bismark Genome Preparation - Step III: Launching the Bowtie 2 indexer
Preparing indexing of CT converted genome in /scratch3/workspace/zdellaert_uri_edu-deep_dive/deep-dive-expression/F-Ptuh/data/Bisulfite_Genome/CT_conversion/
Building a SMALL index
Preparing indexing of GA converted genome in /scratch3/workspace/zdellaert_uri_edu-deep_dive/deep-dive-expression/F-Ptuh/data/Bisulfite_Genome/GA_conversion/
Building a SMALL index
Renaming BS_CT.3.bt2.tmp to BS_CT.3.bt2
Renaming BS_CT.4.bt2.tmp to BS_CT.4.bt2
Renaming BS_CT.1.bt2.tmp to BS_CT.1.bt2
Renaming BS_CT.2.bt2.tmp to BS_CT.2.bt2
Renaming BS_CT.rev.1.bt2.tmp to BS_CT.rev.1.bt2
Renaming BS_CT.rev.2.bt2.tmp to BS_CT.rev.2.bt2
Renaming BS_GA.3.bt2.tmp to BS_GA.3.bt2
Renaming BS_GA.4.bt2.tmp to BS_GA.4.bt2
Renaming BS_GA.1.bt2.tmp to BS_GA.1.bt2
Renaming BS_GA.2.bt2.tmp to BS_GA.2.bt2
Renaming BS_GA.rev.1.bt2.tmp to BS_GA.rev.1.bt2
Renaming BS_GA.rev.2.bt2.tmp to BS_GA.rev.2.bt2
```

### 0.1.2 Compress and generate md5

``` bash
cd ../data
tar -czvf Bisulfite_Genome.tar.gz Bisulfite_Genome
md5sum Bisulfite_Genome.tar.gz | tee Bisulfite_Genome.tar.gz.md5
```

``` bash
1d8399fe951b4f0113ddee438f1bfad1  Bisulfite_Genome.tar.gz
```

## 0.2 Test parameters

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
reads_dir="../output/01.00-F-Ptuh-WGBS-trimming-cutadapt-FastQC-MultiQC/"
genome_folder="../data/"

output_dir="../output/12-Ptuh-WGBS/bismark_paramtest_cutadapt"
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

### 0.2.1 Results from parameter tests:

| Sample                 | Score_Min | Alignment_Rate |
|:-----------------------|:----------|:---------------|
| trimmed_POC-47-TP2_S13 | L0-0.4    | 34.30%         |
| trimmed_POC-47-TP2_S13 | L0-0.6    | 42.00%         |
| trimmed_POC-47-TP2_S13 | L0-0.8    | 48.20%         |
| trimmed_POC-47-TP2_S13 | L0-1.0    | 54.90%         |
| trimmed_POC-47-TP2_S13 | L-1-0.6   | 42.40%         |
| trimmed_POC-48-TP2_S11 | L0-0.4    | 35.20%         |
| trimmed_POC-48-TP2_S11 | L0-0.6    | 44.10%         |
| trimmed_POC-48-TP2_S11 | L0-0.8    | 51.40%         |
| trimmed_POC-48-TP2_S11 | L0-1.0    | 58.80%         |
| trimmed_POC-48-TP2_S11 | L-1-0.6   | 44.60%         |
| trimmed_POC-50-TP2_S14 | L0-0.4    | 32.80%         |
| trimmed_POC-50-TP2_S14 | L0-0.6    | 40.40%         |
| trimmed_POC-50-TP2_S14 | L0-0.8    | 45.70%         |
| trimmed_POC-50-TP2_S14 | L0-1.0    | 51.40%         |
| trimmed_POC-50-TP2_S14 | L-1-0.6   | 40.90%         |
| trimmed_POC-53-TP2_S15 | L0-0.4    | 31.00%         |
| trimmed_POC-53-TP2_S15 | L0-0.6    | 38.60%         |
| trimmed_POC-53-TP2_S15 | L0-0.8    | 45.20%         |
| trimmed_POC-53-TP2_S15 | L0-1.0    | 52.00%         |
| trimmed_POC-53-TP2_S15 | L-1-0.6   | 38.90%         |
| trimmed_POC-57-TP2_S12 | L0-0.4    | 34.70%         |
| trimmed_POC-57-TP2_S12 | L0-0.6    | 43.80%         |
| trimmed_POC-57-TP2_S12 | L0-0.8    | 50.70%         |
| trimmed_POC-57-TP2_S12 | L0-1.0    | 58.10%         |
| trimmed_POC-57-TP2_S12 | L-1-0.6   | 44.20%         |

## 0.3 Align to genome

``` bash
#!/usr/bin/env bash
#SBATCH --ntasks=1 --cpus-per-task=24 #split one task over multiple CPU
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
reads_dir="../output/01.00-F-Ptuh-WGBS-trimming-cutadapt-FastQC-MultiQC/"
genome_folder="../data/"

output_dir="../output/12-Ptuh-WGBS/bismark_cutadapt"
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
        -p 8 \
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

## 0.4 Post-alignment code is based once again on [Steven’s code](https://github.com/urol-e5/timeseries_molecular/blob/main/D-Apul/code/15.5-Apul-bismark.qmd)

### 0.4.1 Deduplication, Sorting, and methylation extraction & calling

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
bismark_dir="../output/12-Ptuh-WGBS/bismark_cutadapt/"
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

This took just over 5 hours with max memory used per node as 249.99GiB.

### 0.4.2 View output

``` bash
head ${bismark_dir}*evidence.cov
```

### 0.4.3 Make summary reports

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

cd ../output/12-Ptuh-WGBS/bismark_cutadapt/

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
