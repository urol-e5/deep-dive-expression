01-Peve-RNA-trimming-FastQC
================
Kathleen Durkin
2024-08-20

- <a href="#1-create-a-bash-variables-file"
  id="toc-1-create-a-bash-variables-file">1 Create a Bash variables
  file</a>
- <a href="#2-raw-reads" id="toc-2-raw-reads">2 Raw reads</a>
  - <a href="#21-download-raw-rna-seq-reads"
    id="toc-21-download-raw-rna-seq-reads">2.1 Download raw RNA-seq
    reads</a>
  - <a href="#22-verify-raw-read-checksums"
    id="toc-22-verify-raw-read-checksums">2.2 Verify raw read checksums</a>
  - <a href="#23-fastqcmultiqc-on-raw-reads"
    id="toc-23-fastqcmultiqc-on-raw-reads">2.3 FastQC/MultiQC on raw
    reads</a>
- <a href="#3-trimmed-reads" id="toc-3-trimmed-reads">3 Trimmed reads</a>
  - <a href="#31-download-trimmed-rna-seq-reads"
    id="toc-31-download-trimmed-rna-seq-reads">3.1 Download trimmed RNA-seq
    reads</a>
  - <a href="#32-verify-raw-read-checksums"
    id="toc-32-verify-raw-read-checksums">3.2 Verify raw read checksums</a>
  - <a href="#33-fastqcmultiqc-on-trimmed-reads"
    id="toc-33-fastqcmultiqc-on-trimmed-reads">3.3 FastQC/MultiQC on trimmed
    reads</a>
- <a href="#4-summary" id="toc-4-summary">4 Summary</a>

Code for trimming and QCing RNAseq data, to be used on *Porites
evermanni*.

For now I’m just going to QC the [raw
reads](https://owl.fish.washington.edu/nightingales/P_evermanni/30-789513166/)
and the [trimmed
reads](https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_evermanni/trimmed/)
generated as a part of the
[deep-dive](https://github.com/urol-e5/deep-dive/tree/main) project. If
additional/different trimming needs to be done for this expression work,
it will be performed here.

Inputs:

- RNA-seq gzipped FastQs (e.g. `*.fastq.gz`)

Outputs:

- [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  HTML reports for raw and trimmed reads.

- [`MultiQC`](https://multiqc.info/) HTML summaries of
  [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  for raw and trimmed reads.

[Sequencing
report](https://github.com/urol-e5/deep-dive/wiki/Azenta_30-789513166_Data_Report.html)
[Trimming
details](https://robertslab.github.io/sams-notebook/posts/2023/2023-05-19-FastQ-QC-and-Trimming---E5-Coral-RNA-seq-Data-for-A.pulchra-P.evermanni-and-P.meandrina-Using-FastQC-fastp-and-MultiQC-on-Mox/)

------------------------------------------------------------------------

# 1 Create a Bash variables file

This allows usage of Bash variables (e.g. paths to common directories)
across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export expression_dir=/home/shared/8TB_HDD_02/shedurkin/expression'
echo 'export output_dir_top=${expression_dir}/E-Peve/output/01-Peve-RNA-trimming-FastQC'
echo 'export raw_fastqc_dir=${output_dir_top}/raw-fastqc'
echo 'export raw_reads_dir=${expression_dir}/E-Peve/data/01-Peve-RNA-trimming-FastQC/raw-reads'
echo 'export raw_reads_url="https://owl.fish.washington.edu/nightingales/P_evermanni/30-789513166/"'
echo 'export trimmed_fastqc_dir=${output_dir_top}/trimmed-fastqc'
echo 'export trimmed_reads_dir=${expression_dir}/E-Peve/data/01-Peve-RNA-trimming-FastQC/trimmed-reads'
echo 'export trimmed_reads_url="https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_evermanni/trimmed/"'
echo ""

echo "# Paths to programs"
echo 'export fastqc=/home/shared/FastQC-0.12.1/fastqc'
echo 'export multiqc=/home/sam/programs/mambaforge/bin/multiqc'
echo 'export flexbar=/home/shared/flexbar-3.5.0-linux/flexbar'
echo ""

echo "# Set FastQ filename patterns"
echo "export fastq_pattern='*.fastq.gz'"
echo "export R1_fastq_pattern='*_R1_*.fastq.gz'"
echo "export R2_fastq_pattern='*_R2_*.fastq.gz'"
echo ""

echo "# Set number of CPUs to use"
echo 'export threads=20'
echo ""

echo "# Input/output files"
echo 'export raw_checksums=checksums.md5'
echo 'export trimmed_checksums=trimmed_fastq_checksums.md5'
echo ""



echo "## Inititalize arrays"
echo 'export fastq_array_R1=()'
echo 'export fastq_array_R2=()'
echo 'export raw_fastqs_array=()'
echo 'export R1_names_array=()'
echo 'export R2_names_array=()'
echo 'export trimmed_fastqs_array=()'
echo ""

echo "# Programs associative array"
echo "declare -A programs_array"
echo "programs_array=("
echo '[fastqc]="${fastqc}" \'
echo '[multiqc]="${multiqc}" \'
echo '[flexbar]="${flexbar}"'
echo ")"
} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # Data directories
    export expression_dir=/home/shared/8TB_HDD_02/shedurkin/expression
    export output_dir_top=${expression_dir}/E-Peve/output/01-Peve-RNA-trimming-FastQC
    export raw_fastqc_dir=${output_dir_top}/raw-fastqc
    export raw_reads_dir=${expression_dir}/E-Peve/data/01-Peve-RNA-trimming-FastQC/raw-reads
    export raw_reads_url="https://owl.fish.washington.edu/nightingales/P_evermanni/30-789513166/"
    export trimmed_fastqc_dir=${output_dir_top}/trimmed-fastqc
    export trimmed_reads_dir=${expression_dir}/E-Peve/data/01-Peve-RNA-trimming-FastQC/trimmed-reads
    export trimmed_reads_url="https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_evermanni/trimmed/"

    # Paths to programs
    export fastqc=/home/shared/FastQC-0.12.1/fastqc
    export multiqc=/home/sam/programs/mambaforge/bin/multiqc
    export flexbar=/home/shared/flexbar-3.5.0-linux/flexbar

    # Set FastQ filename patterns
    export fastq_pattern='*.fastq.gz'
    export R1_fastq_pattern='*_R1_*.fastq.gz'
    export R2_fastq_pattern='*_R2_*.fastq.gz'

    # Set number of CPUs to use
    export threads=20

    # Input/output files
    export raw_checksums=checksums.md5
    export trimmed_checksums=trimmed_fastq_checksums.md5

    ## Inititalize arrays
    export fastq_array_R1=()
    export fastq_array_R2=()
    export raw_fastqs_array=()
    export R1_names_array=()
    export R2_names_array=()
    export trimmed_fastqs_array=()

    # Programs associative array
    declare -A programs_array
    programs_array=(
    [fastqc]="${fastqc}" \
    [multiqc]="${multiqc}" \
    [flexbar]="${flexbar}"
    )

# 2 Raw reads

## 2.1 Download raw RNA-seq reads

Reads are downloaded from:
<https://owl.fish.washington.edu/nightingales/P_evermanni/30-789513166/>

The `--cut-dirs 3` command cuts the preceding directory structure
(i.e. `nightingales/P_evermanni/30-789513166/`) so that we just end up
with the reads.

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${raw_reads_dir} \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 3 \
--no-host-directories \
--no-parent \
--quiet \
--accept ${fastq_pattern} ${raw_reads_url}
```

``` bash
# Load bash variables into memory
source .bashvars

ls -lh "${raw_reads_dir}"
```

    total 69G
    -rw-r--r-- 1 shedurkin labmembers 1.3K May 17  2023 checksums.md5
    -rw-r--r-- 1 shedurkin labmembers 3.5G May 16  2023 POR-71-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.6G May 16  2023 POR-71-S1-TP2_R2_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.4G May 16  2023 POR-73-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.5G May 16  2023 POR-73-S1-TP2_R2_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.4G May 16  2023 POR-76-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.5G May 16  2023 POR-76-S1-TP2_R2_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.3G May 16  2023 POR-79-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.4G May 16  2023 POR-79-S1-TP2_R2_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.4G May 16  2023 POR-82-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.6G May 16  2023 POR-82-S1-TP2_R2_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.5G May 17  2023 RNA-POR-71-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.6G May 17  2023 RNA-POR-71-S1-TP2_R2_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.4G May 17  2023 RNA-POR-73-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.5G May 17  2023 RNA-POR-73-S1-TP2_R2_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.4G May 17  2023 RNA-POR-76-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.5G May 17  2023 RNA-POR-76-S1-TP2_R2_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.3G May 17  2023 RNA-POR-79-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.4G May 17  2023 RNA-POR-79-S1-TP2_R2_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.4G May 17  2023 RNA-POR-82-S1-TP2_R1_001.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 3.6G May 17  2023 RNA-POR-82-S1-TP2_R2_001.fastq.gz

## 2.2 Verify raw read checksums

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${raw_reads_dir} \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 3 \
--no-host-directories \
--no-parent \
--quiet \
--accept checksums.md5 ${raw_reads_url}

cd "${raw_reads_dir}"

md5sum checksums.md5 --check
```

## 2.3 FastQC/MultiQC on raw reads

``` bash
# Load bash variables into memory
source .bashvars

############ RUN FASTQC ############


# Create array of raw FastQs
raw_fastqs_array=(${raw_reads_dir}/${fastq_pattern})

# Pass array contents to new variable as space-delimited list
raw_fastqc_list=$(echo "${raw_fastqs_array[*]}")

echo "Beginning FastQC on raw reads..."
echo ""

# Run FastQC
### NOTE: Do NOT quote raw_fastqc_list
${programs_array[fastqc]} \
--threads ${threads} \
--outdir ${raw_fastqc_dir} \
--quiet \
${raw_fastqc_list}

echo "FastQC on raw reads complete!"
echo ""

############ END FASTQC ############

############ RUN MULTIQC ############
echo "Beginning MultiQC on raw FastQC..."
echo ""

${programs_array[multiqc]} ${raw_fastqc_dir} -o ${raw_fastqc_dir}

echo ""
echo "MultiQC on raw FastQs complete."
echo ""

############ END MULTIQC ############

echo "Removing FastQC zip files."
echo ""
rm ${raw_fastqc_dir}/*.zip
echo "FastQC zip files removed."
echo ""
```

``` bash
# Load bash variables into memory
source .bashvars

# View directory contents
ls -lh ${raw_fastqc_dir}
```

    total 14M
    drwxr-xr-x 2 shedurkin labmembers 4.0K Aug 20 15:08 multiqc_data
    -rw-r--r-- 1 shedurkin labmembers 1.4M Aug 20 15:08 multiqc_report.html
    -rw-r--r-- 1 shedurkin labmembers 603K Aug 20 15:07 POR-71-S1-TP2_R1_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 603K Aug 20 15:08 POR-71-S1-TP2_R2_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 613K Aug 20 15:08 POR-73-S1-TP2_R1_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 610K Aug 20 15:08 POR-73-S1-TP2_R2_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 609K Aug 20 15:07 POR-76-S1-TP2_R1_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 611K Aug 20 15:07 POR-76-S1-TP2_R2_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 617K Aug 20 15:07 POR-79-S1-TP2_R1_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 609K Aug 20 15:07 POR-79-S1-TP2_R2_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 605K Aug 20 15:07 POR-82-S1-TP2_R1_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 600K Aug 20 15:07 POR-82-S1-TP2_R2_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 603K Aug 20 15:08 RNA-POR-71-S1-TP2_R1_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 603K Aug 20 15:08 RNA-POR-71-S1-TP2_R2_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 613K Aug 20 15:08 RNA-POR-73-S1-TP2_R1_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 610K Aug 20 15:08 RNA-POR-73-S1-TP2_R2_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 609K Aug 20 15:07 RNA-POR-76-S1-TP2_R1_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 611K Aug 20 15:08 RNA-POR-76-S1-TP2_R2_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 617K Aug 20 15:07 RNA-POR-79-S1-TP2_R1_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 609K Aug 20 15:08 RNA-POR-79-S1-TP2_R2_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 605K Aug 20 15:07 RNA-POR-82-S1-TP2_R1_001_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 600K Aug 20 15:07 RNA-POR-82-S1-TP2_R2_001_fastqc.html

# 3 Trimmed reads

## 3.1 Download trimmed RNA-seq reads

Reads are downloaded from:
<https://gannet.fish.washington.edu/Atumefaciens/20230519-E5_coral-fastqc-fastp-multiqc-RNAseq/P_evermanni/trimmed/>

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${trimmed_reads_dir} \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 4 \
--no-host-directories \
--no-parent \
--quiet \
--accept ${fastq_pattern} ${trimmed_reads_url}
```

``` bash
# Load bash variables into memory
source .bashvars

ls -lh "${trimmed_reads_dir}"
```

    total 24G
    drwxr-xr-x 2 shedurkin labmembers 4.0K Aug 20 15:09 multiqc_data
    -rw-r--r-- 1 shedurkin labmembers 2.5G May 19  2023 RNA-POR-71-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.6G May 19  2023 RNA-POR-71-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.2G May 19  2023 RNA-POR-73-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.4G May 19  2023 RNA-POR-73-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.4G May 19  2023 RNA-POR-76-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.5G May 19  2023 RNA-POR-76-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.1G May 19  2023 RNA-POR-79-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.3G May 19  2023 RNA-POR-79-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.6G May 19  2023 RNA-POR-82-S1-TP2_R1_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers 2.7G May 19  2023 RNA-POR-82-S1-TP2_R2_001.fastp-trim.20230519.fastq.gz
    -rw-r--r-- 1 shedurkin labmembers  880 May 19  2023 trimmed_fastq_checksums.md5

## 3.2 Verify raw read checksums

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${trimmed_reads_dir} \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 4 \
--no-host-directories \
--no-parent \
--quiet \
--accept trimmed_fastq_checksums.md5 ${trimmed_reads_url}

cd "${trimmed_reads_dir}"

md5sum trimmed_fastq_checksums.md5 --check
```

## 3.3 FastQC/MultiQC on trimmed reads

``` bash
# Load bash variables into memory
source .bashvars

############ RUN FASTQC ############

### NOTE: Do NOT quote raw_fastqc_list
# Create array of trimmed FastQs
trimmed_fastqs_array=(${trimmed_reads_dir}/${fastq_pattern})

# Pass array contents to new variable as space-delimited list
trimmed_fastqc_list=$(echo "${trimmed_fastqs_array[*]}")

echo "Beginning FastQC on raw reads..."
echo ""

# Run FastQC
${programs_array[fastqc]} \
--threads ${threads} \
--outdir ${trimmed_fastqc_dir} \
--quiet \
${trimmed_fastqc_list}

echo "FastQC on trimmed reads complete!"
echo ""

############ END FASTQC ############

############ RUN MULTIQC ############
echo "Beginning MultiQC on raw FastQC..."
echo ""

${programs_array[multiqc]} ${trimmed_fastqc_dir} -o ${trimmed_fastqc_dir}

echo ""
echo "MultiQC on trimmed FastQs complete."
echo ""

############ END MULTIQC ############

echo "Removing FastQC zip files."
echo ""
rm ${trimmed_fastqc_dir}/*.zip
echo "FastQC zip files removed."
echo ""
```

``` bash
# Load bash variables into memory
source .bashvars

# View directory contents
ls -lh ${trimmed_fastqc_dir}
```

    total 8.1M
    drwxr-xr-x 2 shedurkin labmembers 4.0K Aug 20 15:19 multiqc_data
    -rw-r--r-- 1 shedurkin labmembers 1.3M Aug 20 15:19 multiqc_report.html
    -rw-r--r-- 1 shedurkin labmembers 685K Aug 20 15:19 RNA-POR-71-S1-TP2_R1_001.fastp-trim.20230519_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 689K Aug 20 15:19 RNA-POR-71-S1-TP2_R2_001.fastp-trim.20230519_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 701K Aug 20 15:19 RNA-POR-73-S1-TP2_R1_001.fastp-trim.20230519_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 694K Aug 20 15:19 RNA-POR-73-S1-TP2_R2_001.fastp-trim.20230519_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 696K Aug 20 15:19 RNA-POR-76-S1-TP2_R1_001.fastp-trim.20230519_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 695K Aug 20 15:19 RNA-POR-76-S1-TP2_R2_001.fastp-trim.20230519_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 704K Aug 20 15:19 RNA-POR-79-S1-TP2_R1_001.fastp-trim.20230519_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 692K Aug 20 15:19 RNA-POR-79-S1-TP2_R2_001.fastp-trim.20230519_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 687K Aug 20 15:19 RNA-POR-82-S1-TP2_R1_001.fastp-trim.20230519_fastqc.html
    -rw-r--r-- 1 shedurkin labmembers 681K Aug 20 15:19 RNA-POR-82-S1-TP2_R2_001.fastp-trim.20230519_fastqc.html

# 4 Summary
