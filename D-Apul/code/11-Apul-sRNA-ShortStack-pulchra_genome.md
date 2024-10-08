11-Apul-sRNA-ShortStack-pulchra_genome
================
Kathleen Durkin
2024-10-02

- [1 Set R variables](#1-set-r-variables)
- [2 Create a Bash variables file](#2-create-a-bash-variables-file)
- [3 Load ShortStack conda
  environment](#3-load-shortstack-conda-environment)
  - [3.1 A.pulchra genome](#31-apulchra-genome)
  - [3.2 Cnidarian+miRBase database](#32-cnidarianmirbase-database)
  - [3.3 Trimmed sRNA-seq reads](#33-trimmed-srna-seq-reads)
- [4 Run ShortStack](#4-run-shortstack)
  - [4.1 Modify genome filename for ShortStack
    compatability](#41-modify-genome-filename-for-shortstack-compatability)
  - [4.2 Excecute ShortStack command](#42-excecute-shortstack-command)
  - [4.3 Check runtime](#43-check-runtime)
- [5 Results](#5-results)
  - [5.1 ShortStack synopsis](#51-shortstack-synopsis)
  - [5.2 Inspect `Results.txt`](#52-inspect-resultstxt)
    - [5.2.1 Directory tree of all ShortStack
      outputs](#521-directory-tree-of-all-shortstack-outputs)
  - [5.3 Visualize](#53-visualize)
- [6 Citations](#6-citations)

Use [ShortStack](https://github.com/MikeAxtell/ShortStack) ([Axtell
2013](#ref-axtell2013a); [Shahid and Axtell 2014](#ref-shahid2014);
[Johnson et al. 2016](#ref-johnson2016a))to perform alignment of sRNAseq
data and annotation of sRNA-producing genes.

This is the same ShortStack analysis as performed in the `deeep-dive`
project However, we now have a genome for our actual study species, *A.
pulchra*, so we can rerun sRNA discovery with the actual *pulchra*
genome for reference, instead of using *A. millipora*. As in
`deep-dive`, we will also include a customized miRBase database,
utilizing cnidarian miRNAs curated by Jill Ashley, which includes
published cnidarian miRNAs:

- [`cnidarian-mirbase-mature-v22.1.fasta`](../../data/cnidarian-mirbase-mature-v22.1.fasta)

------------------------------------------------------------------------

Inputs:

- Requires trimmed sRNAseq files generated in `deep-dive`
  [08.2-Dapul-sRNAseq-trimming-31bp-fastp-merged.Rmd](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/08.2-Dapul-sRNAseq-trimming-31bp-fastp-merged.Rmd)

  - Filenames formatted: `*fastp-adapters-polyG-31bp-merged.fq.gz`

- *A.pulchra* genome FastA. Not currently publicly available (still
  being annotated by colaborators)

Outputs:

- See [ShortStack outputs
  documentation](https://github.com/MikeAxtell/ShortStack#outputs) for
  full list and detailed descriptions.

Software requirements:

- Utilizes a
  [ShortStack](https://github.com/MikeAxtell/ShortStack#installation)
  Conda/Mamba environment, per the installation instructions.

Replace with name of your ShortStack environment and the path to the
corresponding conda installation (find this *after* you’ve activated the
environment).

E.g.

``` bash
# Activate environment
conda activate ShortStack4_env

# Find conda path
which conda
```

------------------------------------------------------------------------

# 1 Set R variables

``` r
shortstack_conda_env_name <- c("ShortStack-4.1.0_env")
shortstack_cond_path <- c("/home/sam/programs/mambaforge/condabin/conda")
```

# 2 Create a Bash variables file

This allows usage of Bash variables across R Markdown chunks.

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Trimmed FastQ naming pattern"
echo "export trimmed_fastqs_pattern='*fastp-adapters-polyG-31bp-merged.fq.gz'"

echo "# Data directories"
echo 'export expression_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive-expression'
echo 'export expression_data_dir="${expression_dir}/data"'
echo 'export output_dir_top=${expression_dir}/D-Apul/output/11-Apul-sRNA-ShortStack-pulchra_genome'
echo ""

echo "# Input/Output files"
echo 'export genome_fasta_dir=${expression_dir}/D-Apul/data'
echo 'export genome_fasta_name="Apulchra-genome.fa"'
echo 'export shortstack_genome_fasta_name="Apulchra-genome.fa"'
echo 'export trimmed_fastqs_dir="${genome_fasta_dir}/sRNA-trimmed-reads"'

echo 'export mirbase_mature_fasta_version=cnidarian-mirbase-mature-v22.1.fasta'
echo 'export genome_fasta="${genome_fasta_dir}/${shortstack_genome_fasta_name}"'
echo ""

echo "# Set number of CPUs to use"
echo 'export threads=40'
echo ""

echo "# Initialize arrays"
echo 'export trimmed_fastqs_array=()'


} > .bashvars

cat .bashvars
```

    #### Assign Variables ####

    # Trimmed FastQ naming pattern
    export trimmed_fastqs_pattern='*fastp-adapters-polyG-31bp-merged.fq.gz'
    # Data directories
    export expression_dir=/home/shared/8TB_HDD_02/shedurkin/deep-dive-expression
    export expression_data_dir="${expression_dir}/data"
    export output_dir_top=${expression_dir}/D-Apul/output/11-Apul-sRNA-ShortStack-pulchra_genome

    # Input/Output files
    export genome_fasta_dir=${expression_dir}/D-Apul/data
    export genome_fasta_name="Apulchra-genome.fa"
    export shortstack_genome_fasta_name="Apulchra-genome.fa"
    export trimmed_fastqs_dir="${genome_fasta_dir}/sRNA-trimmed-reads"
    export mirbase_mature_fasta_version=cnidarian-mirbase-mature-v22.1.fasta
    export genome_fasta="${genome_fasta_dir}/${shortstack_genome_fasta_name}"

    # Set number of CPUs to use
    export threads=40

    # Initialize arrays
    export trimmed_fastqs_array=()

# 3 Load [ShortStack](https://github.com/MikeAxtell/ShortStack) conda environment

If this is successful, the first line of output should show that the
Python being used is the one in your
$$ShortStack$$(<https://github.com/MikeAxtell/ShortStack> conda
environment path.

E.g.

`python:         /home/sam/programs/mambaforge/envs/mirmachine_env/bin/python`

``` r
use_condaenv(condaenv = shortstack_conda_env_name, conda = shortstack_cond_path)
py_config()
```

    python:         /home/sam/programs/mambaforge/envs/ShortStack-4.1.0_env/bin/python
    libpython:      /home/sam/programs/mambaforge/envs/ShortStack-4.1.0_env/lib/libpython3.12.so
    pythonhome:     /home/sam/programs/mambaforge/envs/ShortStack-4.1.0_env:/home/sam/programs/mambaforge/envs/ShortStack-4.1.0_env
    version:        3.12.7 | packaged by conda-forge | (main, Oct  4 2024, 16:05:46) [GCC 13.3.0]
    numpy:          /home/sam/programs/mambaforge/envs/ShortStack-4.1.0_env/lib/python3.12/site-packages/numpy
    numpy_version:  2.1.1

    NOTE: Python version was forced by use_python() function

Note: I sometimes get an error “failed to initialize requested version
of Python,” which seems to stem from the `reticulate` package default
loading a python environment. I’ve been able to fix this by manually
uninstalling the `reticulate` package, then restarting R and
reinstalling `reticulate` before rerunning this code document. \#
Download reference files

## 3.1 A.pulchra genome

``` bash
# Load bash variables into memory
source .bashvars

wget -O ${genome_fasta_dir}/${shortstack_genome_fasta_name} "https://osf.io/download/kn96u/"
```

## 3.2 Cnidarian+miRBase database

Available n `deep-dive` repo,
[here](https://github.com/urol-e5/deep-dive/blob/main/data/cnidarian-mirbase-mature-v22.1.fasta)

``` bash
# Load bash variables into memory
source .bashvars

wget -O ${expression_data_dir}/"${mirbase_mature_fasta_version}" "https://raw.githubusercontent.com/urol-e5/deep-dive/refs/heads/main/data/cnidarian-mirbase-mature-v22.1.fasta"
```

``` bash
# Load bash variables into memory
source .bashvars

head -5 ${expression_data_dir}/"${mirbase_mature_fasta_version}"
```

    >cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p
    UGAGGUAGUAGGUUGUAUAGUU
    >cel-let-7-3p MIMAT0015091 Caenorhabditis elegans let-7-3p
    CUAUGCAAUUUUCUACCUUACC
    >cel-lin-4-5p MIMAT0000002 Caenorhabditis elegans lin-4-5p

## 3.3 Trimmed sRNA-seq reads

Trimmed in `deep-dive`,
[08.2-Dapul-sRNAseq-trimming-31bp-fastp-merged.Rmd](https://github.com/urol-e5/deep-dive/blob/main/D-Apul/code/08.2-Dapul-sRNAseq-trimming-31bp-fastp-merged.Rmd)

``` bash
# Load bash variables into memory
source .bashvars

wget \
--directory-prefix ${trimmed_fastqs_dir} \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 7 \
--no-host-directories \
--no-parent \
--quiet \
--accept ${trimmed_fastqs_pattern} \
"https://gannet.fish.washington.edu/Atumefaciens/gitrepos/deep-dive/D-Apul/output/08.2-Apul-sRNAseq-trimming-31bp-fastp-merged/trimmed-reads/"
```

# 4 Run ShortStack

## 4.1 Modify genome filename for ShortStack compatability

``` bash
# Load bash variables into memory
source .bashvars

# Check for FastA file first
# Then create rename file if doesn't exist
if [ -f "${genome_fasta_dir}/${shortstack_genome_fasta_name}" ]; then
  echo "${genome_fasta_dir}/${shortstack_genome_fasta_name}"
  echo ""
  echo "Already exists. Nothing to do."
  echo ""
else

  # Copy genome FastA to ShortStack-compatible filename (ending with .fa)
  cp ${genome_fasta_dir}/${genome_fasta_name} ${genome_fasta_dir}/${shortstack_genome_fasta_name}
fi

# Confirm
ls -lh ${genome_fasta_dir}/${shortstack_genome_fasta_name}
```

    /home/shared/8TB_HDD_02/shedurkin/deep-dive-expression/D-Apul/data/Apulchra-genome.fa

    Already exists. Nothing to do.

    -rw-r--r-- 1 shedurkin labmembers 505M Oct  1 13:31 /home/shared/8TB_HDD_02/shedurkin/deep-dive-expression/D-Apul/data/Apulchra-genome.fa

## 4.2 Excecute ShortStack command

Uses the `--dn_mirna` option to identify miRNAs in the genome, without
relying on the `--known_miRNAs`.

This part of the code redirects the output of `time` to the end of
`shortstack.log` file.

- `; } \ 2>> ${output_dir_top}/shortstack.log`

``` bash
# Load bash variables into memory
source .bashvars

# Make output directory, if it doesn't exist
mkdir --parents "${output_dir_top}"

# Create array of trimmed FastQs
trimmed_fastqs_array=(${trimmed_fastqs_dir}/${trimmed_fastqs_pattern})


# Pass array contents to new variable as space-delimited list
trimmed_fastqs_list=$(echo "${trimmed_fastqs_array[*]}")


###### Run ShortStack ######
{ time \
ShortStack \
--genomefile "${genome_fasta}" \
--readfile ${trimmed_fastqs_list} \
--known_miRNAs ${expression_data_dir}/${mirbase_mature_fasta_version} \
--dn_mirna \
--threads ${threads} \
--outdir ${output_dir_top}/ShortStack_out \
&> ${output_dir_top}/shortstack.log ; } \
2>> ${output_dir_top}/shortstack.log
```

## 4.3 Check runtime

``` bash
# Load bash variables into memory
source .bashvars

tail -n 3 ${output_dir_top}/shortstack.log \
| grep "real" \
| awk '{print "ShortStack runtime:" "\t" $2}'
```

    ShortStack runtime: 13m43.312s

# 5 Results

## 5.1 ShortStack synopsis

``` bash
# Load bash variables into memory
source .bashvars

tail -n 25 ${output_dir_top}/shortstack.log
```

    Writing final files

    Found a total of 39 MIRNA loci


    Non-MIRNA loci by DicerCall:
    N 20608
    22 78
    23 27
    21 23
    24 14

    Creating visualizations of microRNA loci with strucVis
    <<< WARNING >>>
    Do not rely on these results alone to annotate new MIRNA loci!
    The false positive rate for de novo MIRNA identification is low, but NOT ZERO
    Insepct each mirna locus, especially the strucVis output, and see
    https://doi.org/10.1105/tpc.17.00851 , https://doi.org/10.1093/nar/gky1141

    Tue 08 Oct 2024 15:54:04 -0700 PDT
    Run Completed!

    real    13m43.312s
    user    170m20.326s
    sys 55m0.571s

ShortStack identified 39 miRNAs. (This is in comparison to running
ShortStack using the *A. millipora* genome, which yielded 38 miRNAs)

## 5.2 Inspect `Results.txt`

``` bash
# Load bash variables into memory
source .bashvars

head ${output_dir_top}/ShortStack_out/Results.txt

echo ""
echo "----------------------------------------------------------"
echo ""

echo "Nummber of potential loci:"
awk '(NR>1)' ${output_dir_top}/ShortStack_out/Results.txt | wc -l
```

    Locus   Name    Chrom   Start   End Length  Reads   DistinctSequences   FracTop Strand  MajorRNA    MajorRNAReads   Short   Long    21  22  23  24  DicerCall   MIRNA   known_miRNAs
    ntLink_7:2229-2654  Cluster_1   ntLink_7    2229    2654    426 177 78  0.023   -   UUCUAAGGUCCGUAUCUCCAGUCCUCUGU   47  2   174 1   0   0   0   N   N   NA
    ntLink_7:3045-3472  Cluster_2   ntLink_7    3045    3472    428 154 66  0.013   -   UGAACGUAUUUUCUGAAGAAACUGCAA 15  2   151 0   0   0   1   N   N   NA
    ntLink_7:9749-10310 Cluster_3   ntLink_7    9749    10310   562 680 184 0.531   .   UCUCAUCAUAUUUAGAAUUUUUGGAAGGA   52  5   651 1   6   11  6   N   N   NA
    ntLink_7:21735-22163    Cluster_4   ntLink_7    21735   22163   429 196 72  0.015   -   UUCUAAGGUCCGUAUCUCCAGUCCUCUGC   59  0   195 0   0   0   1   N   N   NA
    ntLink_7:22554-22980    Cluster_5   ntLink_7    22554   22980   427 146 59  0.062   -   UGAACGUAUUUUCUGAAGAAACUGCAA 15  1   140 0   3   0   2   N   N   NA
    ntLink_7:29261-30185    Cluster_6   ntLink_7    29261   30185   925 808 201 0.532   .   GUCAAGUGCAUCGAUCAAGGAUGGAUCAGG  86  5   777 2   7   9   8   N   N   NA
    ntLink_7:39950-40377    Cluster_7   ntLink_7    39950   40377   428 253 26  0.032   -   UAGAUAUACUUCAAGGGAUAAACCAGGU    216 2   240 0   2   1   8   N   N   NA
    ntLink_7:41332-41775    Cluster_8   ntLink_7    41332   41775   444 910 242 0.053   -   UGCUGGUCUUGUAAUUCUGAACGUAUAUUG  93  7   887 4   1   3   8   N   N   NA
    ntLink_7:42041-42467    Cluster_9   ntLink_7    42041   42467   427 169 60  0.036   -   UCUGAAGAAACUGCAAAGUUCACUGUC 50  3   164 0   1   0   1   N   N   NA

    ----------------------------------------------------------

    Nummber of potential loci:
    20789

Column 20 of the `Results.txt` file identifies if a cluster is a miRNA
or not (`Y` or `N`).

``` bash
# Load bash variables into memory
source .bashvars

echo "Number of loci characterized as miRNA:"
awk '$20=="Y" {print $0}' ${output_dir_top}/ShortStack_out/Results.txt \
| wc -l
echo ""

echo "----------------------------------------------------------"

echo ""
echo "Number of loci _not_ characterized as miRNA:"
awk '$20=="N" {print $0}' ${output_dir_top}/ShortStack_out/Results.txt \
| wc -l
```

    Number of loci characterized as miRNA:
    39

    ----------------------------------------------------------

    Number of loci _not_ characterized as miRNA:
    20750

Column 21 of the `Results.txt` file identifies if a cluster aligned to a
known miRNA (miRBase) or not (`Y` or `NA`).

Since there are no miRNAs, the following code will *not* print any
output.

The `echo` command after the `awk` command is simply there to prove that
the chunk executed.

``` bash
# Load bash variables into memory
source .bashvars

echo "Number of loci matching miRBase miRNAs:"
awk '$21!="NA" {print $0}' ${output_dir_top}/ShortStack_out/Results.txt \
| wc -l
echo ""

echo "----------------------------------------------------------"

echo ""
echo "Number of loci _not_ matching miRBase miRNAs:"
awk '$21=="NA" {print $0}' ${output_dir_top}/ShortStack_out/Results.txt \
| wc -l
```

    Number of loci matching miRBase miRNAs:
    131

    ----------------------------------------------------------

    Number of loci _not_ matching miRBase miRNAs:
    20659

Although there are loci with matches to miRBase miRNAs, ShortStack did
*not* annotate these clusters as miRNAs likely [because they do not
*also* match secondary structure
criteria](https://github.com/MikeAxtell/ShortStack#mirna-annotation).

### 5.2.1 Directory tree of all ShortStack outputs

Many of these are large (by GitHub standards) BAM files, so will not be
added to the repo.

Additionally, it’s unlikely we’ll utilize most of the other files
(bigwig) generated by ShortStack.

``` bash
# Load bash variables into memory
source .bashvars

tree -h ${output_dir_top}/
```

    /home/shared/8TB_HDD_02/shedurkin/deep-dive-expression/D-Apul/output/11-Apul-sRNA-ShortStack-pulchra_genome/
    ├── [4.0K]  figures
    │   ├── [140K]  Apul_ShortStack_dbmatch_histogram.png
    │   ├── [168K]  Apul_ShortStack_dbmatch_histogram_reduced.png
    │   ├── [270K]  Apul_ShortStack_miRNA_histogram.png
    │   ├── [268K]  Apul_ShortStack_miRNA_histogram_reduced.png
    │   └── [198K]  Apul_ShortStack_venn.png
    ├── [9.0K]  shortstack.log
    └── [ 12K]  ShortStack_out
        ├── [ 30K]  alignment_details.tsv
        ├── [1.2M]  Counts.txt
        ├── [110K]  known_miRNAs.gff3
        ├── [1.8M]  known_miRNAs_unaligned.fasta
        ├── [328M]  merged_alignments.bam
        ├── [213K]  merged_alignments.bam.csi
        ├── [ 11K]  mir.fasta
        ├── [2.0M]  Results.gff3
        ├── [2.8M]  Results.txt
        ├── [ 66M]  sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam
        ├── [218K]  sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam.csi
        ├── [255M]  sRNA-ACR-140-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.fa
        ├── [ 73M]  sRNA-ACR-145-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam
        ├── [221K]  sRNA-ACR-145-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam.csi
        ├── [279M]  sRNA-ACR-145-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.fa
        ├── [ 70M]  sRNA-ACR-150-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam
        ├── [220K]  sRNA-ACR-150-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam.csi
        ├── [274M]  sRNA-ACR-150-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.fa
        ├── [ 61M]  sRNA-ACR-173-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam
        ├── [219K]  sRNA-ACR-173-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam.csi
        ├── [235M]  sRNA-ACR-173-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.fa
        ├── [ 58M]  sRNA-ACR-178-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam
        ├── [220K]  sRNA-ACR-178-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam.csi
        ├── [218M]  sRNA-ACR-178-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.fa
        └── [4.0K]  strucVis
            ├── [9.8K]  Cluster_10051.ps.pdf
            ├── [ 22K]  Cluster_10051.txt
            ├── [9.0K]  Cluster_10057.ps.pdf
            ├── [ 42K]  Cluster_10057.txt
            ├── [9.5K]  Cluster_10093.ps.pdf
            ├── [ 30K]  Cluster_10093.txt
            ├── [ 11K]  Cluster_10207.ps.pdf
            ├── [ 40K]  Cluster_10207.txt
            ├── [8.5K]  Cluster_10228.ps.pdf
            ├── [9.9K]  Cluster_10228.txt
            ├── [8.7K]  Cluster_10419.ps.pdf
            ├── [ 15K]  Cluster_10419.txt
            ├── [8.6K]  Cluster_14402.ps.pdf
            ├── [ 10K]  Cluster_14402.txt
            ├── [7.7K]  Cluster_14768.ps.pdf
            ├── [1.2K]  Cluster_14768.txt
            ├── [8.8K]  Cluster_15316.ps.pdf
            ├── [ 15K]  Cluster_15316.txt
            ├── [ 12K]  Cluster_15340.ps.pdf
            ├── [ 29K]  Cluster_15340.txt
            ├── [8.5K]  Cluster_15671.ps.pdf
            ├── [9.6K]  Cluster_15671.txt
            ├── [10.0K]  Cluster_15775.ps.pdf
            ├── [ 46K]  Cluster_15775.txt
            ├── [8.4K]  Cluster_15851.ps.pdf
            ├── [6.2K]  Cluster_15851.txt
            ├── [8.4K]  Cluster_15854.ps.pdf
            ├── [3.4K]  Cluster_15854.txt
            ├── [7.6K]  Cluster_15928.ps.pdf
            ├── [9.1K]  Cluster_15928.txt
            ├── [8.7K]  Cluster_16409.ps.pdf
            ├── [7.3K]  Cluster_16409.txt
            ├── [7.9K]  Cluster_17776.ps.pdf
            ├── [3.6K]  Cluster_17776.txt
            ├── [ 10K]  Cluster_17791.ps.pdf
            ├── [ 52K]  Cluster_17791.txt
            ├── [8.7K]  Cluster_1826.ps.pdf
            ├── [ 36K]  Cluster_1826.txt
            ├── [8.7K]  Cluster_1832.ps.pdf
            ├── [ 22K]  Cluster_1832.txt
            ├── [ 10K]  Cluster_1862.ps.pdf
            ├── [ 24K]  Cluster_1862.txt
            ├── [ 10K]  Cluster_18711.ps.pdf
            ├── [ 42K]  Cluster_18711.txt
            ├── [9.8K]  Cluster_18723.ps.pdf
            ├── [ 11K]  Cluster_18723.txt
            ├── [9.6K]  Cluster_18728.ps.pdf
            ├── [ 63K]  Cluster_18728.txt
            ├── [8.7K]  Cluster_18772.ps.pdf
            ├── [ 14K]  Cluster_18772.txt
            ├── [8.2K]  Cluster_19193.ps.pdf
            ├── [7.3K]  Cluster_19193.txt
            ├── [ 10K]  Cluster_1951.ps.pdf
            ├── [ 20K]  Cluster_1951.txt
            ├── [ 11K]  Cluster_2463.ps.pdf
            ├── [ 33K]  Cluster_2463.txt
            ├── [8.3K]  Cluster_2859.ps.pdf
            ├── [ 21K]  Cluster_2859.txt
            ├── [8.1K]  Cluster_3250.ps.pdf
            ├── [5.4K]  Cluster_3250.txt
            ├── [ 10K]  Cluster_3366.ps.pdf
            ├── [ 20K]  Cluster_3366.txt
            ├── [9.0K]  Cluster_3367.ps.pdf
            ├── [ 24K]  Cluster_3367.txt
            ├── [8.5K]  Cluster_3437.ps.pdf
            ├── [ 24K]  Cluster_3437.txt
            ├── [7.9K]  Cluster_4220.ps.pdf
            ├── [4.9K]  Cluster_4220.txt
            ├── [7.9K]  Cluster_4254.ps.pdf
            ├── [4.5K]  Cluster_4254.txt
            ├── [ 10K]  Cluster_5012.ps.pdf
            ├── [ 59K]  Cluster_5012.txt
            ├── [8.6K]  Cluster_5899.ps.pdf
            ├── [9.3K]  Cluster_5899.txt
            ├── [ 11K]  Cluster_5900.ps.pdf
            ├── [ 44K]  Cluster_5900.txt
            ├── [8.7K]  Cluster_5981.ps.pdf
            └── [ 23K]  Cluster_5981.txt

    3 directories, 108 files

## 5.3 Visualize

We noticed that a) not all of the identified miRNAs have database
matches, and b) some reads have a match in the database but are *not*
classified as miRNAs. Let’s look at this in more depth.

``` r
Apul_shortstack_results <- read.csv("../output/11-Apul-sRNA-ShortStack-pulchra_genome/ShortStack_out/Results.txt", sep="\t")
```

``` r
# Reads identified as miRNAs (but not necessarily known)
Apul_shortstack_results %>% 
  filter(MIRNA == "Y") %>%
  mutate(known_miRNAs = str_sub(known_miRNAs, 1, 40)) %>%
  mutate(Locus = str_sub(Locus, 20, 40)) %>%
  ggplot(aes(x = reorder(Locus, Reads), y = Reads, fill = known_miRNAs)) +
  geom_bar(stat = "identity") +
   geom_text(aes(label = Reads), vjust = 0.5, hjust = 0, color = "black", size = 2.5, angle = 90) +
  labs(x = "miRNA", y = "Read count", 
       title = "Reads identified by ShortStack as miRNAs",
       fill = "Annotation") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
```

![](11-Apul-sRNA-ShortStack-pulchra_genome_files/figure-gfm/generate-plots-1.png)<!-- -->

``` r
ggsave("../output/11-Apul-sRNA-ShortStack-pulchra_genome/figures/Apul_ShortStack_miRNA_histogram.png", width = 12, height = 7, units = "in")


# Reads matched in the reference db (but not necessarily identified as miRNA)
Apul_shortstack_results %>% 
  filter(!is.na(known_miRNAs)) %>%
  mutate(known_miRNAs = str_sub(known_miRNAs, 1, 40)) %>%
  mutate(Locus = str_sub(Locus, 20, 40)) %>%
  ggplot(aes(x = reorder(Locus, Reads), y = Reads, fill = MIRNA)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Reads), vjust = 0.5, hjust = 0, color = "black", size = 2.5, angle = 90) +
  labs(x = "miRNA", y = "Read count", 
       title = "Reads with miRBase+cnidarian database matches",
       fill = "Identified as miRNA?") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
```

![](11-Apul-sRNA-ShortStack-pulchra_genome_files/figure-gfm/generate-plots-2.png)<!-- -->

``` r
ggsave("../output/11-Apul-sRNA-ShortStack-pulchra_genome/figures/Apul_ShortStack_dbmatch_histogram.png", width = 12, height = 7, units = "in")
```

There’s two miRNAs with very high read counts, and it’s making
visualization of the rest difficult. Let’s remove them and retry
visualizing the rest.

``` r
# Reads identified as miRNAs (but not necessarily known)
Apul_shortstack_results %>% 
  filter(MIRNA == "Y") %>%
  filter(Reads < 200000) %>%
  mutate(known_miRNAs = str_sub(known_miRNAs, 1, 40)) %>%
  mutate(Locus = str_sub(Locus, 20, 40)) %>%
  ggplot(aes(x = reorder(Locus, Reads), y = Reads, fill = known_miRNAs)) +
  geom_bar(stat = "identity") +
   geom_text(aes(label = Reads), vjust = 0.5, hjust = 0, color = "black", size = 2.5, angle = 90) +
  labs(x = "miRNA", y = "Read count", 
       title = "Reads identified by ShortStack as miRNAs",
       fill = "Annotation") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
```

![](11-Apul-sRNA-ShortStack-pulchra_genome_files/figure-gfm/generate-plots-reduced-1.png)<!-- -->

``` r
ggsave("../output/11-Apul-sRNA-ShortStack-pulchra_genome/figures/Apul_ShortStack_miRNA_histogram_reduced.png", width = 12, height = 7, units = "in")


# Reads matched in the reference db (but not necessarily identified as miRNA)
Apul_shortstack_results %>% 
  filter(!is.na(known_miRNAs)) %>%
  filter(Reads < 200000) %>%
  mutate(known_miRNAs = str_sub(known_miRNAs, 1, 40)) %>%
  mutate(Locus = str_sub(Locus, 20, 40)) %>%
  ggplot(aes(x = reorder(Locus, Reads), y = Reads, fill = MIRNA)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Reads), vjust = 0.5, hjust = 0, color = "black", size = 2.5, angle = 90) +
  labs(x = "miRNA", y = "Read count", 
       title = "Reads with miRBase+cnidarian database matches",
       fill = "Identified as miRNA?") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
```

![](11-Apul-sRNA-ShortStack-pulchra_genome_files/figure-gfm/generate-plots-reduced-2.png)<!-- -->

``` r
ggsave("../output/11-Apul-sRNA-ShortStack-pulchra_genome/figures/Apul_ShortStack_dbmatch_histogram_reduced.png", width = 12, height = 7, units = "in")
```

``` r
# Make list
mirnas <- Apul_shortstack_results %>% filter(MIRNA == "Y") %>% pull(Locus)
matches <- Apul_shortstack_results %>% filter(!is.na(known_miRNAs)) %>% pull(Locus)

Apul_shortstack_vennlist <- list(
  "Identified as miRNA" = mirnas,
  "Database match" = matches
)

# Make venn diagrams
ggvenn(Apul_shortstack_vennlist)
```

![](11-Apul-sRNA-ShortStack-pulchra_genome_files/figure-gfm/venn-diagram-1.png)<!-- -->

``` r
ggsave("../output/11-Apul-sRNA-ShortStack-pulchra_genome/figures/Apul_ShortStack_venn.png", width = 12, height = 7, units = "in")
```

------------------------------------------------------------------------

# 6 Citations

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-axtell2013a" class="csl-entry">

Axtell, Michael J. 2013. “ShortStack: Comprehensive Annotation and
Quantification of Small RNA Genes.” *RNA* 19 (6): 740–51.
<https://doi.org/10.1261/rna.035279.112>.

</div>

<div id="ref-johnson2016a" class="csl-entry">

Johnson, Nathan R, Jonathan M Yeoh, Ceyda Coruh, and Michael J Axtell.
2016. “Improved Placement of Multi-Mapping Small RNAs.” *G3
Genes\|Genomes\|Genetics* 6 (7): 2103–11.
<https://doi.org/10.1534/g3.116.030452>.

</div>

<div id="ref-shahid2014" class="csl-entry">

Shahid, Saima, and Michael J. Axtell. 2014. “Identification and
Annotation of Small RNA Genes Using ShortStack.” *Methods* 67 (1):
20–27. <https://doi.org/10.1016/j.ymeth.2013.10.004>.

</div>

</div>
