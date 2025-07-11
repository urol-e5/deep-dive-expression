05-Ptuh-sRNA-ShortStack_4.1.0
================
Kathleen Durkin
2024-10-22

- <a href="#1-set-r-variables" id="toc-1-set-r-variables">1 Set R
  variables</a>
- <a href="#2-create-a-bash-variables-file"
  id="toc-2-create-a-bash-variables-file">2 Create a Bash variables
  file</a>
- <a href="#3-load-shortstack-conda-environment"
  id="toc-3-load-shortstack-conda-environment">3 Load ShortStack conda
  environment</a>
  - <a href="#31-pmeandrina-genome" id="toc-31-pmeandrina-genome">3.1
    P.meandrina genome</a>
  - <a href="#32-cnidarianmirbase-database"
    id="toc-32-cnidarianmirbase-database">3.2 Cnidarian+miRBase database</a>
  - <a href="#33-trimmed-srna-seq-reads"
    id="toc-33-trimmed-srna-seq-reads">3.3 Trimmed sRNA-seq reads</a>
- <a href="#4-run-shortstack" id="toc-4-run-shortstack">4 Run
  ShortStack</a>
  - <a href="#41-modify-genome-filename-for-shortstack-compatability"
    id="toc-41-modify-genome-filename-for-shortstack-compatability">4.1
    Modify genome filename for ShortStack compatability</a>
  - <a href="#42-excecute-shortstack-command"
    id="toc-42-excecute-shortstack-command">4.2 Excecute ShortStack
    command</a>
  - <a href="#43-check-runtime" id="toc-43-check-runtime">4.3 Check
    runtime</a>
- <a href="#5-results" id="toc-5-results">5 Results</a>
  - <a href="#51-shortstack-synopsis" id="toc-51-shortstack-synopsis">5.1
    ShortStack synopsis</a>
  - <a href="#52-inspect-resultstxt" id="toc-52-inspect-resultstxt">5.2
    Inspect <code>Results.txt</code></a>
    - <a href="#521-directory-tree-of-all-shortstack-outputs"
      id="toc-521-directory-tree-of-all-shortstack-outputs">5.2.1 Directory
      tree of all ShortStack outputs</a>
  - <a href="#53-visualize" id="toc-53-visualize">5.3 Visualize</a>
- <a href="#6-citations" id="toc-6-citations">6 Citations</a>

Use [ShortStack](https://github.com/MikeAxtell/ShortStack) ([Axtell
2013](#ref-axtell2013a); [Shahid and Axtell 2014](#ref-shahid2014);
[Johnson et al. 2016](#ref-johnson2016a))to perform alignment of sRNAseq
data and annotation of sRNA-producing genes.

This is the same ShortStack analysis as performed in the `deeep-dive`
project. However, ShortStack has an updated version, [ShortStack
4.1.0](https://github.com/MikeAxtell/ShortStack?tab=readme-ov-file#shortstack-version-4-major-changes),
which provides much faster analysis times *and* additional functionality
for visualizing miRNA hairpin structures and generating
genome-browser-ready quantitative coverage tracks of aligned small RNAs.

As in `deep-dive`, we will also include a customized miRBase database,
utilizing cnidarian miRNAs curated by Jill Ashley, which includes
published cnidarian miRNAs:

- [`cnidarian-mirbase-mature-v22.1.fasta`](../../data/cnidarian-mirbase-mature-v22.1.fasta)

------------------------------------------------------------------------

Inputs:

- Requires trimmed sRNAseq files generated by
  [08.2-Pmea-sRNAseq-trimming-31bp-fastp-merged.Rmd](https://github.com/urol-e5/deep-dive/blob/main/F-Pmea/code/08.2-Pmea-sRNAseq-trimming-31bp-fastp-merged.Rmd)

  - Filenames formatted: `*fastp-R1-31bp-auto_adapters-polyG.fq.gz`

- Genome FastA. Stored on `deep-dive` wiki

- Modified MiRBase v22.1 FastA. Includes cnidarian miRNAs provided by
  Jill Ashley.

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
echo 'export output_dir_top=${expression_dir}/F-Ptuh/output/05-Ptuh-sRNA-ShortStack_4.1.0'
echo ""

echo "# Input/Output files"
echo 'export genome_fasta_dir=${expression_dir}/F-Ptuh/data'
echo 'export genome_fasta_name="Pocillopora_meandrina_HIv1.assembly.fa"'
echo 'export shortstack_genome_fasta_name="Pocillopora_meandrina_HIv1.assembly.fa"'
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

Note: I sometimes get an error “failed to initialize requested version
of Python,” which seems to stem from the `reticulate` package default
loading a python environment. I’ve been able to fix this by manually
uninstalling the `reticulate` package, then restarting R and
reinstalling `reticulate` before rerunning this code document. \#
Download reference files

## 3.1 P.meandrina genome

``` bash
# Load bash variables into memory
source .bashvars

wget -O ${genome_fasta_dir}/${shortstack_genome_fasta_name} "https://owl.fish.washington.edu/halfshell/genomic-databank/Pocillopora_meandrina_HIv1.assembly.fasta"
```

``` bash
# Load bash vairables into memory
source .bashvars

head ${genome_fasta_dir}/${shortstack_genome_fasta_name} | cut -c 1-100
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

## 3.3 Trimmed sRNA-seq reads

Trimmed in `deep-dive`,
[06.2-Ptuh-sRNAseq-trimming-31bp-fastp-merged](https://github.com/urol-e5/deep-dive/blob/main/F-Ptuh/code/06.2-Ptuh-sRNAseq-trimming-31bp-fastp-merged.Rmd)

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
"https://gannet.fish.washington.edu/Atumefaciens/gitrepos/deep-dive/F-Pmea/output/08.2-Pmea-sRNAseq-trimming-31bp-fastp-merged/trimmed-reads/"
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

# 5 Results

## 5.1 ShortStack synopsis

``` bash
# Load bash variables into memory
source .bashvars

tail -n 25 ${output_dir_top}/shortstack.log
```

    Writing final files

    Found a total of 37 MIRNA loci


    Non-MIRNA loci by DicerCall:
    N 7053
    23 35
    22 34
    24 15
    21 15

    Creating visualizations of microRNA loci with strucVis
    <<< WARNING >>>
    Do not rely on these results alone to annotate new MIRNA loci!
    The false positive rate for de novo MIRNA identification is low, but NOT ZERO
    Insepct each mirna locus, especially the strucVis output, and see
    https://doi.org/10.1105/tpc.17.00851 , https://doi.org/10.1093/nar/gky1141

    Tue 22 Oct 2024 17:28:40 -0700 PDT
    Run Completed!

    real    12m18.158s
    user    200m22.234s
    sys 25m26.799s

ShortStack identified 37 miRNAs.

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
    Pocillopora_meandrina_HIv1___Sc0000000:9092-9521    Cluster_1   Pocillopora_meandrina_HIv1___Sc0000000  9092    9521    430 10830   384 1.0 +   GGGGGUAUAGCUCAGUGGUAGA  3856    1437    3765    634 4389    177 428 N   N   NA
    Pocillopora_meandrina_HIv1___Sc0000000:53578-53997  Cluster_2   Pocillopora_meandrina_HIv1___Sc0000000  53578   53997   420 285 13  0.996   +   GCCUAAGUUGCUUGGAACA 136 283 2   0   0   0   0   N   N   NA
    Pocillopora_meandrina_HIv1___Sc0000000:150243-150718    Cluster_3   Pocillopora_meandrina_HIv1___Sc0000000  150243  150718  476 2548    246 0.0 -   UGGCUAUGAUGAAAAUGACU    335 848 380 634 376 139 171 N   N   NA
    Pocillopora_meandrina_HIv1___Sc0000000:173728-174150    Cluster_4   Pocillopora_meandrina_HIv1___Sc0000000  173728  174150  423 1257    65  0.997   +   UUUGAUUGCUGUGAUCUGGUUG  432 106 2   39  636 444 30  22  N   apa-mir-2050_Exaiptasia_pallida_Baumgarten_et_al._2017_miR-2050;_Nve;_Spis;_Adi
    Pocillopora_meandrina_HIv1___Sc0000000:187562-188076    Cluster_5   Pocillopora_meandrina_HIv1___Sc0000000  187562  188076  515 185 38  0.427   .   AUAAAUGUCACUACAAGAAACCUGAAAUCGU 25  1   176 1   1   2   4   N   N   NA
    Pocillopora_meandrina_HIv1___Sc0000000:485730-486254    Cluster_6   Pocillopora_meandrina_HIv1___Sc0000000  485730  486254  525 286 127 1.0 +   GAUGGGUGUUAUUACUCCUCAGACAGAC    48  66  183 7   11  3   16  N   N   NA
    Pocillopora_meandrina_HIv1___Sc0000000:496020-496432    Cluster_7   Pocillopora_meandrina_HIv1___Sc0000000  496020  496432  413 72  24  1.0 +   AUGUAGUCGAGCAAAGUCCAUGUGGACGA   27  0   66  2   1   1   2   N   N   NA
    Pocillopora_meandrina_HIv1___Sc0000000:525310-527341    Cluster_8   Pocillopora_meandrina_HIv1___Sc0000000  525310  527341  2032    14822   2990    0.178   -   UUUUCGUCACUUUCUUCAGCCUCAGAGU    966 143 13706   53  107 313 500 N   N   NA
    Pocillopora_meandrina_HIv1___Sc0000000:541262-541723    Cluster_9   Pocillopora_meandrina_HIv1___Sc0000000  541262  541723  462 724 148 0.086   -   UUGGACGAAAUUUCGAGGUUCACACUCGUU  87  0   717 0   4   2   1   N   N   NA

    ----------------------------------------------------------

    Nummber of potential loci:
    7189

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
    37

    ----------------------------------------------------------

    Number of loci _not_ characterized as miRNA:
    7152

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
    117

    ----------------------------------------------------------

    Number of loci _not_ matching miRBase miRNAs:
    7073

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

    /home/shared/8TB_HDD_02/shedurkin/deep-dive-expression/F-Ptuh/output/05-Ptuh-sRNA-ShortStack_4.1.0/
    ├── [4.0K]  figures
    │   ├── [130K]  Ptuh_ShortStack_dbmatch_histogram.png
    │   ├── [138K]  Ptuh_ShortStack_dbmatch_histogram_reduced.png
    │   ├── [203K]  Ptuh_ShortStack_miRNA_histogram.png
    │   ├── [190K]  Ptuh_ShortStack_miRNA_histogram_reduced.png
    │   └── [199K]  Ptuh_ShortStack_venn.png
    ├── [3.7K]  Ptuh_ShortStack_4.1.0_mature.fasta
    ├── [4.7K]  Ptuh_ShortStack_4.1.0_mature.gff3
    ├── [6.1K]  Ptuh_ShortStack_4.1.0_precursor.fasta
    ├── [4.6K]  Ptuh_ShortStack_4.1.0_precursor.gff3
    ├── [3.7K]  Ptuh_ShortStack_4.1.0_star.fasta
    ├── [4.5K]  Ptuh_ShortStack_4.1.0_star.gff3
    ├── [9.1K]  shortstack.log
    └── [ 20K]  ShortStack_out
        ├── [ 29K]  alignment_details.tsv
        ├── [603K]  Counts.txt
        ├── [179K]  known_miRNAs.gff3
        ├── [1.8M]  known_miRNAs_unaligned.fasta
        ├── [184M]  merged_alignments.bam
        ├── [142K]  merged_alignments.bam.csi
        ├── [ 13K]  mir.fasta
        ├── [ 13K]  Ptuh_Results_mature_named_miRNAs.csv
        ├── [915K]  Results.gff3
        ├── [1.4M]  Results.txt
        ├── [ 33M]  sRNA-POC-47-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam
        ├── [162K]  sRNA-POC-47-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam.csi
        ├── [126M]  sRNA-POC-47-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.fa
        ├── [ 30M]  sRNA-POC-48-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam
        ├── [163K]  sRNA-POC-48-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam.csi
        ├── [113M]  sRNA-POC-48-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.fa
        ├── [ 31M]  sRNA-POC-50-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam
        ├── [162K]  sRNA-POC-50-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam.csi
        ├── [116M]  sRNA-POC-50-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.fa
        ├── [ 42M]  sRNA-POC-53-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam
        ├── [162K]  sRNA-POC-53-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam.csi
        ├── [159M]  sRNA-POC-53-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.fa
        ├── [ 45M]  sRNA-POC-57-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam
        ├── [155K]  sRNA-POC-57-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.bam.csi
        ├── [169M]  sRNA-POC-57-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed.fa
        └── [4.0K]  strucVis
            ├── [ 12K]  Cluster_1015.ps.pdf
            ├── [ 35K]  Cluster_1015.txt
            ├── [8.3K]  Cluster_1068.ps.pdf
            ├── [3.5K]  Cluster_1068.txt
            ├── [8.6K]  Cluster_1080.ps.pdf
            ├── [ 10K]  Cluster_1080.txt
            ├── [8.3K]  Cluster_1116.ps.pdf
            ├── [7.3K]  Cluster_1116.txt
            ├── [9.4K]  Cluster_1289.ps.pdf
            ├── [ 25K]  Cluster_1289.txt
            ├── [ 11K]  Cluster_1296.ps.pdf
            ├── [ 64K]  Cluster_1296.txt
            ├── [ 11K]  Cluster_1793.ps.pdf
            ├── [ 41K]  Cluster_1793.txt
            ├── [9.6K]  Cluster_1938.ps.pdf
            ├── [ 35K]  Cluster_1938.txt
            ├── [ 10K]  Cluster_1952.ps.pdf
            ├── [ 39K]  Cluster_1952.txt
            ├── [9.7K]  Cluster_1953.ps.pdf
            ├── [ 30K]  Cluster_1953.txt
            ├── [ 10K]  Cluster_21.ps.pdf
            ├── [ 30K]  Cluster_21.txt
            ├── [9.9K]  Cluster_2793.ps.pdf
            ├── [8.2K]  Cluster_2793.txt
            ├── [9.4K]  Cluster_2837.ps.pdf
            ├── [5.1K]  Cluster_2837.txt
            ├── [8.7K]  Cluster_2839.ps.pdf
            ├── [7.4K]  Cluster_2839.txt
            ├── [9.0K]  Cluster_2859.ps.pdf
            ├── [ 14K]  Cluster_2859.txt
            ├── [8.4K]  Cluster_2973.ps.pdf
            ├── [5.6K]  Cluster_2973.txt
            ├── [7.6K]  Cluster_3392.ps.pdf
            ├── [ 30K]  Cluster_3392.txt
            ├── [ 10K]  Cluster_360.ps.pdf
            ├── [ 40K]  Cluster_360.txt
            ├── [7.8K]  Cluster_3661.ps.pdf
            ├── [ 14K]  Cluster_3661.txt
            ├── [7.1K]  Cluster_36.ps.pdf
            ├── [3.7K]  Cluster_36.txt
            ├── [8.9K]  Cluster_390.ps.pdf
            ├── [ 13K]  Cluster_390.txt
            ├── [9.6K]  Cluster_4039.ps.pdf
            ├── [7.7K]  Cluster_4039.txt
            ├── [9.3K]  Cluster_4040.ps.pdf
            ├── [5.1K]  Cluster_4040.txt
            ├── [9.4K]  Cluster_4118.ps.pdf
            ├── [ 14K]  Cluster_4118.txt
            ├── [9.5K]  Cluster_4435.ps.pdf
            ├── [ 39K]  Cluster_4435.txt
            ├── [ 11K]  Cluster_4437.ps.pdf
            ├── [ 53K]  Cluster_4437.txt
            ├── [9.8K]  Cluster_4438.ps.pdf
            ├── [ 20K]  Cluster_4438.txt
            ├── [ 10K]  Cluster_4439.ps.pdf
            ├── [ 30K]  Cluster_4439.txt
            ├── [10.0K]  Cluster_4572.ps.pdf
            ├── [ 64K]  Cluster_4572.txt
            ├── [8.7K]  Cluster_4754.ps.pdf
            ├── [6.9K]  Cluster_4754.txt
            ├── [9.7K]  Cluster_4823.ps.pdf
            ├── [ 31K]  Cluster_4823.txt
            ├── [7.9K]  Cluster_5253.ps.pdf
            ├── [4.4K]  Cluster_5253.txt
            ├── [ 11K]  Cluster_5612.ps.pdf
            ├── [ 16K]  Cluster_5612.txt
            ├── [8.9K]  Cluster_5740.ps.pdf
            ├── [ 17K]  Cluster_5740.txt
            ├── [9.1K]  Cluster_6382.ps.pdf
            ├── [ 21K]  Cluster_6382.txt
            ├── [8.7K]  Cluster_757.ps.pdf
            ├── [ 27K]  Cluster_757.txt
            ├── [8.6K]  Cluster_925.ps.pdf
            └── [ 29K]  Cluster_925.txt

    3 directories, 111 files

## 5.3 Visualize

We noticed that a) not all of the identified miRNAs have database
matches, and b) some reads have a match in the database but are *not*
classified as miRNAs. Let’s look at this in more depth.

``` r
Ptuh_shortstack_results <- read.csv("../output/05-Ptuh-sRNA-ShortStack_4.1.0/ShortStack_out/Results.txt", sep="\t")
```

``` r
# Reads identified as miRNAs (but not necessarily known)
Ptuh_shortstack_results %>% 
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

![](05-Ptuh-sRNA-ShortStack_4.1.0_files/figure-gfm/generate-plots-1.png)<!-- -->

``` r
ggsave("../output/05-Ptuh-sRNA-ShortStack_4.1.0/figures/Ptuh_ShortStack_miRNA_histogram.png", width = 12, height = 7, units = "in")


# Reads matched in the reference db (but not necessarily identified as miRNA)
Ptuh_shortstack_results %>% 
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

![](05-Ptuh-sRNA-ShortStack_4.1.0_files/figure-gfm/generate-plots-2.png)<!-- -->

``` r
ggsave("../output/05-Ptuh-sRNA-ShortStack_4.1.0/figures/Ptuh_ShortStack_dbmatch_histogram.png", width = 12, height = 7, units = "in")
```

There’s two miRNAs with very high read counts, and it’s making
visualization of the rest difficult. Let’s remove them and retry
visualizing the rest.

``` r
# Reads identified as miRNAs (but not necessarily known)
Ptuh_shortstack_results %>% 
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

![](05-Ptuh-sRNA-ShortStack_4.1.0_files/figure-gfm/generate-plots-reduced-1.png)<!-- -->

``` r
ggsave("../output/05-Ptuh-sRNA-ShortStack_4.1.0/figures/Ptuh_ShortStack_miRNA_histogram_reduced.png", width = 12, height = 7, units = "in")


# Reads matched in the reference db (but not necessarily identified as miRNA)
Ptuh_shortstack_results %>% 
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

![](05-Ptuh-sRNA-ShortStack_4.1.0_files/figure-gfm/generate-plots-reduced-2.png)<!-- -->

``` r
ggsave("../output/05-Ptuh-sRNA-ShortStack_4.1.0/figures/Ptuh_ShortStack_dbmatch_histogram_reduced.png", width = 12, height = 7, units = "in")
```

``` r
# Make list
mirnas <- Ptuh_shortstack_results %>% filter(MIRNA == "Y") %>% pull(Locus)
matches <- Ptuh_shortstack_results %>% filter(!is.na(known_miRNAs)) %>% pull(Locus)

Ptuh_shortstack_vennlist <- list(
  "Identified as miRNA" = mirnas,
  "Database match" = matches
)

# Make venn diagrams
ggvenn(Ptuh_shortstack_vennlist)
```

![](05-Ptuh-sRNA-ShortStack_4.1.0_files/figure-gfm/venn-diagram-1.png)<!-- -->

``` r
ggsave("../output/05-Ptuh-sRNA-ShortStack_4.1.0/figures/Ptuh_ShortStack_venn.png", width = 12, height = 7, units = "in")
```

Generate separate fastas and gffs for the precursor, star, and mature
miRNAs

``` bash
cd ../output/05-Ptuh-sRNA-ShortStack_4.1.0

awk '/^>/{keep=($0 !~ /mature|star/)} keep' ShortStack_out/mir.fasta > Ptuh_ShortStack_4.1.0_precursor.fasta
awk '/^>/{keep=($0 ~ /star/)} keep' ShortStack_out/mir.fasta > Ptuh_ShortStack_4.1.0_star.fasta
awk '/^>/{keep=($0 ~ /mature/)} keep' ShortStack_out/mir.fasta > Ptuh_ShortStack_4.1.0_mature.fasta

awk -F '\t' '$3 == "MIRNA_hairpin"' ShortStack_out/Results.gff3 > Ptuh_ShortStack_4.1.0_precursor.gff3
awk -F '\t' '$3 == "mature_miRNA"' ShortStack_out/Results.gff3 > Ptuh_ShortStack_4.1.0_mature.gff3
awk -F '\t' '$3 == "miRNA-star"' ShortStack_out/Results.gff3 > Ptuh_ShortStack_4.1.0_star.gff3
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
