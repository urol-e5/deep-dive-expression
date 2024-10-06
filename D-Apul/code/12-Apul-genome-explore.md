12-Apul genome-explore
================
Steven Roberts
06 October, 2024

- <a href="#01-genome-analysis" id="toc-01-genome-analysis">0.1 Genome
  Analysis</a>
- <a href="#02-gff-stats" id="toc-02-gff-stats">0.2 GFF Stats</a>

## 0.1 Genome Analysis

``` bash
# Set the fasta file variable
fasta_file="../data/Apulcra-genome.fa"

# Check if input file exists
if [ ! -f "$fasta_file" ]; then
    echo "File not found: $fasta_file"
    exit 1
fi

# Number of sequences
num_seqs=$(grep -c '^>' "$fasta_file")
echo "Number of sequences: $num_seqs"

# Get sequence lengths
awk '
    /^>/ { if (seqlen) { print seqlen; seqlen=0 } next }
    { seqlen += length($0) }
    END { if (seqlen) print seqlen }
' "$fasta_file" > seq_lengths.txt

# Total length
total_length=$(awk '{sum+=$1} END{print sum}' seq_lengths.txt)
printf "Total length: %'d\n" "$total_length"

# Longest sequence
longest=$(sort -nr seq_lengths.txt | head -n1)
printf "Longest sequence length: %'d\n" "$longest"

# Shortest sequence
shortest=$(sort -n seq_lengths.txt | head -n1)
printf "Shortest sequence length: %'d\n" "$shortest"

# N50 calculation
sort -nr seq_lengths.txt > seq_lengths_sorted.txt

n50=$(awk -v total_length="$total_length" '
BEGIN { half_total = total_length / 2; sum = 0 }
{
    sum += $1
    if (sum >= half_total) {
        print $1
        exit
    }
}
' seq_lengths_sorted.txt)
printf "N50: %'d\n" "$n50"

# L50 calculation
l50=$(awk -v total_length="$total_length" '
BEGIN { half_total = total_length / 2; sum = 0; count = 0 }
{
    sum += $1
    count++
    if (sum >= half_total) {
        print count
        exit
    }
}
' seq_lengths_sorted.txt)
echo "L50: $l50"

# GC content and base counts
awk '
    /^>/ { next }
    {
        seq = toupper($0)
        g += gsub(/G/, "", seq)
        c += gsub(/C/, "", seq)
        a += gsub(/A/, "", seq)
        t += gsub(/T/, "", seq)
        n += gsub(/N/, "", seq)
    }
    END {
        total = a + c + g + t + n
        gc = g + c
        printf "Total bases: %'\''d\n", total
        printf "A: %'\''d\n", a
        printf "C: %'\''d\n", c
        printf "G: %'\''d\n", g
        printf "T: %'\''d\n", t
        printf "N: %'\''d\n", n
        printf "GC Content: %.2f%%\n", (gc / total) * 100
    }
' "$fasta_file"

# Cleanup temporary files
rm seq_lengths.txt seq_lengths_sorted.txt
```

## 0.2 GFF Stats

``` bash
# Define the path to the GFF file
gff_file="../data/Apulcra-genome.gff"


# Check if the file exists
if [ ! -f "$gff_file" ]; then
    echo "File does not exist: $gff_file"
    exit 1
fi

echo "Processing file: $gff_file"

# Total number of entries
total_entries=$(grep -vc '^#' $gff_file)
echo "Total entries: $total_entries"

# Number of unique features
unique_features=$(cut -f3 $gff_file | grep -v '^#' | sort | uniq | wc -l)
echo "Unique features: $unique_features"

# Number of entries per feature type
echo "Entries per feature type:"
cut -f3 $gff_file | grep -v '^#' | sort | uniq -c | sort -nr

# Number of unique sources and list each with counts
echo "Unique sources and their counts:"
unique_sources=$(cut -f2 $gff_file | grep -v '^#' | sort | uniq -c | sort -nr)
echo "$unique_sources"
```