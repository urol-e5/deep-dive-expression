#!/bin/bash

# Check if input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 genome.fasta"
    exit 1
fi

fasta_file="$1"

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
echo "Total length: $total_length"

# Longest sequence
longest=$(sort -nr seq_lengths.txt | head -n1)
echo "Longest sequence length: $longest"

# Shortest sequence
shortest=$(sort -n seq_lengths.txt | head -n1)
echo "Shortest sequence length: $shortest"

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
echo "N50: $n50"

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
        printf "Total bases: %d\n", total
        printf "A: %d\n", a
        printf "C: %d\n", c
        printf "G: %d\n", g
        printf "T: %d\n", t
        printf "N: %d\n", n
        printf "GC Content: %.2f%%\n", (gc / total) * 100
    }
' "$fasta_file"

# Cleanup temporary files
rm seq_lengths.txt seq_lengths_sorted.txt
