#!/usr/bin/env bash

sample_id=$1
step_name=$2
bam_file=$3
cpus=$4
out_file=$5

# 1. Total unique reads 
# (cut -f1 extracts the read_id, sort -u deduplicates them)
total_reads=$(samtools view -@ "$cpus" "$bam_file" | cut -f1 | LC_ALL=C sort -u | wc -l)

# 2. Unmapped reads 
# (awk checks if RNAME (col 3) is '*' and POS (col 4) is '0' based on SAM specs)
unmapped_reads=$(samtools view -@ "$cpus" "$bam_file" | awk -F'\t' '$3 == "*" && $4 == 0 {print $1}' | LC_ALL=C sort -u | wc -l)

# 3. Mapped reads (Unique vs Multi-mapped)
# Extract read IDs where RNAME (col 3) is NOT '*', then count how many times each ID appears
samtools view -@ "$cpus" "$bam_file" | awk -F'\t' '$3 != "*" {print $1}' | LC_ALL=C sort | uniq -c > mapped_counts.txt

# If a read ID appears exactly 1 time, it is a unique mapper
unique_reads=$(awk '$1 == 1' mapped_counts.txt | wc -l)

# If a read ID appears > 1 time, it is a multi-mapper (e.g. split/supplementary alignments)
mm_reads=$(awk '$1 > 1' mapped_counts.txt | wc -l)

# Clean up intermediate file
rm mapped_counts.txt

# Output to the exact TSV file requested by Nextflow
echo -e "${sample_id}\t${step_name}\t${total_reads}\t${unmapped_reads}\t${unique_reads}\t${mm_reads}" > "$out_file"