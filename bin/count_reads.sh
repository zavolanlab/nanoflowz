#!/bin/bash

SAMPLE=$1
STEP=$2
BAM=$3
CPUS=$4
OUT=$5

# Extract read name and whether it mapped to a chromosome (RNAME != *)
samtools view -@ $CPUS -F 2048 $BAM | awk '{if($3=="*") print $1"\tUN"; else print $1"\tMA"}' | sort -S 2G | uniq -c > counts.tmp

# Sum up the occurrences
UNMAPPED=$(awk '$3=="UN" {sum++} END {print sum+0}' counts.tmp)
UM=$(awk '$1==1 && $3=="MA" {sum++} END {print sum+0}' counts.tmp)
MM=$(awk '$1>1 && $3=="MA" {sum++} END {print sum+0}' counts.tmp)
TOTAL=$((UNMAPPED + UM + MM))

# Output as TSV row
echo -e "${SAMPLE}\t${STEP}\t${TOTAL}\t${UNMAPPED}\t${UM}\t${MM}" > $OUT

rm counts.tmp