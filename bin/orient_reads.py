#!/usr/bin/env python3

import pysam
import argparse
import csv
import sys
import os

def get_rc(seq):
    """Return the reverse complement of a sequence."""
    trans = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")
    return seq.translate(trans)[::-1]

def main():
    parser = argparse.ArgumentParser(description="Orient unmapped BAM reads to forward strand based on TS tag.")
    parser.add_argument("--input_bam", required=True, help="Input unmapped BAM file")
    parser.add_argument("--output_bam", required=True, help="Output oriented BAM file")
    parser.add_argument("--stats_tsv", required=True, help="Output TSV with TS tag statistics")
    parser.add_argument("--sample_id", required=True, help="Sample ID for stats tracking")
    args = parser.parse_args()

    reads_with_ts = 0
    reads_without_ts = 0
    reads_flipped = 0

    # check_sq=False is critical because unmapped BAMs lack reference contig headers
    with pysam.AlignmentFile(args.input_bam, "rb", check_sq=False) as bam_in, \
         pysam.AlignmentFile(args.output_bam, "wb", template=bam_in) as bam_out:
        
        for read in bam_in:
            try:
                ts = read.get_tag("TS")
                reads_with_ts += 1
                
                if ts == "-":
                    # Reverse complement sequence
                    orig_seq = read.query_sequence
                    read.query_sequence = get_rc(orig_seq)
                    
                    # Reverse quality scores
                    orig_qual = read.query_qualities
                    if orig_qual is not None:
                        read.query_qualities = orig_qual[::-1]
                    
                    # Update TS tag
                    read.set_tag("TS", "+", value_type="A")
                    read.set_tag("rF", 1, value_type="i") # Custom traceability tag
                    reads_flipped += 1
                    
            except KeyError:
                reads_without_ts += 1
                
            bam_out.write(read)

    # Write the chunk statistics to a TSV
    with open(args.stats_tsv, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerow(["sample_id", "chunk_filename", "reads_with_TS", "reads_without_TS", "reads_flipped"])
        chunk_name = os.path.basename(args.input_bam)
        writer.writerow([args.sample_id, chunk_name, reads_with_ts, reads_without_ts, reads_flipped])

if __name__ == "__main__":
    main()