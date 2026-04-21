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
    parser = argparse.ArgumentParser(description="Orient unmapped BAM reads to forward strand based on basecaller tag.")
    parser.add_argument("--input_bam", required=True, help="Input unmapped BAM file")
    parser.add_argument("--output_bam", required=True, help="Output oriented BAM file")
    parser.add_argument("--stats_tsv", required=True, help="Output TSV with tag statistics")
    parser.add_argument("--sample_id", required=True, help="Sample ID for stats tracking")
    
    # NEW ARGUMENTS FOR DYNAMIC TAGGING
    parser.add_argument("--tag_basecaller_ts", required=True, help="Tag storing the basecaller strand (e.g., TS)")
    parser.add_argument("--tag_original_ts", required=True, help="Custom tag to store the original basecaller strand (e.g., ZS)")
    args = parser.parse_args()

    reads_with_ts = 0
    reads_without_ts = 0
    reads_reverse_complemented = 0

    # check_sq=False is critical because unmapped BAMs lack reference contig headers
    with pysam.AlignmentFile(args.input_bam, "rb", check_sq=False) as bam_in, \
         pysam.AlignmentFile(args.output_bam, "wb", template=bam_in) as bam_out:
        
        for read in bam_in:
            try:
                # 1. Fetch the original basecaller tag
                ts = read.get_tag(args.tag_basecaller_ts)
                reads_with_ts += 1
                
                # 2. Store it as a custom tag for traceability
                read.set_tag(args.tag_original_ts, ts, value_type="A")
                
                # 3. Safely delete the original tag using pysam's in-place deletion
                # This bypasses the set_tags() bug with 'B' type binary arrays!
                read.set_tag(args.tag_basecaller_ts, None)
                
                # 4. Flip the sequence and quality scores if Dorado marked it as negative
                if ts == "-":
                    # Reverse complement sequence
                    orig_seq = read.query_sequence
                    read.query_sequence = get_rc(orig_seq)
                    
                    # Reverse quality scores
                    orig_qual = read.query_qualities
                    if orig_qual is not None:
                        read.query_qualities = orig_qual[::-1]
                    
                    reads_reverse_complemented += 1
                    
            except KeyError:
                reads_without_ts += 1
                
            bam_out.write(read)

    # Write the chunk statistics to a TSV
    with open(args.stats_tsv, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        # Dynamically set header based on the provided tag name
        writer.writerow(["sample_id", "chunk_filename", f"reads_with_{args.tag_basecaller_ts}", f"reads_without_{args.tag_basecaller_ts}", "reads_reverse_complemented"])
        chunk_name = os.path.basename(args.input_bam)
        writer.writerow([args.sample_id, chunk_name, reads_with_ts, reads_without_ts, reads_reverse_complemented])

if __name__ == "__main__":
    main()