#!/usr/bin/env python3

import pysam
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Extract 1-based cleavage sites from a BAM custom tag and convert to 0-based BED format.")
    parser.add_argument("--bam_in", required=True, help="Input BAM file with uniquely mapped reads")
    parser.add_argument("--bed_out", required=True, help="Output BED file")
    parser.add_argument("--tag", required=True, help="The SAM tag storing the 1-based cleavage site (e.g., XF or XO)")
    parser.add_argument("--MAPQ_min", type=int, required=False, default=255, help="minimal MAPQ for an alignment to be included")
    args = parser.parse_args()

    with pysam.AlignmentFile(args.bam_in, "rb") as bam, open(args.bed_out, "w") as f:
        for r in bam:
            # MAPQ 255 represents unique alignments in your pipeline
            if r.mapping_quality >= args.MAPQ_min and not r.is_unmapped:
                try:
                    cs_1based = r.get_tag(args.tag)
                    # Convert 1-based tag coordinate to 0-based BED interval (length = 1)
                    start = cs_1based - 1
                    end = cs_1based
                    strand = "-" if r.is_reverse else "+"
                    
                    # Write in BED6 format: chrom, start, end, name, score, strand
                    f.write(f"{r.reference_name}\t{start}\t{end}\t{r.query_name}\t0\t{strand}\n")
                except KeyError:
                    # Skip reads that are missing the designated cleavage site tag
                    pass

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(1)