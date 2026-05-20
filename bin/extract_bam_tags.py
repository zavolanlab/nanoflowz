#!/usr/bin/env python3
import pysam
import argparse
import gzip

def main():
    parser = argparse.ArgumentParser(description="Extract BAM tags to TSV.")
    parser.add_argument("--bam_in", required=True)
    parser.add_argument("--tsv_out", required=True)
    parser.add_argument("--tags", nargs='+', required=True)
    parser.add_argument("--include_multimappers", type=str, default="true")
    args = parser.parse_args()

    include_mm = args.include_multimappers.lower() == 'true'

    with pysam.AlignmentFile(args.bam_in, "rb") as bam, gzip.open(args.tsv_out, "wt") as f:
        # Set up columns: Read ID, Chrom, Strand, [Tags...]
        header = ["read_id", "chrom", "strand"] + args.tags
    
        f.write("\t".join(header) + "\n")
        for read in bam:
            if read.is_unmapped:
                continue
                
            if not include_mm and read.get_tag("NH") > 1:
                continue
                
            chrom = read.reference_name
            strand = '-' if read.is_reverse else '+'
            row = [read.query_name, chrom, strand]
            
            for t in args.tags:
                try:
                    val = read.get_tag(t)
                    row.append(str(val))
                except KeyError:
                    row.append("NA")
            f.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()