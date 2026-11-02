#!/usr/bin/env python3
import pysam
import argparse
import gzip

def get_args():
    parser = argparse.ArgumentParser(description="Extract tags from BAM to TSV")
    parser.add_argument('--bam_in', required=True)
    parser.add_argument('--tsv_out', required=True)
    parser.add_argument('--tags', nargs='+', required=True)
    parser.add_argument('--include_multimappers', action='store_true')
    return parser.parse_args()

def main():
    args = get_args()
    bam_in = pysam.AlignmentFile(args.bam_in, "rb")
    tags = args.tags
    include_mm = args.include_multimappers
    
    # Set up columns: Read ID, Chrom, Strand, Read Length, [Tags...]
    header = ["read_id", "chrom", "strand", "read_length"] + tags
    
    with gzip.open(args.tsv_out, 'wt') as f:
        f.write("\t".join(header) + "\n")
        for read in bam_in:
            if read.is_unmapped:
                continue
                
            if not include_mm:
                try:
                    if read.get_tag("NH") > 1:
                        continue
                except KeyError:
                    pass
                    
            chrom = read.reference_name
            strand = '-' if read.is_reverse else '+'
            
            # Extract read length natively from the PySAM object 
            read_length = read.query_length
            
            row = [read.query_name, chrom, strand, str(read_length)]
            
            for t in tags:
                try:
                    val = read.get_tag(t)
                    row.append(str(val))
                except KeyError:
                    row.append("NA")
            f.write("\t".join(row) + "\n")

    bam_in.close()

if __name__ == "__main__":
    main()