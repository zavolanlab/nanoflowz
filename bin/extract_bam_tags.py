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
        f.write("read_id\tchromosome\t" + "\t".join(args.tags) + "\n")
        
        for r in bam:
            if r.is_unmapped:
                continue
                
            if not include_mm and r.mapping_quality != 255:
                continue

            chrom = r.reference_name if r.reference_name else "NA"
            row = [r.query_name, chrom]
            
            for tag in args.tags:
                try:
                    val = str(r.get_tag(tag))
                except KeyError:
                    val = "NA"
                row.append(val)
            f.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()