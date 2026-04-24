#!/usr/bin/env python3
import pysam
import argparse

def main():
    parser = argparse.ArgumentParser(description="Extract 1-based cleavage sites to BED format with 1/NH weights.")
    parser.add_argument("--bam_in", required=True)
    parser.add_argument("--bed_out", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--include_multimappers", type=str, default="true")
    args = parser.parse_args()

    include_mm = args.include_multimappers.lower() == 'true'

    with pysam.AlignmentFile(args.bam_in, "rb") as bam, open(args.bed_out, "w") as f:
        for r in bam:
            if r.is_unmapped:
                continue
            
            # Filter multi-mappers if requested
            if not include_mm and r.mapping_quality != 255:
                continue

            try:
                cs_1based = r.get_tag(args.tag)
                # Fetch NH tag, defaulting to 1 if it somehow went missing
                nh = r.get_tag("NH") if r.has_tag("NH") else 1
                weight = 1.0 / nh
                
                # Convert 1-based tag coordinate to 0-based BED interval (length = 1)
                start = cs_1based - 1
                end = cs_1based
                strand = "-" if r.is_reverse else "+"
                
                # Write in BED6 format: chrom, start, end, name, score (weight), strand
                f.write(f"{r.reference_name}\t{start}\t{end}\t{r.query_name}\t{weight:.4f}\t{strand}\n")
            except KeyError:
                pass

if __name__ == "__main__":
    main()