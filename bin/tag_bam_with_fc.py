#!/usr/bin/env python3
import pysam
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Tag BAM file using featureCounts CORE output.")
    parser.add_argument("--bam_in", required=True)
    parser.add_argument("--core_txt", required=True)
    parser.add_argument("--bam_out", required=True)
    args = parser.parse_args()

    # 1. Load featureCounts assignments into memory
    # Format: Read_Name \t Status \t Num_Hits \t Target_Features
    assignments = {}
    statuses = {}

    with open(args.core_txt, "r") as f:
        for line in f:
            if line.startswith("#") or line.startswith("Read_Name"):
                continue
            cols = line.strip().split("\t")
            if len(cols) >= 2:
                qname = cols[0]
                status = cols[1]
                statuses[qname] = status
                
                # If assigned, grab the comma-separated list of genes
                if status == "Assigned" and len(cols) >= 4:
                    assignments[qname] = cols[3]

    # 2. Stream BAM and apply tags
    with pysam.AlignmentFile(args.bam_in, "rb") as bam_in, \
         pysam.AlignmentFile(args.bam_out, "wb", header=bam_in.header) as bam_out:
         
        for r in bam_in:
            qname = r.query_name
            
            # Apply Gene Assignment (XT)
            if qname in assignments:
                r.set_tag("XT", assignments[qname], value_type="Z")
            
            # Apply Assignment Status (XS)
            if qname in statuses:
                r.set_tag("XS", statuses[qname], value_type="Z")
                
            bam_out.write(r)

if __name__ == "__main__":
    main()