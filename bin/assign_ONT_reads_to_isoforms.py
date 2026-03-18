#!/usr/bin/env python3

import warnings
warnings.simplefilter('ignore')
import sys
import os
import pandas as pd
import numpy as np
import HTSeq
from argparse import ArgumentParser, RawTextHelpFormatter

def main():
    parser = ArgumentParser(description="ONT Isoform Management", formatter_class=RawTextHelpFormatter)
    parser.add_argument("--mode", dest="mode", choices=['enrich', 'assign'], required=True, 
                        help="enrich: Create master isoforms from all BAMs\nassign: Map alignments to the enriched GTF")
    parser.add_argument("--input_bam_files", dest="input_bam_files", required=True)
    parser.add_argument("--input_gtf_file", dest="input_gtf_file", required=True)
    parser.add_argument("--output_prefix", dest="output_prefix", required=True)
    parser.add_argument("--ThreePrimeEnd_clustering_distance", type=int, default=25)
    parser.add_argument("--gtf_skip_rows", type=int, default=5)

    args = parser.parse_args()

    # Determine BAM files to process
    if args.mode == 'enrich' and os.path.isfile(args.input_bam_files):
        with open(args.input_bam_files, 'r') as f:
            bam_paths = [line.strip() for line in f if line.strip()]
    else:
        bam_paths = [args.input_bam_files]

    # 1. Parse Alignments
    read_coords = []
    for bam_path in bam_paths:
        # Using SAM_Reader as it handles both SAM/BAM correctly in HTSeq and matches the original script
        bam_file = HTSeq.SAM_Reader(bam_path)
        for algn in bam_file:
            if not algn.aligned or algn.supplementary:
                continue
            
            introns_list = []
            cigar_string = ""
            
            # Explicitly parse the CIGAR string to find splice junctions ('N')
            if algn.cigar is not None:
                for cigar_elem in algn.cigar:
                    cigar_string += f"{cigar_elem.size}{cigar_elem.type}"
                    if cigar_elem.type == 'N': 
                        reference_region = cigar_elem.ref_iv
                        # 1-based end coordinate matching original script logic
                        introns_list.append((reference_region.start, reference_region.end + 1))
            
            # Convert to string so pandas can group by it later
            introns_str = str(tuple(introns_list))
            
            # Unique ID: read_id + region + CIGAR string
            uid = f"{algn.read.name}_{algn.iv.chrom}_{algn.iv.start}_{algn.iv.end}_{cigar_string}"
            
            read_coords.append({
                'read_id': algn.read.name,
                'alignment_id': uid,
                'chrom': algn.iv.chrom,
                'start': algn.iv.start,
                'end': algn.iv.end,
                'strand': algn.iv.strand,
                'introns': introns_str
            })
    
    read_df = pd.DataFrame(read_coords)

    # Safety check for empty dataframes
    if read_df.empty:
        print("[WARNING] No valid alignments found in the provided BAM file(s).")
        sys.exit(0)

    # 2. Process based on Mode
    if args.mode == 'enrich':
        # Group identical read structures into unique master isoforms
        isoforms_df = read_df.groupby(['chrom', 'strand', 'introns']).agg({'start': 'min', 'end': 'max'}).reset_index()
        isoforms_df['transcript_id'] = [f"isoform_{i}" for i in range(len(isoforms_df))]
        
        # Save Enriched Transcriptome
        out_path = f"{args.output_prefix}_enriched.tsv"
        isoforms_df.to_csv(out_path, sep='\t', index=False)
        print(f"[SUCCESS] Enriched transcriptome saved to {out_path} with {len(isoforms_df)} unique isoforms.")

    elif args.mode == 'assign':
        # Load the Master Map created in 'enrich' mode
        master_map = pd.read_csv(args.input_gtf_file, sep='\t')
        
        # Map reads to Master Map IDs
        assignments = pd.merge(read_df, master_map, on=['chrom', 'strand', 'introns'], how='left', suffixes=('', '_ref'))
        
        # Format output as requested: Associate every alignment with ONE transcript_id
        final_cols = ['alignment_id', 'read_id', 'transcript_id', 'chrom', 'start', 'end', 'strand']
        out_path = f"{args.output_prefix}_assignments.tsv"
        
        # Handle unassigned reads
        unassigned = assignments[assignments['transcript_id'].isna()]
        if not unassigned.empty:
            unassigned[final_cols].to_csv(f"{args.output_prefix}_unassigned.tsv", sep='\t', index=False)
            
        assigned = assignments.dropna(subset=['transcript_id'])
        if not assigned.empty:
            assigned[final_cols].to_csv(out_path, sep='\t', index=False)
        else:
            # Output an empty file with headers if no assignments match
            pd.DataFrame(columns=final_cols).to_csv(out_path, sep='\t', index=False)
            
        print(f"[SUCCESS] Assigned: {len(assigned)} | Unassigned: {len(unassigned)}. Saved to {out_path}")

if __name__ == "__main__":
    main()