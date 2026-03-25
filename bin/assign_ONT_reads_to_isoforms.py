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
    algn_counter = 0 

    for bam_path in bam_paths:
        bam_file = HTSeq.SAM_Reader(bam_path)
        for algn in bam_file:
            if not algn.aligned or algn.supplementary:
                continue
            
            try:
                pt_val = algn.optional_field("pt")
            except (KeyError, AttributeError):
                pt_val = None
                
            try:
                qs_val = algn.optional_field("qs")
            except (KeyError, AttributeError):
                qs_val = None
            
            introns_list = []
            cigar_string = ""
            
            if algn.cigar is not None:
                for cigar_elem in algn.cigar:
                    cigar_string += f"{cigar_elem.size}{cigar_elem.type}"
                    if cigar_elem.type == 'N': 
                        reference_region = cigar_elem.ref_iv
                        introns_list.append((reference_region.start, reference_region.end + 1))
            
            introns_str = str(tuple(introns_list))
            
            read_coords.append({
                'alignment_id': algn_counter, 
                'read_id': algn.read.name,
                'chrom': algn.iv.chrom,
                'start': algn.iv.start,
                'end': algn.iv.end,
                'strand': algn.iv.strand,
                'introns': introns_str,
                'pt': pt_val,
                'qs': qs_val
            })
            algn_counter += 1
    
    read_df = pd.DataFrame(read_coords)

    if read_df.empty:
        print("[WARNING] No valid alignments found in the provided BAM file(s).")
        sys.exit(0)

    # 2. Process based on Mode
    if args.mode == 'enrich':
        isoforms_df = read_df.groupby(['chrom', 'strand', 'introns']).agg({'start': 'min', 'end': 'max'}).reset_index()
        isoforms_df['transcript_id'] = [f"isoform_{i}" for i in range(len(isoforms_df))]
        
        out_path = f"{args.output_prefix}_enriched.tsv"
        isoforms_df.to_csv(out_path, sep='\t', index=False)
        print(f"[SUCCESS] Enriched transcriptome saved to {out_path} with {len(isoforms_df)} unique isoforms.")

    elif args.mode == 'assign':
        master_map = pd.read_csv(args.input_gtf_file, sep='\t')
        assignments = pd.merge(read_df, master_map, on=['chrom', 'strand', 'introns'], how='left', suffixes=('', '_ref'))
        
        final_cols = ['alignment_id', 'read_id', 'transcript_id', 'chrom', 'start', 'end', 'strand', 'pt', 'qs']
        
        # Added .copy() to avoid SettingWithCopy warnings
        unassigned = assignments[assignments['transcript_id'].isna()].copy()
        unassigned_path = f"{args.output_prefix}_unassigned.tsv"
        unassigned[final_cols].to_csv(unassigned_path, sep='\t', index=False)
            
        assigned = assignments.dropna(subset=['transcript_id']).copy()
        out_path = f"{args.output_prefix}_assignments.tsv"
        assigned[final_cols].to_csv(out_path, sep='\t', index=False)
            
        print(f"[SUCCESS] Assigned: {len(assigned)} | Unassigned: {len(unassigned)}. Saved to {out_path}")

if __name__ == "__main__":
    main()