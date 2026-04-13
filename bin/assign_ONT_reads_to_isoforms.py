#!/usr/bin/env python3

import warnings
warnings.simplefilter('ignore')
import sys
import os
import pandas as pd
import numpy as np
import HTSeq
import ast
import gzip
import logging
from argparse import ArgumentParser, RawTextHelpFormatter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def get_introns_from_exons(exons):
    """Converts a list of (start, end) exons to an intron string."""
    if len(exons) <= 1:
        return "()"
    exons = sorted(exons)
    introns = []
    for i in range(len(exons) - 1):
        # Intron: end of current exon to start of next exon
        introns.append((exons[i][1], exons[i+1][0]))
    return str(tuple(introns))

def main():
    parser = ArgumentParser(description="ONT Isoform Management", formatter_class=RawTextHelpFormatter)
    parser.add_argument("--mode", dest="mode", choices=['enrich', 'assign'], required=True, 
                        help="enrich: Create master isoforms from all BAMs + Reference\nassign: Map alignments to the enriched GTF")
    parser.add_argument("--input_bam_files", dest="input_bam_files", required=True)
    parser.add_argument("--input_gtf_file", dest="input_gtf_file", required=True)
    parser.add_argument("--output_prefix", dest="output_prefix", required=True)
    parser.add_argument("--ThreePrimeEnd_clustering_distance", type=int, default=25)
    parser.add_argument("--gtf_skip_rows", type=int, default=5)

    args = parser.parse_args()

    # Determine BAM files to process
    if os.path.isfile(args.input_bam_files) and not args.input_bam_files.endswith(".bam"):
        with open(args.input_bam_files, 'r') as f:
            bam_paths = [line.strip() for line in f if line.strip()]
    else:
        bam_paths = [args.input_bam_files]

    # --- MODE: ENRICH ---
    if args.mode == 'enrich':
        # 1. Extract Reference Transcriptome Structure
        logging.info(f"Parsing Reference GTF: {args.input_gtf_file}")
        gtf_cols = ['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        ref_gtf_full = pd.read_csv(args.input_gtf_file, sep='\t', comment='#', header=None, names=gtf_cols, skiprows=args.gtf_skip_rows)
        
        # We need exons to build intron strings for the reference
        ref_exons = ref_gtf_full[ref_gtf_full['feature'] == 'exon'].copy()
        ref_exons['transcript_id'] = ref_exons['attribute'].str.extract('transcript_id "([^"]+)"')
        
        ref_data = []
        for tid, group in ref_exons.groupby('transcript_id'):
            chrom = group['chrom'].iloc[0]
            strand = group['strand'].iloc[0]
            # Convert 1-based GTF to 0-based for internal logic consistency
            exons = sorted(list(zip(group['start'] - 1, group['end'])))
            introns_str = get_introns_from_exons(exons)
            
            t_start = group['start'].min() - 1
            t_end = group['end'].max()
            three_prime = t_end if strand == "+" else t_start
            
            ref_data.append({
                'transcript_id': tid, 'chrom': chrom, 'strand': strand,
                'start': t_start, 'end': t_end, 'introns': introns_str,
                'three_prime': three_prime, 'is_ref': True
            })
        ref_df = pd.DataFrame(ref_data)

        # 2. Extract Observed Transcriptome from BAMs
        logging.info("Parsing BAM Alignments...")
        read_coords = []
        algn_counter = 0 
        for bam_path in bam_paths:
            bam_file = HTSeq.SAM_Reader(bam_path)
            for algn in bam_file:
                if not algn.aligned or algn.supplementary:
                    continue
                
                introns_list = []
                if algn.cigar is not None:
                    for cigar_elem in algn.cigar:
                        if cigar_elem.type == 'N': 
                            reference_region = cigar_elem.ref_iv
                            introns_list.append((reference_region.start, reference_region.end + 1))
                
                three_prime = algn.iv.end if algn.iv.strand == "+" else algn.iv.start
                read_coords.append({
                    'transcript_id': None, 'chrom': algn.iv.chrom, 'strand': algn.iv.strand,
                    'start': algn.iv.start, 'end': algn.iv.end, 'introns': str(tuple(introns_list)),
                    'three_prime': three_prime, 'is_ref': False
                })
                algn_counter += 1
        
        read_df = pd.DataFrame(read_coords)

        # 3. Merge and Cluster 3' Ends
        logging.info(f"Clustering 3' ends (distance: {args.ThreePrimeEnd_clustering_distance}bp)...")
        combined_df = pd.concat([ref_df, read_df], ignore_index=True)
        
        # Sort and calculate clusters based on proximity
        all_ends = combined_df[['chrom', 'strand', 'three_prime']].drop_duplicates().sort_values(['chrom', 'strand', 'three_prime'])
        all_ends['diff'] = all_ends.groupby(['chrom', 'strand'])['three_prime'].diff()
        all_ends['new_cluster'] = (all_ends['diff'] > args.ThreePrimeEnd_clustering_distance).astype(int)
        all_ends['cluster_id'] = all_ends.groupby(['chrom', 'strand'])['new_cluster'].cumsum()
        
        combined_df = pd.merge(combined_df, all_ends[['chrom', 'strand', 'three_prime', 'cluster_id']], on=['chrom', 'strand', 'three_prime'])

        # 4. Final Enrichment Logic (Collapse)
        enriched_isoforms = []
        # Group by intron structure and 3' cluster
        for (chrom, strand, introns, cluster_id), group in combined_df.groupby(['chrom', 'strand', 'introns', 'cluster_id']):
            ref_subset = group[group['is_ref'] == True]
            
            if not ref_subset.empty:
                # Use reference metadata if it exists in this structural group
                t_id = ref_subset['transcript_id'].iloc[0]
                t_start = ref_subset['start'].min()
                t_end = ref_subset['end'].max()
            else:
                # Create novel isoform ID if no reference matches
                t_id = f"novel_isoform_{len(enriched_isoforms)}"
                t_start = group['start'].min()
                t_end = group['end'].max()
            
            enriched_isoforms.append({
                'transcript_id': t_id, 'chrom': chrom, 'strand': strand,
                'start': t_start, 'end': t_end, 'introns': introns, 'cluster_id': cluster_id
            })
        
        isoforms_df_final = pd.DataFrame(enriched_isoforms)
        
        # Save master TSV for the 'assign' mode
        out_tsv = f"{args.output_prefix}_enriched.tsv"
        isoforms_df_final.to_csv(out_tsv, sep='\t', index=False)
        logging.info(f"Enriched TSV saved to {out_tsv}")

        # 5. Generate GTF file
        gtf_lines = []
        for _, row in isoforms_df_final.iterrows():
            attr = f'transcript_id "{row["transcript_id"]}"; gene_id "gene_cluster_{row["cluster_id"]}";'
            # Transcript row
            gtf_lines.append([row['chrom'], 'ONT_pipeline', 'transcript', row['start']+1, row['end'], '.', row['strand'], '.', attr])
            
            # Exon rows derived from intron structure
            intron_coords = ast.literal_eval(row['introns'])
            if not intron_coords:
                gtf_lines.append([row['chrom'], 'ONT_pipeline', 'exon', row['start']+1, row['end'], '.', row['strand'], '.', attr])
            else:
                current_start = row['start']
                for i_start, i_end in sorted(intron_coords):
                    gtf_lines.append([row['chrom'], 'ONT_pipeline', 'exon', current_start+1, i_start, '.', row['strand'], '.', attr])
                    current_start = i_end
                gtf_lines.append([row['chrom'], 'ONT_pipeline', 'exon', current_start+1, row['end'], '.', row['strand'], '.', attr])

        out_gtf_path = f"{args.output_prefix}_enriched.gtf"
        pd.DataFrame(gtf_lines).to_csv(out_gtf_path, sep='\t', header=False, index=False, quoting=3)
        logging.info(f"Enriched GTF saved to {out_gtf_path}")

    # --- MODE: ASSIGN ---
    elif args.mode == 'assign':
        logging.info("Loading enriched reference for assignment...")
        # Note: In assign mode, input_gtf_file should be the .tsv produced by 'enrich'
        master_ref = pd.read_csv(args.input_gtf_file, sep='\t')
        isoform_lookup = dict(zip(zip(master_ref.chrom, master_ref.strand, master_ref.introns), master_ref.transcript_id))
        
        assigned_count = 0
        unassigned_count = 0
        algn_counter = 0

        out_assigned = f"{args.output_prefix}_assignments.tsv.gz"
        out_unassigned = f"{args.output_prefix}_unassigned.tsv.gz"
        header = "alignment_id\tread_id\ttranscript_id\tchrom\tstart\tend\tstrand\tpt\tqs\n"

        with gzip.open(out_assigned, 'wt') as f_assign, gzip.open(out_unassigned, 'wt') as f_unassign:
            f_assign.write(header)
            f_unassign.write(header)

            for bam_path in bam_paths:
                bam_file = HTSeq.SAM_Reader(bam_path)
                for algn in bam_file:
                    if not algn.aligned or algn.supplementary:
                        continue
                    
                    try: pt_val = algn.optional_field("pt")
                    except (KeyError, AttributeError): pt_val = None
                    try: qs_val = algn.optional_field("qs")
                    except (KeyError, AttributeError): qs_val = None

                    introns_list = []
                    if algn.cigar is not None:
                        for cigar_elem in algn.cigar:
                            if cigar_elem.type == 'N':
                                reference_region = cigar_elem.ref_iv
                                introns_list.append((reference_region.start, reference_region.end + 1))
                    
                    introns_str = str(tuple(introns_list))
                    query_key = (algn.iv.chrom, algn.iv.strand, introns_str)
                    transcript_id = isoform_lookup.get(query_key, np.nan)

                    row_data = [
                        str(algn_counter), str(algn.read.name), str(transcript_id),
                        str(algn.iv.chrom), str(algn.iv.start), str(algn.iv.end),
                        str(algn.iv.strand), str(pt_val), str(qs_val)
                    ]
                    line = "\t".join(row_data) + "\n"

                    if pd.isna(transcript_id):
                        f_unassign.write(line)
                        unassigned_count += 1
                    else:
                        f_assign.write(line)
                        assigned_count += 1
                    
                    algn_counter += 1

        logging.info(f"Assignment Complete | Assigned: {assigned_count} | Unassigned: {unassigned_count}")
        logging.info(f"Results saved to: {out_assigned}")

if __name__ == "__main__":
    main()