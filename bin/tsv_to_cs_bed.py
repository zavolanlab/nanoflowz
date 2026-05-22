#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Convert read tags TSV to Cleavage Site BED")
    parser.add_argument('--tsv', required=True, help="Input tags.tsv.gz file")
    parser.add_argument('--out_bed', required=True, help="Output .bed.gz file")
    parser.add_argument('--cs_tag', required=True, help="The tag containing the cleavage site (e.g. XF)")
    args = parser.parse_args()

    # Read TSV in chunks to ensure memory efficiency on massive datasets
    chunk_iter = pd.read_csv(args.tsv, sep='\t', chunksize=1000000, dtype={args.cs_tag: str})
    
    counts = None
    for chunk in chunk_iter:
        # Filter out reads where the CS tag is missing or recorded as 'NA'
        chunk = chunk[chunk[args.cs_tag].notna() & (chunk[args.cs_tag] != 'NA')].copy()
        chunk[args.cs_tag] = chunk[args.cs_tag].astype(int)
        
        # Group this chunk by chromosome, strand, and cleavage site
        chunk_counts = chunk.groupby(['chrom', 'strand', args.cs_tag]).size().reset_index(name='count')
        
        if counts is None:
            counts = chunk_counts
        else:
            counts = pd.concat([counts, chunk_counts])
            # Regroup to sum the overlapping counts across concatenated chunks
            counts = counts.groupby(['chrom', 'strand', args.cs_tag])['count'].sum().reset_index()

    if counts is None or counts.empty:
        # Output empty BED if no sites were found
        pd.DataFrame(columns=['chrom', 'start', 'end', 'name', 'count', 'strand']).to_csv(
            args.out_bed, sep='\t', header=False, index=False, compression='gzip'
        )
        return

    # Calculate standard BED coordinates (0-based start, 1-based end)
    counts['end'] = counts[args.cs_tag]
    counts['start'] = counts['end'] - 1
    
    # Create the standard name string (chrom:start:end:strand)
    counts['name'] = counts['chrom'].astype(str) + ':' + \
                     counts['start'].astype(str) + ':' + \
                     counts['end'].astype(str) + ':' + \
                     counts['strand'].astype(str)
    
    # Reorder columns to match standard BED6 format
    bed_df = counts[['chrom', 'start', 'end', 'name', 'count', 'strand']]
    
    # Sort the BED file by chromosome and start position
    bed_df = bed_df.sort_values(['chrom', 'start'])
    
    # Export without header
    bed_df.to_csv(args.out_bed, sep='\t', header=False, index=False, compression='gzip')

if __name__ == "__main__":
    main()