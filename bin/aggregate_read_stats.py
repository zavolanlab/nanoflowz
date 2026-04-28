#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Aggregate read counts from Nextflow chunks")
    parser.add_argument('--input', nargs='+', required=True, help="List of chunk TSV files")
    parser.add_argument('--output', required=True, help="Output aggregated TSV")
    args = parser.parse_args()

    df_list = []
    for f in args.input:
        df = pd.read_csv(f, sep='\t', names=["sample_id", "step", "total_reads", "unmapped_reads", "um_reads", "mm_reads"])
        df_list.append(df)
    
    combined = pd.concat(df_list)
    
    # Group by sample and step, then sum the chunks
    agg = combined.groupby(["sample_id", "step"]).sum().reset_index()
    
    # Sort alphabetically (our step names will start with 01_, 02_, etc.)
    agg = agg.sort_values(["sample_id", "step"])
    
    agg.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()