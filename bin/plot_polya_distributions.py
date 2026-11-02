#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import re

def natural_sort_key(s):
    parsed = re.split('([0-9]+)', str(s))
    return [(0, int(text)) if text.isdigit() else (1, text.lower()) for text in parsed]

def weighted_mean(group):
    weights = 1.0 / group['NH']
    return np.average(group['pt'], weights=weights)

def weighted_median(group):
    df_sorted = group.sort_values('pt')
    weights = 1.0 / df_sorted['NH']
    cumsum = weights.cumsum()
    cutoff = weights.sum() / 2.0
    # Return the first pt value where cumulative weight >= 50%
    return df_sorted[cumsum >= cutoff]['pt'].iloc[0]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', nargs='+', required=True)
    parser.add_argument('--output_prefix', required=True)
    args = parser.parse_args()

    df_list = []
    for f in args.input:
        sample_id = os.path.basename(f).split('.tags.tsv')[0]
        df = pd.read_csv(f, sep='\t', compression='gzip')
        df['sample_id'] = sample_id
        df_list.append(df)

    combined = pd.concat(df_list, ignore_index=True)

    # Clean data
    combined['pt'] = pd.to_numeric(combined['pt'], errors='coerce')
    combined['NH'] = pd.to_numeric(combined['NH'], errors='coerce').fillna(1)
    valid_pt = combined.dropna(subset=['pt'])
    valid_pt = valid_pt[valid_pt['pt'] > 0]

    # Process Chromosomes
    top_chrs = valid_pt['chrom'].value_counts().nlargest(30).index
    unique_top_chrs = [c for c in valid_pt['chrom'].unique() if c in top_chrs]
    sorted_chrs = sorted(unique_top_chrs, key=natural_sort_key)
    
    valid_pt = valid_pt[valid_pt['chrom'].isin(sorted_chrs)]
    valid_pt['chrom'] = pd.Categorical(valid_pt['chrom'], categories=sorted_chrs, ordered=True)

    sns.set_theme(style="whitegrid")

    # ==========================================
    # PLOT 1: All Reads (Deduplicated to avoid MM overcounting)
    # ==========================================
    unique_reads = valid_pt.drop_duplicates(subset=['sample_id', 'read_id'])

    plt.figure(figsize=(10, 6))
    sns.boxplot(data=unique_reads, x='sample_id', y='pt', showfliers=False)
    plt.title('PolyA Tail Length (All Reads)')
    plt.xticks(rotation=45, ha='right'); plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_all_reads_pooled.pdf")
    plt.close()

    plt.figure(figsize=(16, 6))
    sns.boxplot(data=unique_reads, x='sample_id', y='pt', hue='chrom', showfliers=False)
    plt.title('PolyA Tail Length (All Reads, per chrom)')
    plt.xticks(rotation=45, ha='right'); plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_all_reads_per_chr.pdf")
    plt.close()

    # ==========================================
    # Filter for Gene-Assigned Reads
    # ==========================================
    gene_assigned = valid_pt[~valid_pt['XT'].isna() & (valid_pt['XT'] != 'NA')]

    # ==========================================
    # PLOT 2 & 3: Per-Gene Weighted Math
    # ==========================================
    # Calculate Weighted Median
    gene_median = gene_assigned.groupby(['sample_id', 'chrom', 'XT'], observed=True).apply(weighted_median).reset_index(name='pt')
    
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=gene_median, x='sample_id', y='pt', showfliers=False)
    plt.title('PolyA Tail Length (Per-Gene Median)')
    plt.xticks(rotation=45, ha='right'); plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_per_gene_median_pooled.pdf")
    plt.close()

    plt.figure(figsize=(16, 6))
    sns.boxplot(data=gene_median, x='sample_id', y='pt', hue='chrom', showfliers=False)
    plt.title('PolyA Tail Length (Per-Gene Median per chrom)')
    plt.xticks(rotation=45, ha='right'); plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_per_gene_median_per_chr.pdf")
    plt.close()

    # Calculate Weighted Mean
    gene_mean = gene_assigned.groupby(['sample_id', 'chrom', 'XT'], observed=True).apply(weighted_mean).reset_index(name='pt')
    
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=gene_mean, x='sample_id', y='pt', showfliers=False)
    plt.title('PolyA Tail Length (Per-Gene Mean)')
    plt.xticks(rotation=45, ha='right'); plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_per_gene_mean_pooled.pdf")
    plt.close()

    plt.figure(figsize=(16, 6))
    sns.boxplot(data=gene_mean, x='sample_id', y='pt', hue='chrom', showfliers=False)
    plt.title('PolyA Tail Length (Per-Gene Mean per chrom)')
    plt.xticks(rotation=45, ha='right'); plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_per_gene_mean_per_chr.pdf")
    plt.close()

if __name__ == "__main__":
    main()