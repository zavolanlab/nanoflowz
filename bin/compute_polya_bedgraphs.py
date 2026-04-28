#!/usr/bin/env python3
import pysam
import argparse
import sys
import numpy as np
from collections import defaultdict

def weighted_mean(pt_list, nh_list):
    weights = 1.0 / np.array(nh_list)
    return np.average(pt_list, weights=weights)

def weighted_median(pt_list, nh_list):
    # Sort points and weights
    sorted_indices = np.argsort(pt_list)
    pts = np.array(pt_list)[sorted_indices]
    weights = (1.0 / np.array(nh_list))[sorted_indices]
    
    cumsum = np.cumsum(weights)
    cutoff = np.sum(weights) / 2.0
    return pts[cumsum >= cutoff][0]

def main():
    parser = argparse.ArgumentParser(description="Calculate weighted PolyA length bedgraphs.")
    parser.add_argument("--bam_in", required=True)
    parser.add_argument("--tag_cs", required=True)
    parser.add_argument("--metric", choices=['mean', 'median'], required=True)
    parser.add_argument("--include_multimappers", type=str, default="true")
    args = parser.parse_args()

    include_mm = args.include_multimappers.lower() == 'true'

    # Dictionaries to accumulate (pt, NH) lists by genomic coordinate
    # Keys: (chrom, pos)
    plus_data = defaultdict(lambda: ({'pt': [], 'nh': []}))
    minus_data = defaultdict(lambda: ({'pt': [], 'nh': []}))

    with pysam.AlignmentFile(args.bam_in, "rb") as bam:
        for r in bam:
            if r.is_unmapped:
                continue
                
            # Filter multi-mappers if requested
            if not include_mm and r.mapping_quality != 255:
                continue

            try:
                cs_1based = r.get_tag(args.tag_cs)
                pt_len = r.get_tag("pt")
                nh = r.get_tag("NH") if r.has_tag("NH") else 1
                
                if pt_len > 0:
                    cs_0based = cs_1based - 1
                    chrom = r.reference_name
                    
                    if r.is_reverse:
                        minus_data[(chrom, cs_0based)]['pt'].append(pt_len)
                        minus_data[(chrom, cs_0based)]['nh'].append(nh)
                    else:
                        plus_data[(chrom, cs_0based)]['pt'].append(pt_len)
                        plus_data[(chrom, cs_0based)]['nh'].append(nh)
            except KeyError:
                pass

    # Function to write data to bedgraph
    def write_bg(data_dict, filename):
        with open(filename, "w") as f:
            for (chrom, pos), vals in data_dict.items():
                if args.metric == 'mean':
                    val = weighted_mean(vals['pt'], vals['nh'])
                else:
                    val = weighted_median(vals['pt'], vals['nh'])
                
                # BedGraph format: chrom start end value
                f.write(f"{chrom}\t{pos}\t{pos+1}\t{val:.2f}\n")

    write_bg(plus_data, "plus.bg")
    write_bg(minus_data, "minus.bg")

if __name__ == "__main__":
    main()