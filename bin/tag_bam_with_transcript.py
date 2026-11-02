#!/usr/bin/env python3
"""
Assign transcript IDs to BAM reads using an enriched isoform TSV produced by
assign_ONT_reads_to_isoforms.py (enrich mode), and store the result as a BAM
tag.  Reads whose intron chain does not match any isoform are written unmodified.
"""

import logging
import sys
from argparse import ArgumentParser
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def build_introns_str(read: pysam.AlignedSegment) -> str:
    """Return a string representation of splice junctions identical to the one
    produced by assign_ONT_reads_to_isoforms.py so that the lookup key matches."""
    introns = []
    if read.cigartuples:
        ref_pos = read.reference_start
        for op, length in read.cigartuples:
            if op == 3:  # N  – intron / spliced-out region
                introns.append((ref_pos, ref_pos + length))
                ref_pos += length
            elif op in (0, 2, 7, 8):  # M, D, =, X  – consume reference
                ref_pos += length
    return str(tuple(introns))


def main():
    parser = ArgumentParser(description="Tag BAM reads with assigned transcript IDs")
    parser.add_argument("--bam_in",       required=True, help="Input BAM file (coordinate-sorted)")
    parser.add_argument("--enriched_tsv", required=True, help="Enriched isoform TSV from assign_ONT_reads_to_isoforms.py (--mode enrich)")
    parser.add_argument("--bam_out",      required=True, help="Output BAM file")
    parser.add_argument("--tag",          default="YT",  help="BAM tag name for transcript ID (default: YT)")
    args = parser.parse_args()

    # ------------------------------------------------------------------ #
    # 1.  Load enriched isoform table and build a lookup by intron chain  #
    # ------------------------------------------------------------------ #
    logging.info(f"Loading enriched isoform TSV: {args.enriched_tsv}")
    master_ref = pd.read_csv(args.enriched_tsv, sep='\t')

    master_ref['three_prime'] = master_ref.apply(
        lambda r: r['end'] if r['strand'] == '+' else r['start'], axis=1
    )

    # key: (chrom, strand, introns_str)  →  list of (three_prime_coord, transcript_id)
    isoform_candidates: dict = defaultdict(list)
    for _, row in master_ref.iterrows():
        key = (row['chrom'], row['strand'], row['introns'])
        isoform_candidates[key].append((row['three_prime'], row['transcript_id']))

    logging.info(f"Loaded {len(master_ref)} isoforms ({len(isoform_candidates)} unique intron chains)")

    # ------------------------------------------------------------------ #
    # 2.  Iterate over the BAM, tag matched reads, write output           #
    # ------------------------------------------------------------------ #
    assigned = 0
    unassigned = 0

    with pysam.AlignmentFile(args.bam_in, "rb") as bam_in, \
         pysam.AlignmentFile(args.bam_out, "wb", header=bam_in.header) as bam_out:

        for read in bam_in.fetch(until_eof=True):
            if read.is_unmapped or read.is_supplementary:
                bam_out.write(read)
                continue

            chrom  = read.reference_name
            strand = '-' if read.is_reverse else '+'
            introns_str = build_introns_str(read)

            candidate_list = isoform_candidates.get((chrom, strand, introns_str))

            if candidate_list is not None:
                read_three_prime = read.reference_end if strand == '+' else read.reference_start
                transcript_id = min(candidate_list, key=lambda x: abs(x[0] - read_three_prime))[1]
                read.set_tag(args.tag, str(transcript_id), value_type='Z')
                assigned += 1
            else:
                unassigned += 1

            bam_out.write(read)

    logging.info(f"Assignment complete | Tagged: {assigned} | Untagged: {unassigned}")

    # ------------------------------------------------------------------ #
    # 3.  Index the output BAM                                            #
    # ------------------------------------------------------------------ #
    logging.info(f"Indexing {args.bam_out}")
    pysam.index(args.bam_out)
    logging.info("Done")


if __name__ == "__main__":
    main()
