#!/usr/bin/env python3

import warnings
warnings.simplefilter('ignore')

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import numpy as np
import pysam
import csv

def get_rc(seq):
    """Return the reverse complement of a sequence."""
    trans = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")
    return seq.translate(trans)[::-1]

def main():
    parser = ArgumentParser(description="Append polyA tails based on pt and fixed cleavage site tags.")
    parser.add_argument("--input_bam", required=True)
    parser.add_argument("--output_appended_bam", required=True)
    parser.add_argument("--output_skipped_bam", required=True)
    parser.add_argument("--stats_tsv", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--tag_orig_cs", default="XO")
    parser.add_argument("--tag_fixed_cs", default="XF")
    
    args = parser.parse_args()
    
    almnt_file = pysam.AlignmentFile(args.input_bam, "rb")
    bam_appended = pysam.AlignmentFile(args.output_appended_bam, "wb", header=almnt_file.header)
    bam_skipped = pysam.AlignmentFile(args.output_skipped_bam, "wb", header=almnt_file.header)
    
    reads_total = 0
    reads_appended = 0
    reads_skipped_no_pt = 0
    reads_skipped_pt_outside_valid_range = 0
    reads_skipped_error = 0
    
    for almnt in almnt_file:
        reads_total += 1
        
        # 1. Check for valid pt tag
        try:
            pt = almnt.get_tag('pt')
        except KeyError:
            pt = -1
            reads_skipped_no_pt += 1
            bam_skipped.write(almnt)
            continue
            
        if pt <= 0:
            bam_skipped.write(almnt)
            reads_skipped_pt_outside_valid_range += 1
            continue
            
        # 2. Get Cleavage Site Shift Difference
        try:
            OCS = almnt.get_tag(args.tag_orig_cs)
            FCS = almnt.get_tag(args.tag_fixed_cs)
            # The absolute difference is exactly how many bases were "rescued" from the softclip
            difference = abs(OCS - FCS)
        except KeyError:
            bam_skipped.write(almnt)
            reads_skipped_error += 1
            continue

        # 3. Orient to transcript 5' -> 3'
        if almnt.is_forward:
            read_seq = almnt.query_sequence
            read_qualstr = almnt.query_qualities
            cigar = list(almnt.cigar)
        else:
            read_seq = almnt.get_forward_sequence()
            read_qualstr = almnt.query_qualities[::-1] if almnt.query_qualities else None
            cigar = list(almnt.cigar)[::-1]
            
        # 4. Calculate exact truncation
        old_clip_len = cigar[-1][1] if cigar[-1][0] == 4 else 0
        trim_len = old_clip_len - difference
        
        if trim_len < 0:
            bam_skipped.write(almnt)
            reads_skipped_error += 1
            continue

        # 5. Modify Sequence
        new_read_seq = read_seq[:(-trim_len if trim_len > 0 else None)] + "A" * pt

        # 6. Modify Quality String
        if read_qualstr is not None:
            if trim_len > 0:
                trimmed_quals = read_qualstr[-trim_len:]
            else:
                trimmed_quals = read_qualstr[-(min(5, len(read_seq))):] if len(read_seq) > 0 else [30]
            
            quality_val = int(np.round(np.mean(trimmed_quals), 0)) if len(trimmed_quals) > 0 else 30
                
            new_read_qualstr = list(read_qualstr[:(-trim_len if trim_len > 0 else None)])
            new_read_qualstr.extend([quality_val] * pt)
        else:
            new_read_qualstr = None

        # 7. Modify CIGAR
        new_cigar = [[op, length] for op, length in cigar]
            
        if new_cigar[-1][0] == 4:
            new_cigar.pop() # Remove old 3' soft-clip
            
        if difference > 0:
            for idx in range(len(new_cigar)-1, -1, -1):
                if new_cigar[idx][0] == 0: # Add rescued bases to the last Match (M) block
                    new_cigar[idx][1] += difference
                    break
                    
        if pt > 0:
            new_cigar.append([4, pt]) # Append new soft-clip of length pt
            
        # Validation Check
        cigar_query_len = sum(length for op, length in new_cigar if op in [0, 1, 4, 7, 8])
        if len(new_read_seq) != cigar_query_len or (new_read_qualstr and len(new_read_seq) != len(new_read_qualstr)):
            bam_skipped.write(almnt)
            reads_skipped_error += 1
            continue
            
        # 8. Convert back to original BAM orientation
        if almnt.is_forward:
            almnt.query_sequence = new_read_seq
            almnt.query_qualities = new_read_qualstr
            almnt.cigar = new_cigar
        else:
            almnt.query_sequence = get_rc(new_read_seq)
            if new_read_qualstr is not None:
                almnt.query_qualities = new_read_qualstr[::-1]
            almnt.cigar = new_cigar[::-1]
            
        almnt.set_tag('pa', 1, 'i')
        bam_appended.write(almnt)
        reads_appended += 1
        
    almnt_file.close()
    bam_appended.close()
    bam_skipped.close()
    
    # 9. Write TSV Stats
    with open(args.stats_tsv, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerow(["sample_id", "chunk_filename", "total_alignments", "appended_alignments", "skipped_no_pt", "skipped_error", "skipped_pt_outside_valid_range"])
        chunk_name = os.path.basename(args.input_bam)
        writer.writerow([args.sample_id, chunk_name, reads_total, reads_appended, reads_skipped_no_pt, reads_skipped_error, reads_skipped_pt_outside_valid_range])

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!")
        sys.exit(1)