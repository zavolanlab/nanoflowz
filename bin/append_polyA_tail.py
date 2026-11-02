#!/usr/bin/env python3
import pysam
import argparse
import sys

def get_args():
    parser = argparse.ArgumentParser(description="Append polyA tails to reads")
    parser.add_argument('--input_bam', required=True)
    parser.add_argument('--output_appended_bam', required=True)
    parser.add_argument('--output_skipped_bam', required=True)
    parser.add_argument('--stats_tsv', required=True)
    parser.add_argument('--sample_id', required=True)
    parser.add_argument('--tag_orig_cs', default="XO")
    parser.add_argument('--tag_fixed_cs', default="XF")
    return parser.parse_args()

def main():
    args = get_args()
    
    bam_in = pysam.AlignmentFile(args.input_bam, "rb")
    bam_appended = pysam.AlignmentFile(args.output_appended_bam, "wb", header=bam_in.header)
    bam_skipped = pysam.AlignmentFile(args.output_skipped_bam, "wb", header=bam_in.header)
    
    appended_count = 0
    skipped_count = 0
    
    for read in bam_in:
        # Skip unmapped reads
        if read.is_unmapped:
            bam_skipped.write(read)
            skipped_count += 1
            continue
            
        # Extract tags safely
        try:
            pt = int(read.get_tag('pt'))
            if pt <= 0:
                raise ValueError
            orig_cs = int(read.get_tag(args.tag_orig_cs))
            fixed_cs = int(read.get_tag(args.tag_fixed_cs))
        except (KeyError, ValueError):
            bam_skipped.write(read)
            skipped_count += 1
            continue
            
        rev = read.is_reverse
        cigar = read.cigartuples
        seq = read.query_sequence
        qual = read.query_qualities
        
        try:
            if not rev:  # '+' strand
                # Outward shift (rescue): fixed > orig -> diff > 0
                # Inward shift (discard): fixed < orig -> diff < 0
                diff = fixed_cs - orig_cs
                
                sc_len = cigar[-1][1] if cigar[-1][0] == 4 else 0
                cut_right = sc_len - diff
                
                if cut_right > 0:
                    new_seq = seq[:-cut_right] + ("A" * pt)
                    new_qual = list(qual[:-cut_right]) + [30] * pt
                else:
                    new_seq = seq + ("A" * pt)
                    new_qual = list(qual) + [30] * pt
                    
                # Fix CIGAR
                if cigar[-1][0] == 4:
                    new_cigar = cigar[:-1]
                else:
                    new_cigar = list(cigar)
                    
                last_match_idx = len(new_cigar) - 1
                while last_match_idx >= 0 and new_cigar[last_match_idx][0] not in [0, 7, 8]:
                    last_match_idx -= 1
                    
                if last_match_idx >= 0:
                    adj_match = new_cigar[last_match_idx][1] + diff
                    if adj_match > 0:
                        new_cigar[last_match_idx] = (new_cigar[last_match_idx][0], adj_match)
                    else:
                        raise ValueError("Inward shift exceeds match length")
                        
                new_cigar.append((4, pt))
                
            else:  # '-' strand
                # Outward shift: orig > fixed -> diff > 0
                # Inward shift: orig < fixed -> diff < 0
                diff = orig_cs - fixed_cs
                
                sc_len = cigar[0][1] if cigar[0][0] == 4 else 0
                cut_left = sc_len - diff
                
                if cut_left > 0:
                    new_seq = ("T" * pt) + seq[cut_left:]
                    new_qual = [30] * pt + list(qual[cut_left:])
                else:
                    new_seq = ("T" * pt) + seq
                    new_qual = [30] * pt + list(qual)
                    
                # Fix CIGAR
                if cigar[0][0] == 4:
                    new_cigar = cigar[1:]
                else:
                    new_cigar = list(cigar)
                    
                first_match_idx = 0
                while first_match_idx < len(new_cigar) and new_cigar[first_match_idx][0] not in [0, 7, 8]:
                    first_match_idx += 1
                    
                if first_match_idx < len(new_cigar):
                    adj_match = new_cigar[first_match_idx][1] + diff
                    if adj_match > 0:
                        new_cigar[first_match_idx] = (new_cigar[first_match_idx][0], adj_match)
                    else:
                        raise ValueError("Inward shift exceeds match length")
                        
                new_cigar = [(4, pt)] + new_cigar
                
                # Calculate the physical POS shift for the minus strand ---
                # If diff is negative (inward shift), new_start increases.
                # If diff is positive (outward shift), new_start decreases.
                new_start = read.reference_start - diff
            
            # Apply adjustments to the PySAM object
            read.query_sequence = new_seq
            read.query_qualities = new_qual
            read.cigartuples = new_cigar
            
            # Apply the POS shift for the minus strand
            if rev:
                read.reference_start = new_start
            
            bam_appended.write(read)
            appended_count += 1
            
        except Exception as e:
            # Safely skip reads that fail extreme boundary math
            bam_skipped.write(read)
            skipped_count += 1

    bam_in.close()
    bam_appended.close()
    bam_skipped.close()
    
    # Write stats tracker
    with open(args.stats_tsv, "w") as f:
        f.write("sample_id\tappended_reads\tskipped_reads\n")
        f.write(f"{args.sample_id}\t{appended_count}\t{skipped_count}\n")

if __name__ == "__main__":
    main()