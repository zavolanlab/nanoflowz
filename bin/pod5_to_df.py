#!/usr/bin/env python3
import pod5 as p5
import pandas as pd
import numpy as np
import pysam
import argparse
import logging
import sys

# Configure logging to stderr
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stderr
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pod5', required=True)
    parser.add_argument('--bam', required=True)
    parser.add_argument('--sample_id', required=True)
    args = parser.parse_args()

    try:
        logging.info(f"Starting processing for sample: {args.sample_id}")
        logging.info(f"Loading BAM file: {args.bam}")
        
        # 1. Load BAM data including Poly-A estimation tags
        bam_data = {}
        try:
            with pysam.AlignmentFile(args.bam, "rb", check_sq=False) as bam:
                for read in bam.fetch(until_eof=True):
                    bam_data[read.query_name] = {
                        'seq': read.query_sequence,
                        'mv': read.get_tag("mv") if read.has_tag("mv") else None,
                        'ts': read.get_tag("ts") if read.has_tag("ts") else 0,
                        'ns': read.get_tag("ns") if read.has_tag("ns") else 0,
                        'pt': read.get_tag("pt") if read.has_tag("pt") else 0,  # Poly-A length
                        'pa': read.get_tag("pa") if read.has_tag("pa") else None # Signal anchors
                    }
            logging.info(f"Successfully loaded {len(bam_data)} reads from BAM file")
        except Exception as e:
            logging.error(f"Failed to load BAM file {args.bam}: {e}")
            raise

        logging.info(f"Processing POD5 file: {args.pod5}")
        processed_reads = 0
        
        try:
            with p5.Reader(args.pod5) as reader:
                for read_record in reader.reads():
                    read_id = str(read_record.read_id)
                    if read_id not in bam_data:
                        continue

                    try:
                        # Signal extraction
                        signal = read_record.signal
                        sample_rate = read_record.run_info.sample_rate
                        time = np.arange(len(signal)) / sample_rate
                        raw_signal_df = pd.DataFrame({'time': time, 'signal': signal})

                        info = bam_data[read_id]
                        
                        # Initialize polyA column (default to -1 for non-tail region)
                        raw_signal_df['polyA'] = -1

                        # 2. Map Poly-A Region if tags exist
                        if info['pa'] is not None:
                            # pa array indices: 1 = start of polyA, 2 = end of polyA
                            pa_start = info['pa'][1]
                            pa_end = info['pa'][2]
                            
                            # Fill the polyA column with the estimated length (pt) for that region
                            # This makes it easy to filter/color the tail in plots
                            raw_signal_df.loc[pa_start:pa_end, 'polyA'] = info['pt']

                        # 3. Handle Moves and Base Mapping
                        if info['mv'] is not None:
                            stride = info['mv'][0]
                            moves = info['mv'][1:]
                            sequence = info['seq']
                            
                            # Build 'ann' (moves) array
                            a = info['ts'] * [-1]
                            for elem in moves:
                                a += [elem] * stride
                            a += (len(raw_signal_df) - info['ns']) * [-1]
                            a = [-2] * max((len(raw_signal_df) - len(a)), 0) + a
                            
                            # Trim or pad 'a' to match signal length exactly
                            if len(a) > len(raw_signal_df):
                                a = a[:len(raw_signal_df)]
                            else:
                                a += [-1] * (len(raw_signal_df) - len(a))
                            
                            raw_signal_df['ann'] = a

                            # 4. Map Nucleotides to Signal
                            base_labels = ['N'] * len(raw_signal_df)
                            seq_idx = 0
                            seq_len = len(sequence)

                            for i, val in enumerate(a):
                                if val == 1 and seq_idx < seq_len:
                                    current_base = sequence[seq_idx]
                                    base_labels[i] = current_base
                                    seq_idx += 1
                                elif val == 0 and seq_idx > 0:
                                    base_labels[i] = sequence[seq_idx - 1]

                            raw_signal_df['base'] = base_labels

                        # Save annotated CSV
                        output_name = f"{args.sample_id}_{read_id}_mapped.csv"
                        raw_signal_df.to_csv(output_name, index=False)
                        processed_reads += 1
                        logging.info(f"Successfully processed read {read_id} -> {output_name}")
                    except Exception as e:
                        logging.error(f"Error processing read {read_id}: {e}")
                        raise
            
            logging.info(f"Successfully processed {processed_reads} reads from POD5 file")
        except Exception as e:
            logging.error(f"Failed to process POD5 file {args.pod5}: {e}")
            raise
        
        logging.info(f"Completed processing for sample: {args.sample_id}")
    except Exception as e:
        logging.error(f"Fatal error in main processing: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()