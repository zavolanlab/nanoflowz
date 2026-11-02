#!/usr/bin/env python3
import argparse
import pysam
import re

def get_rc(seq):
    return seq.translate(str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn"))[::-1]

def get_args():
    parser = argparse.ArgumentParser(description="Annotate BAM with polyA motif presence.")
    parser.add_argument('--bam_in', required=True)
    parser.add_argument('--bam_out', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--motifs', nargs='+', required=True)
    parser.add_argument('--window_up', type=int, default=40)
    parser.add_argument('--window_down', type=int, default=20)
    parser.add_argument('--tag_cs', default="XF")
    parser.add_argument('--tag_motif', default="YF")
    return parser.parse_args()

def main():
    args = get_args()
    fasta = pysam.FastaFile(args.fasta)
    bam_in = pysam.AlignmentFile(args.bam_in, "rb")
    bam_out = pysam.AlignmentFile(args.bam_out, "wb", header=bam_in.header)

    motif_regexes = {m: re.compile(f'(?=({m.upper().replace("U", "T")}))') for m in args.motifs}

    # CACHE: Key: (chrom, strand, cs_pos) -> Value: "AAUAAA:(-21,-10);AUUAAA:(-35)"
    cs_cache = {}

    for read in bam_in:
        if read.is_unmapped:
            bam_out.write(read)
            continue

        try:
            # cleavage site tag is expected to be 1-based
            cs_pos = int(read.get_tag(args.tag_cs))
        except KeyError:
            bam_out.write(read)
            continue

        chrom = read.reference_name
        strand = '-' if read.is_reverse else '+'
        cache_key = (chrom, strand, cs_pos)

        if cache_key in cs_cache:
            motif_str = cs_cache[cache_key]
        else:
            try:
                chrom_len = fasta.get_reference_length(chrom)
            except KeyError:
                cs_cache[cache_key] = ""
                bam_out.write(read)
                continue

            # Convert to 0-based for precise FASTA fetching
            cs_0based = cs_pos - 1

            # Fetch strand-aware FASTA coordinates
            if strand == '+':
                fetch_start = cs_0based - args.window_up
                fetch_end = cs_0based + args.window_down + 1
            else:
                fetch_start = cs_0based - args.window_down
                fetch_end = cs_0based + args.window_up + 1

            pad_left = max(0, -fetch_start)
            pad_right = max(0, fetch_end - chrom_len)
            safe_start = max(0, fetch_start)
            safe_end = min(chrom_len, fetch_end)

            try:
                seq = fasta.fetch(chrom, safe_start, safe_end).upper()
                seq = ("N" * pad_left) + seq + ("N" * pad_right)
                if strand == '-':
                    seq = get_rc(seq)

                # Find Motifs
                motif_dict = {}
                for m, regex in motif_regexes.items():
                    matches = [match.start() for match in regex.finditer(seq)]
                    if matches:
                        # Because of symmetric padding, the CS is always exactly at index args.window_up
                        rel_positions = [str(idx - args.window_up) for idx in matches]
                        motif_dict[m] = f"({','.join(rel_positions)})"

                if motif_dict:
                    # Format: AAUAAA:(-38,-10);AUUAAA:(-35)
                    motif_str = ";".join([f"{k}:{v}" for k, v in motif_dict.items()])
                else:
                    motif_str = ""

            except (KeyError, ValueError):
                motif_str = ""

            cs_cache[cache_key] = motif_str

        # Apply tag and write
        if motif_str:
            read.set_tag(args.tag_motif, motif_str, value_type='Z')

        bam_out.write(read)

    bam_in.close()
    bam_out.close()

if __name__ == "__main__":
    main()