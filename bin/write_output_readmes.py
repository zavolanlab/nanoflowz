#!/usr/bin/env python3
"""
Writes README.md files that explain the nanoflowz outputs: one in the results directory (overview, samples,
BAM tags), one in every output subfolder (files, formats, how they were made) and one in the run directory.
Only the Python standard library is used, so the script also runs outside the pipeline, e.g. to document
the results of an older run:

    write_output_readmes.py --outdir <rundir>/results --params <rundir>/run_params.json \
        --rundir_readme <rundir>/README.md --only_existing
"""

import argparse
import csv
import datetime
import json
import os
import re
from collections import OrderedDict

# Defaults of nextflow.config; values from --params files override them
DEFAULTS = {
    "tsv": None, "rundir": "./run_01", "outdir": None, "workdir_param": None, "conda_envs_dir": None,
    "dorado": None, "model": None, "polyA": None, "ref": None, "reference_gtf": None,
    "viz_gene_list": "", "num_reads": 10,
    "tag_basecaller_ts": "TS", "tag_original_ts": "ZS", "tag_orig_cs": "XO", "tag_fixed_cs": "XF",
    "shift_ambiguous_cs": False, "tag_quantification_cs": "XF",
    "min_phred_score": 0, "tag_phred_mapped": "ZM", "tag_phred_softclipped": "ZC",
    "polya_tail_metric": "median", "include_multimappers": True,
    "motif_list": "AAUAAA AUUAAA UAUAAA AGUAAA AAUACA CAUAAA AAUAUA GAUAAA AAUGAA AAGAAA ACUAAA AAUAGA",
    "motif_window_up": 40, "motif_window_down": 20, "motif_anchor": "start", "motif_weighting": "uniform",
    "motif_score_bins": 10, "tag_motif_info": "YF", "tag_transcript_id": "YT",
    "visualize_polyA_tail_length_distribution_flag": False, "visualize_cleavage_site_motifs_flag": False,
    "bam_export_tags": "pt XT YT NH XF YF", "merge_non_polyA": True, "merge_skipped_pt": True,
}

# Values hard-coded in main.nf (keep in sync)
FIXED = {"polyA_min_A_percent": 80, "polyA_min_length": 8, "umi_target_len": 26, "umi_tag": "RX",
         "three_prime_cluster_nt": 25}


def as_bool(v):
    return str(v).strip().lower() in ("true", "1", "yes")


class Fill(dict):
    """format_map dictionary that leaves unknown placeholders untouched."""
    def __missing__(self, key):
        return "{" + key + "}"


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join([":---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x).replace("|", "\\|") for x in r) + " |")
    return "\n".join(out)


# ----------------------------------------------------------------------------------------------------------
# Output folders: (path relative to outdir, condition on params, one-line summary, README body).
# Placeholders {name} are filled from the run parameters and FIXED.
# ----------------------------------------------------------------------------------------------------------
FOLDERS = [
    ("basecalling", lambda p: True,
     "poly(A) reads of each sample, aligned, merged over POD5 chunks, **before** UMI deduplication",
     """
**Produced by:** process `merge_bams` (`samtools merge`).

| file | content |
|:---|:---|
| `<sample>.dorado.sorted.bam` (+ `.bai`) | coordinate-sorted alignments of one sample |

Despite the folder name, these are not raw basecalls. Each POD5 file was processed separately, and the chunks were then merged per sample:
1. `dorado basecaller` with `--estimate-poly-a` (tag `pt`);
2. reads reverse-complemented to the mRNA orientation (`{tag_original_ts}`);
3. aligned with `minimap2 -ax splice:hq --secondary=no -Y`, supplementary alignments removed;
4. UMIs in `{umi_tag}` padded or trimmed to {umi_target_len} nt;
5. alignment 3' ends corrected (`{tag_orig_cs}` → `{tag_fixed_cs}`);
6. only poly(A) reads kept;
7. the poly(A) tail appended.

So these files hold the same reads as the later BAMs, but with PCR duplicates still present. See the README in the results folder for the steps and tags.

**Use:** read numbers before deduplication, or repeating the deduplication with other settings. For anything else use the later BAMs (`../assignments/`).
"""),
    ("map_genome_merged_UMIdedup", lambda p: True,
     "UMI-deduplicated alignments with corrected NH tag and MAPQ: the basis of all quantifications",
     """
**Produced by:** `umi_tools dedup` (process `umi_tools_dedup`), then `redefine_qual_and_NHtag` from zavolab_pyutils (process `redefine_nh_tags`).

| file | content |
|:---|:---|
| `<sample>.redefined_NH.sorted.bam` (+ `.bai`) | coordinate-sorted, deduplicated alignments |

- **Deduplication:** `umi_tools dedup --extract-umi-method=tag --umi-tag={umi_tag} --method unique`. Reads with the same alignment position, strand and UMI count once; the one with the best mapping quality is kept.
- **NH and MAPQ:** recomputed from the alignments left after deduplication. NH is the number of alignments; MAPQ is 255 for unique reads (as in STAR) and 0 for multi-mappers.
- minimap2 is run with `--secondary=no` and supplementary alignments are removed, so normally every read has one alignment (NH = 1).

**Use:** input of the cleavage-site bigWigs (`../cleavage_sites_bigwigs/`) and of featureCounts (`../gene_assignments/`). The later BAMs contain the same reads with extra tags.
"""),
    ("gene_assignments", lambda p: True,
     "gene-level read counts (featureCounts) and BAMs with the gene of every read",
     """
**Produced by:** process `assign_alignments_to_genes`: `featureCounts -L -M -O -s 1 -R CORE` on `../map_genome_merged_UMIdedup/`, with `{reference_gtf}`; tags added by `bin/tag_bam_with_fc.py`.

| file | content |
|:---|:---|
| `<sample>.featureCounts.txt` | gene × read-count table: `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length` and the count (last column) |
| `<sample>.featureCounts.txt.summary` | reads per assignment status (`Assigned`, `Unassigned_NoFeatures`, `Unassigned_Ambiguity`, ...) |
| `<sample>.gene_assigned.sorted.bam` (+ `.bai`) | the deduplicated alignments with tags `XT` (gene ID) and `XS` (featureCounts status) |

- `-L`: long-read mode. `-s 1`: stranded, because reads are in the mRNA orientation.
- `-M -O`: multi-mapping reads, and reads overlapping several genes, are counted for every gene they hit (no fractional counting). For such reads `XT` lists all genes, comma-separated.
- `XT` is missing for unassigned reads; `XS` gives the reason.
- Gene IDs are those of the GTF (with version suffix). `Length` is the summed exon length of the gene; it is not needed for 3'-end counts.

**Use:** the count tables are the gene expression matrix (one column per sample, deduplicated reads).
"""),
    ("motif_annotated_bams", lambda p: True,
     "the gene-assigned BAMs plus the poly(A)-signal motifs found near each read's cleavage site",
     """
**Produced by:** process `annotate_reads_with_motifs` (`bin/annotate_motifs_in_bam.py`).

| file | content |
|:---|:---|
| `<sample>.motif_annotated.bam` (+ `.bai`) | alignments of `../gene_assignments/` plus tag `{tag_motif_info}` |

`{tag_motif_info}` lists the motifs found in the genome, on the read's strand, in a window from {motif_window_up} nt upstream to {motif_window_down} nt downstream of the cleavage site (tag `{tag_quantification_cs}`). Motifs searched: `{motif_list}`.

The format is `MOTIF:(pos,pos);MOTIF:(pos)`, e.g. `AAUAAA:(-23,13)`. `pos` is the start of the motif relative to the cleavage site; negative values are upstream. The tag is empty when no motif was found (`NA` in `../read_tags/`).

**Use:** input of the isoform steps. For analysis use `../assignments/`, which has all tags.
"""),
    ("transcriptome", lambda p: True,
     "the \"enriched\" transcript annotation: reference transcripts plus novel isoforms observed in the reads of all samples",
     """
**Produced by:** process `transcriptome_annotation_enrichment` (`bin/assign_ONT_reads_to_isoforms.py --mode enrich`), from the BAMs of all samples in `../motif_annotated_bams/` and `{reference_gtf}`.

| file | content |
|:---|:---|
| `master_enriched.tsv` | one row per isoform (reference and novel), with the columns below |
| `master_enriched.gtf` | the reference GTF plus one `transcript` and its `exon` lines per novel isoform (`gene_id "novel_gene_<n>"`) |

Columns of `master_enriched.tsv`:
- `transcript_id`: reference ID or `novel_isoform_<n>`;
- `chrom`, `strand`, `start`, `end`;
- `introns`: the intron chain as a tuple of (start, end) pairs;
- `three_prime`: the 3' end;
- `is_ref`;
- `cluster_id`: the cluster of 3' ends.

**How novel isoforms are defined.** Each read and each reference transcript is described by its intron chain and its 3' end. 3' ends on the same chromosome and strand are clustered: neighbours ≤ {three_prime_cluster_nt} nt apart join a cluster. Every combination of intron chain and 3'-end cluster that contains reads but no reference transcript becomes a novel isoform, spanning the reads' minimal start to maximal end.

Several kinds of reads therefore also give novel isoforms:
- unspliced reads;
- reads whose 3' end lies away from annotated transcript ends;
- 5'-truncated reads, whose intron chain is only the 3' part of a reference chain;
- reads with a splice junction shifted by alignment errors.

Many novel isoforms are supported by a single read; filter them by read support before use. `novel_gene_<n>` is not linked to the reference gene; use `XT` in the BAMs or `../read_tags/` for the gene.

**Use:** the isoform catalogue behind tag `{tag_transcript_id}` in `../assignments/`; the GTF can be loaded in IGV.
"""),
    ("assignments", lambda p: True,
     "**final BAMs**: all reads with all tags, including the transcript isoform of every read",
     """
**Produced by:** process `read_to_transcript_assignment` (`bin/tag_bam_with_transcript.py`), from `../motif_annotated_bams/` and `../transcriptome/master_enriched.tsv`.

| file | content |
|:---|:---|
| `<sample>.transcript_assigned.bam` (+ `.bai`) | coordinate-sorted, UMI-deduplicated poly(A) reads with all tags |

`{tag_transcript_id}` holds the isoform, matched by:
1. the same chromosome, strand and exact intron chain;
2. among those, the isoform whose 3' end is closest to the read's 3' end.

Reads without a matching intron chain carry no `{tag_transcript_id}`. Since the enriched annotation is built from the same reads, nearly every read gets one. Many of them are novel isoforms (`novel_isoform_<n>`; see `../transcriptome/README.md`).

All tags are explained in the README of the results folder.

**Use:** the most complete per-read data. The light-weight table version is in `../read_tags/`. For IGV, load these BAMs together with `../cleavage_sites_bigwigs/` and `../transcriptome/master_enriched.gtf`.
"""),
    ("read_tags", lambda p: True,
     "one table per sample with one row per alignment and its key tags (the light-weight version of the final BAMs)",
     """
**Produced by:** process `extract_read_tags_tsv` (`bin/extract_bam_tags.py`), from `../assignments/`.

| file | content |
|:---|:---|
| `<sample>.tags.tsv.gz` | tab-separated, gzipped; header line |

Columns:
- `read_id`;
- `chrom`;
- `strand`: strand of the alignment, which is the strand of the RNA (reads are in the mRNA orientation);
- `read_length`: length of the aligned read sequence, including soft clips and the appended poly(A) tail;
- the tags `{bam_export_tags}`.

The tags are explained in the README of the results folder. A missing tag is written as `NA`.

Multi-mapping reads are {multimapper_rows}.

**Use:** per-read analyses in Python or R without reading BAMs, e.g. counts per gene (`XT`), per isoform (`YT`) or per cleavage site (`XF` with `chrom` and `strand`), and poly(A) tail lengths (`pt`).
"""),
    ("cleavage_sites_bigwigs", lambda p: True,
     "read support of every cleavage site, as strand-specific bigWigs and as a BED table",
     """
**Produced by:** process `make_bigwig_for_cleavage_sites` (`bin/extract_cs_bigwig.py`, bedtools, `bedGraphToBigWig`), from `../map_genome_merged_UMIdedup/`.

| file | content |
|:---|:---|
| `<sample>.plus.bigwig`, `<sample>.minus.bigwig` | per genomic position: the number of reads whose cleavage site is there, on the + / − strand |
| `<sample>.cs.bed.gz` | the same sites as BED6: `chrom`, `start` (0-based), `end` (= `start` + 1), `name` (`chrom:start:end:strand`), `score` (number of reads), `strand` |

- The cleavage site of a read is tag `{tag_quantification_cs}`: the 3' end of the alignment after the soft-clip correction, i.e. the last nucleotide before the poly(A) tail.
- Each read counts with weight 1/NH. Multi-mappers are {multimapper_included}.
- Counts are deduplicated reads. They are not normalised and not clustered: nearby sites are separate entries.
- A bigWig file is empty (0 bytes) when a sample has no reads on that strand.

**Use:**
- IGV tracks of 3'-end usage.
- The BED table is the input for cleavage-site clustering and for PAQR-AI.
"""),
    ("polya_length_bigwigs", lambda p: True,
     "poly(A) tail length at every cleavage site, as strand-specific bigWigs",
     """
**Produced by:** process `make_bigwig_for_polya_length` (`bin/compute_polya_bedgraphs.py`, `bedGraphToBigWig`), from `../assignments/`.

| file | content |
|:---|:---|
| `<sample>.polya_{polya_tail_metric}.plus.bigwig`, `<sample>.polya_{polya_tail_metric}.minus.bigwig` | per cleavage site (tag `{tag_quantification_cs}`): the {polya_tail_metric} poly(A) tail length (nt) of the reads ending there |

- Tail lengths are Dorado's estimates (tag `pt`), from reads with `pt` > 0.
- Reads are weighted by 1/NH. Multi-mappers are {multimapper_included}.
- The number of reads behind each value is in `../cleavage_sites_bigwigs/`; values from single reads are noisy.

**Use:** IGV tracks; poly(A) length per 3' end. For statistics, use the per-read `pt` in `../read_tags/`.
"""),
    ("QC", lambda p: True,
     "tables that track the reads through the pipeline, plus BAMs of the reads that were left out",
     """
**Produced by:** processes `merge_ts_stats`, `merge_polyA_stats`, `merge_pt_stats` and `aggregate_read_stats` (`bin/count_reads.sh`, `bin/aggregate_read_stats.py`).

| file | one row per | columns |
|:---|:---|:---|
| `master_TS_orientation_stats.tsv` | POD5 chunk | `reads_with_TS`, `reads_without_TS` (no strand from Dorado), `reads_reverse_complemented` (reads Dorado placed on the − strand, flipped to the mRNA orientation) |
| `master_polyA_stats.tsv` | POD5 chunk | `polyA_alignments` (kept), `non_polyA_alignments` (left out, see `non_polyA_alignments/`), `low_q_polyA_alignments` (poly(A) reads below the base-quality threshold, left out) |
| `master_pt_append_stats.tsv` | POD5 chunk (no chunk name) | `appended_reads` (poly(A) tail appended, kept), `skipped_reads` (no valid Dorado `pt`, left out, see `skipped_pt_alignments/`) |
| `master_read_tracking_stats.tsv` | sample × step | `total_reads` (distinct read IDs), `unmapped_reads`, `um_reads` (one alignment), `mm_reads` (> 1 alignment) |

Steps in `master_read_tracking_stats.tsv`:
- `01_dorado_basecall` to `04_normalized_umi_lengths`: all reads; aligned from step 03.
- `05`: unmapped reads removed.
- `06_extracted_polyA_reads`: poly(A) reads only.
- `07_appended_polyA_tails`: reads with a valid `pt`.
- `08_umi_deduped`: after UMI deduplication.
- `09` to `11`: tags added, no reads removed.

Step 08 against step 07 gives the duplication rate; step 06 against step 05, the share of poly(A) reads.

**Poly(A) read criterion** (SCINPAS `get_polyA`). After the soft-clip correction, the 3' soft clip of the read must have:
- ≥ {polyA_min_length} nt;
- ≥ {polyA_min_A_percent}% A;
- a median base quality ≥ {min_phred_score}.
"""),
    ("QC/non_polyA_alignments", lambda p: as_bool(p["merge_non_polyA"]),
     "aligned reads without a poly(A) tail at their 3' end; not used downstream",
     """
**Produced by:** process `merge_non_polyA_bams` (`samtools merge`).

| file | content |
|:---|:---|
| `<sample>.non_polyA.sorted.bam` (+ `.bai`) | coordinate-sorted, **not** deduplicated |

These are reads whose 3' soft clip does not pass the poly(A) criterion (see `../README.md`): for example truncated reads, reads from internal priming that align through genomic A stretches, or reads with adapter-only ends.

**Use:** inspection in IGV; checking what was lost.
"""),
    ("QC/skipped_pt_alignments", lambda p: as_bool(p["merge_skipped_pt"]),
     "poly(A) reads without a valid Dorado tail-length estimate; not used downstream",
     """
**Produced by:** process `merge_skipped_pt_bams` (`samtools merge`).

| file | content |
|:---|:---|
| `<sample>.skipped_pt.sorted.bam` (+ `.bai`) | coordinate-sorted, **not** deduplicated |

These are poly(A) reads (soft-clip criterion passed) for which Dorado gave no tail length (`pt` ≤ 0 or missing) or the cleavage-site tags were missing, so no poly(A) tail could be appended. Their number per chunk is in `../master_pt_append_stats.tsv`.

**Use:** inspection; checking that the dropped reads are not biased toward particular genes.
"""),
    ("analysis_figures", lambda p: as_bool(p["visualize_cleavage_site_motifs_flag"]) or as_bool(p["visualize_polyA_tail_length_distribution_flag"]),
     "overview figures (PDF)",
     """
| folder | content |
|:---|:---|
| `cleavage_site_motifs/<sample>/` | poly(A)-signal motif meta-plots around the cleavage sites (`visualize_cleavage_site_motifs_flag`) |
| `polyA_tail_lengths/` | poly(A) tail length distributions of all samples (`visualize_polyA_tail_length_distribution_flag`) |

See the README in each folder.
"""),
    ("analysis_figures/cleavage_site_motifs", lambda p: as_bool(p["visualize_cleavage_site_motifs_flag"]),
     "meta-plots of poly(A)-signal motifs around the cleavage sites, one folder per sample",
     """
**Produced by:** process `visualize_cleavage_site_motifs` (`plot_cs_motifs` from zavolab_pyutils), from the bigWigs in `../../cleavage_sites_bigwigs/`.

The cleavage sites of a sample are split into up to {motif_score_bins} quantile bins of read support. Ties can merge bins, so the files may say e.g. `bin1of5`; bin 1 has the lowest support. For each bin, the genomic window from {motif_window_up} nt upstream to {motif_window_down} nt downstream is scanned for `{motif_list}`, with motif positions anchored at the motif {motif_anchor} and sites weighted `{motif_weighting}`.

| file in `<sample>/` | content |
|:---|:---|
| `<sample>_motif_metaplot_bin<i>of<n>_counts.pdf` | share of sites with 0, 1, 2, ... occurrences of each motif in the window |
| `<sample>_motif_metaplot_bin<i>of<n>_pos_all_sites.pdf` | position of each motif relative to the cleavage site, as a share of all sites of the bin |
| `<sample>_motif_metaplot_bin<i>of<n>_pos_conditional.pdf` | the same, as a share of the sites that contain the motif |

These plots count sites, not reads, and do not consider relative usage. They serve as a rough check: true cleavage sites show AAUAAA / AUUAAA peaking about 20 nt upstream, and more so in the better-supported bins.
"""),
    ("analysis_figures/polyA_tail_lengths", lambda p: as_bool(p["visualize_polyA_tail_length_distribution_flag"]),
     "poly(A) tail length distributions of all samples (PDF)",
     """
**Produced by:** process `visualize_polyA_tail_length_distribution` (`bin/plot_polya_distributions.py`), from `../../read_tags/`.

| file | content |
|:---|:---|
| `polya_distribution_all_reads_pooled.pdf`, `..._all_reads_per_chr.pdf` | `pt` of all reads per sample, pooled or per chromosome |
| `polya_distribution_per_gene_median_pooled.pdf`, `..._per_chr.pdf` | per-gene median `pt` (weighted by 1/NH) |
| `polya_distribution_per_gene_mean_pooled.pdf`, `..._per_chr.pdf` | per-gene mean `pt` |

Only reads with `pt` > 0 are used, on the 30 chromosomes with the most reads.
"""),
    ("subsampled_pod5", lambda p: True,
     "raw signal (POD5) of the few reads chosen for the signal plots",
     """
**Produced by:** process `filter_pod5_combined` (`pod5 filter`).

| file | content |
|:---|:---|
| `<sample>.subset.pod5` | the raw current signal of the selected reads of one sample |

Read selection (process `extract_read_ids`), among the gene-assigned reads of `../read_tags/`: {read_selection}. A sample with no gene-assigned reads has no file.

**Use:** input of `../moves_bam/` and `../annotated_data/`; can be opened with the `pod5` tools.
"""),
    ("moves_bam", lambda p: True,
     "the selected reads basecalled again, with move tables",
     """
**Produced by:** process `dorado_emit_moves`: `dorado basecaller --emit-moves --estimate-poly-a`, aligned by Dorado (`--mm2-opts "-x splice -Y"`), on `../subsampled_pod5/`.

| file | content |
|:---|:---|
| `<sample>.subset.moves.bam` | unsorted, not indexed; tag `mv` (move table), `pt`/`pa` (poly(A) length and signal anchors), `ts`/`ns` (trimmed and total signal samples) |

Reads here are **not** oriented or processed as in the main pipeline: this BAM only links bases to the raw signal.

**Use:** input of `../annotated_data/`.
"""),
    ("annotated_data", lambda p: True,
     "raw current signal of each selected read with the basecalled bases and the poly(A) region, one CSV per read",
     """
**Produced by:** process `generate_signal_df` (`bin/pod5_to_df.py`), from `../subsampled_pod5/` and `../moves_bam/`.

| file | content |
|:---|:---|
| `<sample>_<read_id>_mapped.csv` | one row per signal sample |

Columns:
- `time`: seconds from the start of the read;
- `signal`: raw current, in uncalibrated digitiser units;
- `polyA`: Dorado's tail length `pt` inside the poly(A) segment of the signal, −1 elsewhere;
- `ann`: move-table value (1 = a new base starts, 0 = the same base continues, −1 / −2 = trimmed adapter or signal outside the called region);
- `base`: the base the signal sample belongs to, `N` outside the called region.

**Use:** input of the plots in `../QC_plots/raw_signal_annotation/`; custom plots of the signal.
"""),
    ("QC_plots", lambda p: True,
     "QC plots",
     """
| folder | content |
|:---|:---|
| `raw_signal_annotation/` | raw current of single reads with bases and poly(A) segment marked |
"""),
    ("QC_plots/raw_signal_annotation", lambda p: True,
     "raw current of single selected reads, with the called bases and the poly(A) segment marked (PDF per read)",
     """
**Produced by:** process `visualize_signal` (`bin/plot_signal.py`), from `../../annotated_data/`.

| path | content |
|:---|:---|
| `<sample>/<sample>.subset/<sample>_<read_id>.pdf` | one read; the title gives the read ID, gene (`XT`) and poly(A) length (`pt`) from `../../assignments/` |

**Use:** visual check of basecalling and of the poly(A) tail estimate, e.g. whether the marked poly(A) segment matches the flat low-variance stretch of the signal and whether `pt` is plausible.
"""),
]

TAGS = [
    ("pt", "Dorado", "estimated poly(A) tail length (nt); −1 if none was found"),
    ("pa", "Dorado", "poly(A) anchor positions in the raw signal"),
    ("RX", "Dorado", "UMI sequence; padded or trimmed to {umi_target_len} nt for `umi_tools`"),
    ("qs, du, ns, ts, ch, st, rn, fn, ...", "Dorado", "mean base quality, duration, signal samples, trimmed samples, channel, start time, read number, source POD5 file, ..."),
    ("{tag_original_ts}", "`orient_reads.py`", "read strand from Dorado (`{tag_basecaller_ts}`); reads with `-` were reverse-complemented, so all reads are in the mRNA orientation"),
    ("{tag_orig_cs}", "SCINPAS `fix_softclipped`", "original cleavage site: 3' end of the alignment, 1-based"),
    ("{tag_fixed_cs}", "SCINPAS `fix_softclipped`", "corrected cleavage site, 1-based: the 3' end extended over soft-clipped bases that match the genome{shift_note}"),
    ("{tag_phred_mapped}, {tag_phred_softclipped}", "SCINPAS `get_polyA`", "median base quality of the aligned part and of the 3' soft clip"),
    ("NH", "`redefine_qual_and_NHtag`", "number of alignments of the read after deduplication (MAPQ 255 if 1, else 0)"),
    ("XT", "featureCounts", "gene ID(s) of the reference GTF; comma-separated if several; missing if unassigned"),
    ("XS", "featureCounts", "assignment status (`Assigned`, `Unassigned_NoFeatures`, `Unassigned_Ambiguity`, ...)"),
    ("{tag_motif_info}", "`annotate_motifs_in_bam.py`", "poly(A)-signal motifs near `{tag_quantification_cs}`, e.g. `AAUAAA:(-23,13)` (positions relative to the cleavage site)"),
    ("{tag_transcript_id}", "`tag_bam_with_transcript.py`", "isoform from `transcriptome/master_enriched.tsv` (reference or `novel_isoform_<n>`)"),
    ("ts (Z), tp, cm, s1, s2, de, rl, NM, AS, ms", "minimap2", "standard minimap2 tags (`ts` with value `+`/`-` is minimap2's transcript strand, not Dorado's `ts` (trimmed samples, integer))"),
]


def read_params(files):
    p = dict(DEFAULTS)
    for f in files:
        with open(f) as fh:
            p.update({k: v for k, v in json.load(fh).items() if v is not None})
    if not p.get("outdir"):
        p["outdir"] = os.path.join(str(p["rundir"]), "results")
    p.update(FIXED)
    multi = as_bool(p["include_multimappers"])
    p["multimapper_rows"] = ("included, one row per alignment; weight them by 1/NH" if multi
                             else "left out (only reads with NH = 1)")
    p["multimapper_included"] = "included" if multi else "left out"
    genes = str(p.get("viz_gene_list") or "").strip()
    p["read_selection"] = (f"up to {p['num_reads']} random reads of each of the genes `{genes}`" if genes else
                           f"{p['num_reads']} random genes, and up to {p['num_reads']} random reads of each")
    p["shift_note"] = ("; ambiguous sites in A/T stretches are shifted to the stretch boundary"
                       if as_bool(p["shift_ambiguous_cs"]) else "")
    return p


def sample_table(p, annotation_file):
    """Samples of the sample sheet with their number of POD5 files, joined with an optional annotation table."""
    counts = OrderedDict()
    tsv = p.get("tsv")
    if tsv and os.path.exists(str(tsv)):
        with open(tsv) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                counts[row["sample_id"]] = counts.get(row["sample_id"], 0) + 1
    ann, ann_cols = {}, []
    if annotation_file:
        with open(annotation_file) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            ann_cols = [c for c in reader.fieldnames if c != "sample_id"]
            for row in reader:
                ann[row["sample_id"]] = row
                counts.setdefault(row["sample_id"], "")
    if not counts:
        return "*(sample sheet not found: `{}`)*".format(tsv)
    natural = lambda s: [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]
    rows = [[s, counts[s]] + [ann.get(s, {}).get(c, "") for c in ann_cols] for s in sorted(counts, key=natural)]
    return md_table(["sample_id", "POD5 files"] + ann_cols, rows)


def write(path, text, overwrite):
    if os.path.exists(path) and not overwrite:
        print(f"exists, not overwritten: {path}")
        return
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text.strip("\n") + "\n")
    print(f"written: {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", required=True, help="results directory to write the READMEs into")
    parser.add_argument("--params", nargs="*", default=[], help="JSON files with run parameters (later files override earlier ones)")
    parser.add_argument("--rundir_readme", default=None, help="also write a README for the run directory to this path")
    parser.add_argument("--sample_annotation", default=None, help="TSV with column sample_id and any annotation columns, shown in the sample table")
    parser.add_argument("--notes", default=None, help="Markdown file with run-specific notes, added to the results README")
    parser.add_argument("--nanoflowz_version", default="", help="e.g. git commit of nanoflowz")
    parser.add_argument("--run_name", default="", help="Nextflow run name")
    parser.add_argument("--only_existing", action="store_true", help="write READMEs only into folders that exist (for finished runs)")
    parser.add_argument("--no_overwrite", action="store_true", help="keep existing README.md files")
    args = parser.parse_args()

    p = read_params(args.params)
    fill = Fill({k: v for k, v in p.items()})
    overwrite = not args.no_overwrite

    written = []
    for rel, cond, summary, body in FOLDERS:
        folder = os.path.join(args.outdir, rel)
        if args.only_existing:
            if not os.path.isdir(folder):
                continue
        elif not cond(p):
            continue
        os.makedirs(folder, exist_ok=True)
        text = f"# `{rel}/`\n\n{summary.format_map(fill)}\n" + body.format_map(fill)
        text += f"\n\n---\n*Written by nanoflowz `bin/write_output_readmes.py`; overview of all outputs: `{os.path.relpath(args.outdir, folder)}/README.md`.*\n"
        write(os.path.join(folder, "README.md"), text, overwrite)
        written.append((rel, summary.format_map(fill)))

    # results README
    run_rows = [("run directory", f"`{p['rundir']}`"), ("results directory", f"`{p['outdir']}`"),
                ("nanoflowz version", args.nanoflowz_version or "not recorded"),
                ("Nextflow run", args.run_name or "not recorded"),
                ("README written", datetime.date.today().isoformat()),
                ("sample sheet", f"`{p['tsv']}`"), ("Dorado", f"`{p['dorado']}`"), ("basecalling model", f"`{p['model']}`"),
                ("poly(A) config", f"`{p['polyA']}`"), ("genome", f"`{p['ref']}`"), ("annotation", f"`{p['reference_gtf']}`")]
    folder_rows = [(f"[`{rel}/`]({rel}/README.md)", s) for rel, s in written]
    notes = ""
    if args.notes:
        with open(args.notes) as fh:
            notes = "\n## Notes on this run\n\n" + fh.read().strip("\n") + "\n"
    top = f"""
# nanoflowz results

Results of the nanoflowz pipeline (Oxford Nanopore cDNA with UMIs: basecalling, poly(A) reads, UMI deduplication, gene, cleavage-site, isoform and poly(A)-length quantification). Every folder has its own README with the files, their formats and how they were made.

{md_table(["", ""], run_rows)}
{notes}
## Samples

{sample_table(p, args.sample_annotation)}

## Where is what

{md_table(["folder", "content"], folder_rows)}

**What to use for what:**

| question | file |
|:---|:---|
| gene expression | `gene_assignments/<sample>.featureCounts.txt` |
| cleavage-site (3' end) usage | `cleavage_sites_bigwigs/<sample>.cs.bed.gz` (table), `*.plus/minus.bigwig` (IGV) |
| per-read data (gene, isoform, cleavage site, poly(A) length, motifs) | `read_tags/<sample>.tags.tsv.gz` |
| poly(A) tail length | column `pt` of `read_tags/`; `polya_length_bigwigs/` per cleavage site |
| isoforms | tag `{p['tag_transcript_id']}` in `assignments/` and `read_tags/`; catalogue in `transcriptome/` |
| alignments for IGV | `assignments/<sample>.transcript_assigned.bam` |
| read numbers through the pipeline | `QC/master_read_tracking_stats.tsv` |

**The BAM folders are successive stages of the same reads.**
1. `basecalling/`: poly(A) reads, before deduplication.
2. `map_genome_merged_UMIdedup/`: deduplicated.
3. `gene_assignments/`: + `XT`, `XS`.
4. `motif_annotated_bams/`: + `{p['tag_motif_info']}`.
5. `assignments/`: + `{p['tag_transcript_id']}`.

From step 2 on, the reads are the same and only tags are added, so `assignments/` alone holds everything.

## Processing steps

1. `dorado basecaller` (GPU) on every POD5 file, with poly(A) tail estimation (`pt`).
2. Reads reverse-complemented to the mRNA orientation, using Dorado's strand tag.
3. `minimap2 -ax splice:hq --secondary=no -Y` to the genome; supplementary alignments removed.
4. UMIs (`{p['umi_tag']}`) normalised to {p['umi_target_len']} nt.
5. SCINPAS: cleavage site = 3' end of the alignment, corrected for soft-clipped bases (`{p['tag_orig_cs']}` → `{p['tag_fixed_cs']}`). Unmapped reads are removed.
6. SCINPAS: poly(A) reads kept. Their 3' soft clip must have ≥ {p['polyA_min_length']} nt with ≥ {p['polyA_min_A_percent']}% A. Others go to `QC/non_polyA_alignments/`.
7. The 3' soft clip is replaced by a poly(A) tail of Dorado's length `pt`. Reads without a valid `pt` go to `QC/skipped_pt_alignments/`.
8. Chunks merged per sample (`basecalling/`); UMI deduplication with `umi_tools`; NH and MAPQ recomputed (`map_genome_merged_UMIdedup/`).
9. Cleavage-site counts (`cleavage_sites_bigwigs/`); featureCounts gene assignment (`gene_assignments/`); motifs (`motif_annotated_bams/`).
10. Enriched transcript annotation from all samples (`transcriptome/`); isoform of every read (`assignments/`); tag tables (`read_tags/`); poly(A)-length bigWigs (`polya_length_bigwigs/`).
11. QC: read numbers per step (`QC/`); raw-signal plots of a few reads (`subsampled_pod5/`, `moves_bam/`, `annotated_data/`, `QC_plots/`); optional figures (`analysis_figures/`).

## BAM tags

{md_table(["tag", "set by", "meaning"], [(t.format_map(fill), s, m.format_map(fill)) for t, s, m in TAGS])}

Coordinates in tags are 1-based genomic positions. Reads are in the mRNA orientation, so the strand of an alignment is the strand of the RNA.

---
*Written by nanoflowz `bin/write_output_readmes.py`.*
"""
    write(os.path.join(args.outdir, "README.md"), top, overwrite)

    if args.rundir_readme:
        outdir_rel = os.path.relpath(str(p["outdir"]), str(p["rundir"]))
        rundir = f"""
# nanoflowz run

| path | content |
|:---|:---|
| [`{outdir_rel}/`]({outdir_rel}/README.md) | **the results**: see its README |
| `work/` | Nextflow working directories, one per task: intermediate files, `.command.sh` (the exact command), `.command.log` / `.command.err`. Not needed for the results; can be deleted after the run, but then `-resume` restarts from scratch |
| `conda_envs/` | the conda environments of the run |
| `run_params.json` | the run parameters (`-params-file`); unset parameters take the defaults of nanoflowz `nextflow.config` |
| sample sheet | `{p['tsv']}`: `sample_id` and the POD5 file(s) of each sample |
| poly(A) config | `{p['polyA']}`: Dorado poly(A) estimation settings |

The Nextflow log (`.nextflow.log`) and the run history are in the directory from which `nextflow run` was started.

---
*Written by nanoflowz `bin/write_output_readmes.py`.*
"""
        write(args.rundir_readme, rundir, overwrite)


if __name__ == "__main__":
    main()
