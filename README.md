# Nanoflowz: Nanopore Signal Processing Pipeline

Project homepage: [Nanoflowz on GitHub](https://github.com/zavolanlab/nanoflowz)

Nanoflowz is a Nextflow DSL2 pipeline for processing Oxford Nanopore sequencing data from POD5 files to:

- Dorado-basecalled and aligned reads  
- Subsampled POD5 signal  
- Annotated signal-level CSVs  
- Per-read PDF signal plots

The main workflow is defined in `main.nf`.

---

## Overview

Given a TSV file describing samples and their POD5 input files, the pipeline:

1. **Basecalls reads with Dorado**
   - Aligns reads to a reference with minimap2 via Dorado
   - Sorts and indexes the resulting BAM

2. **Extracts random read IDs**
   - Samples a subset of read IDs from the basecalled BAM using `params.num_reads`

3. **Subsamples POD5 files**
   - Filters the original POD5 to only include the sampled read IDs

4. **Re-runs Dorado with moves**
   - Generates a BAM with movement (`--emit-moves`) information on the subsampled POD5

5. **Generates annotated signal data**
   - Runs `pod5_to_df.py` to produce CSVs with signal-level annotations

6. **Visualizes signal**
   - Runs `plot_signal.py` to create per-read PDF plots

---

## Input TSV format

The pipeline expects a tab-delimited file with at least the following columns:

- `sample_id` – unique identifier for the sample  
- `pod5` – path to the input POD5 file

Example:

```text
sample_id    pod5
sampleA      /path/to/sampleA.pod5
sampleB      /path/to/sampleB.pod5
```

You must pass this file with `--tsv`.

---

## Key processes

- **`dorado_basecall`**: Dorado basecalling + alignment, outputs `${sample_id}.dorado.sup.sorted.bam` and `.bai`.
- **`extract_read_ids`**: Extracts unique read IDs from BAM and randomly selects `params.num_reads`.
- **`filter_pod5`**: Creates `${sample_id}.subset.pod5` containing only selected reads.
- **`dorado_emit_moves`**: Runs Dorado with `--emit-moves` on the subset POD5, outputting `${sample_id}.moves.bam`.
- **`generate_signal_df`**: Calls `pod5_to_df.py` to create annotated signal-level CSVs.
- **`visualize_signal`**: Calls `plot_signal.py` to generate per-read PDF plots.

---

## Requirements

- **Nextflow** (DSL2 enabled)
- **Dorado** (GPU-capable, CUDA)
- **samtools**
- **pod5** CLI tools
- **Python** environment with dependencies for:
  - `pod5_to_df.py`
  - `plot_signal.py` (e.g. `pandas`, `matplotlib`, etc.)
- Access to a **reference genome** FASTA (for mapping)

---

## Parameters (commonly used)

Configured via CLI or `nextflow.config`:

- `--tsv` – path to the input TSV file **(required)**
- `--outdir` – output directory for results (e.g. `results/`)
- `--dorado` – path or name of the Dorado executable
- `--model` – Dorado model name/path
- `--polyA` – poly-A configuration file for Dorado
- `--ref` – path to reference genome FASTA
- `--num_reads` – number of read IDs to sample per sample
- `--figwidth` – plot width for `plot_signal.py`
- `--figheight` – plot height for `plot_signal.py`

See `main.nf` for the authoritative parameter list.

---

## Running the pipeline

Basic example (parameters from `nextflow.config`):

```bash
nextflow run main.nf --tsv samples.tsv
```

All other parameters (e.g. `outdir`, `dorado`, `model`, `polyA`, `ref`, `num_reads`, plotting settings) should be set in your `nextflow.config`.

---

## Outputs

Under `--outdir` you should see:

- `basecalling/`
  - `${sample_id}.dorado.sup.sorted.bam`
  - `${sample_id}.dorado.sup.sorted.bam.bai`
- `subsampled_pod5/`
  - `${sample_id}.subset.pod5`
- `moves_bam/`
  - `${sample_id}.moves.bam`
- `annotated_data/`
  - `*.csv` annotated signal tables
- `plots/`
  - `*.pdf` signal plots per read
