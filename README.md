# Nanoflowz
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/) [![Conda](https://img.shields.io/badge/conda-supported-lightgrey.svg)](https://docs.conda.io/en/latest/)

**Nanoflowz** is a bioinformatics analysis pipeline for Oxford Nanopore Technologies (ONT) RNA and cDNA sequencing data. 
It is tailored for GPU-accelerated basecalling, producing quality control plots, _de novo_ transcriptome assembly and assigning alignments to transcripts. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Conda environments to automatically manage and isolate software dependencies.

## Pipeline Summary
1. **GPU Basecalling (`dorado`)**: Performs basecalling and poly-A tail estimation natively on NVIDIA GPUs.
2. **Alignment & Merging (`minimap2`)**: Maps reads to the reference genome using minimap2 and merges BAMs (e.g. technical replicates from samples).
3. **De novo transcriptome assembly**: Enriches input transcriptome annotations based on observed read alignments
4. **Isoform Analysis (`custom python`)**: Assigns alignments to specific transcript isoforms.
5. **Annotating raw current data from a subsample of reads**: 
5.1 **Signal Extraction (`pod5`)**: Subsamples reads of interest and extracts corresponding raw signal chunks from `.pod5` files.
5.2 **Move-table Emission (`dorado`)**: Re-processes subsetted reads to emit basecaller move-tables. 
5.3 **Dataframe Generation & QC Visualization (`seaborn`)**: Synchronizes sequence strings with raw signal variations and renders PDF plots for individual reads.

## Quick Start

### 1. Prepare your Input Data (`samples.tsv`)
You will need to create a tab-separated file (referred to as `tsv` below in step 4) containing your sample IDs and the **absolute paths** to your raw `.pod5` files. 
Each `sample_id` may correspond to multiple `.pod5` files which will be basecalled and aligned independently and then merged to `sample_id`-level `.bam` files.
*(Note: Do not use relative symlinks).*

```tsv
sample_id	pod5
barcode01	/absolute/path/to/data/barcode01/file_1.pod5
barcode01	/absolute/path/to/data/barcode01/file_2.pod5
barcode02	/absolute/path/to/data/barcode02/file_1.pod5
barcode02	/absolute/path/to/data/barcode02/file_2.pod5
```

### 2. Install Dorado
Currently, the only officially supported way is to install Dorado from their [github repository](https://github.com/nanoporetech/dorado). 
Follow the instructions there.
Get the path to `dorado` executor (e.g. `$HOME/packages/dorado-1.4.0-linux-x64/bin/dorado`).

### 3. Download basecalling model(s) for Dorado
The core of the Dorado is a neural network that transforms original raw current data (current vs time) for every read into a nucleotide read sequence.

For that, ONT trained several neural network models optimized for different protocols and molecule types (DNA/cDNA or RNA).

Currently, the list of available models is available [here](https://software-docs.nanoporetech.com/dorado/latest/models/list/).

By default, we recommend to use "super-accurate" (`_sup`) models which at least in particular circumstances can give drastically more accurate read sequences as an output (internal research, not published). This is in contrast to the [current recommendation of ONT](https://software-docs.nanoporetech.com/dorado/latest/models/models/#understanding-model-names) to stick to "high-accuracy" `_hac` models by default.

### 4. Define your Parameters (`run_params.json`)
Instead of modifying the core `nextflow.config` file, Nanoflowz accepts a JSON file containing your run-specific inputs, reference genomes, and basecaller model paths. Here is the example of the specification for a human cell line data:

```json
{
    "tsv": "/path/to/samples.tsv",
    "rundir": "/path/to/run_output_directory",
    "dorado": "/path/to/dorado-1.4.0-linux-x64/bin/dorado",
    "model": "/path/to/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0",
    "polyA": "/path/to/polyA_config.toml",
    "ref": "/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    "reference_gtf": "/path/to/gencode.v42.annotation.gtf",
    "conda_envs_dir": "/path/to/where_nanoflowz_should_create_conda_envs", # default is rundir/conda_envs
}
```

Output directory can be specified by setting `outdir` parameter, by default set to `${params.rundir}/results`.

Nextflow working directory with all the intermediate files can be specified by setting `workdir_param` parameter, by default set to `${params.rundir}/work`.

### 5. Run the Pipeline
Execute the pipeline using the `-profile conda` flag. This ensures Nextflow handles all Python/Samtools dependencies automatically by automatically creating conda environments from `.yml` files stored in [envs](envs/) subdirectory.

Look [here](https://docs.seqera.io/nextflow/cli) to see various ways to execute nextflow from CLI.

Currently, we recommend to stick to **standard Local Execution** method. Clone the reposity first to your local machine:
```bash
git clone https://github.com/zavolanlab/nanoflowz.git
cd nanoflowz
echo "$(pwd)/main.nf" # this will print you the absolute path to main nanoflowz executor script
```
Then you can use the `run_params.json` file created at step 4 as an argument:
```bash
nextflow run <put the printed path to main.nf here> -params-file <put the path to your created run_params.json file> -profile conda -resume
```

## Configuring in Python / Jupyter notebook
Nanoflowz is designed to be easily wrapped by Python scripts or Jupyter Notebooks. You can dynamically generate the required parameter files and trigger the pipeline programmatically. See the example in the jupyter notebook [ONT_analysis.ipynb](https://github.com/zavolanlab/APA_localization/blob/9ddaf03675811d69191fa9e6aa0efb2211728d58/ONT_analysis.ipynb) from [another repository](https://github.com/zavolanlab/APA_localization/tree/ont_analysis) of Zavolan Lab.

```python
import json
import subprocess
from pathlib import Path

# 1. Define your parameters dictionary
json_params_content = {
    "tsv": str(Path('samples.tsv').resolve()),
    "rundir": str(Path('./wf_runs/run_01').resolve()),
    "dorado": "/absolute/path/to/dorado",
    "model": "/absolute/path/to/model",
    "polyA": "/absolute/path/to/polyA_config.toml",
    "ref": "/absolute/path/to/genome.fa",
    "reference_gtf": "/absolute/path/to/annotation.gtf"
}

# 2. Write the JSON file
params_file = Path('run_params.json')
with open(params_file, "w") as f:
    json.dump(json_params_content, f, indent=4)

# 3. Trigger Nextflow
cmd = [
    "nextflow", "run", "main.nf", # assuming we are in the nanoflowz directory
    "-params-file", str(params_file.resolve()),
    "-profile", "conda"
]
subprocess.run(cmd, check=True)
```

## Output Structure
Nanoflowz uses a "Sibling Architecture" inside your designated `rundir` to keep your final biological results perfectly separated from heavy temporary files.

```tree
rundir/
├── results/              # ✨ Pristine final outputs (BAMs, CSVs, PDFs)
│   ├── basecalling/
│   ├── transcriptome/
│   ├── QC_plots/
│   └── ...
├── work/                 # 🗑️ Nextflow temporary execution directories (Safe to delete post-run)
└── conda_envs/           # 📦 Isolated software environments for this specific run
```

## Credits
Nanoflowz was originally written by the Zavolan Lab.

**CURRENT STATUS**: please cite this github repository if you use **nanoflowz** in your research.