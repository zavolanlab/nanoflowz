# Nanoflowz

Automated Dorado basecalling.

## Setup
1. Ensure your `samples.tsv` is tab-separated with headers `sample_id` and `pod5`.

## Running the Pipeline
To start a new run:
```bash
nextflow run main.nf --tsv samples.tsv
```