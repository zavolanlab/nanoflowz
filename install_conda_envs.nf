#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    NANOFLOWZ CONDA ENVIRONMENT INSTALLER
    Pre-builds every conda environment from envs/*.yml into params.conda_envs_dir.
    Run it with the same -params-file as main.nf on a node WITH internet access
    (login / vscode node) before launching main.nf on a compute node without it.
    Nextflow names env directories by hashing the .yml content, so the envs built
    here are exactly the ones main.nf will find and reuse.
========================================================================================
*/

workflow {
    channel.fromPath("${projectDir}/envs/*.yml")
        | build_conda_env
        | view { it.trim() }
}

process build_conda_env {
    tag "${yml.baseName}"
    // always run on the current (internet-connected) node, never submit to SLURM
    executor 'local'
    conda { yml.toString() }

    input:
    val yml

    output:
    stdout

    script:
    """
    echo "[ready] ${yml.name} -> \$CONDA_PREFIX"
    """
}
