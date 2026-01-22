#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Check if TSV is provided
if( !params.tsv ) { exit 1, "Please provide the input TSV file with --tsv" }

workflow {

    // Using lowercase 'channel' avoids the ConfigObject error
    samples_ch = channel
        .fromPath(params.tsv)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            if (!row.sample_id || !row.pod5) {
                error "TSV missing required columns 'sample_id' or 'pod5'"
            }
            tuple(row.sample_id, file(row.pod5)) 
        }

    dorado_basecall(samples_ch)
}

process dorado_basecall {
    tag "${sample_id}"
    
    // This moves the final files to your results folder
    publishDir "${params.outdir}/basecalling", mode: 'copy'

    input:
    tuple val(sample_id), path(pod5_input)

    output:
    tuple val(sample_id), path("${sample_id}.dorado.sup.sorted.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.dorado.sup.sorted.bam.bai"), emit: bai

    script:
    """
    OUT_BAM="${sample_id}.dorado.sup.sorted.bam"

    # Run Dorado and pipe directly to Samtools for sorting
    ${params.dorado} basecaller \\
        ${params.model} \\
        ${pod5_input} \\
        --estimate-poly-a \\
        --poly-a-config ${params.polyA} \\
        --mm2-opts "-x splice -Y" \\
        --reference ${params.ref} \\
        --device "cuda:all" \\
        | samtools sort -@ ${task.cpus} -o \$OUT_BAM -
    
    # Index the resulting BAM
    samtools index -@ ${task.cpus} \$OUT_BAM
    """
}