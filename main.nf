#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Check if TSV is provided
if( !params.tsv ) { exit 1, "Please provide the input TSV file with --tsv" }

workflow {
    samples_ch = channel
        .fromPath(params.tsv)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_id, file(row.pod5)) }

    // 1. Initial basecall
    dorado_basecall(samples_ch)

    // 2. Extract random IDs
    extract_read_ids(dorado_basecall.out.bam)

    // 3. Create subset POD5
    filter_input_ch = samples_ch.join(extract_read_ids.out.ids_file)
    
    filter_pod5(filter_input_ch)

    // 4. Re-run Dorado for moves (Needs the POD5 subset)
    dorado_emit_moves(filter_pod5.out)

    // Join POD5 and BAM before generating CSV
    // filter_pod5.out is [sample_id, pod5]
    // dorado_emit_moves.out.bam is [sample_id, bam]
    final_input_ch = filter_pod5.out.join(dorado_emit_moves.out.bam)
    
    generate_signal_df(final_input_ch)

    visualize_signal(generate_signal_df.out.flatten())
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

process extract_read_ids {
    tag "${sample_id}"
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_read_ids.txt"), emit: ids_file

    script:
    """
    samtools view ${bam} | cut -f1 | sort -u | shuf -n ${params.num_reads} > ${sample_id}_read_ids.txt
    """
}

process filter_pod5 {
    tag "${sample_id}"
    publishDir "${params.outdir}/subsampled_pod5", mode: 'copy'

    input:
    // This matches the output of the .join()
    tuple val(sample_id), path(original_pod5), path(read_ids_txt)

    output:
    tuple val(sample_id), path("${sample_id}.subset.pod5")

    script:
    """
    pod5 filter ${original_pod5} --output ${sample_id}.subset.pod5 --ids ${read_ids_txt} --force-overwrite
    """
}

process dorado_emit_moves {
    tag "${sample_id}"
    publishDir "${params.outdir}/moves_bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(subset_pod5)

    output:
    tuple val(sample_id), path("${sample_id}.moves.bam"), emit: bam

    script:
    """
    # Pipe Dorado output to samtools to ensure a valid, compressed BAM with header
    ${params.dorado} basecaller \\
        ${params.model} \\
        ${subset_pod5} \\
        --emit-moves \\
        --estimate-poly-a \\
        --poly-a-config ${params.polyA} \\
        --mm2-opts "-x splice -Y" \\
        --reference ${params.ref} \\
        --device "cuda:all" | samtools view -bS - > ${sample_id}.moves.bam
    """
}

process generate_signal_df {
    tag "${sample_id}"
    publishDir "${params.outdir}/annotated_data", mode: 'copy'

    input:
    tuple val(sample_id), path(subset_pod5), path(moves_bam)

    output:
    path "*.csv"

    script:
    """
    pod5_to_df.py --pod5 ${subset_pod5} --bam ${moves_bam} --sample_id ${sample_id}
    """
}

process visualize_signal {
    tag "${csv.baseName}"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path csv

    output:
    path "*.pdf"

    script:
    def read_id = csv.baseName.replace("_mapped", "")
    """
    plot_signal.py \\
        --csv ${csv} \\
        --output ${read_id}.pdf \\
        --title "Read: ${read_id}" \\
        --figwidth ${params.figwidth} \\
        --figheight ${params.figheight}
    """
}

