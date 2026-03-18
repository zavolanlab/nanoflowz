#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    NANAFLOWZ MAIN PIPELINE
    Optimized for consolidated subsampling, signal visualization, and isoform analysis.
========================================================================================
*/

// Check if TSV is provided
if( !params.tsv ) { exit 1, "Please provide the input TSV file with --tsv" }
if( !params.reference_gtf ) { exit 1, "Please provide the reference GTF file with --reference_gtf" }

// Default base directory for QC plots
params.qc_base_dir = "${params.outdir}/QC_plots/raw_signal_annotation"

workflow {
    // 1. Initial Data Ingest: Create a stream of (sample_id, pod5_file)
    samples_ch = channel
        .fromPath(params.tsv)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_id, file(row.pod5)) }

    // 2. Initial basecall: Parallel processing per POD5 chunk
    dorado_basecall(samples_ch)

    // 3. Merge: Combine small BAMs into a single sample-level BAM
    merge_input_ch = dorado_basecall.out.bam.groupTuple()
    merge_bams(merge_input_ch)

    // ==============================================================================
    // ISOFORM ANALYSIS
    // ==============================================================================

    // A. Collect all merged BAM files to build the master enriched transcriptome
    all_bams_ch = merge_bams.out.bam.map { it[1] }.collect()
    
    transcriptome_annotation_enrichment(
        all_bams_ch, 
        params.reference_gtf
    )

    // B. Assign individual alignments for each sample to the enriched transcriptome
    read_to_transcript_assignment(
        merge_bams.out.bam, 
        transcriptome_annotation_enrichment.out.gtf
    )

    // ==============================================================================

    // 4. Selection: Pick the random Read IDs to investigate
    extract_read_ids(merge_bams.out.bam)

    // 5. We take all original POD5 paths and group them by sample_id
    // This allows one process to search all files at once.
    all_pod5s_per_sample = samples_ch.map { id, pod5 -> [id, pod5] }.groupTuple()
    
    filter_input_ch = extract_read_ids.out.ids_file.join(all_pod5s_per_sample)
    filter_pod5_combined(filter_input_ch)

    // 6. Deep Dive: Re-run Dorado for moves on the subsampled POD5
    dorado_emit_moves(filter_pod5_combined.out.pod5)

    // 7. Prepare Data: Join the subsampled POD5 and its Move-BAM to generate CSVs
    final_input_ch = filter_pod5_combined.out.pod5.join(dorado_emit_moves.out.bam)
    generate_signal_df(final_input_ch)

    // 8. Visualize: Plot the squiggles for each individual read
    visualize_input = generate_signal_df.out.results.transpose()
    visualize_signal(visualize_input)
}

/*
========================================================================================
    PROCESS DEFINITIONS
========================================================================================
*/

process dorado_basecall {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(pod5_input)

    output:
    tuple val(sample_id), path("chunk_${pod5_input.baseName}.bam"), emit: bam

    script:
    """
    OUT_BAM="chunk_${pod5_input.baseName}.bam"
    ${params.dorado} basecaller \\
        ${params.model} \\
        ${pod5_input} \\
        --estimate-poly-a \\
        --poly-a-config ${params.polyA} \\
        --mm2-opts "-x splice -Y" \\
        --reference ${params.ref} \\
        --device "cuda:all" \\
        | samtools sort -@ ${task.cpus} -o \$OUT_BAM -
    """
}

process merge_bams {
    tag "${sample_id}"
    publishDir "${params.outdir}/basecalling", mode: 'copy'

    input:
    tuple val(sample_id), path(bams)

    output:
    tuple val(sample_id), path("${sample_id}.dorado.sup.sorted.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.dorado.sup.sorted.bam.bai"), emit: bai

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_id}.dorado.sup.sorted.bam ${bams}
    samtools index -@ ${task.cpus} ${sample_id}.dorado.sup.sorted.bam
    """
}

process transcriptome_annotation_enrichment {
    label 'process_high'
    publishDir "${params.outdir}/transcriptome", mode: 'copy'

    input:
    path bams
    path reference_gtf

    output:
    path "master_enriched.tsv", emit: gtf

    script:
    """
    # Create a manifest file containing the paths of all staged BAM files
    ls *.bam > bam_list.txt
    
    assign_ONT_reads_to_isoforms.py \\
        --mode enrich \\
        --input_bam_files bam_list.txt \\
        --input_gtf_file ${reference_gtf} \\
        --output_prefix master
    """
}

process read_to_transcript_assignment {
    tag "${sample_id}"
    publishDir "${params.outdir}/assignments", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path enriched_gtf

    output:
    tuple val(sample_id), path("${sample_id}_assignments.tsv.gz"), emit: tsv
    path "${sample_id}_unassigned.tsv", optional: true

    script:
    """
    assign_ONT_reads_to_isoforms.py \\
        --mode assign \\
        --input_bam_files ${bam} \\
        --input_gtf_file ${enriched_gtf} \\
        --output_prefix ${sample_id}
    
    gzip ${sample_id}_assignments.tsv
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

process filter_pod5_combined {
    tag "${sample_id}"
    publishDir "${params.outdir}/subsampled_pod5", mode: 'copy'

    input:
    tuple val(sample_id), path(read_ids_txt), path(all_pod5_chunks)

    output:
    tuple val(sample_id), path("${sample_id}.subset.pod5"), emit: pod5

    script:
    """
    # Scans all chunks in one pass to find the target IDs
    pod5 filter ${all_pod5_chunks} \\
        --output ${sample_id}.subset.pod5 \\
        --ids ${read_ids_txt} \\
        --missing-ok \\
        --force-overwrite
    """
}

process dorado_emit_moves {
    tag "${sample_id}"
    publishDir "${params.outdir}/moves_bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(subset_pod5)

    output:
    tuple val(sample_id), path("${subset_pod5.baseName}.moves.bam"), emit: bam

    script:
    """
    ${params.dorado} basecaller \\
        ${params.model} \\
        ${subset_pod5} \\
        --emit-moves \\
        --estimate-poly-a \\
        --poly-a-config ${params.polyA} \\
        --mm2-opts "-x splice -Y" \\
        --reference ${params.ref} \\
        --device "cuda:all" | samtools view -bS - > ${subset_pod5.baseName}.moves.bam
    """
}

process generate_signal_df {
    tag "${sample_id}"
    publishDir "${params.outdir}/annotated_data", mode: 'copy'

    input:
    tuple val(sample_id), path(subset_pod5), path(moves_bam)

    output:
    tuple val(sample_id), val(subset_pod5.baseName), path("*.csv"), emit: results

    script:
    """
    pod5_to_df.py --pod5 ${subset_pod5} --bam ${moves_bam} --sample_id ${sample_id}
    """
}

process visualize_signal {
    tag "${sample_id} - ${pod5_name}"
    publishDir "${params.qc_base_dir}/${sample_id}/${pod5_name}", mode: 'copy'

    input:
    tuple val(sample_id), val(pod5_name), path(csv)

    output:
    path "*.pdf"

    script:
    def read_id = csv.baseName.replace("_mapped", "")
    def circle_flag = params.emit_circles ? "--emit_circles" : ""
    """
    plot_signal.py \\
        --csv ${csv} \\
        --output ${read_id}.pdf \\
        --title "Read: ${read_id}" \\
        --figwidth ${params.figwidth} \\
        --figheight ${params.figheight} \\
        ${circle_flag}
    """
}