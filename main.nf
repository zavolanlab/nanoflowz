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
    
    // 0. Build the minimap2 index once from the reference FASTA
    build_minimap2_index(params.ref)

    // 1. Initial Data Ingest: Create a stream of (sample_id, pod5_file)
    samples_ch = channel
        .fromPath(params.tsv)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_id, file(row.pod5)) }

    // 2.a Initial basecall: Parallel processing per POD5 chunk
    dorado_basecall(samples_ch)

    // 2a. Orient Strands (Emits both the oriented ubam and the stats tsv)
    orient_strands(dorado_basecall.out.ubam)

    // Collect all chunk TSVs and compile them into a master QC report
    merge_ts_stats(orient_strands.out.stats.collect())

    // 2.c Alignment: Pass both the unmapped BAMs and the compiled index
    minimap2_align(orient_strands.out.ubam, build_minimap2_index.out.mmi)

    // ==============================================================================
    // CHUNK-LEVEL PROCESSING (Highly Parallel)
    // ==============================================================================

    // 2c. Fix Softclipped Alignments
    scinpas_fix_softclipped(minimap2_align.out.bam, params.ref)

    // 2d. Extract PolyA reads
    // scinpas_get_polyA(scinpas_fix_softclipped.out.bam, params.ref)

    // 2e. Append PolyA tails to the isolated PolyA reads
    //append_polyA_tails(scinpas_get_polyA.out.polyA_bam)

    // 2f. Normalize UMIs on the appended chunks
    normalize_umi_lengths(scinpas_fix_softclipped.out.bam)
    
    // ==============================================================================
    // SAMPLE-LEVEL PROCESSING (Merged Data)
    // ==============================================================================
    
    // 3. Merge: Combine small BAMs into a single sample-level BAM
    merge_input_ch = normalize_umi_lengths.out.bam.groupTuple()
    merge_bams(merge_input_ch)

    // 4.a UMI Deduplication
    umi_tools_dedup(merge_bams.out.bam.join(merge_bams.out.bai))

    // 4.b Redefine NH Tags to correct for alignment filtering after UMI deduplication
    redefine_nh_tags(umi_tools_dedup.out.bam)

    // ==============================================================================
    // ISOFORM ANALYSIS
    // ==============================================================================
    
    // 5.a Collect all redefined BAM files
    // all_bams_ch = redefine_nh_tags.out.bam.map { it[1] }.collect()
    
    // 5.b build "enriched" transcriptome annotation using all aligned reads across all input samples
    // transcriptome_annotation_enrichment(
    //     all_bams_ch, 
    //     params.reference_gtf
    // )

    // 5.c Assign individual alignments for each sample to the transcript isoforms in enriched transcriptome
    // read_to_transcript_assignment(
    //     redefine_nh_tags.out.bam, 
    //     transcriptome_annotation_enrichment.out.tsv
    // )

    // ==============================================================================
    // QC: Raw Current Signal visualization and annotation
    // ==============================================================================

    // a. Selection: Pick the random Read IDs to investigate
    // extract_read_ids(redefine_nh_tags.out.bam)

    // // b. We take all original POD5 paths and group them by sample_id
    // // This allows one process to search all files at once.
    // all_pod5s_per_sample = samples_ch.map { id, pod5 -> [id, pod5] }.groupTuple()
    
    // filter_input_ch = extract_read_ids.out.ids_file.join(all_pod5s_per_sample)
    // filter_pod5_combined(filter_input_ch)

    // // c. Deep Dive: Re-run Dorado for moves on the subsampled POD5
    // dorado_emit_moves(filter_pod5_combined.out.pod5)

    // // d. Prepare Data: Join the subsampled POD5 and its Move-BAM to generate CSVs
    // final_input_ch = filter_pod5_combined.out.pod5.join(dorado_emit_moves.out.bam)
    // generate_signal_df(final_input_ch)

    // // e. Visualize: Plot the squiggles for each individual read
    // visualize_input = generate_signal_df.out.results.transpose()
    // visualize_signal(visualize_input)
}

/*
========================================================================================
    PROCESS DEFINITIONS
========================================================================================
*/

process build_minimap2_index {
    label 'process_medium'
    label 'env_samtools'

    input:
    path fasta

    output:
    path "${fasta.baseName}.mmi", emit: mmi

    script:
    """
    minimap2 -x splice -t ${task.cpus} -d ${fasta.baseName}.mmi ${fasta}
    """
}

process dorado_basecall {
    tag "${sample_id}"
    label 'process_gpu'
    label 'env_samtools'  
    
    input:
    tuple val(sample_id), path(pod5_input)

    output:
    // Notice we emit 'ubam' (unmapped BAM) here
    tuple val(sample_id), path("chunk_${pod5_input.baseName}.ubam"), emit: ubam

    script:
    """
    OUT_UBAM="chunk_${pod5_input.baseName}.ubam"
    
    # Run dorado without alignment flags, outputting directly to an unmapped BAM
    ${params.dorado} basecaller \\
        ${params.model} \\
        ${pod5_input} \\
        --estimate-poly-a \\
        --poly-a-config ${params.polyA} \\
        --device "cuda:all" > \$OUT_UBAM
    """
}

process orient_strands {
    tag "${sample_id} - chunk"
    label 'process_low'
    label 'env_bam_processing_with_python' 

    input:
    tuple val(sample_id), path(ubam)

    output:
    tuple val(sample_id), path("${ubam.baseName}.oriented.ubam"), emit: ubam
    path "${ubam.baseName}.ts_stats.tsv", emit: stats

    script:
    """
    orient_reads.py \\
        --input_bam ${ubam} \\
        --output_bam ${ubam.baseName}.oriented.ubam \\
        --stats_tsv ${ubam.baseName}.ts_stats.tsv \\
        --sample_id ${sample_id}
    """
}

process merge_ts_stats {
    publishDir "${params.qc_base_dir}", mode: 'copy'
    label 'process_single'
    
    input:
    path tsv_files
    
    output:
    path "master_TS_orientation_stats.tsv"
    
    script:
    """
    # 1. Extract the header from the first TSV file
    head -n 1 \$(ls ${tsv_files} | head -n 1) > master_TS_orientation_stats.tsv
    
    # 2. Append the data from all chunk TSVs (skipping the header line in each)
    for file in ${tsv_files}; do
        tail -n +2 \$file >> master_TS_orientation_stats.tsv
    done
    """
}

process minimap2_align {
    tag "${sample_id}"
    label 'process_medium'
    label 'env_samtools'

    input:
    tuple val(sample_id), path(ubam)
    path mmi_index

    output:
    tuple val(sample_id), path("chunk_${ubam.baseName}.bam"), emit: bam

    script:
    """
    OUT_BAM="chunk_${ubam.baseName}.bam"
    
    # 1. samtools fastq -T "*" extracts the fastq AND appends all BAM tags to the header.
    # 2. minimap2 -y reads those tags and securely copies them into the aligned BAM output.
    samtools fastq -@ ${task.cpus} -T "*" ${ubam} \\
        | minimap2 -y -ax splice -Y -t ${task.cpus} ${mmi_index} - \\
        | samtools sort -m 2G -@ ${task.cpus} -o \$OUT_BAM -
    """
}

process scinpas_fix_softclipped {
    tag "${sample_id} - chunk"
    label 'process_medium_low_cpu'
    label 'env_bam_processing_with_python'

    input:
    tuple val(sample_id), path(bam)
    path fasta

    output:
    tuple val(sample_id), path("${bam.baseName}.fixed.bam"), emit: bam
    path "${bam.baseName}.fix_stats.csv", emit: csv

    script:
    """
    # Create index for pysam
    samtools index -@ ${task.cpus} ${bam}
    
    scinpas_fix_softclipped \\
        --bam_file ${bam} \\
        --fasta ${fasta} \\
        --bam_out unsorted_fixed.bam \\
        --csv_out ${bam.baseName}.fix_stats.csv \\
        --exact_out \\
        --one_based_tags
        
    samtools sort -@ ${task.cpus} -m 2G unsorted_fixed.bam > ${bam.baseName}.fixed.bam
    rm unsorted_fixed.bam
    """
}

process scinpas_get_polyA {
    tag "${sample_id} - chunk"
    label 'process_medium_low_cpu'
    label 'env_bam_processing_with_python'

    input:
    tuple val(sample_id), path(bam)
    path fasta

    output:
    tuple val(sample_id), path("${bam.baseName}.polyA.bam"), emit: polyA_bam
    tuple val(sample_id), path("${bam.baseName}.non_polyA.bam"), emit: non_polyA_bam

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    
    scinpas_get_polyA \\
        --bam_input ${bam} \\
        --o_polyA ${bam.baseName}.polyA.bam \\
        --o_nonpolyA ${bam.baseName}.non_polyA.bam \\
        --o_low_q_polyA ${bam.baseName}.lowQ_polyA.bam \\
        --fasta ${fasta} \\
        --percentage_threshold 80 \\
        --length_threshold 8 \\
        --use_fc 1 \\
        --exact_out
    """
}

process scinpas_get_unique_cs {
    tag "${sample_id} - chunk"
    label 'process_low'
    label 'env_bam_processing_with_python'

    input:
    tuple val(sample_id), path(polyA_bam)

    output:
    tuple val(sample_id), path("${polyA_bam.baseName}.cleavage_sites.bed"), emit: bed

    script:
    """
    samtools index -@ ${task.cpus} ${polyA_bam}
    
    scinpas_get_unique_cs \\
        --bam ${polyA_bam} \\
        --bed_out ${polyA_bam.baseName}.cleavage_sites.bed \\
        --use_fc 1 \\
        --sample_name ${sample_id} \\
        --exact_out
    """
}

process append_polyA_tails {
    tag "${sample_id} - chunk"
    label 'process_medium'
    label 'env_bam_processing_with_python'

    input:
    tuple val(sample_id), path(polyA_bam)

    output:
    tuple val(sample_id), path("${polyA_bam.baseName}.pA_appended.bam"), emit: bam

    script:
    """
    # Runs the custom script placed in nanoflowz/bin/
    append_polyA_tail.py \\
        --input_bam_file ${polyA_bam} \\
        --output_bam_file ${polyA_bam.baseName}.pA_appended.bam
    """
}

process normalize_umi_lengths {
    tag "${sample_id} - chunk"
    label 'process_low'
    label 'env_bam_processing_with_python'

    input:
    tuple val(sample_id), path(bam)

    output:
    // Emits the normalized chunk
    tuple val(sample_id), path("${bam.baseName}.normalized.bam"), emit: bam

    script:
    """
    # Call the CLI tool from zavolab_pyutils to normalize UMI lengths
    # This is a crucial step to prevent umi_tools from crashing due to UMIs of variable lengths.
    normalize_umi_lengths \\
        --input_bam ${bam} \\
        --output_bam ${bam.baseName}.normalized.bam \\
        --target_len 26
        
    # We skip samtools index here because samtools merge doesn't need it
    """
}

process merge_bams {
    tag "${sample_id}"
    publishDir "${params.outdir}/basecalling", mode: 'copy'
    label 'process_medium'
    label 'env_samtools'

    input:
    tuple val(sample_id), path(bams)

    output:
    tuple val(sample_id), path("${sample_id}.dorado.sorted.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.dorado.sorted.bam.bai"), emit: bai

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_id}.dorado.sorted.bam ${bams}
    samtools index -@ ${task.cpus} ${sample_id}.dorado.sorted.bam
    """
}

process umi_tools_dedup {
    tag "${sample_id}"
    // Deduplication holds UMIs in RAM, requires high memory
    label 'process_high_memory_low_cpu' 
    label 'env_umi_tools'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.name_sorted.bam"), emit: bam

    script:
    """
    umi_tools dedup \\
        --extract-umi-method=tag \\
        --umi-tag=RX \\
        --method unique \\
        -I ${bam} \\
        -S unsorted_dedup.bam
    
    # After deduplication, we sort the BAM by NAME to prepare for downstream processing
    samtools sort -n -@ ${task.cpus} -m 2G unsorted_dedup.bam > ${sample_id}.dedup.name_sorted.bam
    
    rm unsorted_dedup.bam
    """
}

process redefine_nh_tags {
    tag "${sample_id}"
    label 'process_high_memory_low_cpu'
    label 'env_bam_processing_with_python' 
    publishDir "${params.outdir}/map_genome_merged_UMIdedup", mode: 'copy'
    
    // we assume that input BAM files were name-sorted in the previous step.
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.redefined_NH.sorted.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.redefined_NH.sorted.bam.bai"), emit: bai

    script:
    """
    # 1. Use custom Python script from zavolab_pyutils to correct NH tags and assign MAPQ=255 for unique alignments (as STAR aligner does).
    # For MM reads, MAPQ=0 is assigned.
    redefine_qual_and_NHtag \\
        --input_bam_file ${bam} \\
        --out_bam_file unsorted.bam
        
    # 2. Sort and index
    samtools sort -@ ${task.cpus} -m 2G unsorted.bam > ${sample_id}.redefined_NH.sorted.bam
    samtools index -@ ${task.cpus} ${sample_id}.redefined_NH.sorted.bam
    
    # 3. Clean up intermediate file
    rm unsorted.bam
    """
}

process transcriptome_annotation_enrichment {
    publishDir "${params.outdir}/transcriptome", mode: 'copy'
    label 'process_high_memory_low_cpu'
    label 'env_isoform'

    input:
    path bams
    path reference_gtf

    output:
    path "master_enriched.tsv", emit: tsv
    path "master_enriched.gtf", emit: gtf

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
    label 'process_dynamic_memory'
    label 'env_isoform'

    input:
    tuple val(sample_id), path(bam)
    path enriched_tsv

    output:
    tuple val(sample_id), path("${sample_id}_assignments.tsv.gz"), emit: tsv
    path "${sample_id}_unassigned.tsv.gz", optional: true

    script:
    """
    assign_ONT_reads_to_isoforms.py \\
        --mode assign \\
        --input_bam_files ${bam} \\
        --input_gtf_file ${enriched_tsv} \\
        --output_prefix ${sample_id}    
    """
}

process extract_read_ids {
    tag "${sample_id}"
    label 'process_single'
    label 'env_samtools'
    
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
    label 'process_single'
    label 'env_pod5'

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
    label 'process_gpu'
    label 'env_samtools'    
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
    label 'process_single'
    label 'env_pod5'

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
    label 'process_low'
    label 'env_plot'

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