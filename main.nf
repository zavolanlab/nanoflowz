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

    // 0a. Optional: verify that the GTF and genome FASTA are compatible before doing any real work
    if (params.check_gtf_compatibility) {
        check_gtf_genome_compatibility(params.ref, params.reference_gtf)
    }

    // 0b. Build the minimap2 index once from the reference FASTA
    build_minimap2_index(params.ref)

    // 1. Initial Data Ingest: Create a stream of (sample_id, pod5_file)
    samples_ch = channel
        .fromPath(params.tsv)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample_id, file(row.pod5)) }

    // 2.a Initial basecall: Parallel processing per POD5 chunk
    dorado_basecall(samples_ch)

    // 2.b. Orient Strands (Emits both the oriented ubam and the stats tsv)
    orient_strands(dorado_basecall.out.ubam)

    // Collect all chunk TSVs and compile them into a master QC report
    merge_ts_stats(orient_strands.out.stats.collect())

    // 2.c Alignment: Pass both the unmapped BAMs and the compiled index
    minimap2_align(orient_strands.out.ubam, build_minimap2_index.out.mmi)
    
    // ==============================================================================
    // CHUNK-LEVEL PROCESSING (Highly Parallel)
    // ==============================================================================

    // 2.d. Normalize UMIs on the appended chunks
    normalize_umi_lengths(minimap2_align.out.bam)

    // 2.e. Fix Softclipped Alignments
    scinpas_fix_softclipped(normalize_umi_lengths.out.bam, params.ref)
    
    // 2.f. Extract PolyA reads
    scinpas_get_polyA(scinpas_fix_softclipped.out.bam, params.ref)

    // Merge the PolyA stats into a master QC report
    merge_polyA_stats(scinpas_get_polyA.out.stats.collect())

    // 2.g. Append PolyA tails based on `pt` tags
    append_polyA_tails(scinpas_get_polyA.out.polyA_bam)
    
    // Merge the append stats into a master QC report
    merge_pt_stats(append_polyA_tails.out.stats.collect())    
    
    // ==============================================================================
    // SAMPLE-LEVEL PROCESSING (Merged Data)
    // ==============================================================================
    
    // 3. Merge: Combine small BAMs into a single sample-level BAM
    merge_input_ch = append_polyA_tails.out.bam.groupTuple()
    merge_bams(merge_input_ch)
    
    // Optional: Also merge the non-polyA BAMs for inspection
    if (params.merge_non_polyA) {
        merge_non_polyA_input_ch = scinpas_get_polyA.out.non_polyA_bam.groupTuple()
        merge_non_polyA_bams(merge_non_polyA_input_ch)
    }

    // Optional merge of skipped pt reads
    if (params.merge_skipped_pt) {
        merge_skipped_pt_input_ch = append_polyA_tails.out.skipped_bam.groupTuple()
        merge_skipped_pt_bams(merge_skipped_pt_input_ch)
    }

    // 4.a UMI Deduplication
    umi_tools_dedup(merge_bams.out.bam.join(merge_bams.out.bai))

    // 4.b Redefine NH Tags to correct for alignment filtering after UMI deduplication
    redefine_nh_tags(umi_tools_dedup.out.bam)
    
    // 4.c Generate BigWigs for Cleavage Sites using the custom scripts
    make_bigwig_for_cleavage_sites(redefine_nh_tags.out.bam.join(redefine_nh_tags.out.bai), params.ref)

    // ==============================================================================
    // GENE ASSIGNMENT
    // ==============================================================================

    // Assign reads to genes via featureCounts
    assign_alignments_to_genes(
        redefine_nh_tags.out.bam.join(redefine_nh_tags.out.bai), 
        params.reference_gtf
    )

    // ==============================================================================
    // polyA tail length bigwig generation, collection of tabular data, and visualization
    // ==============================================================================

    // Generate PolyA tail length BigWigs (Mean/Median depending on configuration)
    make_bigwig_for_polya_length(
        assign_alignments_to_genes.out.bam.join(assign_alignments_to_genes.out.bai), params.ref)

    // Extract Read ID, pt, and XT tags to a TSV Table
    extract_read_tags_tsv(
        assign_alignments_to_genes.out.bam.join(assign_alignments_to_genes.out.bai)
    )

    // Generate Comparative Boxplots across all samples
    // Extract just the file (index 1 of the tuple) and collect them into a list
    tsv_list_ch = extract_read_tags_tsv.out.tsv.map { sample_id, tsv_file -> tsv_file }.collect()
    
    visualize_polyA_tail_length_distribution(tsv_list_ch)

    // ==============================================================================
    // QC: MAPPING STATISTICS
    // ==============================================================================
    
    ch_to_count = dorado_basecall.out.ubam.map{ id, bam -> [id, bam, "01_dorado_basecall"] }
        .mix(
            orient_strands.out.ubam.map{ id, bam -> [id, bam, "02_reverse_complemented_backward_oriented_reads"] },
            minimap2_align.out.bam.map{ id, bam -> [id, bam, "03_aligned_with_minimap2"] },
            normalize_umi_lengths.out.bam.map{ id, bam -> [id, bam, "04_normalized_umi_lengths"] },
            scinpas_fix_softclipped.out.bam.map{ id, bam -> [id, bam, "05_fixed_softclipped_alignments"] },
            scinpas_get_polyA.out.polyA_bam.map{ id, bam -> [id, bam, "06_extracted_polyA_reads"] },
            append_polyA_tails.out.bam.map{ id, bam -> [id, bam, "07_appended_polyA_tails"] },
            umi_tools_dedup.out.bam.map{ id, bam -> [id, bam, "08_umi_deduped"] },
            redefine_nh_tags.out.bam.map{ id, bam -> [id, bam, "09_nh_tags_and_MAPQ_redefined"] },
            assign_alignments_to_genes.out.bam.map{ id, bam -> [id, bam, "10_geneID_assigned"] }
        )

    // Call count_reads
    all_read_counts = count_reads(ch_to_count)

    aggregate_read_stats(all_read_counts.collect())

    // ==============================================================================
    // ISOFORM ANALYSIS
    // ==============================================================================
    
    // 5.a Collect all redefined BAM files
    all_bams_ch = redefine_nh_tags.out.bam.map { it[1] }.collect()
    
    // 5.b build "enriched" transcriptome annotation using all aligned reads across all input samples
    transcriptome_annotation_enrichment(
        all_bams_ch, 
        params.reference_gtf
    )

    // 5.c Assign individual alignments for each sample to the transcript isoforms in enriched transcriptome
    read_to_transcript_assignment(
        redefine_nh_tags.out.bam, 
        transcriptome_annotation_enrichment.out.tsv
    )

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
    label 'process_medium_low_cpu'
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
        --sample_id ${sample_id} \\
        --tag_basecaller_ts ${params.tag_basecaller_ts} \\
        --tag_original_ts ${params.tag_original_ts}
    """
}

process merge_ts_stats {
    publishDir "${params.outdir}/QC", mode: 'copy'
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
        | minimap2 -y -ax splice:hq --secondary=no -Y -t ${task.cpus} ${mmi_index} - \\
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
    def shift_flag = params.shift_ambiguous_cs ? "--shift_ambiguous_cs" : ""
    """
    samtools index -@ ${task.cpus} ${bam}
    
    scinpas_fix_softclipped \\
        --bam_file ${bam} \\
        --fasta ${fasta} \\
        --bam_out unsorted_fixed.bam \\
        --csv_out ${bam.baseName}.fix_stats.csv \\
        --exact_out \\
        --one_based_tags \\
        --tag_orig_cs ${params.tag_orig_cs} \\
        --tag_fixed_cs ${params.tag_fixed_cs} \\
        ${shift_flag}

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
    path "${bam.baseName}.polyA_stats.tsv", emit: stats

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
        --exact_out \\
        --stats_tsv ${bam.baseName}.polyA_stats.tsv \\
        --sample_id ${sample_id} \\
        --min_phred ${params.min_phred_score} \\
        --tag_phred_mapped ${params.tag_phred_mapped} \\
        --tag_phred_softclipped ${params.tag_phred_softclipped} \\
        --tag_orig_cs ${params.tag_orig_cs} \\
        --tag_fixed_cs ${params.tag_fixed_cs}
    """
}

process merge_polyA_stats {
    publishDir "${params.outdir}/QC", mode: 'copy'
    label 'process_single'
    
    input:
    path tsv_files
    
    output:
    path "master_polyA_stats.tsv"
    
    script:
    """
    # Extract header from the first file
    head -n 1 \$(ls ${tsv_files} | head -n 1) > master_polyA_stats.tsv
    
    # Append the body of all chunk TSVs
    for file in ${tsv_files}; do
        tail -n +2 \$file >> master_polyA_stats.tsv
    done
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
        --exact_out \\
        --tag_orig_cs ${params.tag_orig_cs} \\
        --tag_fixed_cs ${params.tag_fixed_cs}
    """
}

process append_polyA_tails {
    tag "${sample_id} - chunk"
    label 'process_medium_low_cpu'
    label 'env_bam_processing_with_python'

    input:
    tuple val(sample_id), path(polyA_bam)

    output:
    tuple val(sample_id), path("${polyA_bam.baseName}.pA_appended.bam"), emit: bam
    tuple val(sample_id), path("${polyA_bam.baseName}.pA_skipped.bam"), emit: skipped_bam
    path "${polyA_bam.baseName}.pt_stats.tsv", emit: stats

    script:
    """
    append_polyA_tail.py \\
        --input_bam ${polyA_bam} \\
        --output_appended_bam ${polyA_bam.baseName}.pA_appended.bam \\
        --output_skipped_bam ${polyA_bam.baseName}.pA_skipped.bam \\
        --stats_tsv ${polyA_bam.baseName}.pt_stats.tsv \\
        --sample_id ${sample_id} \\
        --tag_orig_cs ${params.tag_orig_cs} \\
        --tag_fixed_cs ${params.tag_fixed_cs}
    """
}

process merge_pt_stats {
    publishDir "${params.outdir}/QC", mode: 'copy'
    label 'process_single'
    
    input:
    path tsv_files
    
    output:
    path "master_pt_append_stats.tsv"
    
    script:
    """
    head -n 1 \$(ls ${tsv_files} | head -n 1) > master_pt_append_stats.tsv
    for file in ${tsv_files}; do
        tail -n +2 \$file >> master_pt_append_stats.tsv
    done
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

process merge_non_polyA_bams {
    tag "${sample_id}"
    publishDir "${params.outdir}/QC/non_polyA_alignments", mode: 'copy'
    label 'process_medium'
    label 'env_samtools'

    input:
    tuple val(sample_id), path(bams)

    output:
    tuple val(sample_id), path("${sample_id}.non_polyA.sorted.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.non_polyA.sorted.bam.bai"), emit: bai

    script:
    """
    # Merge the chunked non-polyA BAMs into a sample-level BAM
    samtools merge -@ ${task.cpus} ${sample_id}.non_polyA.sorted.bam ${bams}
    
    # Index it so it can be immediately loaded into IGV
    samtools index -@ ${task.cpus} ${sample_id}.non_polyA.sorted.bam
    """
}

process merge_skipped_pt_bams {
    tag "${sample_id}"
    publishDir "${params.outdir}/QC/skipped_pt_alignments", mode: 'copy'
    label 'process_medium'
    label 'env_samtools'

    input:
    tuple val(sample_id), path(bams)

    output:
    tuple val(sample_id), path("${sample_id}.skipped_pt.sorted.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.skipped_pt.sorted.bam.bai"), emit: bai

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_id}.skipped_pt.sorted.bam ${bams}
    samtools index -@ ${task.cpus} ${sample_id}.skipped_pt.sorted.bam
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

// for QC - collecting mapping stats at each step of the pipeline

process count_reads {
    tag "${sample_id} - ${step_name}"
    label 'process_medium'
    label 'env_samtools'

    input:
    tuple val(sample_id), path(bam), val(step_name)

    output:
    path "${bam.baseName}.${step_name}.tsv", emit: tsv

    script:
    """
    count_reads.sh ${sample_id} ${step_name} ${bam} ${task.cpus} ${bam.baseName}.${step_name}.tsv
    """
}

process aggregate_read_stats {
    publishDir "${params.outdir}/QC", mode: 'copy'
    label 'process_single'
    label 'env_bam_processing_with_python'
    
    input:
    path tsv_files
    
    output:
    path "master_read_tracking_stats.tsv"
    
    script:
    """
    aggregate_read_stats.py --input ${tsv_files} --output master_read_tracking_stats.tsv
    """
}

process assign_alignments_to_genes {
    tag "${sample_id}"
    publishDir "${params.outdir}/gene_assignments", mode: 'copy'
    label 'process_medium'
    label 'env_featureCounts' 

    input:
    tuple val(sample_id), path(bam), path(bai)
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}.gene_assigned.sorted.bam"), emit: bam
    tuple val(sample_id), path("${sample_id}.gene_assigned.sorted.bam.bai"), emit: bai
    path "${sample_id}.featureCounts.txt", emit: summary
    path "${sample_id}.featureCounts.txt.summary", emit: feature_stats

    script:
    """
    # 1. Run featureCounts with -R CORE to generate a text map of assignments
    featureCounts \\
        -T ${task.cpus} \\
        -L \\
        -M \\
        -O \\
        -s 1 \\
        -a ${gtf} \\
        -o ${sample_id}.featureCounts.txt \\
        -R CORE \\
        ${bam}

    # featureCounts automatically writes the text assignments to "<input_bam>.featureCounts"
    
    # 2. Inject the XT and XS tags directly into a new BAM
    tag_bam_with_fc.py \\
        --bam_in ${bam} \\
        --core_txt ${bam}.featureCounts \\
        --bam_out unsorted_tagged.bam

    # 3. Sort and index the tagged BAM
    samtools sort -@ ${task.cpus} -m 2G unsorted_tagged.bam > ${sample_id}.gene_assigned.sorted.bam
    samtools index -@ ${task.cpus} ${sample_id}.gene_assigned.sorted.bam

    # 4. Clean up
    rm unsorted_tagged.bam ${bam}.featureCounts
    """
}

process make_bigwig_for_cleavage_sites {
    tag "${sample_id}"
    publishDir "${params.outdir}/cleavage_sites_bigwigs", mode: 'copy'
    label 'process_medium_low_cpu'
    label 'env_bam_processing_with_python' 

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta

    output:
    tuple val(sample_id), path("${sample_id}.plus.bigwig"), emit: bw_plus
    tuple val(sample_id), path("${sample_id}.minus.bigwig"), emit: bw_minus
    tuple val(sample_id), path("${sample_id}.read_sum.tsv"), emit: tsv

    script:
    """
    # 1. Create genome index if not already present
    samtools faidx ${fasta}
    
    # 2. Extract Cleavage Sites with 1/NH weights
    extract_cs_bigwig.py \\
        --bam_in ${bam} \\
        --bed_out cs.bed \\
        --tag ${params.tag_quantification_cs} \\
        --include_multimappers ${params.include_multimappers}

    # 3. Sort the extracted BED
    sort -k1,1 -k2,2n cs.bed > cs.sorted.bed

    # 4. Generate Weighted BedGraphs per strand
    # Since intervals are exactly 1bp long, we just group by coordinate and sum the weights (col 5).
    # This natively outputs the exact 4-column format required for BedGraphs!
    awk '\$6 == "+"' cs.sorted.bed | bedtools groupby -g 1,2,3 -c 5 -o sum > plus.bg
    awk '\$6 == "-"' cs.sorted.bed | bedtools groupby -g 1,2,3 -c 5 -o sum > minus.bg

    # 5. Sort BedGraphs (bedGraphToBigWig strictly requires coordinate-sorted input)
    sort -k1,1 -k2,2n plus.bg > plus.sorted.bg
    sort -k1,1 -k2,2n minus.bg > minus.sorted.bg

    # 6. Convert to BigWig
    [ -s plus.sorted.bg ]  && bedGraphToBigWig plus.sorted.bg  ${fasta}.fai ${sample_id}.plus.bigwig  || touch ${sample_id}.plus.bigwig
    [ -s minus.sorted.bg ] && bedGraphToBigWig minus.sorted.bg ${fasta}.fai ${sample_id}.minus.bigwig || touch ${sample_id}.minus.bigwig

    # 7. Generate Summary TSV (chr, start, end, strand, weighted_count)
    bedtools groupby -i cs.sorted.bed -g 1,2,3,6 -c 5 -o sum > ${sample_id}.read_sum.tsv
    
    # Clean up intermediate large files
    rm cs.bed plus.bg minus.bg plus.sorted.bg minus.sorted.bg
    """
}

process check_gtf_genome_compatibility {
    label 'process_single'
    label 'env_samtools'

    input:
    path fasta
    path gtf

    output:
    path "compatibility_check.txt", emit: report

    script:
    """
    # Index the genome FASTA and extract chromosome names from column 1 of the .fai
    samtools faidx ${fasta}
    cut -f1 ${fasta}.fai | sort > genome_chroms.txt
    TOTAL_GENOME=\$(wc -l < genome_chroms.txt)

    # Extract unique chromosome names from column 1 of the GTF (skip comment lines)
    grep -v '^#' ${gtf} | cut -f1 | sort -u > gtf_chroms.txt

    # Count how many genome chromosomes appear in the GTF
    FOUND_IN_GTF=\$(comm -12 genome_chroms.txt gtf_chroms.txt | wc -l)

    {
        echo "Genome chromosomes (from .fai):    \$TOTAL_GENOME"
        echo "Genome chromosomes found in GTF:   \$FOUND_IN_GTF"
    } | tee compatibility_check.txt

    if [ "\$TOTAL_GENOME" -eq 0 ]; then
        echo "ERROR: No chromosomes found in genome FASTA index." >&2
        exit 1
    fi

    # Require at least 50% of genome chromosomes to be present in the GTF.
    # Uses integer arithmetic: found*100 < total*50  ↔  found/total < 0.5
    PCT=\$(( FOUND_IN_GTF * 100 / TOTAL_GENOME ))
    if [ \$(( FOUND_IN_GTF * 100 )) -lt \$(( TOTAL_GENOME * 50 )) ]; then
        echo "ERROR: Only \${FOUND_IN_GTF}/\${TOTAL_GENOME} genome chromosomes (\${PCT}%) are present in the GTF." >&2
        echo "At least 50% of genome chromosomes must appear in the GTF." >&2
        echo "Please verify that the genome FASTA and GTF correspond to the same reference assembly," >&2
        echo "or set '--check_gtf_compatibility false' to skip this check." >&2
        exit 1
    fi

    echo "OK: \${FOUND_IN_GTF}/\${TOTAL_GENOME} genome chromosomes (\${PCT}%) present in GTF." | tee -a compatibility_check.txt
    """
}
process make_bigwig_for_polya_length {
    tag "${sample_id}"
    publishDir "${params.outdir}/polya_length_bigwigs", mode: 'copy'
    label 'process_medium'
    label 'env_bam_processing_with_python' 

    input:
    tuple val(sample_id), path(bam), path(bai)
    path fasta

    output:
    tuple val(sample_id), path("${sample_id}.polya_${params.polya_tail_metric}.plus.bigwig"), emit: bw_plus
    tuple val(sample_id), path("${sample_id}.polya_${params.polya_tail_metric}.minus.bigwig"), emit: bw_minus

    script:
    """
    samtools faidx ${fasta}
    
    # Python script natively computes the weighted mean/median bedgraphs!
    compute_polya_bedgraphs.py \\
        --bam_in ${bam} \\
        --tag_cs ${params.tag_quantification_cs} \\
        --metric ${params.polya_tail_metric} \\
        --include_multimappers ${params.include_multimappers}

    # Sort BedGraphs
    sort -k1,1 -k2,2n plus.bg > plus.sorted.bg
    sort -k1,1 -k2,2n minus.bg > minus.sorted.bg

    # Convert to BigWig
    bedGraphToBigWig plus.sorted.bg ${fasta}.fai ${sample_id}.polya_${params.polya_tail_metric}.plus.bigwig
    bedGraphToBigWig minus.sorted.bg ${fasta}.fai ${sample_id}.polya_${params.polya_tail_metric}.minus.bigwig
    
    rm plus.bg minus.bg plus.sorted.bg minus.sorted.bg
    """
}

process extract_read_tags_tsv {
    tag "${sample_id}"
    publishDir "${params.outdir}/read_tag_tables", mode: 'copy'
    label 'process_medium_low_cpu'
    label 'env_bam_processing_with_python'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.tags.tsv.gz"), emit: tsv

    script:
    """
    extract_bam_tags.py \\
        --bam_in ${bam} \\
        --tsv_out ${sample_id}.tags.tsv.gz \\
        --tags ${params.bam_export_tags} \\
        --include_multimappers ${params.include_multimappers}
    """
}

process visualize_polyA_tail_length_distribution {
    publishDir "${params.outdir}/analysis_figures/polyA_tail_lengths", mode: 'copy'
    label 'process_medium'
    label 'env_plot' 

    input:
    path tsv_files

    output:
    path "*.pdf"

    script:
    """
    plot_polya_distributions.py \\
        --input ${tsv_files} \\
        --output_prefix polya_distribution
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