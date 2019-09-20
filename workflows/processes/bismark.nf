

process makeBismarkIndex {
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
               saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file fasta

    output:
    file "BismarkIndex"

    script:
    aligner = params.aligner == 'bismark_hisat' ? '--hisat2' : '--bowtie2'
    slam = params.slamseq ? '--slam' : ''
    """
    mkdir BismarkIndex
    cp $fasta BismarkIndex/
    bismark_genome_preparation $aligner $slam BismarkIndex
    """
}

process bismark_align {
    tag "$name"
    publishDir "${params.outdir}/bismark_alignments", mode: 'copy',
        saveAs: {filename ->
            if( filename.indexOf(".fq.gz") > 0 ) "unmapped/$filename"
            else if( filename.indexOf("report.txt") > 0 ) "logs/$filename"
            else if( (!params.saveAlignedIntermediates && !params.nodedup && !params.rrbs).every() && filename == "where_are_my_files.txt" ) filename
            else if( (params.saveAlignedIntermediates || params.nodedup || params.rrbs).any() && filename != "where_are_my_files.txt" ) filename
            else null
        }

    input:
    tuple val(name), path(reads)
    file index
    file wherearemyfiles
    file knownsplices

    output:
    tuple val(name), path("*.bam"), emit: ch_bam
    tuple val(name), path("*report.txt"), emit: ch_bismark_align_log
    file "*.fq.gz" optional true
    file "where_are_my_files.txt"

    script:
    aligner = params.aligner == "bismark_hisat" ? "--hisat2" : "--bowtie2"
    splicesites = params.aligner == "bismark_hisat" && knownsplices.name != 'null' ? "--known-splicesite-infile <(hisat2_extract_splice_sites.py ${knownsplices})" : ''
    pbat = params.pbat ? "--pbat" : ''
    non_directional = params.single_cell || params.zymo || params.non_directional ? "--non_directional" : ''
    unmapped = params.unmapped ? "--unmapped" : ''
    mismatches = params.relaxMismatches ? "--score_min L,0,-${params.numMismatches}" : ''
    multicore = ''
    if( task.cpus ){
        // Numbers based on recommendation by Felix for a typical mouse genome
        if( params.single_cell || params.zymo || params.non_directional ){
            cpu_per_multicore = 5
            mem_per_multicore = (18.GB).toBytes()
        } else {
            cpu_per_multicore = 3
            mem_per_multicore = (13.GB).toBytes()
        }
        // How many multicore splits can we afford with the cpus we have?
        ccore = ((task.cpus as int) / cpu_per_multicore) as int
        // Check that we have enough memory, assuming 13GB memory per instance (typical for mouse alignment)
        try {
            tmem = (task.memory as nextflow.util.MemoryUnit).toBytes()
            mcore = (tmem / mem_per_multicore) as int
            ccore = Math.min(ccore, mcore)
        } catch (all) {
            log.debug "Not able to define bismark align multicore based on available memory"
        }
        if( ccore > 1 ){
          multicore = "--multicore $ccore"
        }
    }
    if( params.singleEnd ) {
        """
        bismark $aligner \\
            --bam $pbat $non_directional $unmapped $mismatches $multicore \\
            --genome $index \\
            $reads \\
            $splicesites
        """
    } else {
        """
        bismark $aligner \\
            --bam $pbat $non_directional $unmapped $mismatches $multicore \\
            --genome $index \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            $splicesites
        """
    }
}


process bismark_deduplicate {
    tag "$name"
    publishDir "${params.outdir}/bismark_deduplicated", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename"}

    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("*.deduplicated.bam"), emit: ch_bam_dedup
    tuple val(name), path("*.deduplication_report.txt"), emit: ch_bismark_dedup_log

    script:
    if( params.singleEnd ) {
        """
        deduplicate_bismark -s --bam $bam
        """
    } else {
        """
        deduplicate_bismark -p --bam $bam
        """
    }
}

process bismark_methXtract {
    tag "$name"
    publishDir "${params.outdir}/bismark_methylation_calls", mode: 'copy',
        saveAs: {filename ->
            if( filename.indexOf("splitting_report.txt" ) > 0 ) "logs/$filename"
            else if( filename.indexOf("M-bias" ) > 0) "m-bias/$filename"
            else if( filename.indexOf(".cov" ) > 0 ) "methylation_coverage/$filename"
            else if( filename.indexOf("bedGraph" ) > 0 ) "bedGraph/$filename"
            else "methylation_calls/$filename"
        }

    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("*splitting_report.txt"), emit: ch_bismark_splitting_report
    tuple val(name), path("*.M-bias.txt"), emit: ch_bismark_mbias
    file '*.{png,gz}'

    script:
    comprehensive = params.comprehensive ? '--comprehensive --merge_non_CpG' : ''
    meth_cutoff = params.meth_cutoff ? "--cutoff ${params.meth_cutoff}" : ''
    multicore = ''
    if( task.cpus ){
        // Numbers based on Bismark docs
        ccore = ((task.cpus as int) / 10) as int
        if( ccore > 1 ){
          multicore = "--multicore $ccore"
        }
    }
    buffer = ''
    if( task.memory ){
        mbuffer = (task.memory as nextflow.util.MemoryUnit) - 2.GB
        // only set if we have more than 6GB available
        if( mbuffer.compareTo(4.GB) == 1 ){
          buffer = "--buffer_size ${mbuffer.toGiga()}G"
        }
    }
    if(params.singleEnd) {
        """
        bismark_methylation_extractor $comprehensive $meth_cutoff \\
            $multicore $buffer \\
            --bedGraph \\
            --counts \\
            --gzip \\
            -s \\
            --report \\
            $bam
        """
    } else {
        """
        bismark_methylation_extractor $comprehensive $meth_cutoff \\
            $multicore $buffer \\
            --ignore_r2 2 \\
            --ignore_3prime_r2 2 \\
            --bedGraph \\
            --counts \\
            --gzip \\
            -p \\
            --no_overlap \\
            --report \\
            $bam
        """
    }
}


process bismark_report {
    tag "$name"
    publishDir "${params.outdir}/bismark_reports", mode: 'copy'

    input:
    tuple val(name), path(align_log), path(dedup_log), path(splitting_report), path(mbias) from ch_bismark_logs_for_bismark_report

    output:
    file '*{html,txt}' into ch_bismark_reports_results_for_multiqc

    script:
    """
    bismark2report \\
        --alignment_report $align_log \\
        --dedup_report $dedup_log \\
        --splitting_report $splitting_report \\
        --mbias_report $mbias
    """
}


process bismark_summary {
    publishDir "${params.outdir}/bismark_summary", mode: 'copy'

    input:
    file ('*') from ch_bam_for_bismark_summary.collect()
    file ('*') from ch_bismark_align_log_for_bismark_summary.collect()
    file ('*') from ch_bismark_dedup_log_for_bismark_summary.collect()
    file ('*') from ch_bismark_splitting_report_for_bismark_summary.collect()
    file ('*') from ch_bismark_mbias_for_bismark_summary.collect()

    output:
    file '*{html,txt}' into ch_bismark_summary_results_for_multiqc

    script:
    """
    bismark2summary
    """
}
