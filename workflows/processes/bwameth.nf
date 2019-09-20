
process makeBwaMemIndex {
    tag "$fasta"
    publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta}*" into ch_bwa_meth_indices_for_bwamem_align

    script:
    """
    bwameth.py index $fasta
    """
}

process makeFastaIndex {
    tag "$fasta"
    publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta}.fai" into ch_fasta_index_for_methyldackel

    script:
    """
    samtools faidx $fasta
    """
}

process bwamem_align {
    tag "$name"
    publishDir "${params.outdir}/bwa-mem_alignments", mode: 'copy',
        saveAs: {filename ->
            if( !params.saveAlignedIntermediates && filename == "where_are_my_files.txt" ) filename
            else if( params.saveAlignedIntermediates && filename != "where_are_my_files.txt" ) filename
            else null
        }

    input:
    tuple val(name), path(reads) from ch_trimmed_reads_for_alignment
    file bwa_meth_indices from ch_bwa_meth_indices_for_bwamem_align.collect()
    file wherearemyfiles from ch_wherearemyfiles.collect()

    output:
    tuple val(name), path('*.bam') into ch_bam_for_samtools_sort_index_flagstat, ch_bam_for_preseq
    file "where_are_my_files.txt"

    script:
    fasta = bwa_meth_indices[0].toString() - '.bwameth' - '.c2t' - '.amb' - '.ann' - '.bwt' - '.pac' - '.sa'
    prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    bwameth.py \\
        --threads ${task.cpus} \\
        --reference $fasta \\
        $reads | samtools view -bS - > ${prefix}.bam
    """
}

process samtools_sort_index_flagstat {
    tag "$name"
    publishDir "${params.outdir}/bwa-mem_alignments", mode: 'copy',
        saveAs: {filename ->
            if(filename.indexOf("report.txt") > 0) "logs/$filename"
            else if( (!params.saveAlignedIntermediates && !params.nodedup && !params.rrbs).every() && filename == "where_are_my_files.txt") filename
            else if( (params.saveAlignedIntermediates || params.nodedup || params.rrbs).any() && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    tuple val(name), path(bam) from ch_bam_for_samtools_sort_index_flagstat
    file wherearemyfiles from ch_wherearemyfiles.collect()

    output:
    tuple val(name), path("${bam.baseName}.sorted.bam") into ch_bam_sorted_for_markDuplicates
    file "${bam.baseName}.sorted.bam.bai" into ch_bam_index
    file "${bam.baseName}_flagstat_report.txt" into ch_flagstat_results_for_multiqc
    file "${bam.baseName}_stats_report.txt" into ch_samtools_stats_results_for_multiqc
    file "where_are_my_files.txt"

    script:
    def avail_mem = task.memory ? ((task.memory.toBytes() - 6000000000) / task.cpus) : false
    def sort_mem = avail_mem && avail_mem > 2000000000 ? "-m $avail_mem" : ''
    """
    samtools sort $bam \\
        -@ ${task.cpus} $sort_mem \\
        -o ${bam.baseName}.sorted.bam
    samtools index ${bam.baseName}.sorted.bam
    samtools flagstat ${bam.baseName}.sorted.bam > ${bam.baseName}_flagstat_report.txt
    samtools stats ${bam.baseName}.sorted.bam > ${bam.baseName}_stats_report.txt
    """
}


process markDuplicates {
    tag "$name"
    publishDir "${params.outdir}/bwa-mem_markDuplicates", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".bam") == -1 ? "logs/$filename" : "$filename"}

    input:
    tuple val(name), path(bam) from ch_bam_sorted_for_markDuplicates

    output:
    tuple val(name), path("${bam.baseName}.markDups.bam") into ch_bam_dedup_for_methyldackel, ch_bam_dedup_for_qualimap
    file "${bam.baseName}.markDups.bam.bai" into ch_bam_index_for_methyldackel //ToDo check if this correctly overrides the original channel
    file "${bam.baseName}.markDups_metrics.txt" into ch_markDups_results_for_multiqc

    script:
    if( !task.memory ){
        log.info "[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this."
        avail_mem = 3
    } else {
        avail_mem = task.memory.toGiga()
    }
    """
    picard -Xmx${avail_mem}g MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${bam.baseName}.markDups.bam \\
        METRICS_FILE=${bam.baseName}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT
    samtools index ${bam.baseName}.markDups.bam
    """
}


process methyldackel {
    tag "$name"
    publishDir "${params.outdir}/MethylDackel", mode: 'copy'

    input:
    tuple val(name), path(bam) from ch_bam_dedup_for_methyldackel
    file bam_index from ch_bam_index_for_methyldackel
    file fasta from ch_fasta
    file fasta_index from ch_fasta_index_for_methyldackel

    output:
    file "${bam.baseName}*" into ch_methyldackel_results_for_multiqc

    script:
    allcontexts = params.comprehensive ? '--CHG --CHH' : ''
    mindepth = params.mindepth > 0 ? "--minDepth ${params.mindepth}" : ''
    ignoreFlags = params.ignoreFlags ? "--ignoreFlags" : ''
    methylKit = params.methylKit ? "--methylKit" : ''
    """
    MethylDackel extract $allcontexts $ignoreFlags $methylKit $mindepth $fasta $bam
    MethylDackel mbias $allcontexts $ignoreFlags $fasta $bam ${bam.baseName} --txt > ${bam.baseName}_methyldackel.txt
    """
}
