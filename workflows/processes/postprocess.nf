process qualimap {
    tag "$name"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    input:
    tuple val(name), path(bam) from ch_bam_dedup_for_qualimap

    output:
    file "${bam.baseName}_qualimap" into ch_qualimap_results_for_multiqc

    script:
    gcref = params.genome == 'GRCh37' ? '-gd HUMAN' : ''
    gcref = params.genome == 'GRCm38' ? '-gd MOUSE' : ''
    def avail_mem = task.memory ? ((task.memory.toBytes() - 6000000000) / task.cpus) : false
    def sort_mem = avail_mem && avail_mem > 2000000000 ? "-m $avail_mem" : ''
    """
    samtools sort $bam \\
        -@ ${task.cpus} $sort_mem \\
        -o ${bam.baseName}.sorted.bam
    qualimap bamqc $gcref \\
        -bam ${bam.baseName}.sorted.bam \\
        -outdir ${bam.baseName}_qualimap \\
        --collect-overlap-pairs \\
        --java-mem-size=${task.memory.toGiga()}G \\
        -nt ${task.cpus}
    """
}


process preseq {
    tag "$name"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    input:
    tuple val(name), path(bam) from ch_bam_for_preseq

    output:
    file "${bam.baseName}.ccurve.txt" into preseq_results

    script:
    def avail_mem = task.memory ? ((task.memory.toBytes() - 6000000000) / task.cpus) : false
    def sort_mem = avail_mem && avail_mem > 2000000000 ? "-m $avail_mem" : ''
    """
    samtools sort $bam \\
        -@ ${task.cpus} $sort_mem \\
        -o ${bam.baseName}.sorted.bam
    preseq lc_extrap -v -B ${bam.baseName}.sorted.bam -o ${bam.baseName}.ccurve.txt
    """
}



process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    path multiqc_config from ch_multiqc_config
    path ('fastqc/*') from ch_fastqc_results_for_multiqc.collect().ifEmpty([])
    path ('trimgalore/*') from ch_trim_galore_results_for_multiqc.collect().ifEmpty([])
    path ('bismark/*') from ch_bismark_align_log_for_multiqc.collect().ifEmpty([])
    path ('bismark/*') from ch_bismark_dedup_log_for_multiqc.collect().ifEmpty([])
    path ('bismark/*') from ch_bismark_splitting_report_for_multiqc.collect().ifEmpty([])
    path ('bismark/*') from ch_bismark_mbias_for_multiqc.collect().ifEmpty([])
    path ('bismark/*') from ch_bismark_reports_results_for_multiqc.collect().ifEmpty([])
    path ('bismark/*') from ch_bismark_summary_results_for_multiqc.collect().ifEmpty([])
    path ('samtools/*') from ch_flagstat_results_for_multiqc.flatten().collect().ifEmpty([])
    path ('samtools/*') from ch_samtools_stats_results_for_multiqc.flatten().collect().ifEmpty([])
    path ('picard/*') from ch_markDups_results_for_multiqc.flatten().collect().ifEmpty([])
    path ('methyldackel/*') from ch_methyldackel_results_for_multiqc.flatten().collect().ifEmpty([])
    path ('qualimap/*') from ch_qualimap_results_for_multiqc.collect().ifEmpty([])
    path ('preseq/*') from preseq_results.collect().ifEmpty([])
    path ('software_versions/*') from ch_software_versions_yaml_for_multiqc.collect()
    path workflow_summary from create_workflow_summary(summary)

    output:
    path "*multiqc_report.html" into ch_multiqc_report
    path "*_data"
    path "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config . \\
        -m custom_content -m picard -m qualimap -m bismark -m samtools -m preseq -m cutadapt -m fastqc
    """
}


process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}
