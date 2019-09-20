
include './processes/preprocess.nf' params(params)

fastqc()

if( params.notrim ){
    ch_trimmed_reads_for_alignment = ch_read_files
    ch_trim_galore_results_for_multiqc = Channel.from(false)
} else {
    trim_galore()
}
