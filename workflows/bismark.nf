
include './processes/bismark.nf' params(params)

workflow bismark {

    // Setup
    ch_splicesites_for_bismark_hisat_align = params.known_splices ? Channel.fromPath("${params.known_splices}", checkIfExists: true).collect() : file('null')

    // Build reference genome index
    assert params.bismark_index || params.fasta : "No reference genome index or fasta file specified"

    if( params.bismark_index ){
        // We already have an index - create the channel
        Channel
            .fromPath(params.bismark_index, checkIfExists: true)
            .ifEmpty { exit 1, "Bismark index file not found: ${params.bismark_index}" }
            .set { ch_bismark_index }
    } else {
        // No bismark index - build from fasta
        Channel
            .fromPath(params.fasta, checkIfExists: true)
            .ifEmpty { exit 1, "fasta file not found : ${params.fasta}" }
            .set { ch_fasta }

        makeBismarkIndex(ch_fasta).set { ch_bismark_index }
    }

    bismark_align(
        ch_trimmed_reads_for_alignment,
        ch_bismark_index,
        ch_wherearemyfiles,
        ch_splicesites_for_bismark_hisat_align
    )
    bismark_deduplicate(
        bismark_align.out.ch_bam
    )
    bismark_methXtract(
        bismark_deduplicate.out.ch_bam_dedup
    )
    bismark_report()
    bismark_summary()

}
