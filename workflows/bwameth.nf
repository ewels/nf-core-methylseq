
include './processes/bwameth.nf' params(params)

/*
 * PREPROCESSING - Build bwa-mem index
 */
if( !params.bwa_meth_index ){
    makeBwaMemIndex(ch_fasta)
}

/*
 * PREPROCESSING - Index Fasta file
 */
if( !params.fasta_index ){
    makeFastaIndex()
}
