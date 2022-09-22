#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VARSCAN } from '../../modules/local/varscan' addParams( options: [:] )

workflow test_varscan {

    input = [ [ id:'test' ], // meta map
              file(params.test_data['homo_sapiens']['illumina']['mpileup']['test.mpileup.gz'], checkIfExists: true)
            
              

    VARSCAN ( input,
              file(params.test_data['homo_sapiens']['genome']['genome_fasta'], checkIfExists: true),
              "hg19",
              "100" )
}