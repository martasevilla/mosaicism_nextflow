#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

workflow test_input_check {

  // Define each part in a variable to simplify future modifications
  input = [ [ id:'test' ], // meta map
    file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam'], checkIfExists: true),
    file(params.test_data['sarscov2']['illumina']['test_single_end_sorted_bam_bai'], checkIfExists: true) ]
  fasta = file(params.test_data['sarscov2']['genome']['genome_fasta'], checkIfExists: true)
  genome   = "hg19"
  binSize  = "100"
  gapsRef  = file(params.test_data['sarscov2']['genome']['test_bed'], checkIfExists: true)

  RUN_CNVNATOR (input, fasta, genome, binSize, gapsRef)

}