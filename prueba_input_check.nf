#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

ch_input = Channel.fromPath( './assets/samplesheet.csv' )
ch_fasta = Channel.fromPath( './testdata/hg19.fa')
ch_fasta_fai = Channel.fromPath( './testdata/hg19.fa.fai')

include { INPUT_CHECK } from './subworkflows/local/input_check'
include { SAMTOOLS_SORT } from './modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_MPILEUP } from './modules/nf-core/modules/samtools/mpileup/main'
include { VARDICTJAVA } from './modules/nf-core/modules/vardictjava/main'
include { VARSCAN } from './modules/local/varscan'

workflow {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //INPUT_CHECK.out.reads.flatten().first().concat(INPUT_CHECK.out.reads.flatten().filter(~/.*.bam/).toList()).toList().view()

    //INPUT_CHECK.out.reads_bam.view()

    SAMTOOLS_SORT (
    	INPUT_CHECK.out.reads_bam
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions.first())

    //
    // MODULE: Run Samtools mpileup
    //

    ch_mpileup = SAMTOOLS_SORT.out.bam.flatten().first()
    .concat(SAMTOOLS_SORT.out.bam.flatten().filter(~/.*.bam/)
    .concat(INPUT_CHECK.out.reads.flatten().filter(~/.*.bed/))).toList()

    SAMTOOLS_MPILEUP (
    	ch_mpileup, ch_fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions.first())

    //
    // MODULE: Run VarScan
    //

    VARSCAN (
    	SAMTOOLS_MPILEUP.out.mpileup
    )

    ch_versions = ch_versions.mix(VARSCAN.out.versions.first())

    //
    // MODULE: Run Vardictjava
    //

    VARDICTJAVA (
    	INPUT_CHECK.out.reads, ch_fasta, ch_fasta_fai
    )

    ch_versions = ch_versions.mix(VARDICTJAVA.out.versions.first())

}
