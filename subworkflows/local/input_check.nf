//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK_BAM {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_bam_channel(it) }
        .set { bam }

    emit:
    bam                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_bam_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    // add path(s) of the bam and bai file(s) to the meta map
    def bam_meta = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam file does not exist!\n${row.bam}"
    }
   if (!file(row.bed).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bed file does not exist!\n${row.bed}"
    }
        bam_meta = [ meta, [file(row.bam)] ]

    return bam_meta


}

workflow INPUT_CHECK_BED {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_bed_channel(it) }
        .set { bam }

    emit:
    bed                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_bed_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    // add path(s) of the bam and bai file(s) to the meta map
    def bem_meta = []
   if (!file(row.bed).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bed file does not exist!\n${row.bed}"
    }
        bem_meta = [ meta, [ file(row.bed) ] ]

    return bed_meta
