//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

//
// Check input samplesheet and get read channels
//
workflow INPUT_CHECK_BAM {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_bam_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ bam, bai, bed] ]
def create_bam_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample

    // add path(s) of the bam/bai files to the meta map
    def bam_meta = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam file does not exist!\n${row.bam}"
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
        .set { bed }

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
        bed_meta = [ meta, [ file(row.bed) ] ]

    return bed_meta
}

workflow INPUT_CHECK_VARDICTJAVA {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_vardictjava_channel(it) }
        .set { vardictjava }

    emit:
    vardictjava                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_vardictjava_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    // add path(s) of the bam and bai file(s) to the meta map
    def vardictjava_meta = []
        vardictjava_meta = [ meta, file(row.bam),file(row.bai), file(row.bed) ]

    return vardictjava_meta
}

