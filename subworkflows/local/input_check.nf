//
// Check input samplesheet and get read, sample, and case channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow CHECK_INPUT {
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:
        SAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true, sep:',' )
            .set { sheet }
        reads        = sheet.map { create_bam_channel(it) }
        bed          = sheet.map{create_bed_channel (it)}
        vardictjava  = sheet.map{create_vardictjava_channel(it)}

    emit:
        reads           // channel: [ val(meta), [ reads ] ]
        bed
        vardictjava         // channel: [ sample_id, sex, phenotype, paternal_id, maternal_id, case_id ]
        versions  = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

def create_bed_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    // add path(s) of the bam and bai file(s) to the meta map
    def bed_meta = []
   if (!file(row.bed).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bed file does not exist!\n${row.bed}"
    }
        bed_meta = [ meta, [ file(row.bed) ] ]

    return bed_meta
}

def create_bam_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample

    // add path(s) of the bam/bai files to the meta map
    def bam_meta = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam file does not exist!\n${row.bam}"
    }
    
    if (!file(row.bai).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam.bai file does not exist!\n${row.bai}"
    }
    
    if (!file(row.bed).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bed file does not exist!\n${row.bed}"
    }
    
    bam_meta = [ meta, [file(row.bam)] ]
    //bam_meta = [ meta, [ file(row.bam)] ]
    
    return bam_meta
}

def create_vardictjava_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample

    // add path(s) of the bam/bai files to the meta map
    def vardictjava_meta = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam file does not exist!\n${row.bam}"
    }
    
    if (!file(row.bai).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bam.bai file does not exist!\n${row.bai}"
    }
    
    if (!file(row.bed).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bed file does not exist!\n${row.bed}"
    }
    
    vardictjava_meta = [ meta, file(row.bam), file(row.bai), file(row.bed)]
    //bam_meta = [ meta, [ file(row.bam)] ]
    
    return vardictjava_meta
}
