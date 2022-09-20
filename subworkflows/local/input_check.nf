//
// Check input samplesheet and get read, sample, and case channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

<<<<<<< HEAD
workflow INPUT_CHECK {
=======
workflow CHECK_INPUT {
>>>>>>> main
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:
<<<<<<< HEAD
    canal_sp = SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' ) // canal donde cada elemento es una fila del csv

    canal_sp
        .map {create_bam_bai_bed_channel(it)}
        .set { reads_bam_bai_bed }

    canal_sp
        .map {create_bam_channel(it)}
        .set { reads_bam }

    canal_sp
        .map {create_bed_channel(it)}
        .set { reads_bed }

    //reads_bam = reads.flatten().first().concat(reads.flatten().filter(~/.*.bam/).toList()).toList()

    emit:
    reads_bam_bai_bed // channel: [ meta, [ bam, bai, bed] ]
    reads_bam // channel: [ meta, [ bam ] ]
    reads_bed // channel: [ meta, [ bed ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]

}

//Function to get list of [ meta, [ bam, bai, bed] ]

def create_bam_bai_bed_channel(LinkedHashMap row) {
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

    bam_meta = [ meta, file(row.bam), file(row.bai), file(row.bed) ]
    //bam_meta = [ meta, [ file(row.bam)] ]

    return bam_meta
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

    bam_meta = [ meta, file(row.bam) ]
    //bam_meta = [ meta, [ file(row.bam)] ]

    return bam_meta
}

def create_bed_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample

    // add path(s) of the bam/bai files to the meta map

    if (!file(row.bed).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Bed file does not exist!\n${row.bed}"
    }

    bam_meta = [ meta, file(row.bed) ]
    //bam_meta = [ meta, [ file(row.bam)] ]

    return bam_meta
}
=======
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
>>>>>>> main
