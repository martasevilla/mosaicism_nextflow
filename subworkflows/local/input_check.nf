//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map {create_bam_bai_bed_channel(it)}
        .set { reads }

    reads_bam = reads.flatten().first().concat(reads.flatten().filter(~/.*.bam/).toList()).toList()

    emit:
    reads // channel: [ meta, [ bam, bai, bed] ]
    reads_bam // channel: [ meta, [ bam ] ]
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
