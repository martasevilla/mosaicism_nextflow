//
// Run GATK mutect2 in tumor only mode, getepileupsummaries, calculatecontamination and filtermutectcalls
//

include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/gatk4/createsequencedictionary'
include { GATK4_MUTECT2                as MUTECT2 }                  from '../../modules/nf-core/gatk4/mutect2/main'
// include { GATK4_GETPILEUPSUMMARIES     as GETPILEUPSUMMARIES }       from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
// include { GATK4_CALCULATECONTAMINATION as CALCULATECONTAMINATION }   from '../../../modules/nf-core/gatk4/calculatecontamination/main'
// include { GATK4_FILTERMUTECTCALLS      as FILTERMUTECTCALLS }        from '../../../modules/nf-core/gatk4/filtermutectcalls/main'

workflow BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK {
    take:
    ch_input                 // channel: [ val(meta), [ input ], [ input_index ], [] ]
    ch_fasta                 // channel: /path/to/reference/fasta
    ch_fai                   // channel: /path/to/reference/fasta/index
    ch_germline_resource     // channel: /path/to/germline/resource
    ch_germline_resource_tbi // channel: /path/to/germline/index
    ch_interval_file         // channel: /path/to/interval/file
    // ch_dict                  // channel: /path/to/reference/fasta/dictionary
    // ch_panel_of_normals      // channel: /path/to/panel/of/normals
    // ch_panel_of_normals_tbi  // channel: /path/to/panel/of/normals/index
    

    main:
    ch_versions = Channel.empty()


    // Channel for FASTA, FAI and DICT with META
    ch_fasta
        .map{ create_fasta_input(it) }
        .set{ ch_fasta_meta }

    ch_fai
        .map{ create_fai_input(it) }
        .set{ ch_fai_meta }

    GATK4_CREATESEQUENCEDICTIONARY (
        ch_fasta_meta
    )

    ch_dict_meta = GATK4_CREATESEQUENCEDICTIONARY.out.dict

    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    //
    //Perform variant calling using mutect2 module in tumor single mode.
    //

    ch_input
        .combine(ch_interval_file, by: 0)
        .set{ch_input_mutect2}

    MUTECT2 (
        ch_input_mutect2,
        ch_fasta_meta,
        ch_fai_meta,
        ch_dict_meta,
        ch_germline_resource,
        ch_germline_resource_tbi,
        // ch_panel_of_normals,
        // ch_panel_of_normals_tbi
    )

    ch_versions = ch_versions.mix(MUTECT2.out.versions)
/*
    //
    //Generate pileup summary table using getpileupsummaries.
    //

    ch_pileup_input = ch_input.combine(ch_interval_file).map {
        meta, input_file, input_index, which_norm, intervals ->
        [meta, input_file, input_index, intervals]
    }

    GETPILEUPSUMMARIES (
        ch_pileup_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi
    )

    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES.out.versions)

    //
    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
    //
    ch_pileup = GETPILEUPSUMMARIES.out.table.collect()

    //[] is a placeholder for the optional input where the matched normal sample would be passed in for tumor-normal samples, which is not necessary for this workflow.
    ch_pileup.add([])

    CALCULATECONTAMINATION ( ch_pileup )

    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)

    //
    //Mutect2 calls filtered by filtermutectcalls using the contamination and segmentation tables.
    //

    ch_vcf = MUTECT2.out.vcf.collect()
    ch_tbi = MUTECT2.out.tbi.collect()
    ch_stats = MUTECT2.out.stats.collect()

    ch_segment = CALCULATECONTAMINATION.out.segmentation.collect()
    ch_contamination = CALCULATECONTAMINATION.out.contamination.collect()

    ch_filtermutect_in = ch_vcf
        .combine(ch_tbi, by: 0)
        .combine(ch_stats, by: 0)
        .combine(ch_segment, by: 0)
        .combine(ch_contamination, by: 0)
    // Adding [] as a placeholder for the optional input file artifact priors, which is only used for tumor-normal samples and therefor isn't needed in this workflow.
    // and [] as a placeholder for entering a contamination estimate value, which is not needed as this workflow uses the contamination table instead.
        .map{ meta, vcf, tbi, stats, segment, contamination -> [meta, vcf, tbi, stats, [], segment, contamination, [] ] }

    ch_filtermutect_in.view()

    FILTERMUTECTCALLS ( ch_filtermutect_in, ch_fasta, ch_fai, ch_dict )
    ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)
*/
    emit:
    mutect2_vcf         = MUTECT2.out.vcf.collect()                          // channel: [ val(meta), [ vcf ] ]
    mutect2_index       = MUTECT2.out.tbi.collect()                          // channel: [ val(meta), [ tbi ] ]
    mutect2_stats       = MUTECT2.out.stats.collect()                        // channel: [ val(meta), [ stats ] ]
    versions            = ch_versions                                        // channel: [ versions.yml ]
/*
    pileup_table        = GETPILEUPSUMMARIES.out.table.collect()             // channel: [ val(meta), [ table ] ]

    contamination_table = CALCULATECONTAMINATION.out.contamination.collect() // channel: [ val(meta), [ contamination ] ]
    segmentation_table  = CALCULATECONTAMINATION.out.segmentation.collect()  // channel: [ val(meta), [ segmentation ] ]

    filtered_vcf        = FILTERMUTECTCALLS.out.vcf.collect()                // channel: [ val(meta), [ vcf ] ]
    filtered_index      = FILTERMUTECTCALLS.out.tbi.collect()                // channel: [ val(meta), [ tbi ] ]
    filtered_stats      = FILTERMUTECTCALLS.out.stats.collect()              // channel: [ val(meta), [ stats ] ]
*/
}

def create_fasta_input (fasta) {
    // create meta map
    def meta = [:]
    meta.id = fasta.baseName

    // add path of the fasta file to the meta map
    def fasta_meta = []
    if (!file(fasta).exists()) {
        exit 1, "ERROR: Please check input parameters -> Fasta file does not exist!\n${fasta}"
    }

    fasta_meta = [ meta, file(fasta) ]

    return fasta_meta

}

def create_fai_input (fai) {
    // create meta map
    def meta = [:]
    meta.id = fai.baseName

    // add path of the fai file to the meta map
    def fai_meta = []
    if (!file(fai).exists()) {
        exit 1, "ERROR: Please check input parameters -> Fai file does not exist!\n${fai}"
    }

    fai_meta = [ meta, file(fai) ]

    return fai_meta
}

