/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMosaicism.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta , params.fai]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Fasta file not specified!' }
if (params.fai) { ch_fasta_fai = file(params.fai) } else { exit 1, 'Fai file not specified!' }
if (params.chrom_sizes) { ch_chrom_sizes = file(params.chrom_sizes) } else { exit 1, 'chrom.sizes file not specified!' }
if (params.germline_resource) { ch_germline_resource = file(params.germline_resource) } else { exit 1, 'germline_resource VCF file not specified!' }
if (params.germline_resource_tbi) { ch_germline_resource_tbi = file(params.germline_resource_tbi) } else { exit 1, 'germline_resource TBI file not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { VARSCAN_WF } from '../subworkflows/local/varscan_workflow'
include { BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK } from '../subworkflows/local/bam_tumor_only_somatic_variant_calling_gatk'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { VARDICTJAVA } from '../modules/nf-core/vardictjava/main'
include { BEDTOOLS_INTERSECT } from '../modules/nf-core/bedtools/intersect/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MOSAICISM {

  ch_versions = Channel.empty()

  //
  // SUBWORKFLOW: Read in samplesheet, validate and stage input files
  //
  INPUT_CHECK (
      ch_input
  )
  ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

  //
  // SUBWORKFLOW: Run Samtools_sort, Samtools_mpileup and Varscan
  //

  VARSCAN_WF(

    INPUT_CHECK.out.reads_bam, INPUT_CHECK.out.reads_bed, ch_fasta

    )

  ch_versions = ch_versions.mix(VARSCAN_WF.out.ch_versions)
  
  //
  // SUBWORKFLOW: Run Mutect2 and related software
  //
  BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK(
    INPUT_CHECK.out.reads_bam_bai, Channel.fromPath(ch_fasta).first(), Channel.fromPath(ch_fasta_fai).first(), Channel.fromPath(ch_germline_resource).first(), Channel.fromPath(ch_germline_resource_tbi).first(), INPUT_CHECK.out.reads_bed
  )
  
  ch_versions = ch_versions.mix(BAM_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_GATK.out.versions.first())

  //
  // MODULE: Run Vardictjava
  //
  
  Channel.fromPath(ch_fasta)
    .map{ create_fa_meta(it) }
    .set{ ch_fasta_meta }
  
  Channel.fromPath(ch_fasta_fai)
    .map{ create_fai_meta(it) }
    .set{ ch_fai_meta }

  VARDICTJAVA (
    INPUT_CHECK.out.reads_bam_bai_bed, ch_fasta_meta, ch_fai_meta
  )

  ch_versions = ch_versions.mix(VARDICTJAVA.out.versions.first())

  //ch_versions.mix(VARDICTJAVA.out.versions.first().view())
  //ch_versions.mix(VARDICTJAVA.out.versions.view())

  ch_bedtools = VARDICTJAVA.out.vcf.join(VARSCAN_WF.out.varscan_out)
  
  Channel.fromPath(ch_chrom_sizes)
    .map{ create_chrom_sizes_input(it) }
    .set{ ch_chrom_sizes_meta }

  BEDTOOLS_INTERSECT (
    ch_bedtools, ch_chrom_sizes_meta
  )

  ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions.first())

  CUSTOM_DUMPSOFTWAREVERSIONS (
      ch_versions.unique().collectFile(name: 'collated_versions.yml')
  )

}

// Add meta for fasta input for VARDICTJAVA

def create_fa_meta (fasta) {
    // create meta map
    def meta = [:]
    meta.id = fasta.baseName
    
    // add path of fasta file to the meta map
    def chrom_sizes_meta = []
    if (!file(fasta).exists()) {
        exit 1, "ERROR: Please check input parameters -> fasta file does not exist!\n${fasta}"
    }

    fasta_meta = [ meta, file(fasta) ]

    return fasta_meta

}

// Add meta for fai input for VARDICTJAVA

def create_fai_meta (fai) {
    // create meta map
    def meta = [:]
    meta.id = fai.baseName
    
    // add path of fasta file to the meta map
    def fai_meta = []
    if (!file(fai).exists()) {
        exit 1, "ERROR: Please check input parameters -> fai file does not exist!\n${fai}"
    }

    fai_meta = [ meta, file(fai) ]

    return fai_meta

}

// Add meta for chrom.sizes input for BEDTOOLS_INTERSECT

def create_chrom_sizes_input (chrom_sizes) {
    // create meta map
    def meta = [:]
    meta.id = chrom_sizes.baseName
    
    // add path of the chrom.sizes file to the meta map
    def chrom_sizes_meta = []
    if (!file(chrom_sizes).exists()) {
        exit 1, "ERROR: Please check input parameters -> chrom_sizes file does not exist!\n${chrom_sizes}"
    }

    chrom_sizes_meta = [ meta, file(chrom_sizes) ]

    return chrom_sizes_meta

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
