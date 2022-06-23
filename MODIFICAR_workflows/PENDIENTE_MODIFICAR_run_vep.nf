/* 
 * Workflow to run VEP on chromosome based VCF files
 *
 * This workflow relies on Nextflow (see https://www.nextflow.io/tags/workflow.html)
 *
 */

nextflow.enable.dsl=2

 // params default
params.help = false
params.cpus = 1
params.outdir = "outdir"


 // print usage
if (params.help) {
  log.info ''
  log.info 'Pipeline to intersect'
  log.info '-------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info '  nextflow -C nf_config/nextflow.config run workflows/run_vep.nf --vcf <path-to-vcf> --chros 1,2 --vep_config'
  log.info ''
  log.info 'Options:'
  log.info '  --vcf VCF                 VCF that will be intersect. Currently supports sorted and bgzipped file'
  log.info '  --outdir DIRNAME          Name of output dir. Default: outdir'
  log.info '  --cpus INT	        Number of CPUs to use. Default 1.'
  exit 1
}

###############################falta esto!

// Input validation

if( !params.vcf) {
  exit 1, "Undefined --vcf parameter. Please provide the path to a VCF file"
}

vcfFile = file(params.vcf)
if( !vcfFile.exists() ) {
  exit 1, "The specified VCF file does not exist: ${params.vcf}"
}

##############################esto se podria quedar????? 
check_bgzipped = "bgzip -t $params.vcf".execute()
check_bgzipped.waitFor()
if(check_bgzipped.exitValue()){
  exit 1, "The specified VCF file is not bgzipped: ${params.vcf}"
}


############3esto!!!!!???
//def sout = new StringBuilder(), serr = new StringBuilder()
//check_parsing = "$params.singularity_dir/vep.sif tabix -p vcf -f $params.vcf".execute()
//check_parsing.consumeProcessOutput(sout, serr)
//check_parsing.waitFor()
//if( serr ){
//  exit 1, "The specified VCF file has issues in parsing: $serr"
//}
vcf_index = "${params.vcf}.tbi"


log.info 'Starting workflow.....'


##########################FALTA
workflow run_vep {
  take: vep_config
  main:
    chr_str = params.chros.toString()
    chr = Channel.of(chr_str.split(',')).view()
    splitVCF(chr, file( params.vcf ), file( vcf_index ))
    chrosVEP(splitVCF.out, vep_config)
    mergeVCF(chrosVEP.out.vcfFile.collect(), chrosVEP.out.indexFile.collect())
  emit:
    vcfFile = mergeVCF.out.vcfFile
    indexFile = mergeVCF.out.indexFile
}

workflow {
  run_vep( file( params.vep_config ) )
}
