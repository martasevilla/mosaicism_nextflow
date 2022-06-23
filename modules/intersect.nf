#!/usr/bin/env nextflow

/* 
 * Script to intersect the two vcf outputs from vardict and varsome
 */

nextflow.enable.dsl=2

// defaults
prefix = "out"
mergedVCF = "merged-file"
params.outdir = ""
params.cpus = 1

process intersectVCF {
  /*
  Function to intersect

  Returns
  -------
  Returns 2 files:
      1) A VCF format file 
      2) A tabix index for that VCF
  */

  publishDir "${params.outdir}", 
    enabled: "${params.outdir}" as Boolean,
    mode:'copy'
    
  cpus params.cpus
  #ESTO NO NO? ------> container "quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0"
   
  input:
  path(vcfFiles)
  path(indexFiles)

  output:
  path("${ mergedVCF }.vcf.gz"), emit: vcfFile                                     
  path("${ mergedVCF }.vcf.gz.tbi"), emit: indexFile

  script: 
  """
  #bcftools concat ${ vcfFiles } -Oz -o temp-${ mergedVCF}.vcf.gz
  #bcftools sort -Oz temp-${ mergedVCF}.vcf.gz -o ${ mergedVCF}.vcf.gz 
  #bcftools  index -t ${ mergedVCF}.vcf.gz


  bcftools isec -p path("${ mergedVCF }.vcf.gz") -n=2 -w1 ${ vcfFiles } -o ${ mergedVCF}.vcf.gz -Oz


  """
}


