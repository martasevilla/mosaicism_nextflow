//
// Ejecuta todos los pasos necesarios para usar Varscan (samtools_sort, samtools_mpileup, unzip, varscan)
//



include { BCFTOOLS_ANNOTATE      } from '../../modules/nf-core/modules/bcftools/annotate/main'
include { BCFTOOLS_QUERY      } from '../../modules/nf-core/modules/bcftools/query/main'
include { TABIX_TABIX      } from '../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_BGZIP      } from '../../modules/nf-core/modules/tabix/bgzip/main'
include { POSTPROCESSING      } from '../../modules/local/postprocessing.nf'

workflow NORMALIZE {
    take:
    
    varscan_out // file: [ [ meta ], bam ] file, fasta

    main:

    ch_versions = Channel.empty()
    
    
    TABIX_BGZIP(

        varscan_out
    )
    
    ch_varscan_gz = TABIX_BGZIP.out.output
    
    
    TABIX_TABIX (
    	ch_varscan_gz
    )
    
    ch_varscan_tbi = TABIX_TABIX.out.tbi
    
    
    BCFTOOLS_QUERY (
    	ch_varscan_gz.join(ch_varscan_tbi), 
        [],
        [],
        []
        
    )

    ch_query=BCFTOOLS_QUERY.out.output
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())


    POSTPROCESSING (
    	ch_query
    )

    ch_post=POSTPROCESSING.out.postprocessing

    //TABIX_BGZIP(

      //  ch_post
    //)

    //ch_post_bgzip = TABIX_BGZIP.out.output
    
    //TABIX_TABIX (
    //	ch_post_bgzip
    //)
    //ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions.first())

 
 
 //path(input), path(index), path(annotations), path(annotations_index), path(header_lines)
 
 
    //ch_vcf
     //       .join(ch_index)
       //     .join(CADD.out.tsv)
         //   .join(TABIX_CADD.out.tbi)
           // .combine(ch_header)
            //.set { ch_annotate_in }




//bcftools annotate -s Sample1 -a af.txt.gz -h hdr.txt -c CHROM,POS,FORMAT/AF varscan.vcf.gz


    //BCFTOOLS_ANNOTATE (
    //	VARSCAN.out.vcf_varscan,
    //)
    //ch_versions = ch_versions.mix(SAMTOOLS_MPILEUP.out.versions.first())

}
