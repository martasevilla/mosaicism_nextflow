include { TABIX_BGZIP as TABIX_BGZIP} from "../../modules/nf-core/tabix/bgzip/main"
include { TABIX_BGZIP as TABIX_BGZIP_2 } from "../../modules/nf-core/tabix/bgzip/main"
include { TABIX_TABIX as TABIX_TABIX } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TABIX_2 } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TABIX_3 } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_ANNOTATE } from '../../modules/nf-core/bcftools/annotate/main'
include { FORMAT_AF_FILE } from '../../modules/local/formataf/main'
include { CREATE_HEADER } from '../../modules/local/create_header/main'
include { BCFTOOLS_QUERY } from '../../modules/nf-core/bcftools/query/main'


workflow PROCESSING_VARSCAN {

    take:
    varscan_out
    

    main:

    ch_versions = Channel.empty()


    TABIX_BGZIP (
        varscan_out
    )

   varscan_gz = TABIX_BGZIP.out.output
   ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions.first())



   // MODULE: Run tabix

   TABIX_TABIX (
       TABIX_BGZIP.out.output
   )


   ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())
   varscan_tabix=TABIX_TABIX.out.tbi
   

   ch_bcftools_query = TABIX_BGZIP.out.output.join(TABIX_TABIX.out.tbi)


   // MODULE: Run BCFTOOLS_QUERY



   BCFTOOLS_QUERY (
         ch_bcftools_query,
         [],
         [],
         []

    )
  
  query_out = BCFTOOLS_QUERY.out.output
  //ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())


    FORMAT_AF_FILE (
        query_out
    )

    //
    // MODULE: Run bgzip in af_file
    //   


    TABIX_BGZIP_2 (
        FORMAT_AF_FILE.out.af_txt
    )

   ch_versions = ch_versions.mix(TABIX_BGZIP_2.out.versions.first())


   // MODULE: Run tabix


   TABIX_TABIX_2 (
       TABIX_BGZIP_2.out.output
   )


   ch_versions = ch_versions.mix(TABIX_TABIX_2.out.versions.first())

   CREATE_HEADER(
   )

    ch_vcf=varscan_gz
    ch_index=varscan_tabix
    ch_anno=TABIX_BGZIP_2.out.output
    ch_anno_index=TABIX_TABIX_2.out.tbi
    ch_header=CREATE_HEADER.out.hdr

    ch_vcf.join(ch_index).join(ch_anno).join(ch_anno_index).combine(ch_header).set { ch_annotate_in }




  BCFTOOLS_ANNOTATE (
         ch_annotate_in
    )


   TABIX_TABIX_3 (
       BCFTOOLS_ANNOTATE.out.vcf
   )


  //ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())
  annotation_out = BCFTOOLS_ANNOTATE.out.vcf
  annotation_tabix=TABIX_TABIX_3.out.tbi

    emit:
    annotation_out
    annotation_tabix
    ch_versions


}