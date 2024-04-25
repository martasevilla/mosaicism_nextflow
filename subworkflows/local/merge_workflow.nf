include { BCFTOOLS_MERGE  } from '../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_REHEADER  } from '../../modules/nf-core/bcftools/reheader/main'
include { RENAME  } from '../../modules/local/rename/main'

workflow MERGE_WORKFLOW {

    take:
    ch_vardictjava
    ch_varscan
    ch_mutect2

    main:


    ch_versions = Channel.empty()


      RENAME(
   )

   ch_samples=RENAME.out.rename
   ch_versions = ch_versions.mix(RENAME.out.versions)




    BCFTOOLS_MERGE (
   ch_vardictjava.mix(ch_varscan).mix(ch_mutect2).groupTuple(),
    tuple([],[]), 
    tuple([],[]),
    []
  )

  ch_merge=BCFTOOLS_MERGE.out.merged_variants
  ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)





BCFTOOLS_MERGE.out.merged_variants
        .combine(ch_samples)
        .map { meta, vcf, samples ->
            [ meta, vcf, [], samples ]
        }
        .set { ch_reheader_input }




    BCFTOOLS_REHEADER(
    ch_reheader_input,
    tuple([],[])
    )


    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)



emit:
    versions            = ch_versions                                        // channel: [ versions.yml ]

}

