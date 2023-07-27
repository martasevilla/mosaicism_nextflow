process POSTPROCESSING {
    tag "$meta.id"
    label 'process_low'
    conda "bioconda::bcftools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'quay.io/biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta), path(query)

    output:
   
    tuple val(meta), path("*.txt"), emit: postprocessing
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    //awk {'print $1"\t"$2"\t"$3/$4'} ${query} |sed 's/\,/./' > ${prefix}_af.txt
    """
    
    awk $args ${query} | $args2 > ${prefix}_af.txt
    
    cat <<-END_VERSIONS > versions.yml
    

    END_VERSIONS
    """
}
