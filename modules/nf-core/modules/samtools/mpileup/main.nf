process SAMTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"
    input:
    tuple val(meta), path(sort)
    //tuple val(meta), path(reads), path(intervals)
    //path  fasta

    output:
    tuple val(meta), path("*.mpileup"), emit: mpileup
    //tuple val(meta), path("*.mpileup"), emit: mpileup
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // esto de momento no lo pongo def intervals = intervals ? "-l ${intervals}" : ""
    """
    samtools mpileup \\
        --no-BAQ --fasta-ref /home/msevilla/Escritorio/mosaicism_nextflow-main/hg19.fa \\
        --output ${prefix}.mpileup \\
        ${sort}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
    """
    

}


    /*
    //bgzip ${prefix}.mpileup
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    */