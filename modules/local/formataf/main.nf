process FORMAT_AF_FILE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(txt)

    output:
    tuple val(meta), path('*.txt'), emit: af_txt
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
     cat ${txt}| awk -v OFS='\t' {'print \$1, \$2, \$3/\$4'} > af_format.txt
     sed 's/,/./g' af_format.txt > af_format.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*python version://; s/ *\$//')
    END_VERSIONS
    """
}