process HDR {
    label 'process_medium'

    output:
    tuple path('*.txt'), emit: hdr
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when
    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def args = task.ext.args ?: ''
    """
    echo -e '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">' >> hdr.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*python version://; s/ *\$//')
    END_VERSIONS
    """
}