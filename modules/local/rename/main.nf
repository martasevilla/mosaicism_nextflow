process RENAME {
    label 'process_medium'

    output:
    tuple path('*.txt'), emit: rename
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when
    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def args = task.ext.args ?: ''
    """
    echo -e 'Vardictjava\nVarscan\nMutect2' > rename.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/^.*python version://; s/ *\$//')
    END_VERSIONS
    """
}