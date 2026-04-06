process RASTAIR_MBIAS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/mstinn_rastair:2.0.0':
        'blancojmskcc/mstinn_rastair:2.0.0' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fasta_fai)

    output:
    tuple val(meta), path("*.mbias"), emit: mbias
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp ${bam} input.bam
    cp ${bai} input.bam.bai

    rastair \\
        mbias \\
        --fasta-file ${fasta} \\
        ${bam} \\
        >> ${prefix}_rastair_output.mbias

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rastair_output.mbias

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
}
