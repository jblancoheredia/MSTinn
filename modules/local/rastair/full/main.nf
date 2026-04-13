process RASTAIR_FULL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/mstinn_rastair:2.0.0':
        'blancojmskcc/mstinn_rastair:2.0.0' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fasta_fai)
    tuple val(meta3), path(intervals)
    val(plot_type)
    val(plot_ax_x)
    val(plot_ax_y)

    output:
    tuple val(meta), path("*.pdf")  , optional:true , emit: pdf
    tuple val(meta), path("*.png")  , optional:true , emit: png
    tuple val(meta), path("*.tsv")  , optional:true , emit: tsv
    tuple val(meta), path("*.svg")  , optional:true , emit: svg
    tuple val(meta), path("*.jpeg") , optional:true , emit: jpeg
    tuple val(meta), path("*.mods")                 , emit: mods
    tuple val(meta), path("*.mbias")                , emit: mbias
    tuple val(meta), path("*cutoffs.tsv")           , emit: cutoffs
    tuple val(meta), path("*_output.summary")       , emit: summary
    tuple val(meta), path("*.targeted.summary")     , emit: targeted_summary
    tuple val(meta), path("*_output_perread.mods")  , emit: perread
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$(mktemp -d)

    rastair \\
        mbias \\
        --fasta-file ${fasta} \\
        ${bam} \\
        >> ${prefix}_rastair_output.mbias

    cutoffs=\$(PyMbias \\
        -p ${prefix} \\
        -s ${plot_type} \\
        -x ${plot_ax_x} \\
        -y ${plot_ax_y} \\
        -i ${prefix}_rastair_output.mbias)

    eval "\$cutoffs"

    rastair \\
        call \\
        -@ ${task.cpus} \\
        --nOT \$not \\
        --nOB \$nob \\
        --fasta-file ${fasta} \\
        ${bam} \\
        >> ${prefix}_rastair_output.mods

    Summarize \\
        --input-file ${prefix}_rastair_output.mods \\
        --output-dir .

    PyGCbias \\
        -m ${prefix}_rastair_output.mods \\
        -r ${fasta} \\
        -o . \\
        -p ${prefix} \\
        --export-beta-table

    PyCpGbias \\
        -m ${prefix}_rastair_output.mods \\
        -r ${fasta} \\
        -o . \\
        -p ${prefix} \\
        --export-beta-table

    bedtools \\
        intersect \\
        -header \\
        -u \\
        -a ${prefix}_rastair_output.mods \\
        -b ${intervals} \\
        > ${prefix}_rastair_output.targeted.mods

    Summarize \\
        --input-file ${prefix}_rastair_output.targeted.mods \\
        --output-dir .

    rastair \\
        per-read \\
        -@ ${task.cpus} \\
        --fasta-file ${fasta} \\
        ${bam} \\
        >> ${prefix}_rastair_output_perread.mods

    PyVbias \\
        -i ${prefix}_rastair_output_perread.mods \\
        -o . \\
        -p ${prefix} \\
        --export-vdata

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}_cutoffs.tsv
    touch ${prefix}_rastair_output.mods
    touch ${prefix}_rastair_output.mbias
    touch ${prefix}_rastair_output.mods.summary
    touch ${prefix}_rastair_output_perread.mods
    touch ${prefix}_rastair_output.targeted.mods
    touch ${prefix}_rastair_output.targeted.mods.summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
}
