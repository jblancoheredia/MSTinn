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
    tuple val(meta), path("*.tsv")                  , emit: cutoffs
    tuple val(meta), path("*.pdf")  , optional:true , emit: pdf
    tuple val(meta), path("*.svg")  , optional:true , emit: svg
    tuple val(meta), path("*.jpeg") , optional:true , emit: jpeg
    tuple val(meta), path("*.mods")                 , emit: mods
    tuple val(meta), path("*.mbias")                , emit: mbias
    tuple val(meta), path("*_output.mods.summary")  , emit: summary
    tuple val(meta), path("*.targeted.mods.summary"), emit: targeted_summary
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
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
        ${prefix}_rastair_output.mbias)

    eval "\$cutoffs"

    rastair \\
        call \\
        --nOT \$not \\
        --nOB \$nob \\
        --fasta-file ${fasta} \\
        ${bam} \\
        >> ${prefix}_rastair_output.mods

    Summarize \\
        --rastair_in  ${prefix}_rastair_output.mods \\
        --summary_out ${prefix}_rastair_output.mods.summary

    bedtools \\
        intersect \\
        -header \\
        -u \\
        -a ${prefix}_rastair_output.mods \\
        -b ${intervals} \\
        > ${prefix}_rastair_output.targeted.mods

    Summarize \\
        --rastair_in  ${prefix}_rastair_output.targeted.mods \\
        --summary_out ${prefix}_rastair_output.targeted.mods.summary

    rastair \\
        per-read \\
        --nOT \$not \\
        --nOB \$nob \\
        --fasta-file ${fasta} \\
        ${bam} \\
        >> ${prefix}_rastair_output_perread.mods

    Summarize \\
        --rastair_in  ${prefix}_rastair_output_perread.mods \\
        --summary_out ${prefix}_rastair_output_perread.mods.summary

    bedtools \\
        intersect \\
        -header \\
        -u \\
        -a ${prefix}_rastair_output_perread.mods \\
        -b ${intervals} \\
        > ${prefix}_rastair_output_perread.targeted.mods

    Summarize \\
        --rastair_in  ${prefix}_rastair_output_perread.targeted.mods \\
        --summary_out ${prefix}_rastair_output_perread.targeted.mods.summary

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
    touch ${prefix}_rastair_output.targeted.mods
    touch ${prefix}_rastair_output.targeted.mods.summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rastair: \$(rastair --version)
    END_VERSIONS
    """
}
