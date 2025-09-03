process FGBIO_CLIPBAM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*_clipped.bam"), path("*_clipped.bai")   , emit: bam
    tuple val(meta), path("*.metrics.txt")                          , emit: metrics
    path "versions.yml"                                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio ClipBam] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }

    """
    samtools \\
        sort \\
        -n \\
        -l 0 \\
        ${bam} |\\
    fgbio \\
        -Xmx${mem_gb}g \\
        ClipBam \\
        -i /dev/stdin \\
        -o ${prefix}_clipped.bam \\
        -r ${fasta} \\
        -m ${prefix}.metrics.txt \\
        -S coordinate \\
        --clip-overlapping-reads=true \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}_clipped.bam
    touch  ${prefix}.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
