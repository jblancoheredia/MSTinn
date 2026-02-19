process FGBIO_CORRECTUMIS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fe/fe9479adc5e6e0a1c125d346fdfa0dd313834249e9c55c40e8d44ec3a48c6559/data' :
        'community.wave.seqera.io/library/fgbio:3.1.1--6c9a88faf1d62b6c' }"

    input:
    tuple val(meta), path(bam)
    val(revcomp)
    val(max_mismatch)
    val(min_distance)
    val(min_corrected)

    output:
    tuple val(meta), path("*.mapped.fixedUMI.bam"), emit: bam
    tuple val(meta), path("*.metrics.txt"),         emit: metrics
    path "versions.yml",                            emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio CorrectUmis] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        CorrectUmis \\
        -i ${bam} \\
        --revcomp=${revcomp} \\
        -o ${prefix}.mapped.fixedUMI.bam \\
        --max-mismatches=${max_mismatch} \\
        --min-distance=${min_distance} \\
        --min-corrected=${min_corrected} \\
        -M ${prefix}.metrics.txt \\
        -r ${prefix}.rejected.bam \\
        -t RX \\
        -u \\
        GAGACGAT TTCCAAGG CGCATGAT ACGGAACA \\
        CGGCTAAT GCTATCCT TGGACTCT ATCCAGAG \\
        CTTAGGAC GTGCCATA TCGCTGTT TTCGTTGG \\
        AAGCACTG GTCGAAGA ACCACGAT GATTACCG \\
        GCACAACT GCGTCATT GAAGGAAG ACTGAGGT \\
        TGAAGACG GTTACGCA AGCGTGTT GATCGAGT \\
        TTGCGAAG CTGTTGAC GATGTGTG ACGTTCAG \\
        TTGCAGAC CAATGTGG ACGACTTG ACTAGGAG \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.mapped.fixedUMI.bam
    touch  ${prefix}.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
