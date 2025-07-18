process SAMBLASTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:60ebac4ad9c6530c0d7bf6844f52ec6916e1e0b1-0' :
        'quay.io/biocontainers/mulled-v2-19fa9f1a5c3966b63a24166365e81da35738c5ab:60ebac4ad9c6530c0d7bf6844f52ec6916e1e0b1-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    path "versions.yml"                         , emit: versions
    tuple val(meta), path("*tagged.bam")        , emit: bam
    tuple val(meta), path("*.split.bam")        , emit: split_bam
    tuple val(meta), path('*.split.fastq.gz')   , emit: split_reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$bam" == "${prefix}.bam" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools view -h $args2 $bam | \\
    samblaster $args -s ${prefix}_split.sam | \\
    samtools view $args3 -Sb - > ${prefix}.tmp.bam

    rm ${prefix}.tmp.bam

    samtools view -h $args2 $bam | \\
    samblaster $args | \\
    samtools view $args3 -Sb - > ${prefix}.bam

    samtools view ${prefix}_split.sam | \\
        awk '{
            name=\$1;
            sub(/_[12]\$/, "", name);
            sequence=\$10;
            quality=\$11;

            if (NR % 2 == 1) {
                print "@" name >> "${prefix}_R1.split.fastq";
                print sequence >> "${prefix}_R1.split.fastq";
                print "+" >> "${prefix}_R1.split.fastq";
                print quality >> "${prefix}_R1.split.fastq"
            }
            else {
                print "@" name >> "${prefix}_R2.split.fastq";
                print sequence >> "${prefix}_R2.split.fastq";
                print "+" >> "${prefix}_R2.split.fastq";
                print quality >> "${prefix}_R2.split.fastq"
            }
        }'

    gzip ${prefix}_R1.split.fastq
    gzip ${prefix}_R2.split.fastq

    samtools view -b ${prefix}_split.sam -o ${prefix}.split.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samblaster: \$( samblaster -h 2>&1 | head -n 1 | sed 's/^samblaster: Version //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$bam" == "${prefix}.bam" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.bam
    touch ${prefix}_R1.split.fastq.gz
    touch ${prefix}_R2.split.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samblaster: \$( samblaster -h 2>&1 | head -n 1 | sed 's/^samblaster: Version //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
