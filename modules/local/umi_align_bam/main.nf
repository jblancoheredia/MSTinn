process ALIGN_BAM_RAW {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta),  path(unmapped_bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(bwa_dict)
    tuple val(meta5), path(bwa_dir)
    val sort_type

    output:
    tuple val(meta), path("*.mapped.bam")       , emit: bam
    tuple val(meta), path("*.sorted.bam")       , emit: sort_bam
    tuple val(meta), path("*.sorted.bam.bai")   , emit: sort_bai
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samtools_fastq_args = task.ext.samtools_fastq_args ?: ''
    def samtools_sort_args = task.ext.samtools_sort_args ?: ''
    def bwa_args = task.ext.bwa_args ?: ''
    def fgbio_args = task.ext.fgbio_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fgbio_mem_gb = 4
    def extra_command = ""

    if (!task.memory) {
        log.info '[fgbio ZipperBams] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else if (fgbio_mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            fgbio_mem_gb = 1
        } else {
            fgbio_mem_gb = task.memory.giga - 1
        }
    }

    fgbio_zipper_bams_output = "/dev/stdout"
    fgbio_zipper_bams_compression = 0
    extra_command = " | samtools sort -n"
    extra_command += " --threads "+ task.cpus
    extra_command += " -o " + prefix + ".mapped.bam "
    extra_command += " -"

    """
    # The real path to the BWA index prefix`
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools fastq ${samtools_fastq_args} ${unmapped_bam} \\
        | bwa-mem2 mem ${bwa_args} -t ${task.cpus} -p -B 3 -K 100000000 -Y -M -R ${meta.read_group} \$BWA_INDEX_PREFIX - \\
        | fgbio -Xmx${fgbio_mem_gb}g \\
            --compression ${fgbio_zipper_bams_compression} \\
            --async-io=true \\
            ZipperBams \\
            --unmapped ${unmapped_bam} \\
            --ref ${fasta} \\
            --output ${fgbio_zipper_bams_output} \\
            --tags-to-reverse Consensus \\
            --tags-to-revcomp Consensus \\
            ${fgbio_args} \\
            ${extra_command};

    samtools sort  -@ ${task.cpus} ${prefix}.mapped.bam -o ${prefix}.sorted.bam
    samtools index -@ ${task.cpus} ${prefix}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mapped.bam
    touch ${prefix}.sorted.bam
    touch ${prefix}.sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

process ALIGN_BAM_CON {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::fgbio=2.4.0 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://blancojmskcc/umi_aligner:1.0.0' :
        'blancojmskcc/umi_aligner:1.0.0' }"

    input:
    tuple val(meta),  path(bam), path(bai), path(duplex_bam), path(duplex_bai), path(simplex_bam), path(simplex_bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(bwa_dir)

    output:
    tuple val(meta), path("*.mapped.bam"),          emit: bam
    tuple val(meta), path("*.mapped.duplex.bam"),   emit: duplex_bam
    tuple val(meta), path("*.mapped.simplex.bam"),  emit: simplex_bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samtools_fastq_args = task.ext.samtools_fastq_args ?: ''
    def samtools_sort_args = task.ext.samtools_sort_args ?: ''
    def bwa_args = task.ext.bwa_args ?: ''
    def fgbio_args = task.ext.fgbio_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fgbio_mem_gb = 4
    def extra_command = ""

    if (!task.memory) {
        log.info '[fgbio ZipperBams] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else if (fgbio_mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            fgbio_mem_gb = 1
        } else {
            fgbio_mem_gb = task.memory.giga - 1
        }
    }
    """
    # The real path to the BWA index prefix`
    BWA_INDEX_PREFIX=\$(find -L ./ -name "*.amb" | sed 's/.amb//')

    samtools fastq ${samtools_fastq_args} ${bam} \\
        | bwa-mem2 mem ${bwa_args} -t ${task.cpus} -p -K 150000000 -Y \$BWA_INDEX_PREFIX - \\
        | fgbio -Xmx${fgbio_mem_gb}g \\
            --compression 1 \\
            --async-io=true \\
            ZipperBams \\
            --unmapped ${bam} \\
            --ref ${fasta} \\
            --output ${prefix}.mapped.bam \\
            --tags-to-reverse Consensus \\
            --tags-to-revcomp Consensus \\
            --sort TemplateCoordinate

    samtools fastq ${samtools_fastq_args} ${duplex_bam} \\
        | bwa-mem2 mem ${bwa_args} -t ${task.cpus} -p -K 150000000 -Y \$BWA_INDEX_PREFIX - \\
        | fgbio -Xmx${fgbio_mem_gb}g \\
            --compression 1 \\
            --async-io=true \\
            ZipperBams \\
            --unmapped ${duplex_bam} \\
            --ref ${fasta} \\
            --output ${prefix}.mapped.duplex.bam \\
            --tags-to-reverse Consensus \\
            --tags-to-revcomp Consensus \\
            --sort TemplateCoordinate

    samtools fastq ${samtools_fastq_args} ${simplex_bam} \\
        | bwa-mem2 mem ${bwa_args} -t ${task.cpus} -p -K 150000000 -Y \$BWA_INDEX_PREFIX - \\
        | fgbio -Xmx${fgbio_mem_gb}g \\
            --compression 1 \\
            --async-io=true \\
            ZipperBams \\
            --unmapped ${simplex_bam} \\
            --ref ${fasta} \\
            --output ${prefix}.mapped.simplex.bam \\
            --tags-to-reverse Consensus \\
            --tags-to-revcomp Consensus \\
            --sort TemplateCoordinate

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mapped.bam
    touch ${prefix}.mapped.duplex.bam
    touch ${prefix}.mapped.simplex.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | tail -n 1)
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}