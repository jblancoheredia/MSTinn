




    //
    // MODULE: Run BedTools with Intersect to subset bam file to the target region
    //
    PREP_BEDTOOLS_INTERSECT(ch_bam_dedup, ch_metbed, 'targeted')
    ch_versions = ch_versions.mix(PREP_BEDTOOLS_INTERSECT.out.versions.first())
    ch_bam_mapped_targeted = PREP_BEDTOOLS_INTERSECT.out.bam

    //
    // MODULE: Run SamTools Index
    //
    SAMTOOLS_INDEX(ch_bam_mapped_targeted)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    ch_bam_mapped_targeted_indexed = PREP_BEDTOOLS_INTERSECT.out.bam.join(SAMTOOLS_INDEX.out.bai)

    //
    // MODULE: Run rasTair in per-read mode
    //
    RASTAIR_MBIAS(ch_bam_mapped_targeted_indexed, ch_metref, ch_metfai)
    ch_versions = ch_versions.mix(RASTAIR_MBIAS.out.versions.first())

    //
    // MODULE: Run PyMbias
    //
    PYMBIAS(RASTAIR_MBIAS.out.mbias, params.pymbias_plot_type, params.pymbias_plot_ax_x, params.pymbias_plot_ax_y)
    ch_versions = ch_versions.mix(PYMBIAS.out.versions.first())
    PYMBIAS.out.cutoffs
        .map { meta, cutoff_file -> getMbiasParams(cutoff_file) }
        .set { ch_cutoffs }

    //
    // MODULE: Run rasTair
    //
    RASTAIR(ch_bam_mapped_targeted_indexed, ch_metref, ch_metfai, ch_cutoffs)
    ch_versions = ch_versions.mix(RASTAIR.out.versions.first())
    ch_rastair_mods = RASTAIR.out.mods

    //
    // MODULE: Run BedTools with Intersect to subset bam file to the target region
    //
    RASTAIR_BEDTOOLS_INTERSECT(ch_rastair_mods, ch_metbed, '_rastair_output_targeted')
    ch_versions = ch_versions.mix(RASTAIR_BEDTOOLS_INTERSECT.out.versions.first())
    ch_rastair_mods_targeted = RASTAIR_BEDTOOLS_INTERSECT.out.mods

    //
    // MODULE: Run rasTair summary for the targeted mods
    //
    RASTAIR_SUMMARY(ch_rastair_mods_targeted)
    ch_versions = ch_versions.mix(RASTAIR_BEDTOOLS_INTERSECT.out.versions.first())

    //
    // MODULE: Run astair
    //
    ASTAIR(ch_bam_mapped_targeted_indexed, ch_metref, ch_metfai, ch_metdir, params.read_length)
    ch_versions = ch_versions.mix(ASTAIR.out.versions.first())
    ch_astair_mods = ASTAIR.out.mods

    //
    // MODULE: Run BedTools with Intersect to subset bam file to the target region
    //
    ASTAIR_BEDTOOLS_INTERSECT(ch_astair_mods, ch_metbed, '_astair_output_targeted')
    ch_versions = ch_versions.mix(ASTAIR_BEDTOOLS_INTERSECT.out.versions.first())
    ch_astair_mods_targeted = ASTAIR_BEDTOOLS_INTERSECT.out.mods
