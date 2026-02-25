/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT PLUGINS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                              } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                                                        } from '../modules/nf-core/fastqc/main'
include { SPADES                                                        } from '../modules/local/spades/main'
include { ASTAIR                        	                            } from '../modules/local/astair/main'
include { MULTIQC                                                       } from '../modules/nf-core/multiqc/main'
include { BWAMEM2                                                       } from '../modules/local/bwamem2/main'
include { BWA_METH                 	                                    } from '../modules/local/bwameth/main'
include { MOSDEPTH                                                      } from '../modules/local/mosdepth/main'
include { CAT_FASTQ                     	                            } from '../modules/nf-core/cat/fastq/main'
include { SAMBLASTER                                                    } from '../modules/local/samblaster/main'
include { MOSDEPTH_DUP                                                  } from '../modules/local/mosdepth/main'
include { MOSDEPTH_CON                                                  } from '../modules/local/mosdepth/main'
include { MOSDEPTH_RAW                                                  } from '../modules/local/mosdepth/main'
include { MOSDEPTH_SIM                                                  } from '../modules/local/mosdepth/main'
include { RASTAIR_FULL                 	                                } from '../modules/local/rastair/full/main'
include { ALIGN_BAM_CON                                                 } from '../modules/local/umi_align_bam/main'
include { ALIGN_BAM_RAW                                                 } from '../modules/local/umi_align_bam/main'
include { FGBIO_CLIPBAM                                                 } from '../modules/local/fgbio/clipbam/main'
include { PRESEQ_CCURVE                                                 } from '../modules/local/preseq/ccurve/main'
include { FILTER_CONTIGS                                                } from '../modules/local/filter_contigs/main'
include { SAMTOOLS_INDEX                	                            } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS                                                } from '../modules/local/samtools/stats/main'
include { PRESEQ_LCEXTRAP                                               } from '../modules/local/preseq/lcextrap/main'
include { SENTIEON_BWAMEM                                               } from '../modules/nf-core/sentieon/bwamem/main'
include { UMI_READ_COUNTS                                               } from '../modules/local/umi_read_counts/main'
include { COLLECTHSMETRICS                                              } from '../modules/local/picard/collecthsmetrics/main'
include { FASTQC_CONSENSUS                                              } from '../modules/local/fastqc_consensus/main'
include { FGBIO_FASTQTOBAM                                              } from '../modules/nf-core/fgbio/fastqtobam/main'
include { FGBIO_SORTCONBAM                                              } from '../modules/local/fgbio/sortconbam/main.nf'
include { MSISENSORPRO_CON                                              } from '../modules/local/msisensorpro/pro/main'
include { MSISENSORPRO_RAW                                              } from '../modules/local/msisensorpro/pro/main'
include { FGBIO_CORRECTUMIS                                             } from '../modules/local/fgbio/correctumis/main'
include { COLLECT_UMI_METRICS                                           } from '../modules/local/collect_umi_metrics/main'
include { DOWNSAMPLINGS_COUNT                                           } from '../modules/local/downsamplings/count'
include { DOWNSAMPLINGS_SEQTK                                           } from '../modules/local/downsamplings/seqtk'
include { COLLECTHSMETRICS_CON                                          } from '../modules/local/picard/collecthsmetrics/main'
include { COLLECTHSMETRICS_DUP                                          } from '../modules/local/picard/collecthsmetrics/main'
include { COLLECTHSMETRICS_RAW                                          } from '../modules/local/picard/collecthsmetrics/main'
include { COLLECTHSMETRICS_SIM                                          } from '../modules/local/picard/collecthsmetrics/main'
include { GATK4_MARKDUPLICATES          	                            } from '../modules/local/gatk4/markduplicates/main'
include { FGBIO_GROUPREADSBYUMI                                         } from '../modules/local/fgbio/groupreadsbyumi/main'
include { GATK4_HAPLOTYPECALLER                                         } from '../modules/local/gatk4/haplotypecaller/main'
include { SAMTOOLS_COLLATEFASTQ                                         } from '../modules/nf-core/samtools/collatefastq/main'
include { GATK4_VARIANTFILTRATION                                       } from '../modules/local/gatk4/variantfiltration/main'
include { PREP_BEDTOOLS_INTERSECT       	                            } from '../modules/local/bedtools/prep_bedtools_intersect' 
include { SAMTOOLS_SORT_INDEX_CON                                       } from '../modules/local/samtools/sort_index/main'
include { SAMTOOLS_SORT_INDEX_RAW                                       } from '../modules/local/samtools/sort_index/main'
include { ASTAIR_BEDTOOLS_INTERSECT     	                            } from '../modules/local/bedtools/astair_bedtools_intersect' 
include { FGBIO_FILTERCONSENSUSREADS                                    } from '../modules/local/fgbio/filterconsensusreads/main'
include { FGBIO_COLLECTDUPLEXSEQMETRICS                                 } from '../modules/local/fgbio/collectduplexseqmetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS                                 } from '../modules/local/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS as COLLECTMULTIPLEMETRICS       } from '../modules/local/picard/collectmultiplemetrics/main'
include { FGBIO_CALLDUPLEXCONSENSUSREADS                                } from '../modules/local/fgbio/callduplexconsensusreads/main'
include { FGBIO_ERRORRATEBYREADPOSITION_CON                             } from '../modules/local/fgbio/errorratebyreadposition/main'
include { FGBIO_ERRORRATEBYREADPOSITION_RAW                             } from '../modules/local/fgbio/errorratebyreadposition/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                 IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMultiqc                                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                        } from '../subworkflows/local/utils_nfcore_mstinn_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                            CREATE CUSTOM CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_msi                                              = Channel.fromPath(params.msi).map                          { it -> [[id:it.Name], it] }.collect()
ch_metbed                                           = Channel.fromPath(params.metbed).map                       { it -> [[id:it.Name], it] }.collect()
ch_metdct                                           = Channel.fromPath(params.metdct).map                       { it -> [[id:it.Name], it] }.collect()
ch_metdir                                           = Channel.fromPath(params.metdir).map                       { it -> [[id:it.Name], it] }.collect()
ch_metref                                           = Channel.fromPath(params.metref).map                       { it -> [[id:it.Name], it] }.collect()
ch_metfai                                           = Channel.fromPath(params.metfai).map                       { it -> [[id:it.Name], it] }.collect()
ch_bwadct                                           = Channel.fromPath(params.bwadct).map                       { it -> [[id:it.Name], it] }.collect()
ch_bwadir                                           = Channel.fromPath(params.bwadir).map                       { it -> [[id:it.Name], it] }.collect()
ch_bwaref                                           = Channel.fromPath(params.bwaref).map                       { it -> [[id:it.Name], it] }.collect()
ch_bwafai                                           = Channel.fromPath(params.bwafai).map                       { it -> [[id:it.Name], it] }.collect()
ch_intervals                                        = Channel.fromPath(params.intervals).map                    { it -> [[id:it.Name], it] }.collect()
ch_known_sites                                      = Channel.fromPath(params.known_sites).map                  { it -> [[id:it.Name], it] }.collect()
ch_metdictgatk                                      = Channel.fromPath(params.metdictgatk).map                  { it -> [[id:it.Name], it] }.collect()
ch_known_sites_tbi                                  = Channel.fromPath(params.known_sites_tbi).map              { it -> [[id:it.Name], it] }.collect()
ch_intervals_gunzip                                 = Channel.fromPath(params.intervals_bed_gunzip).map         { it -> [[id:it.Name], it] }.collect()
ch_intervals_gunzip_index                           = Channel.fromPath(params.intervals_bed_gunzip_index).map   { it -> [[id:it.Name], it] }.collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                            RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2
workflow MSTINN {

    take:
    ch_samplesheet
    
    main:
    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC(ch_samplesheet)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // Create channel branch from validated samplesheet
    //
    ch_samplesheet
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }
        .set { ch_fastq }
    ch_fastq_single     = ch_fastq.single.map   {meta, fastqs -> addReadgroupToMeta(meta, fastqs)}
    ch_fastq_multiple   = ch_fastq.multiple.map {meta, fastqs -> addReadgroupToMeta(meta, fastqs)}

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (ch_fastq_multiple)
    .reads
    .mix(ch_fastq_single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    //
    // MODULE: Run in-house script for counting reads
    //
    DOWNSAMPLINGS_COUNT(ch_cat_fastq)
    ch_versions = ch_versions.mix(DOWNSAMPLINGS_COUNT.out.versions)
    ch_global_min_reads = DOWNSAMPLINGS_COUNT.out.total_reads
        .map { file -> 
            def count = file.text.trim()
            count.toInteger()
        }
        .collect()
        .map { counts -> 
            def min_count = counts.min()
            min_count
        }

    if (params.run_downsamplings) {

        if (params.downsampling_total_reads) {

            //
            // MODULE: Run Downsampling with seqtk
            //
            DOWNSAMPLINGS_SEQTK(ch_cat_fastq, params.downsampling_total_reads)
            ch_versions = ch_versions.mix(DOWNSAMPLINGS_SEQTK.out.versions)
            ch_downsampled_reads = DOWNSAMPLINGS_SEQTK.out.downsampled_reads

            ch_fastqs = ch_downsampled_reads

        } else {

            //
            // MODULE: Run Downsampling with seqtk
            //
            DOWNSAMPLINGS_SEQTK(ch_cat_fastq, ch_global_min_reads)
            ch_versions = ch_versions.mix(DOWNSAMPLINGS_SEQTK.out.versions)
            ch_downsampled_reads = DOWNSAMPLINGS_SEQTK.out.downsampled_reads

            ch_fastqs = ch_downsampled_reads

        }

    } else {

        ch_fastqs = ch_cat_fastq

    }

    if (params.run_umiprocessing) {
        
        //
        // MODULE: Run fgbio FastqToBam
        //
        FGBIO_FASTQTOBAM(ch_fastqs)
        ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions.first())
        ch_ubam = FGBIO_FASTQTOBAM.out.bam

        //
        // MODULE: Align with bwa mem but avoid sort
        //
        sort = false
        ALIGN_BAM_RAW(ch_ubam, ch_metref, ch_metdir, ch_metdct)
        ch_versions = ch_versions.mix(ALIGN_BAM_RAW.out.versions.first())
        ch_raw_bam = ALIGN_BAM_RAW.out.bam
        ch_raw_sort_bam = ALIGN_BAM_RAW.out.sort_bam
        ch_raw_sort_bai = ALIGN_BAM_RAW.out.sort_bai

        //
        // MODULE: Run fgbio correctumis
        //
        FGBIO_CORRECTUMIS(ch_raw_bam, params.correct_revcomp, params.correct_max_mismatch, params.correct_min_distance, params.correct_min_corrected)
        ch_multiqc_files = ch_multiqc_files.mix(FGBIO_CORRECTUMIS.out.metrics.map{it[1]}.collect())
        ch_versions = ch_versions.mix(FGBIO_CORRECTUMIS.out.versions.first())
        ch_bam_fcu = FGBIO_CORRECTUMIS.out.bam

        //
        // MODULE: Run SamToools Sort & Index
        //
        SAMTOOLS_SORT_INDEX_RAW(ch_bam_fcu, ch_metref, params.metfai)
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX_RAW.out.versions.first())
        ch_bam_fcu_sort = SAMTOOLS_SORT_INDEX_RAW.out.bam
        ch_bam_fcu_indx = SAMTOOLS_SORT_INDEX_RAW.out.bai
        ch_bam_fcu_stix = SAMTOOLS_SORT_INDEX_RAW.out.bam_bai

        //
        // MODULE: Run Picard Tool CollectMultipleMetrics
        //
        COLLECTMULTIPLEMETRICS(ch_bam_fcu_stix, ch_metref, ch_metfai)
        ch_multiqc_files = ch_multiqc_files.mix(COLLECTMULTIPLEMETRICS.out.metrics.map{it[1]}.collect())
        ch_versions = ch_versions.mix(COLLECTMULTIPLEMETRICS.out.versions.first())

        //
        // MODULE: Run ErrorRateByReadPosition 
        //
        FGBIO_ERRORRATEBYREADPOSITION_RAW(ch_bam_fcu_sort, ch_metref, ch_metfai, ch_metdct, params.known_sites, params.known_sites_tbi, params.interval_list)
        ch_multiqc_files = ch_multiqc_files.mix(FGBIO_ERRORRATEBYREADPOSITION_RAW.out.metrics.map{it[1]}.collect())
        ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION_RAW.out.versions.first())

        //
        // MODULE: Run MosDepth
        //
        MOSDEPTH_RAW(ch_bam_fcu_stix, ch_metref, ch_metfai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
        ch_versions = ch_versions.mix(MOSDEPTH_RAW.out.versions.first())

        //
        // MODULE: Run Picard's Collect HS Metrics for raw BAM files
        //
        COLLECTHSMETRICS_RAW(ch_bam_fcu_stix, ch_metref, ch_metfai, ch_metdct, params.hsmetrics_baits, params.hsmetrics_trgts, params.seq_library)
        ch_versions = ch_versions.mix(COLLECTHSMETRICS_RAW.out.versions.first())
        ch_coverage_raw  = COLLECTHSMETRICS_RAW.out.coverage
        ch_hsmetrics_raw = COLLECTHSMETRICS_RAW.out.hsmetrics

        //
        // MODULE: Run MSI Sensor PRO
        ///
        MSISENSORPRO_RAW(ch_bam_fcu_stix, ch_msi)
        ch_versions = ch_versions.mix(MSISENSORPRO_RAW.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_RAW.out.summary.map{it[1]}.collect())
        ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_RAW.out.msi_uns.map{it[1]}.collect())
        ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_RAW.out.msi_all.map{it[1]}.collect())

        //
        // MODULE: Run SamBlaster
        //
        SAMBLASTER(ch_bam_fcu)
        ch_versions = ch_versions.mix(SAMBLASTER.out.versions.first())
        ch_split_reads = SAMBLASTER.out.split_reads

        //
        // MODULE: Run SPAdes
        //
        SPADES(ch_split_reads)
        ch_versions = ch_versions.mix(SPADES.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(SPADES.out.log.map{it[1]}.collect())

        //
        // MODULE: Run FilterContigs custom script
        //
        FILTER_CONTIGS(SPADES.out.contigs, '100', '1')
        ch_versions = ch_versions.mix(FILTER_CONTIGS.out.versions.first())
        ch_split_contigs = FILTER_CONTIGS.out.fasta

        //
        // MODULE: Run fgbio GroupReadsByUmi
        //
        FGBIO_GROUPREADSBYUMI(SAMBLASTER.out.bam, params.group_strategy, params.group_edits, params.group_include_secondary, params.group_allow_inter_contig, params.group_include_supplementary, params.group_min_map_q, params.group_include_non_pf_reads, params.group_mark_duplicates)
        ch_multiqc_files = ch_multiqc_files.mix(FGBIO_GROUPREADSBYUMI.out.histogram.map{it[1]}.collect())
        ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions.first())
        ch_grouped_family_sizes = FGBIO_GROUPREADSBYUMI.out.histogram
        ch_bam_grouped = FGBIO_GROUPREADSBYUMI.out.bam

        //
        // MODULE: Run fgbio CollectDuplexSeqMetrics
        //
        FGBIO_COLLECTDUPLEXSEQMETRICS(ch_bam_grouped, params.interval_list)
        ch_multiqc_files = ch_multiqc_files.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.metrics.map{it[1]}.collect())
        ch_multiqc_files = ch_multiqc_files.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.pdf.map{it[1]}.collect())   
        ch_versions = ch_versions.mix(FGBIO_COLLECTDUPLEXSEQMETRICS.out.versions.first())

        //
        // MODULE: Run fgbio CallDuplexConsensusReads
        //
        FGBIO_CALLDUPLEXCONSENSUSREADS(ch_bam_grouped, params.call_min_reads, params.call_min_baseq)
        ch_versions = ch_versions.mix(FGBIO_CALLDUPLEXCONSENSUSREADS.out.versions.first())
        FGBIO_CALLDUPLEXCONSENSUSREADS.out.bam.set { ch_consensus_bam }

        //
        // MODULE: Run fgbio SortBam
        //
        FGBIO_SORTCONBAM(ch_consensus_bam)
        ch_versions = ch_versions.mix(FGBIO_SORTCONBAM.out.versions.first())
        ch_consensus_bam_sorted = FGBIO_SORTCONBAM.out.bam

        //
        // MODULE: Run FgBIO FilterConsensusReads to produce the "Consensus", "Duplex" & "Simplex" BAM files
        //
        FGBIO_FILTERCONSENSUSREADS(ch_consensus_bam_sorted, params.bwaref, params.bwafai, params.filter_min_reads, params.filter_min_base_quality, params.filter_max_base_error_rate, params.filter_max_read_error_rate, params.filter_max_no_call_fraction)
        ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions.first())
        ch_bam_bai_con_fil = FGBIO_FILTERCONSENSUSREADS.out.suplex_bam_bai
        ch_bam_bai_duplex_fil = FGBIO_FILTERCONSENSUSREADS.out.duplex_bam_bai
        ch_bam_bai_simplex_fil = FGBIO_FILTERCONSENSUSREADS.out.simplex_bam_bai

        // Combine BAM fils by meta data
	    ch_align_bam_con_in = ch_bam_bai_con_fil
	        .join(ch_bam_bai_duplex_fil)
	        .join(ch_bam_bai_simplex_fil)

        //
        // MODULE: Align with BWA mem
        //
        ALIGN_BAM_CON(ch_align_bam_con_in, ch_metref, ch_metfai, ch_metdir, ch_metdct)
        ch_versions = ch_versions.mix(ALIGN_BAM_CON.out.versions.first())
        ch_bam_con = ALIGN_BAM_CON.out.bam
        ch_bam_duplex = ALIGN_BAM_CON.out.duplex_bam
        ch_bam_simplex = ALIGN_BAM_CON.out.simplex_bam

        // Combine BAM fils by meta data
	    ch_sort_index_in = ch_bam_con
	        .join(ch_bam_duplex)
	        .join(ch_bam_simplex)

        //
        // MODULE: Run SamToools Sort & Index
        //
        SAMTOOLS_SORT_INDEX_CON(ch_sort_index_in, ch_metref, ch_metfai)
        ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX_CON.out.versions)
        ch_bam_con_sort = SAMTOOLS_SORT_INDEX_CON.out.bam
        ch_bam_con_indx = SAMTOOLS_SORT_INDEX_CON.out.bai
        ch_bam_con_stix = SAMTOOLS_SORT_INDEX_CON.out.bam_bai
        ch_bam_dup_stix = SAMTOOLS_SORT_INDEX_CON.out.bam_duplex
        ch_bam_sim_stix = SAMTOOLS_SORT_INDEX_CON.out.bam_simplex

        // Combine BAM fils by meta data
    	ch_umi_metrics_in = ch_bam_con_stix
    	    .join(ch_bam_dup_stix)
    	    .join(ch_bam_sim_stix)

        //
        // MODULE: Run SamTools View to count reads accross the BAM files
        //
        COLLECT_UMI_METRICS(ch_umi_metrics_in)
        ch_versions = ch_versions.mix(COLLECT_UMI_METRICS.out.versions.first())
        ch_cons_family_sizes = COLLECT_UMI_METRICS.out.cons_family_sizes

        // Combine BAM fils by meta data
    	ch_umi_read_counts_in = ch_ubam
    	    .join(ch_bam_fcu)
    	    .join(ch_bam_grouped)
    	    .join(ch_consensus_bam)
    	    .join(ch_bam_bai_con_fil)
    	    .join(ch_bam_con_stix)
    	    .join(ch_bam_dup_stix)
    	    .join(ch_bam_sim_stix)

        //
        // MODULE: Run SamTools View to count reads accross the BAM files
        //
        UMI_READ_COUNTS(ch_umi_read_counts_in)
        ch_versions = ch_versions.mix(UMI_READ_COUNTS.out.versions.first())

        //
        // MODULE: Run Preseq CCurve
        //
        PRESEQ_CCURVE(ch_cons_family_sizes)
        ch_versions = ch_versions.mix(PRESEQ_CCURVE.out.versions.first())

        //
        // MODULE: Run Preseq LCExtrap
        //
        PRESEQ_LCEXTRAP(ch_grouped_family_sizes)
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())

        //
        // MODULE: Run MosDepth
        //
        MOSDEPTH_CON(ch_bam_con_stix, ch_metref, params.metfai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
        ch_versions = ch_versions.mix(MOSDEPTH_CON.out.versions.first())

        //
        // MODULE: Run MosDepth
        //
        MOSDEPTH_DUP(ch_bam_dup_stix, ch_metref, params.metfai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
        ch_versions = ch_versions.mix(MOSDEPTH_DUP.out.versions.first())

        //
        // MODULE: Run MosDepth
        //
        MOSDEPTH_SIM(ch_bam_sim_stix, ch_metref, params.metfai, params.intervals_bed_gunzip, params.intervals_bed_gunzip_index)
        ch_versions = ch_versions.mix(MOSDEPTH_SIM.out.versions.first())

        //
        // MODULE: Run ErrorRateByReadPosition in Final BAM
        //
        FGBIO_ERRORRATEBYREADPOSITION_CON(ch_bam_con_sort, ch_metref, ch_metfai, ch_metdct, params.known_sites, params.known_sites_tbi, params.interval_list)
        ch_versions = ch_versions.mix(FGBIO_ERRORRATEBYREADPOSITION_CON.out.versions)

        //
        // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
        //
        COLLECTHSMETRICS_CON(ch_bam_con_stix, ch_metref, ch_metfai, ch_metdct, params.hsmetrics_baits, params.hsmetrics_trgts, params.seq_library)
        ch_versions = ch_versions.mix(COLLECTHSMETRICS_CON.out.versions.first())
        ch_coverage_con  = COLLECTHSMETRICS_CON.out.coverage
        ch_hsmetrics_con = COLLECTHSMETRICS_CON.out.hsmetrics

        //
        // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
        //
        COLLECTHSMETRICS_DUP(ch_bam_dup_stix, ch_metref, ch_metfai, ch_metdct, params.hsmetrics_baits, params.hsmetrics_trgts, params.seq_library)
        ch_versions = ch_versions.mix(COLLECTHSMETRICS_DUP.out.versions.first())
        ch_coverage_con  = COLLECTHSMETRICS_DUP.out.coverage
        ch_hsmetrics_con = COLLECTHSMETRICS_DUP.out.hsmetrics

        //
        // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
        //
        COLLECTHSMETRICS_SIM(ch_bam_sim_stix, ch_metref, ch_metfai, ch_metdct, params.hsmetrics_baits, params.hsmetrics_trgts, params.seq_library)
        ch_versions = ch_versions.mix(COLLECTHSMETRICS_SIM.out.versions.first())
        ch_coverage_con  = COLLECTHSMETRICS_SIM.out.coverage
        ch_hsmetrics_con = COLLECTHSMETRICS_SIM.out.hsmetrics

        //
        // MODULE: Run MSI Sensor PRO
        ///
        MSISENSORPRO_CON(ch_bam_con_stix, ch_msi)
        ch_versions = ch_versions.mix(MSISENSORPRO_CON.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_CON.out.summary.map{it[1]}.collect())
        ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_CON.out.msi_uns.map{it[1]}.collect())
        ch_multiqc_files = ch_multiqc_files.mix(MSISENSORPRO_CON.out.msi_all.map{it[1]}.collect())

        //
        // MODULE: Extract FastQ reads from BAM
        //
        SAMTOOLS_COLLATEFASTQ(ch_bam_bai_duplex_fil, ch_metref, [])
        ch_versions = ch_versions.mix(SAMTOOLS_COLLATEFASTQ.out.versions)
        ch_consensus_reads = SAMTOOLS_COLLATEFASTQ.out.fastq

        //
        // MODULE: Run FastQC
        //
        FASTQC_CONSENSUS(ch_consensus_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_CONSENSUS.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC_CONSENSUS.out.versions)

        // User aligner selection
        if (params.aligner == 'bwa-meth') {
            //
            // MODULE: Run BWA-METH
            //
            BWA_METH(ch_consensus_reads, ch_metdir, ch_metref)
            ch_versions = ch_versions.mix(BWA_METH.out.versions.first())
            ch_bam_dedup = BWA_METH.out.bam
        } else if (params.aligner == 'bwa-mem2') {
            //
            // MODULE: Run BWA-MEM2
            //
            sort_bam = 'sort'
            BWAMEM2(ch_consensus_reads, ch_bwadir, ch_bwaref, ch_bwafai, sort_bam)
            ch_versions = ch_versions.mix(BWAMEM2.out.versions.first())
            ch_bam_dedup = BWAMEM2.out.bam
        } else {
            error "Invalid aligner selected: ${params.aligner}. Please choose either 'bwa-meth' or 'bwa-mem2'"
        }

    } else {

            // User aligner selection
        if (params.aligner == 'bwa-meth') {
            //
            // MODULE: Run BWA-METH
            //
            BWA_METH(ch_fastqs, ch_metdir, ch_metref)
            ch_versions = ch_versions.mix(BWA_METH.out.versions.first())
            ch_bam_mapped = BWA_METH.out.bam
        } else if (params.aligner == 'bwa-mem2') {
            //
            // MODULE: Run BWA-MEM2
            //
            sort_bam = 'sort'
            BWAMEM2(ch_fastqs, ch_bwadir, ch_bwaref, ch_bwafai, sort_bam)
            ch_versions = ch_versions.mix(BWAMEM2.out.versions.first())
            ch_bam_mapped = BWAMEM2.out.bam
        } else {
            error "Invalid aligner selected: ${params.aligner}. Please choose either 'bwa-meth' or 'bwa-mem2'"
        }

        //
        // MODULE: Run GATK4 MarkDuplicates
        //
        GATK4_MARKDUPLICATES(ch_bam_mapped, params.metref, params.metfai)
        ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.complex_metrics.map{it[1]}.collect())
        ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.map{it[1]}.collect())
        ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())
        ch_bam_dedup = GATK4_MARKDUPLICATES.out.bam

    }

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
    // MODULE: Run SAMtools Stats
    //
    SAMTOOLS_STATS(ch_bam_mapped_targeted_indexed, ch_fasta)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    //
    // MODULE: Run MosDepth
    //
    ch_bam_raw_mosdepth = ch_bam_mapped_targeted_indexed.map { meta, bam, bai -> tuple(meta, bam) }
    ch_bai_raw_mosdepth = ch_bam_mapped_targeted_indexed.map { meta, bam, bai -> tuple(meta, bai) }
    MOSDEPTH(ch_bam_raw_mosdepth, ch_bai_raw_mosdepth, ch_fasta, params.fai, params.mosdepth_canonical_exomes)
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt)

    //
    // MODULE: Run Picard Tool CollectMultipleMetrics
    //
    PICARD_COLLECTMULTIPLEMETRICS(ch_bam_mapped_targeted_indexed, ch_metref, ch_metfai)
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics.map{it[1]}.collect())
    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())

    //
    // MODULE: Run Picard's Collect HS Metrics for consensus BAM files
    //
    COLLECTHSMETRICS(ch_bam_mapped_targeted_indexed, ch_metref, ch_metfai, ch_metdct, params.hsmetrics_baits, params.hsmetrics_trgts, params.seq_library)
    ch_versions = ch_versions.mix(COLLECTHSMETRICS.out.versions.first())
    ch_coverage_con  = COLLECTHSMETRICS.out.coverage
    ch_hsmetrics_con = COLLECTHSMETRICS.out.hsmetrics

    //
    // MODULE: Run rasTair
    //
    RASTAIR_FULL(ch_bam_mapped_targeted_indexed, ch_metref, ch_metfai, ch_metbed, params.pymbias_plot_type, params.pymbias_plot_ax_x, params.pymbias_plot_ax_y)
    ch_versions = ch_versions.mix(RASTAIR_FULL.out.versions.first())
    ch_rastair_mods = RASTAIR_FULL.out.mods
    
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
    
    //
    // MODULE: Run FgBio ClipBAM 
    //
    FGBIO_CLIPBAM(ch_bam_mapped_targeted_indexed, ch_metref, ch_metfai)
    ch_versions = ch_versions.mix(FGBIO_CLIPBAM.out.versions)
    ch_bam_clipped = FGBIO_CLIPBAM.out.bam

    //
    // MODULE: HaplotypeCaller from GATK4 (Calls germline SNPs and indels via local re-assembly of haplotypes.)
    //

    GATK4_HAPLOTYPECALLER(ch_metfai, ch_metdictgatk, ch_known_sites, ch_metref, ch_known_sites_tbi, ch_intervals, ch_bam_clipped)
    ch_versions  = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.ifEmpty(null))
    ch_haplotypecaller_raw = GATK4_HAPLOTYPECALLER.out.vcf
    ch_haplotypecaller_tbi = GATK4_HAPLOTYPECALLER.out.tbi
    ch_haplotypecaller_raw_tbi_combined = ch_haplotypecaller_raw.join(ch_haplotypecaller_tbi)
    ch_haplotypecaller_vcf_tbi = ch_haplotypecaller_raw_tbi_combined.map { meta, vcf, tbi ->
        [meta, vcf, tbi]
    }

    //
    // MODULE: VariantFiltration from GATK4 (Filter variant calls based on certain criteria.)
    // 
    GATK4_VARIANTFILTRATION(ch_haplotypecaller_raw, ch_haplotypecaller_tbi, ch_metref, ch_metfai, ch_metdictgatk)
    ch_versions  = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions.ifEmpty(null))
    ch_final_vcf = GATK4_VARIANTFILTRATION.out.vcf

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:

    versions            = ch_versions
    multiqc_report      = MULTIQC.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                      FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Add readgroup to meta and remove lane
def addReadgroupToMeta(meta, files) {
    def flowcell = flowcellLaneFromFastq(files[0])
    def CN = params.seq_center ? "CN:${params.seq_center}\\t" : ''
    def read_group = "\"@RG\\tID:${meta.id}\\t${CN}PU:${flowcell}\\tSM:${meta.patient}_${meta.id}\\tLB:${meta.id}\\tPL:${params.seq_platform}\""
    meta  = meta + [read_group: read_group.toString()]

    return [ meta, files ]
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}

// Get the mbias parameters from the TSV file
def getMbiasParams(tsv_file) {
    def values = [:]
    def content = file(tsv_file).readLines()
    def lines = content[1..-1]
    def ot_values = []
    def ob_values = []
    
    lines.each { line ->
        def fields = line.split('\t')
        def strand = fields[1]
        def left = fields[3]
        def right = fields[4]
        
        if (strand == "OT") {
            ot_values.add(left)
            ot_values.add(right)
        } else if (strand == "OB") {
            ob_values.add(left)
            ob_values.add(right)
        }
    }
    
    values.not = ot_values.join(',')
    values.nob = ob_values.join(',')
    
    return values
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                       THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
