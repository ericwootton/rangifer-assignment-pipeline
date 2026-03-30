#!/usr/bin/env nextflow

/*
========================================================================================
    RANGIFER ASSIGNMENT PIPELINE
========================================================================================
    Hierarchical probabilistic assignment of Rangifer tarandus samples to subspecies,
    ecotype, and herd using multi-method ensemble classification.

    Four modes (--mode):
      'snp'       - SNP chip data (default)
      'mito'      - Mitochondrial variants (haploid-specific QC)
      'wgs'       - Whole-genome data (SV calling, phasing, Fst outliers, IBD, fineSTRUCTURE)
      'herd_diff' - Herd-level differentiation analysis

    Author: E. Wootton
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

version = '2.0.0'

/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/

def helpMessage() {
    log.info"""

    RANGIFER ASSIGNMENT PIPELINE v${version}
    ========================================

    Usage:
        nextflow run main.nf --wgs_vcf <vcf> --metadata <csv> --mode <snp|mito|wgs> [options]

    Required:
        --wgs_vcf           Path to QC'd VCF (chip/mito/wgs variants)
        --metadata          CSV with sample metadata (Sample, Subspecies, Ecotype_Cleaned, Herd_Cleaned)

    Mode selection:
        --mode              Analysis mode: 'snp' (default), 'mito', 'wgs', or 'herd_diff'

    Herd differentiation mode inputs:
        --genotype_matrix   Pre-QC'd genotype matrix TSV (from a previous pipeline run)

    SNP mode inputs:
        --finalreport       Path to Illumina FinalReport CSV
        --fr_samplesheet    Path to Illumina sample sheet CSV
        --nwt_genotypes     Path to NWT genotype matrix CSV
        --chip_manifest     Path to chip manifest (required with chip data)
        --fasta             Reference FASTA (required with chip data)

    WGS mode inputs:
        --bam_dir           Path to BAM directory (for SV calling)
        --fasta             Reference FASTA (for SV calling)

    QC parameters:
        --min_call_rate     Minimum sample call rate [default: 0.70]
        --max_het           Maximum heterozygosity [default: 0.50]
        --min_het           Minimum heterozygosity [default: 0.15]
        --snp_call_rate     Minimum SNP call rate [default: 0.70]
        --max_snp_het       Maximum SNP heterozygosity [default: 0.49]
        --ibs_threshold     IBS duplicate detection threshold [default: 0.95]
        --skip_strand_filter Skip strand-ambiguous SNP filter [default: false]

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (params.mode == 'herd_diff') {
    if (!params.genotype_matrix) {
        error "ERROR: --genotype_matrix is required for herd_diff mode"
    }
} else if (!params.wgs_vcf) {
    error "ERROR: --wgs_vcf is required"
}
if (!params.metadata) {
    error "ERROR: --metadata is required"
}

def valid_modes = ['snp', 'mito', 'wgs', 'herd_diff']
if (!(params.mode in valid_modes)) {
    error "ERROR: --mode must be one of: ${valid_modes.join(', ')} (got '${params.mode}')"
}

if (params.mode == 'snp' && (params.finalreport || params.nwt_genotypes)) {
    if (!params.chip_manifest) {
        error "ERROR: --chip_manifest is required when chip data is provided"
    }
    if (!params.fasta) {
        error "ERROR: --fasta is required when chip data is provided"
    }
}

if (params.mode == 'wgs') {
    if (!params.bam_dir) {
        error "ERROR: --bam_dir is required for WGS mode"
    }
    if (!params.fasta) {
        error "ERROR: --fasta is required for WGS mode"
    }
}

/*
========================================================================================
    INCLUDE MODULES
========================================================================================
*/

// Phase 1: Data Integration
include { BUILD_SNP_MAP } from './modules/build_snp_map'
include { CONVERT_FINALREPORT } from './modules/convert_finalreport'
include { CONVERT_NWT_GENOTYPES } from './modules/convert_nwt_genotypes'
include { EXTRACT_VCF_GENOTYPES } from './modules/extract_vcf_genotypes'
include { MERGE_GENOTYPE_SOURCES } from './modules/merge_genotype_sources'

// Phase 2: QC
include { SAMPLE_QC } from './modules/sample_qc'
include { SNP_QC } from './modules/snp_qc'
include { DUPLICATE_DETECTION } from './modules/duplicate_detection'
include { BATCH_CORRECT } from './modules/batch_correct'

// Phase 3: Structure Discovery
include { RUN_PCA } from './modules/run_pca'
include { DAPC_FIND_CLUSTERS } from './modules/dapc_find_clusters'
include { CONVERT_TO_BED } from './modules/convert_to_bed'
include { RUN_ADMIXTURE } from './modules/run_admixture'
include { SYNTHESIZE_STRUCTURE } from './modules/synthesize_structure'

// Phase 4: Hierarchical Classification (shared)
include { PREPARE_NODE_DATA } from './modules/prepare_node_data'
include { CONVERT_NODE_BED } from './modules/convert_node_bed'
include { RUN_NODE_ADMIXTURE } from './modules/run_node_admixture'
include { RUN_DAPC } from './modules/run_dapc'
include { RUN_ASSIGNPOP } from './modules/run_assignpop'
include { RUN_POPFINDER } from './modules/run_popfinder'
include { ENSEMBLE_DECISION } from './modules/ensemble_decision'
include { GENERATE_NODE_PLOTS } from './modules/generate_node_plots'

// SNP mode: Fst node filter
include { FST_NODE_FILTER } from './modules/fst_node_filter'

// WGS mode: SV calling, phasing, per-node analyses
include { CALL_SVS_DISCOVER } from './modules/call_svs'
include { MERGE_SV_SITES } from './modules/call_svs'
include { GENOTYPE_SVS } from './modules/call_svs'
include { MERGE_SV_GENOTYPES } from './modules/merge_sv_genotypes'
include { PHASE_VCF } from './modules/phase_vcf'
include { COMPUTE_FST_DXY } from './modules/compute_fst_dxy'
include { OUTFLANK_OUTLIERS } from './modules/outflank_outliers'
include { PREPARE_WGS_NODE_DATA } from './modules/prepare_wgs_node_data'
include { IBD_NETWORK } from './modules/ibd_network'
include { CHROMOPAINTER_FINESTRUCTURE } from './modules/chromopainter_finestructure'

// Phase 5: Output
include { COMPILE_ASSIGNMENTS } from './modules/compile_assignments'
include { FEATURE_IMPORTANCE } from './modules/feature_importance'
include { GENERATE_REPORT } from './modules/generate_report'

// Herd differentiation mode
include { PREPARE_HERD_DATA } from './modules/prepare_herd_data'
include { PAIRWISE_FST } from './modules/pairwise_fst'
include { HERD_DIFF_SUMMARY } from './modules/herd_diff_summary'

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {

    log.info """
    =========================================
    RANGIFER ASSIGNMENT PIPELINE v${version}
    Mode: ${params.mode}
    =========================================
    """.stripIndent()

    ch_metadata = Channel.fromPath(params.metadata, checkIfExists: true)

  if (params.mode == 'herd_diff') {
    // ==========================================
    // HERD DIFFERENTIATION MODE
    // Takes a pre-QC'd genotype matrix and runs all differentiation
    // analyses with herds as groups (no hierarchical nodes).
    // ==========================================

    ch_matrix = Channel.fromPath(params.genotype_matrix, checkIfExists: true)

    PREPARE_HERD_DATA(ch_matrix, ch_metadata)

    RUN_PCA(PREPARE_HERD_DATA.out.herd_matrix, ch_metadata)
    DAPC_FIND_CLUSTERS(PREPARE_HERD_DATA.out.herd_matrix, ch_metadata)
    CONVERT_TO_BED(PREPARE_HERD_DATA.out.herd_matrix)
    RUN_ADMIXTURE(CONVERT_TO_BED.out.plink_bed)

    ch_herd_node = PREPARE_HERD_DATA.out.herd_node_dir
        .map { dir ->
            def meta_file = file("${dir}/node_meta.json")
            [ dir, meta_file ]
        }

    PAIRWISE_FST(ch_herd_node)

    CONVERT_NODE_BED(ch_herd_node)
    RUN_NODE_ADMIXTURE(CONVERT_NODE_BED.out.plink_bed)
    RUN_DAPC(ch_herd_node)
    RUN_ASSIGNPOP(ch_herd_node)

    // Only run popfinder if enough herds have sufficient samples
    ch_popfinder_node = ch_herd_node
        .filter { node_dir, meta_file ->
            def meta = new groovy.json.JsonSlurper().parse(meta_file)
            def counts = meta.sample_counts instanceof Map ? meta.sample_counts.values() : []
            def n_sufficient = counts.count { it >= params.min_samples_popfinder }
            def n_groups = counts.size()
            n_groups >= 2 && n_sufficient >= 2 && n_sufficient >= (n_groups / 2.0)
        }
    RUN_POPFINDER(ch_popfinder_node)

    HERD_DIFF_SUMMARY(
        PAIRWISE_FST.out.fst_matrix,
        PAIRWISE_FST.out.fst_long,
        RUN_PCA.out.scores,
        RUN_PCA.out.eigenvalues,
        RUN_ADMIXTURE.out.structure_results,
        RUN_ADMIXTURE.out.choosek,
        RUN_DAPC.out.results.map { it[1] },
        RUN_ASSIGNPOP.out.results.map { it[1] },
        RUN_NODE_ADMIXTURE.out.results.map { it[1] },
        ch_metadata
    )

  } else {
    // ==========================================
    // HIERARCHICAL ASSIGNMENT MODE (snp / mito / wgs)
    // ==========================================

    // ==========================================
    // PHASE 1: DATA INTEGRATION
    // ==========================================

    ch_wgs_vcf = Channel.fromPath(params.wgs_vcf, checkIfExists: true)
        .map { vcf ->
            def tbi = file("${vcf}.tbi")
            if (!tbi.exists()) error "Index file not found: ${tbi}"
            [ vcf, tbi ]
        }

    if (params.mode == 'snp') {
        // SNP mode: full chip data integration
        ch_finalreport_geno = Channel.empty()
        ch_nwt_geno = Channel.empty()

        if (params.finalreport || params.nwt_genotypes) {
            ch_manifest = Channel.fromPath(params.chip_manifest, checkIfExists: true)
            ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
            BUILD_SNP_MAP(ch_manifest, ch_fasta)

            if (params.finalreport) {
                ch_fr = Channel.fromPath(params.finalreport, checkIfExists: true)
                ch_fr_sheet = params.fr_samplesheet
                    ? Channel.fromPath(params.fr_samplesheet, checkIfExists: true)
                    : Channel.of(file('NO_SAMPLESHEET'))
                CONVERT_FINALREPORT(ch_fr, ch_fr_sheet, BUILD_SNP_MAP.out.snp_map)
                ch_finalreport_geno = CONVERT_FINALREPORT.out.allele_matrix
            }

            if (params.nwt_genotypes) {
                ch_nwt = Channel.fromPath(params.nwt_genotypes, checkIfExists: true)
                CONVERT_NWT_GENOTYPES(ch_nwt, BUILD_SNP_MAP.out.snp_map)
                ch_nwt_geno = CONVERT_NWT_GENOTYPES.out.allele_matrix
            }
        }

        EXTRACT_VCF_GENOTYPES(ch_wgs_vcf)

        ch_chip_sources = ch_finalreport_geno
            .mix(ch_nwt_geno)
            .collect()
            .ifEmpty(file('NO_CHIP_DATA'))

        MERGE_GENOTYPE_SOURCES(
            EXTRACT_VCF_GENOTYPES.out.genotype_tsv,
            EXTRACT_VCF_GENOTYPES.out.snp_info,
            ch_chip_sources
        )

        ch_genotype_matrix = MERGE_GENOTYPE_SOURCES.out.genotype_matrix
        ch_snp_info = MERGE_GENOTYPE_SOURCES.out.snp_info
        ch_data_source = MERGE_GENOTYPE_SOURCES.out.data_source

    } else {
        // Mito and WGS modes: VCF-only, no chip data integration
        EXTRACT_VCF_GENOTYPES(ch_wgs_vcf)

        ch_chip_sources = Channel.of(file('NO_CHIP_DATA'))

        MERGE_GENOTYPE_SOURCES(
            EXTRACT_VCF_GENOTYPES.out.genotype_tsv,
            EXTRACT_VCF_GENOTYPES.out.snp_info,
            ch_chip_sources
        )

        ch_genotype_matrix = MERGE_GENOTYPE_SOURCES.out.genotype_matrix
        ch_snp_info = MERGE_GENOTYPE_SOURCES.out.snp_info
        ch_data_source = MERGE_GENOTYPE_SOURCES.out.data_source
    }

    // ==========================================
    // PHASE 2: QUALITY CONTROL
    // ==========================================

    SNP_QC(
        ch_genotype_matrix,
        ch_snp_info
    )

    SAMPLE_QC(
        SNP_QC.out.filtered_matrix,
        ch_data_source
    )

    DUPLICATE_DETECTION(
        SAMPLE_QC.out.filtered_matrix
    )

    BATCH_CORRECT(
        DUPLICATE_DETECTION.out.filtered_matrix,
        ch_data_source
    )

    // ==========================================
    // PHASE 3: STRUCTURE DISCOVERY
    // ==========================================

    RUN_PCA(
        BATCH_CORRECT.out.corrected_matrix,
        ch_metadata
    )

    DAPC_FIND_CLUSTERS(
        BATCH_CORRECT.out.corrected_matrix,
        ch_metadata
    )

    CONVERT_TO_BED(
        BATCH_CORRECT.out.corrected_matrix
    )

    RUN_ADMIXTURE(
        CONVERT_TO_BED.out.plink_bed
    )

    SYNTHESIZE_STRUCTURE(
        RUN_PCA.out.pca_results,
        DAPC_FIND_CLUSTERS.out.cluster_results,
        RUN_ADMIXTURE.out.structure_results,
        ch_metadata
    )

    // ==========================================
    // PHASE 4: WGS-SPECIFIC GENOME-WIDE PREPROCESSING
    // (only in WGS mode, done once before per-node analysis)
    // ==========================================

    if (params.mode == 'wgs') {
        // SV calling (Delly)
        ch_bams = Channel
            .fromPath("${params.bam_dir}/*.bam")
            .map { bam ->
                def sample = bam.baseName.replaceAll(/\.sorted\.markdup$/, '')
                def bai = file("${bam}.bai")
                tuple(sample, bam, bai)
            }

        ch_fasta_wgs = Channel.fromPath(params.fasta, checkIfExists: true)
        ch_fai_wgs = Channel.fromPath("${params.fasta}.fai", checkIfExists: true)

        CALL_SVS_DISCOVER(ch_bams, ch_fasta_wgs.first(), ch_fai_wgs.first())

        ch_all_bcfs = CALL_SVS_DISCOVER.out.bcf.map { it[1] }.collect()
        ch_all_csis = CALL_SVS_DISCOVER.out.bcf.map { it[2] }.collect()

        MERGE_SV_SITES(ch_all_bcfs, ch_all_csis)

        GENOTYPE_SVS(
            ch_bams,
            ch_fasta_wgs.first(),
            ch_fai_wgs.first(),
            MERGE_SV_SITES.out.sites_bcf,
            MERGE_SV_SITES.out.sites_csi
        )

        ch_geno_bcfs = GENOTYPE_SVS.out.bcf.map { it[1] }.collect()
        ch_geno_csis = GENOTYPE_SVS.out.bcf.map { it[2] }.collect()

        MERGE_SV_GENOTYPES(ch_geno_bcfs, ch_geno_csis)

        // Phasing (Eagle)
        ch_qc_vcf = Channel.fromPath(params.wgs_vcf, checkIfExists: true)
            .map { vcf -> [ vcf, file("${vcf}.tbi") ] }
        PHASE_VCF(ch_qc_vcf.map { it[0] }, ch_qc_vcf.map { it[1] })
    }

    // ==========================================
    // PHASE 4/5: HIERARCHICAL CLASSIFICATION
    // ==========================================

    PREPARE_NODE_DATA(
        BATCH_CORRECT.out.corrected_matrix,
        ch_metadata,
        SYNTHESIZE_STRUCTURE.out.hierarchy_definition
    )

    ch_nodes = PREPARE_NODE_DATA.out.node_data
        .flatten()
        .filter { it.isDirectory() }
        .map { node_dir ->
            def meta_file = file("${node_dir}/node_meta.json")
            def node_name = node_dir.name
            [ node_name, node_dir, meta_file ]
        }

    // Unfiltered node data for classifiers that use all SNPs
    ch_unfiltered_nodes = ch_nodes
        .map { node_name, node_dir, meta_file ->
            [ node_dir, meta_file ]
        }

    if (params.mode == 'snp') {
        // SNP mode: Fst filter to top 5000 SNPs for assignPOP
        FST_NODE_FILTER(ch_nodes)

        ch_fst_filtered_nodes = FST_NODE_FILTER.out.filtered_node
            .map { node_name, filtered_dir ->
                def meta_file = file("${filtered_dir}/node_meta.json")
                [ filtered_dir, meta_file ]
            }

        // DAPC, popfinder, and ADMIXTURE use all SNPs
        ch_class_nodes = ch_unfiltered_nodes
        // assignPOP uses Fst-filtered SNPs
        ch_assignpop_nodes = ch_fst_filtered_nodes

    } else if (params.mode == 'wgs') {
        // WGS mode: OutFLANK outliers + SV genotypes per node

        COMPUTE_FST_DXY(
            ch_nodes.map { name, dir, meta -> [ name, dir, meta ] },
            PHASE_VCF.out.phased_vcf
        )

        OUTFLANK_OUTLIERS(ch_nodes)

        PREPARE_WGS_NODE_DATA(
            OUTFLANK_OUTLIERS.out.outlier_data,
            MERGE_SV_GENOTYPES.out.sv_matrix
        )

        IBD_NETWORK(
            ch_nodes,
            PHASE_VCF.out.phased_vcf
        )

        CHROMOPAINTER_FINESTRUCTURE(
            ch_nodes,
            PHASE_VCF.out.phased_vcf
        )

        // Use WGS-prepared node data for classification
        ch_class_nodes = PREPARE_WGS_NODE_DATA.out.wgs_node_data
            .map { node_name, wgs_dir ->
                def meta_file = file("${wgs_dir}/node_meta.json")
                [ wgs_dir, meta_file ]
            }
        ch_assignpop_nodes = ch_class_nodes

    } else {
        // Mito mode: use node data as-is
        ch_class_nodes = ch_nodes
            .map { node_name, node_dir, meta_file ->
                [ node_dir, meta_file ]
            }
        ch_assignpop_nodes = ch_class_nodes
    }

    // Run classifiers
    CONVERT_NODE_BED(ch_class_nodes)
    RUN_NODE_ADMIXTURE(CONVERT_NODE_BED.out.plink_bed)

    RUN_DAPC(ch_class_nodes)

    RUN_ASSIGNPOP(ch_assignpop_nodes)

    // popfinder at nodes where majority of classes have sufficient samples
    ch_popfinder_nodes = ch_class_nodes
        .filter { node_dir, meta_file ->
            def meta = new groovy.json.JsonSlurper().parse(meta_file)
            def counts = meta.sample_counts instanceof Map ? meta.sample_counts.values() : []
            def n_sufficient = counts.count { it >= params.min_samples_popfinder }
            def n_groups = counts.size()
            n_groups >= 2 && n_sufficient >= 2 && n_sufficient >= (n_groups / 2.0)
        }
    RUN_POPFINDER(ch_popfinder_nodes)

    // Join classifier results per node
    ch_dapc_results = RUN_DAPC.out.results
    // In SNP mode, assignPOP runs on Fst-filtered directories; strip the suffix to match keys
    ch_assignpop_results = RUN_ASSIGNPOP.out.results
        .map { node_key, result_dir ->
            def clean_key = node_key.replaceAll(/_fst_filtered$/, '')
            [ clean_key, result_dir ]
        }
    ch_popfinder_results = RUN_POPFINDER.out.results

    // Build popfinder lookup as a Groovy Map to handle nodes without popfinder results
    ch_popfinder_map = ch_popfinder_results
        .map { key, dir -> [ key, dir ] }
        .toList()
        .map { entries ->
            def m = [:]
            entries.each { m[it[0]] = it[1] }
            m
        }

    // Placeholder directory for nodes without popfinder results
    def no_popfinder_dir = file("${workflow.workDir}/NO_POPFINDER")
    no_popfinder_dir.mkdirs()

    ch_ensemble_input = ch_dapc_results
        .join(ch_assignpop_results, by: 0)
        .combine(ch_popfinder_map)
        .map { key, dapc_dir, ap_dir, pf_map ->
            def pf_dir = pf_map.containsKey(key) ? pf_map[key] : no_popfinder_dir
            [ key, dapc_dir, ap_dir, pf_dir ]
        }

    ENSEMBLE_DECISION(ch_ensemble_input)

    // ==========================================
    // PHASE 5: PER-NODE PLOTS
    // ==========================================

    ch_admixture_results = RUN_NODE_ADMIXTURE.out.results
    ch_admixture_map = ch_admixture_results
        .map { key, dir -> [ key, dir ] }
        .toList()
        .map { entries ->
            def m = [:]
            entries.each { m[it[0]] = it[1] }
            m
        }

    def no_admixture_dir = file("${workflow.workDir}/NO_ADMIXTURE")
    no_admixture_dir.mkdirs()

    ch_ensemble_keyed = ENSEMBLE_DECISION.out.decisions
        .map { f ->
            def name = f.name.replaceAll(/_ensemble\.tsv$/, '')
            [ name, f ]
        }

    ch_node_dir_keyed = ch_nodes
        .map { name, dir, meta -> [ name, dir ] }

    ch_plot_input = ch_node_dir_keyed
        .join(ch_dapc_results, by: 0)
        .join(ch_assignpop_results, by: 0)
        .combine(ch_popfinder_map)
        .combine(ch_admixture_map)
        .map { key, node_dir, dapc_dir, ap_dir, pf_map, adm_map ->
            def pf_dir = pf_map.containsKey(key) ? pf_map[key] : no_popfinder_dir
            def adm_dir = adm_map.containsKey(key) ? adm_map[key] : no_admixture_dir
            [ key, node_dir, dapc_dir, ap_dir, pf_dir, adm_dir ]
        }
        .join(ch_ensemble_keyed, by: 0)

    GENERATE_NODE_PLOTS(ch_plot_input, ch_data_source)

    // ==========================================
    // PHASE 6: COMPILE AND REPORT
    // ==========================================

    ch_all_ensemble = ENSEMBLE_DECISION.out.decisions.collect()

    COMPILE_ASSIGNMENTS(
        ch_all_ensemble,
        ch_metadata,
        SAMPLE_QC.out.sample_stats
    )

    FEATURE_IMPORTANCE(
        ch_dapc_results.map { it[1] }.collect(),
        ch_assignpop_results.map { it[1] }.collect()
    )

    GENERATE_REPORT(
        COMPILE_ASSIGNMENTS.out.final_assignments,
        FEATURE_IMPORTANCE.out.importance_summary,
        RUN_PCA.out.pca_results,
        SYNTHESIZE_STRUCTURE.out.hierarchy_definition
    )

  } // end else (hierarchical mode)
}
