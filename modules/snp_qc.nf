/*
========================================================================================
    SNP-LEVEL QUALITY CONTROL
========================================================================================
    1. Call rate: remove SNPs genotyped in < 70% of retained samples
    2. Excess heterozygosity: remove SNPs with obs het > 0.49 (paralogy)
    3. Known faulty probes: remove if list provided
    4. Strand-ambiguous SNPs: remove A/T and C/G pairs (batch effect mitigation)
    NO MAF filter (plan explicitly says not to apply one)
    NO LD pruning (plan explicitly says not to)
*/

process SNP_QC {
    tag "snp_qc"
    label 'process_medium'

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path genotype_matrix
    path snp_info

    output:
    path "snp_qc_filtered.tsv", emit: filtered_matrix
    path "snp_info_qc.tsv",     emit: snp_info
    path "snp_qc_report.txt",   emit: report

    script:
    def faulty_cmd = params.faulty_probes ? "faulty <- readLines('${params.faulty_probes}')" : "faulty <- character(0)"
    """
    #!/usr/bin/env Rscript

    library(data.table)

    sink("snp_qc_report.txt")
    cat("=== SNP-Level QC Report ===\\n\\n")

    geno <- fread("${genotype_matrix}", header = TRUE)
    snp_info <- fread("${snp_info}", header = TRUE)
    cat("Input:", nrow(geno), "samples x", ncol(geno) - 1, "SNPs\\n\\n")

    # Load faulty probe list if provided
    ${faulty_cmd}
    cat("Faulty probes to exclude:", length(faulty), "\\n")

    snp_cols <- colnames(geno)[-1]
    geno_mat <- as.matrix(geno[, -1])
    n_samples <- nrow(geno_mat)

    # Per-SNP statistics
    snp_call_rate <- colSums(!is.na(geno_mat)) / n_samples
    snp_het_rate  <- colSums(geno_mat == 1, na.rm = TRUE) / colSums(!is.na(geno_mat))
    snp_het_rate[is.nan(snp_het_rate)] <- 0

    # Filters
    pass_callrate <- snp_call_rate >= ${params.snp_call_rate}
    pass_het      <- snp_het_rate <= ${params.max_snp_het}
    pass_faulty   <- !(snp_cols %in% faulty)

    # Strand-ambiguous SNP filter (A/T and C/G pairs)
    # These cannot be reliably oriented across platforms and are a major
    # source of batch effects in multi-platform genotype analyses.
    # Skipped for single-platform (WGS-only) runs where strand is unambiguous.
    pass_strand <- rep(TRUE, length(snp_cols))
    skip_strand <- as.logical("${params.skip_strand_filter}")
    if (!skip_strand && "ref" %in% colnames(snp_info) && "alt" %in% colnames(snp_info)) {
        # Build lookup by snp_id
        strand_lookup <- snp_info[, .(snp_id, ref, alt)]
        setkey(strand_lookup, snp_id)
        for (idx in seq_along(snp_cols)) {
            si <- strand_lookup[snp_cols[idx]]
            if (nrow(si) == 1 && !is.na(si\$ref) && !is.na(si\$alt)) {
                pair <- paste0(si\$ref, si\$alt)
                if (pair %in% c("AT", "TA", "CG", "GC")) {
                    pass_strand[idx] <- FALSE
                }
            }
        }
    }

    cat("SNPs failing call rate (<", ${params.snp_call_rate}, "):", sum(!pass_callrate), "\\n")
    cat("SNPs failing het filter (>", ${params.max_snp_het}, "):", sum(!pass_het), "\\n")
    cat("SNPs in faulty probe list:", sum(!pass_faulty), "\\n")
    cat("Strand-ambiguous SNPs (A/T, C/G) removed:", sum(!pass_strand), "\\n")

    keep <- pass_callrate & pass_het & pass_faulty & pass_strand
    cat("SNPs passing all filters:", sum(keep), "\\n\\n")

    # Filter
    kept_snps <- snp_cols[keep]
    geno_filtered <- geno[, c("sample_id", kept_snps), with = FALSE]

    cat("Output:", nrow(geno_filtered), "samples x", ncol(geno_filtered) - 1, "SNPs\\n")
    sink()

    # Update SNP info
    snp_info_filtered <- snp_info[snp_id %in% kept_snps]

    fwrite(geno_filtered, "snp_qc_filtered.tsv", sep = "\\t")
    fwrite(snp_info_filtered, "snp_info_qc.tsv", sep = "\\t")
    """
}
