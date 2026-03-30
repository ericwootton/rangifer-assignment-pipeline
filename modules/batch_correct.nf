/*
========================================================================================
    PLATFORM BATCH CORRECTION
========================================================================================
    Removes SNPs with significant allele frequency or missingness differences between
    genotyping platforms (WGS, FinalReport, NWT chip). Uses per-SNP Fisher's exact
    tests with BH multiple testing correction.

    Two tests per SNP:
      1. Differential missingness: [called, missing] x platform
      2. Differential allele frequency: [ref_alleles, alt_alleles] x platform

    SNPs failing either test (BH-adjusted p < alpha) are removed.
*/

process BATCH_CORRECT {
    tag "batch_correct"
    label 'process_medium'

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path genotype_matrix
    path data_source

    output:
    path "batch_corrected_matrix.tsv", emit: corrected_matrix
    path "batch_correction_report.txt", emit: report

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    sink("batch_correction_report.txt")
    cat("=== Platform Batch Correction Report ===\\n\\n")

    geno <- fread("${genotype_matrix}", header = TRUE)
    ds <- fread("${data_source}", header = TRUE)

    sample_ids <- geno\$sample_id
    snp_cols <- colnames(geno)[-1]
    geno_mat <- as.matrix(geno[, -1])
    n_snps <- length(snp_cols)

    cat("Input:", nrow(geno), "samples x", n_snps, "SNPs\\n\\n")

    # Map samples to platforms (only samples present in the matrix)
    platform <- ds\$source[match(sample_ids, ds\$sample_id)]
    platforms <- sort(unique(platform[!is.na(platform)]))

    cat("Platforms detected:", paste(platforms, collapse = ", "), "\\n")
    for (p in platforms) {
        cat("  ", p, ":", sum(platform == p, na.rm = TRUE), "samples\\n")
    }
    cat("\\n")

    if (length(platforms) < 2) {
        cat("Only one platform detected -- no batch correction needed.\\n")
        cat("Output:", nrow(geno), "samples x", n_snps, "SNPs (unchanged)\\n")
        sink()
        fwrite(geno, "batch_corrected_matrix.tsv", sep = "\\t")
        quit(save = "no")
    }

    n_plat <- length(platforms)

    # Vectorize: compute per-platform counts for all SNPs at once
    ref_counts   <- matrix(0L, nrow = n_plat, ncol = n_snps)
    alt_counts   <- matrix(0L, nrow = n_plat, ncol = n_snps)
    call_counts  <- matrix(0L, nrow = n_plat, ncol = n_snps)
    miss_counts  <- matrix(0L, nrow = n_plat, ncol = n_snps)
    n_per_plat   <- integer(n_plat)

    for (p_idx in seq_along(platforms)) {
        mask <- which(platform == platforms[p_idx])
        n_per_plat[p_idx] <- length(mask)
        p_mat <- geno_mat[mask, , drop = FALSE]
        call_counts[p_idx, ]  <- as.integer(colSums(!is.na(p_mat)))
        miss_counts[p_idx, ]  <- as.integer(colSums(is.na(p_mat)))
        ref_counts[p_idx, ]   <- as.integer(2L * colSums(p_mat == 0L, na.rm = TRUE) +
                                             colSums(p_mat == 1L, na.rm = TRUE))
        alt_counts[p_idx, ]   <- as.integer(2L * colSums(p_mat == 2L, na.rm = TRUE) +
                                             colSums(p_mat == 1L, na.rm = TRUE))
    }

    cat("--- Test 1: Differential Missingness ---\\n")
    miss_pvals <- rep(1.0, n_snps)

    for (i in seq_len(n_snps)) {
        tab <- cbind(call_counts[, i], miss_counts[, i])
        # Only test platforms with samples present
        has_samples <- rowSums(tab) > 0
        if (sum(has_samples) < 2) next
        tab <- tab[has_samples, , drop = FALSE]
        # Only test if missingness varies (not all-called or all-missing across platforms)
        if (all(tab[, 2] == 0) || all(tab[, 1] == 0)) next
        miss_pvals[i] <- tryCatch(fisher.test(tab)\$p.value, error = function(e) 1.0)
    }

    miss_padj <- p.adjust(miss_pvals, method = "BH")
    fail_miss <- miss_padj < ${params.batch_alpha}
    cat("SNPs with significant differential missingness (BH < ${params.batch_alpha}):",
        sum(fail_miss), "\\n\\n")

    cat("--- Test 2: Differential Allele Frequency ---\\n")
    af_pvals <- rep(1.0, n_snps)

    for (i in seq_len(n_snps)) {
        tab <- cbind(ref_counts[, i], alt_counts[, i])
        # Only test platforms with genotyped samples
        has_data <- rowSums(tab) > 0
        if (sum(has_data) < 2) next
        tab <- tab[has_data, , drop = FALSE]
        # Need both alleles observed somewhere (otherwise monomorphic across all platforms)
        if (any(colSums(tab) == 0)) next
        af_pvals[i] <- tryCatch(fisher.test(tab)\$p.value, error = function(e) 1.0)
    }

    af_padj <- p.adjust(af_pvals, method = "BH")
    fail_af <- af_padj < ${params.batch_alpha}
    cat("SNPs with significant differential allele frequency (BH < ${params.batch_alpha}):",
        sum(fail_af), "\\n\\n")

    # Combined filter
    fail_batch <- fail_miss | fail_af
    n_fail <- sum(fail_batch)

    cat("--- Summary ---\\n")
    cat("SNPs flagged (missingness OR allele freq):", n_fail, "\\n")
    cat("  Missingness only:", sum(fail_miss & !fail_af), "\\n")
    cat("  Allele frequency only:", sum(fail_af & !fail_miss), "\\n")
    cat("  Both:", sum(fail_miss & fail_af), "\\n")
    cat("SNPs retained:", sum(!fail_batch), "of", n_snps,
        "(", round(100 * sum(!fail_batch) / n_snps, 1), "%)\\n\\n")

    # Report per-platform allele frequencies for the most divergent removed SNPs
    if (n_fail > 0 && n_fail <= 50) {
        cat("--- Removed SNPs Detail ---\\n")
        fail_idx <- which(fail_batch)
        for (fi in fail_idx) {
            afs <- ref_counts[, fi] / (ref_counts[, fi] + alt_counts[, fi])
            afs[is.nan(afs)] <- NA
            af_str <- paste(platforms, "=",
                            ifelse(is.na(afs), "NA", sprintf("%.3f", afs)),
                            collapse = ", ")
            miss_str <- paste(platforms, "=",
                              sprintf("%d/%d", miss_counts[, fi], n_per_plat),
                              collapse = ", ")
            cat(snp_cols[fi], "\\n")
            cat("  Ref AF:", af_str, "\\n")
            cat("  Missing:", miss_str, "\\n")
            cat("  p_miss=", sprintf("%.2e", miss_padj[fi]),
                " p_af=", sprintf("%.2e", af_padj[fi]), "\\n")
        }
    } else if (n_fail > 50) {
        cat("--- Top 50 Most Divergent Removed SNPs (by AF p-value) ---\\n")
        fail_idx <- which(fail_batch)
        ord <- fail_idx[order(af_pvals[fail_idx])]
        for (fi in head(ord, 50)) {
            afs <- ref_counts[, fi] / (ref_counts[, fi] + alt_counts[, fi])
            afs[is.nan(afs)] <- NA
            af_str <- paste(platforms, "=",
                            ifelse(is.na(afs), "NA", sprintf("%.3f", afs)),
                            collapse = ", ")
            cat(snp_cols[fi], ": ", af_str,
                " (p_miss=", sprintf("%.2e", miss_padj[fi]),
                ", p_af=", sprintf("%.2e", af_padj[fi]), ")\\n")
        }
    }

    cat("\\nOutput:", nrow(geno), "samples x", sum(!fail_batch), "SNPs\\n")
    sink()

    # Filter
    kept <- snp_cols[!fail_batch]
    geno_corrected <- geno[, c("sample_id", kept), with = FALSE]
    fwrite(geno_corrected, "batch_corrected_matrix.tsv", sep = "\\t")
    """
}
