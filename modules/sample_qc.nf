/*
========================================================================================
    SAMPLE-LEVEL QUALITY CONTROL
========================================================================================
    1. Call rate filter: remove samples < 70% genotyped (platform-adjusted)
    2. Heterozygosity filter: flag samples > max_het (possible mixed/contaminated)
    3. Species confirmation: flag samples < min_het (possible non-Rangifer)

    Platform-aware call rates: for multi-platform datasets, the call rate denominator
    is the number of SNPs genotypeable by that sample's platform, not total SNPs.
    This prevents penalizing chip samples for SNPs absent from their chip design.
*/

process SAMPLE_QC {
    tag "sample_qc"
    label 'process_medium'

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path genotype_matrix
    path data_source

    output:
    path "sample_qc_filtered.tsv",  emit: filtered_matrix
    path "sample_stats.tsv",        emit: sample_stats
    path "sample_qc_report.txt",    emit: report

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    sink("sample_qc_report.txt")
    cat("=== Sample-Level QC Report ===\\n\\n")

    geno <- fread("${genotype_matrix}", header = TRUE)
    cat("Input:", nrow(geno), "samples x", ncol(geno) - 1, "SNPs\\n\\n")

    sample_ids <- geno\$sample_id
    geno_mat <- as.matrix(geno[, -1])

    n_snps <- ncol(geno_mat)

    # Load data source mapping
    ds <- fread("${data_source}", header = TRUE)

    # Calculate platform-adjusted call rates
    # For each source group, detect platform-absent SNPs (NA in ALL samples of that source)
    # These are SNPs not on the chip, not quality failures
    n_genotypeable <- rep(n_snps, length(sample_ids))
    names(n_genotypeable) <- sample_ids

    for (src in unique(ds\$source)) {
        src_ids <- ds[source == src]\$sample_id
        src_ids <- intersect(src_ids, sample_ids)
        if (length(src_ids) == 0) next

        src_rows <- which(sample_ids %in% src_ids)
        if (length(src_rows) == 0) next

        src_mat <- geno_mat[src_rows, , drop = FALSE]
        # SNPs where ALL samples from this source are NA = platform-absent
        platform_absent <- colSums(!is.na(src_mat)) == 0
        n_geno <- n_snps - sum(platform_absent)

        cat("Source:", src, "- samples:", length(src_ids),
            "- genotypeable SNPs:", n_geno, "of", n_snps,
            "(", sum(platform_absent), "platform-absent)\\n")

        for (sid in src_ids) {
            n_genotypeable[sid] <- n_geno
        }
    }

    # Per-sample statistics
    n_called <- rowSums(!is.na(geno_mat))
    raw_call_rate <- n_called / n_snps
    adj_call_rate <- n_called / n_genotypeable[sample_ids]
    het_rate <- rowSums(geno_mat == 1, na.rm = TRUE) / n_called

    stats <- data.table(
        sample_id       = sample_ids,
        source          = ds\$source[match(sample_ids, ds\$sample_id)],
        call_rate       = round(adj_call_rate, 4),
        raw_call_rate   = round(raw_call_rate, 4),
        het_rate        = round(het_rate, 4),
        n_genotyped     = n_called,
        n_genotypeable  = n_genotypeable[sample_ids],
        n_missing       = as.integer(n_genotypeable[sample_ids] - n_called)
    )

    # Apply filters using platform-adjusted call rate
    stats[, pass_call_rate := call_rate >= ${params.min_call_rate}]
    stats[, pass_het_high  := het_rate <= ${params.max_het}]
    stats[, pass_het_low   := het_rate >= ${params.min_het}]
    stats[, flag_mixed     := het_rate > ${params.max_het}]
    stats[, flag_non_rangifer := het_rate < ${params.min_het} & call_rate >= 0.5]
    stats[, pass_all := pass_call_rate & pass_het_high & pass_het_low]

    cat("\\nCall rate filter (<", ${params.min_call_rate}, ", platform-adjusted):", sum(!stats\$pass_call_rate), "removed\\n")
    cat("High het filter (>", ${params.max_het}, "):", sum(!stats\$pass_het_high), "flagged/removed\\n")
    cat("Low het filter (<", ${params.min_het}, "):", sum(!stats\$pass_het_low), "flagged/removed\\n")
    cat("Samples passing all filters:", sum(stats\$pass_all), "\\n\\n")

    # Report samples that are rescued by platform-aware call rate
    rescued <- stats[raw_call_rate < ${params.min_call_rate} & call_rate >= ${params.min_call_rate}]
    if (nrow(rescued) > 0) {
        cat("Samples rescued by platform-adjusted call rate:", nrow(rescued), "\\n")
        cat("  (raw call rate below threshold, but platform-adjusted rate passes)\\n")
        for (r in 1:nrow(rescued)) {
            cat("  ", rescued\$sample_id[r], ": raw=", rescued\$raw_call_rate[r],
                "adj=", rescued\$call_rate[r], "(", rescued\$source[r], ")\\n")
        }
        cat("\\n")
    }

    # Report flagged samples
    if (any(stats\$flag_mixed)) {
        cat("Flagged as possible mixed samples (het >", ${params.max_het}, "):\\n")
        cat(paste(stats[flag_mixed == TRUE]\$sample_id, collapse = "\\n"), "\\n\\n")
    }
    if (any(stats\$flag_non_rangifer)) {
        cat("Flagged as possible non-Rangifer (het <", ${params.min_het}, "):\\n")
        cat(paste(stats[flag_non_rangifer == TRUE]\$sample_id, collapse = "\\n"), "\\n\\n")
    }

    # Report removed samples
    removed <- stats[pass_all == FALSE]
    if (nrow(removed) > 0) {
        cat("Removed samples:\\n")
        for (r in 1:nrow(removed)) {
            reasons <- c()
            if (!removed\$pass_call_rate[r]) reasons <- c(reasons, paste0("call_rate=", removed\$call_rate[r]))
            if (!removed\$pass_het_high[r]) reasons <- c(reasons, paste0("het=", removed\$het_rate[r]))
            if (!removed\$pass_het_low[r]) reasons <- c(reasons, paste0("low_het=", removed\$het_rate[r]))
            cat("  ", removed\$sample_id[r], "(", removed\$source[r], "):", paste(reasons, collapse = ", "), "\\n")
        }
        cat("\\n")
    }

    # Filter
    keep <- stats[pass_all == TRUE]\$sample_id
    geno_filtered <- geno[sample_id %in% keep]

    cat("Output:", nrow(geno_filtered), "samples retained\\n")
    sink()

    # Write outputs
    fwrite(geno_filtered, "sample_qc_filtered.tsv", sep = "\\t")
    fwrite(stats, "sample_stats.tsv", sep = "\\t")
    """
}
