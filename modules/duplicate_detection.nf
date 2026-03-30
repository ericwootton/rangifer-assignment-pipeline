/*
========================================================================================
    DUPLICATE DETECTION
========================================================================================
    Calculate pairwise IBS. Pairs with IBS > threshold are likely the same individual.
    Retain the sample with the higher call rate.
*/

process DUPLICATE_DETECTION {
    tag "duplicate_detection"
    label 'process_medium'

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path genotype_matrix

    output:
    path "dedup_genotype_matrix.tsv", emit: filtered_matrix
    path "snp_info_dedup.tsv",        emit: snp_info
    path "duplicate_pairs.tsv",       emit: duplicate_pairs
    path "duplicate_report.txt",      emit: report

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    sink("duplicate_report.txt")
    cat("=== Duplicate Detection Report ===\\n\\n")

    geno <- fread("${genotype_matrix}", header = TRUE)
    cat("Input:", nrow(geno), "samples\\n")

    sample_ids <- geno\$sample_id
    geno_mat <- as.matrix(geno[, -1])
    n <- nrow(geno_mat)

    ibs_threshold <- ${params.ibs_threshold}

    # Short-circuit: skip pairwise IBS when threshold >= 1.0 (disabled)
    if (ibs_threshold >= 1.0) {
        cat("IBS duplicate detection disabled (threshold >= 1.0)\\n")
        cat("\\nDuplicate pairs found: 0\\n")
        cat("Samples to remove: 0\\n")
        cat("\\nOutput:", n, "samples\\n")
        sink()
        fwrite(geno, "dedup_genotype_matrix.tsv", sep = "\\t")
        fwrite(data.table(sample1 = character(), sample2 = character(), ibs = numeric()),
               "duplicate_pairs.tsv", sep = "\\t")
        snp_dt <- data.table(snp_id = colnames(geno)[-1])
        fwrite(snp_dt, "snp_info_dedup.tsv", sep = "\\t")
        quit(save = "no")
    }

    # Calculate call rates for tiebreaking
    call_rates <- rowSums(!is.na(geno_mat)) / ncol(geno_mat)

    # Calculate pairwise IBS (proportion of matching genotypes among shared non-missing)
    # For large matrices, compute in blocks
    dup_pairs <- data.table(sample1 = character(), sample2 = character(), ibs = numeric())
    remove_ids <- character()

    if (n > 1) {
        for (i in 1:(n-1)) {
            for (j in (i+1):n) {
                shared <- !is.na(geno_mat[i,]) & !is.na(geno_mat[j,])
                n_shared <- sum(shared)
                if (n_shared < 100) next  # need minimum shared SNPs
                concordance <- sum(geno_mat[i, shared] == geno_mat[j, shared]) / n_shared
                if (concordance > ${params.ibs_threshold}) {
                    dup_pairs <- rbind(dup_pairs, data.table(
                        sample1 = sample_ids[i],
                        sample2 = sample_ids[j],
                        ibs = round(concordance, 6)
                    ))
                    # Remove the one with lower call rate
                    if (call_rates[i] >= call_rates[j]) {
                        remove_ids <- c(remove_ids, sample_ids[j])
                    } else {
                        remove_ids <- c(remove_ids, sample_ids[i])
                    }
                }
            }
            if (i %% 50 == 0) cat("Processed", i, "of", n-1, "samples...\\n")
        }
    }

    remove_ids <- unique(remove_ids)
    cat("\\nDuplicate pairs found:", nrow(dup_pairs), "\\n")
    cat("Samples to remove:", length(remove_ids), "\\n")

    if (nrow(dup_pairs) > 0) {
        cat("\\nDuplicate pairs:\\n")
        print(dup_pairs)
    }

    # Filter
    geno_filtered <- geno[!sample_id %in% remove_ids]
    cat("\\nOutput:", nrow(geno_filtered), "samples\\n")
    sink()

    fwrite(geno_filtered, "dedup_genotype_matrix.tsv", sep = "\\t")
    fwrite(dup_pairs, "duplicate_pairs.tsv", sep = "\\t")

    # Pass through snp_info (just copy column names as IDs)
    snp_dt <- data.table(snp_id = colnames(geno_filtered)[-1])
    fwrite(snp_dt, "snp_info_dedup.tsv", sep = "\\t")
    """
}
