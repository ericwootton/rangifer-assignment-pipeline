/*
========================================================================================
    FST NODE FILTER (SNP mode)
========================================================================================
    At each hierarchy node, compute per-SNP Fst (Weir & Cockerham) between groups,
    rank by Fst, and select top 5,000 differentiating SNPs for classification.
    For multi-group nodes, uses average pairwise Fst per SNP.
*/

process FST_NODE_FILTER {
    tag "${node_name}"

    input:
    tuple val(node_name), path(node_dir), path(node_meta)

    output:
    tuple val(node_name), path("${node_name}_fst_filtered"), emit: filtered_node

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(jsonlite)

    meta <- fromJSON("${node_meta}")
    node_name <- meta[["node_name"]]

    # Read node data
    ref_geno <- fread(file.path("${node_dir}", "reference_genotypes.tsv"), header = TRUE)
    ref_labels <- fread(file.path("${node_dir}", "reference_labels.tsv"), header = TRUE)

    # Merge labels
    ref_labels_map <- setNames(ref_labels[["group"]], ref_labels[["sample_id"]])
    groups <- ref_labels_map[ref_geno[["sample_id"]]]

    # Extract numeric genotype matrix (samples x SNPs)
    geno_mat <- as.matrix(ref_geno[, -1])
    snp_names <- colnames(ref_geno)[-1]
    unique_groups <- unique(groups)
    n_groups <- length(unique_groups)

    cat("Node:", node_name, "\\n")
    cat("  Samples:", nrow(geno_mat), "Groups:", n_groups, "SNPs:", ncol(geno_mat), "\\n")

    # Weir & Cockerham Fst per SNP
    # For each pair of groups, compute per-SNP Fst
    compute_pairwise_fst <- function(g1_mat, g2_mat) {
        n1 <- nrow(g1_mat)
        n2 <- nrow(g2_mat)
        n_bar <- (n1 + n2) / 2
        n_c <- n_bar  # simplified for 2 pops

        fst_vals <- numeric(ncol(g1_mat))
        for (j in seq_len(ncol(g1_mat))) {
            p1 <- mean(g1_mat[, j], na.rm = TRUE) / 2
            p2 <- mean(g2_mat[, j], na.rm = TRUE) / 2
            p_bar <- (n1 * p1 + n2 * p2) / (n1 + n2)

            # Between-group variance
            s2 <- (n1 * (p1 - p_bar)^2 + n2 * (p2 - p_bar)^2) / n_bar

            # Expected heterozygosity
            h_bar <- (n1 * 2 * p1 * (1 - p1) + n2 * 2 * p2 * (1 - p2)) / (2 * (n1 + n2))

            # Fst components (simplified Weir & Cockerham)
            num <- s2 - (1 / (n_bar - 1)) * (p_bar * (1 - p_bar) - s2 / 2 - h_bar / 4)
            den <- p_bar * (1 - p_bar)

            fst_vals[j] <- ifelse(den > 0, num / den, 0)
        }
        fst_vals[is.na(fst_vals)] <- 0
        fst_vals[fst_vals < 0] <- 0
        return(fst_vals)
    }

    # Compute average pairwise Fst across all group pairs
    fst_sum <- numeric(ncol(geno_mat))
    n_pairs <- 0

    if (n_groups >= 2) {
        for (i in seq_len(n_groups - 1)) {
            for (j in (i + 1):n_groups) {
                g1_idx <- which(groups == unique_groups[i])
                g2_idx <- which(groups == unique_groups[j])
                if (length(g1_idx) >= 2 && length(g2_idx) >= 2) {
                    pair_fst <- compute_pairwise_fst(geno_mat[g1_idx, , drop = FALSE],
                                                     geno_mat[g2_idx, , drop = FALSE])
                    fst_sum <- fst_sum + pair_fst
                    n_pairs <- n_pairs + 1
                }
            }
        }
    }

    if (n_pairs > 0) {
        avg_fst <- fst_sum / n_pairs
    } else {
        avg_fst <- rep(0, ncol(geno_mat))
    }

    # Rank and select top 5000
    n_select <- min(5000, length(avg_fst))
    top_idx <- order(avg_fst, decreasing = TRUE)[seq_len(n_select)]
    selected_snps <- snp_names[top_idx]

    cat("  Average Fst range: [", round(min(avg_fst), 4), ",", round(max(avg_fst), 4), "]\\n")
    cat("  Selected top", n_select, "SNPs by Fst\\n")

    # Create output directory with filtered data
    out_dir <- paste0(node_name, "_fst_filtered")
    dir.create(out_dir, showWarnings = FALSE)

    # Write filtered reference genotypes
    ref_out <- ref_geno[, c("sample_id", selected_snps), with = FALSE]
    fwrite(ref_out, file.path(out_dir, "reference_genotypes.tsv"), sep = "\\t")

    # Write filtered unknown genotypes if they exist
    unk_file <- file.path("${node_dir}", "unknown_genotypes.tsv")
    if (file.exists(unk_file)) {
        unk_geno <- fread(unk_file, header = TRUE)
        unk_out <- unk_geno[, c("sample_id", intersect(selected_snps, colnames(unk_geno))), with = FALSE]
        fwrite(unk_out, file.path(out_dir, "unknown_genotypes.tsv"), sep = "\\t")
    }

    # Copy labels and update node_meta
    file.copy(file.path("${node_dir}", "reference_labels.tsv"), file.path(out_dir, "reference_labels.tsv"))
    meta[["n_snps"]] <- n_select
    meta[["fst_filtered"]] <- TRUE
    meta[["fst_range"]] <- c(round(min(avg_fst), 6), round(max(avg_fst), 6))
    write_json(meta, file.path(out_dir, "node_meta.json"), pretty = TRUE, auto_unbox = TRUE)

    # Write Fst values for diagnostics
    fst_df <- data.table(snp = snp_names, fst = round(avg_fst, 6))
    fst_df <- fst_df[order(-fst)]
    fwrite(fst_df, file.path(out_dir, "snp_fst_values.tsv"), sep = "\\t")
    """
}
