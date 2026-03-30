/*
========================================================================================
    PAIRWISE FST
========================================================================================
    Computes pairwise Weir & Cockerham Fst between all herds.
    Outputs a pairwise matrix, long-format table, heatmap, and per-SNP Fst per pair.
*/

process PAIRWISE_FST {
    tag "pairwise_fst"
    label 'process_high'

    publishDir "${params.outdir}/herd_differentiation/pairwise_fst", mode: 'copy'

    input:
    tuple path(node_dir), path(node_meta)

    output:
    path "pairwise_fst_matrix.tsv", emit: fst_matrix
    path "pairwise_fst_long.tsv",   emit: fst_long
    path "pairwise_fst_heatmap.pdf", emit: fst_heatmap
    path "per_snp_fst",              emit: per_snp_fst, type: 'dir'

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(jsonlite)

    ref_geno <- fread(file.path("${node_dir}", "reference_genotypes.tsv"), header = TRUE)
    ref_labels <- fread(file.path("${node_dir}", "reference_labels.tsv"), header = TRUE)

    labels_map <- setNames(ref_labels[["group"]], ref_labels[["sample_id"]])
    groups <- labels_map[ref_geno[["sample_id"]]]
    unique_groups <- sort(unique(groups))
    n_groups <- length(unique_groups)

    geno_mat <- as.matrix(ref_geno[, -1])
    snp_names <- colnames(ref_geno)[-1]

    cat("Computing pairwise Fst between", n_groups, "herds\\n")
    cat("Samples:", nrow(geno_mat), " SNPs:", ncol(geno_mat), "\\n\\n")

    # Weir & Cockerham Fst per SNP for a pair of populations
    compute_pair_fst <- function(g1_mat, g2_mat) {
        n1 <- nrow(g1_mat)
        n2 <- nrow(g2_mat)
        n_bar <- (n1 + n2) / 2

        fst_vals <- numeric(ncol(g1_mat))
        for (j in seq_len(ncol(g1_mat))) {
            p1 <- mean(g1_mat[, j], na.rm = TRUE) / 2
            p2 <- mean(g2_mat[, j], na.rm = TRUE) / 2
            p_bar <- (n1 * p1 + n2 * p2) / (n1 + n2)

            s2 <- (n1 * (p1 - p_bar)^2 + n2 * (p2 - p_bar)^2) / n_bar
            h_bar <- (n1 * 2 * p1 * (1 - p1) + n2 * 2 * p2 * (1 - p2)) / (2 * (n1 + n2))

            num <- s2 - (1 / (n_bar - 1)) * (p_bar * (1 - p_bar) - s2 / 2 - h_bar / 4)
            den <- p_bar * (1 - p_bar)

            fst_vals[j] <- ifelse(den > 0, num / den, 0)
        }
        fst_vals[is.na(fst_vals)] <- 0
        return(fst_vals)
    }

    fst_matrix <- matrix(0, n_groups, n_groups,
                         dimnames = list(unique_groups, unique_groups))
    dir.create("per_snp_fst", showWarnings = FALSE)
    results_long <- list()

    for (i in seq_len(n_groups - 1)) {
        for (j in (i + 1):n_groups) {
            g1 <- unique_groups[i]
            g2 <- unique_groups[j]
            g1_idx <- which(groups == g1)
            g2_idx <- which(groups == g2)

            cat("  ", g1, "(n=", length(g1_idx), ") vs ", g2, "(n=", length(g2_idx), ")")

            if (length(g1_idx) < 2 || length(g2_idx) < 2) {
                cat(" -- SKIPPED (insufficient samples)\\n")
                fst_matrix[i, j] <- NA
                fst_matrix[j, i] <- NA
                next
            }

            snp_fst <- compute_pair_fst(
                geno_mat[g1_idx, , drop = FALSE],
                geno_mat[g2_idx, , drop = FALSE]
            )

            snp_fst_clean <- pmax(snp_fst, 0)
            genome_fst <- mean(snp_fst_clean)

            fst_matrix[i, j] <- genome_fst
            fst_matrix[j, i] <- genome_fst

            cat(" Fst =", round(genome_fst, 4), "\\n")

            results_long[[length(results_long) + 1]] <- data.table(
                pop1 = g1, pop2 = g2,
                n1 = length(g1_idx), n2 = length(g2_idx),
                fst = round(genome_fst, 6),
                median_snp_fst = round(median(snp_fst_clean), 6),
                max_snp_fst = round(max(snp_fst_clean), 6),
                n_snps = length(snp_fst)
            )

            # Per-SNP Fst for this pair
            pair_name <- paste0(gsub(" ", "_", g1), "_vs_", gsub(" ", "_", g2))
            snp_dt <- data.table(snp = snp_names, fst = round(snp_fst, 6))
            snp_dt <- snp_dt[order(-fst)]
            fwrite(snp_dt, file.path("per_snp_fst", paste0(pair_name, ".tsv")),
                   sep = "\\t")
        }
    }

    # Write outputs
    fst_df <- as.data.table(fst_matrix, keep.rownames = "Herd")
    fwrite(fst_df, "pairwise_fst_matrix.tsv", sep = "\\t")

    results_dt <- rbindlist(results_long)
    fwrite(results_dt, "pairwise_fst_long.tsv", sep = "\\t")

    # ---- Heatmap ----
    pdf("pairwise_fst_heatmap.pdf", width = 10, height = 9)

    # Cluster rows/cols by Fst distance
    plot_mat <- fst_matrix
    plot_mat[is.na(plot_mat)] <- 0
    if (n_groups >= 3) {
        hc <- hclust(as.dist(plot_mat), method = "average")
        ord <- hc[["order"]]
    } else {
        ord <- seq_len(n_groups)
    }
    plot_mat <- plot_mat[ord, ord]
    n <- nrow(plot_mat)

    cols <- colorRampPalette(c("white", "#FFF7BC", "#FEC44F", "#D95F0E", "#7F0000"))(100)
    max_fst <- max(plot_mat[upper.tri(plot_mat)], na.rm = TRUE)
    if (max_fst == 0) max_fst <- 0.01

    par(mar = c(10, 10, 3, 2))
    image(1:n, 1:n, t(plot_mat)[, n:1], col = cols, zlim = c(0, max_fst),
          axes = FALSE, xlab = "", ylab = "",
          main = "Pairwise Fst Between Herds")
    axis(1, at = 1:n, labels = rownames(plot_mat), las = 2, cex.axis = 0.75)
    axis(2, at = 1:n, labels = rev(colnames(plot_mat)), las = 2, cex.axis = 0.75)

    for (ii in 1:n) {
        for (jj in 1:n) {
            if (ii != jj) {
                val <- plot_mat[ii, jj]
                text(ii, n - jj + 1, sprintf("%.3f", val), cex = 0.55)
            }
        }
    }
    box()
    dev.off()

    cat("\\nPairwise Fst computation complete\\n")
    """
}
