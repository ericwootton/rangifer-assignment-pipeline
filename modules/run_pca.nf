/*
========================================================================================
    PCA - STRUCTURE DISCOVERY
========================================================================================
    Run PCA on all reference samples, plot PC1 vs PC2, PC1 vs PC3, PC2 vs PC3,
    colored by labeled subspecies, ecotype, and herd.
*/

process RUN_PCA {
    tag "pca"
    label 'process_medium'

    publishDir "${params.outdir}/structure_discovery/pca", mode: 'copy'

    input:
    path genotype_matrix
    path metadata

    output:
    path "pca_results.rds",    emit: pca_results
    path "pca_eigenvalues.tsv", emit: eigenvalues
    path "pca_scores.tsv",      emit: scores
    path "pca_plots.pdf",       emit: plots

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(adegenet)

    # Load data
    geno <- fread("${genotype_matrix}", header = TRUE)
    meta <- fread("${metadata}", header = TRUE)

    # Normalize metadata column names
    if ("Sample" %in% colnames(meta) && !"sample_id" %in% colnames(meta)) setnames(meta, "Sample", "sample_id")
    if ("Subspecies" %in% colnames(meta) && !"subspecies" %in% colnames(meta)) setnames(meta, "Subspecies", "subspecies")
    if ("Ecotype_Cleaned" %in% colnames(meta) && !"ecotype" %in% colnames(meta)) setnames(meta, "Ecotype_Cleaned", "ecotype")
    if ("Herd_Cleaned" %in% colnames(meta) && !"herd" %in% colnames(meta)) setnames(meta, "Herd_Cleaned", "herd")

    cat("Genotype matrix:", nrow(geno), "samples x", ncol(geno) - 1, "SNPs\\n")

    # Convert to genlight object
    geno_mat <- as.matrix(geno[, -1])
    rownames(geno_mat) <- geno\$sample_id

    # Mean-impute missing values for PCA
    for (j in 1:ncol(geno_mat)) {
        col_mean <- mean(geno_mat[, j], na.rm = TRUE)
        geno_mat[is.na(geno_mat[, j]), j] <- round(col_mean)
    }

    # Run PCA
    pca <- prcomp(geno_mat, center = TRUE, scale. = FALSE)

    # Extract scores
    scores <- data.table(
        sample_id = geno\$sample_id,
        as.data.table(pca\$x[, 1:min(20, ncol(pca\$x))])
    )

    # Merge with metadata for plotting
    # Auto-detect ID column
    id_col <- intersect(c("sample_id", "ID", "Sample_ID", "SampleID"), colnames(meta))[1]
    if (!is.na(id_col) && id_col != "sample_id") {
        setnames(meta, id_col, "sample_id")
    }
    scores_meta <- merge(scores, meta, by = "sample_id", all.x = TRUE)

    # Eigenvalues
    var_explained <- (pca\$sdev^2) / sum(pca\$sdev^2) * 100
    eigen_dt <- data.table(
        PC = paste0("PC", 1:length(var_explained)),
        eigenvalue = pca\$sdev^2,
        var_percent = round(var_explained, 2),
        cumulative_var = round(cumsum(var_explained), 2)
    )

    # Save outputs
    fwrite(scores, "pca_scores.tsv", sep = "\\t")
    fwrite(eigen_dt, "pca_eigenvalues.tsv", sep = "\\t")
    saveRDS(list(pca = pca, scores_meta = scores_meta, eigenvalues = eigen_dt), "pca_results.rds")

    # Plots
    pdf("pca_plots.pdf", width = 12, height = 10)

    # Detect grouping columns
    group_cols <- intersect(c("subspecies", "ecotype", "herd", "population",
                              "lineage", "Source_ID", "group"), colnames(scores_meta))

    for (gcol in group_cols) {
        groups <- scores_meta[[gcol]]
        groups[is.na(groups)] <- "Unknown"
        n_groups <- length(unique(groups))
        cols <- rainbow(n_groups)
        names(cols) <- unique(groups)

        # PC1 vs PC2
        par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

        plot(scores_meta\$PC1, scores_meta\$PC2, col = cols[groups], pch = 19, cex = 0.8,
             xlab = paste0("PC1 (", eigen_dt\$var_percent[1], "%)"),
             ylab = paste0("PC2 (", eigen_dt\$var_percent[2], "%)"),
             main = paste0("PCA colored by ", gcol))
        legend("topright", legend = names(cols), col = cols, pch = 19, cex = 0.6, ncol = ceiling(n_groups/10))

        # PC1 vs PC3
        plot(scores_meta\$PC1, scores_meta\$PC3, col = cols[groups], pch = 19, cex = 0.8,
             xlab = paste0("PC1 (", eigen_dt\$var_percent[1], "%)"),
             ylab = paste0("PC3 (", eigen_dt\$var_percent[3], "%)"),
             main = paste0("PC1 vs PC3 by ", gcol))

        # PC2 vs PC3
        plot(scores_meta\$PC2, scores_meta\$PC3, col = cols[groups], pch = 19, cex = 0.8,
             xlab = paste0("PC2 (", eigen_dt\$var_percent[2], "%)"),
             ylab = paste0("PC3 (", eigen_dt\$var_percent[3], "%)"),
             main = paste0("PC2 vs PC3 by ", gcol))

        # Scree plot
        barplot(eigen_dt\$var_percent[1:20], names.arg = 1:20,
                xlab = "PC", ylab = "Variance Explained (%)",
                main = "Scree Plot", col = "steelblue")
    }

    dev.off()
    cat("PCA complete\\n")
    """
}
