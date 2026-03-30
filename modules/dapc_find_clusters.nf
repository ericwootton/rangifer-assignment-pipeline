/*
========================================================================================
    DAPC CLUSTER DISCOVERY
========================================================================================
    Use find.clusters() to identify optimal K without predefined labels, then
    compare discovered clusters to labeled groups.
*/

process DAPC_FIND_CLUSTERS {
    tag "dapc_clusters"
    label 'process_medium'

    publishDir "${params.outdir}/structure_discovery/dapc_clusters", mode: 'copy'

    input:
    path genotype_matrix
    path metadata

    output:
    path "cluster_results.rds",  emit: cluster_results
    path "bic_values.tsv",       emit: bic_values
    path "cluster_membership.tsv", emit: membership
    path "cluster_plots.pdf",    emit: plots

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(adegenet)

    geno <- fread("${genotype_matrix}", header = TRUE)
    meta <- fread("${metadata}", header = TRUE)

    # Normalize metadata column names
    if ("Sample" %in% colnames(meta) && !"sample_id" %in% colnames(meta)) setnames(meta, "Sample", "sample_id")
    if ("Subspecies" %in% colnames(meta) && !"subspecies" %in% colnames(meta)) setnames(meta, "Subspecies", "subspecies")
    if ("Ecotype_Cleaned" %in% colnames(meta) && !"ecotype" %in% colnames(meta)) setnames(meta, "Ecotype_Cleaned", "ecotype")
    if ("Herd_Cleaned" %in% colnames(meta) && !"herd" %in% colnames(meta)) setnames(meta, "Herd_Cleaned", "herd")

    cat("Running DAPC cluster discovery on", nrow(geno), "samples\\n")

    # Convert to matrix
    geno_mat <- as.matrix(geno[, -1])
    rownames(geno_mat) <- geno\$sample_id

    # Mean-impute for DAPC
    for (j in 1:ncol(geno_mat)) {
        col_mean <- mean(geno_mat[, j], na.rm = TRUE)
        geno_mat[is.na(geno_mat[, j]), j] <- round(col_mean)
    }

    # Create genlight object
    gl <- new("genlight", geno_mat)
    indNames(gl) <- geno\$sample_id

    # Find clusters: test K from 1 to max_k
    max_k <- min(${params.max_k_dapc}, nrow(geno) - 1)
    n_pca <- min(nrow(geno) - 1, ncol(geno_mat) - 1, 300)

    set.seed(42)
    grp <- find.clusters(gl, max.n.clust = max_k, n.pca = n_pca,
                         choose.n.clust = FALSE, criterion = "min")

    # BIC values
    bic_dt <- data.table(K = 1:length(grp\$Kstat), BIC = grp\$Kstat)
    optimal_k <- bic_dt[which.min(BIC)]\$K
    cat("Optimal K by BIC:", optimal_k, "\\n")

    # Cluster membership
    membership <- data.table(
        sample_id = geno\$sample_id,
        cluster = as.integer(grp\$grp)
    )

    # Compare to metadata labels
    id_col <- intersect(c("sample_id", "ID", "Sample_ID"), colnames(meta))[1]
    if (!is.na(id_col) && id_col != "sample_id") setnames(meta, id_col, "sample_id")
    membership <- merge(membership, meta, by = "sample_id", all.x = TRUE)

    # Save
    fwrite(bic_dt, "bic_values.tsv", sep = "\\t")
    fwrite(membership, "cluster_membership.tsv", sep = "\\t")
    saveRDS(list(grp = grp, bic = bic_dt, membership = membership, optimal_k = optimal_k),
            "cluster_results.rds")

    # Plots
    pdf("cluster_plots.pdf", width = 12, height = 8)

    # BIC plot
    plot(bic_dt\$K, bic_dt\$BIC, type = "b", pch = 19,
         xlab = "Number of Clusters (K)", ylab = "BIC",
         main = "DAPC Cluster Selection (BIC)")
    abline(v = optimal_k, col = "red", lty = 2)
    text(optimal_k, max(bic_dt\$BIC, na.rm = TRUE), paste("K =", optimal_k), col = "red", pos = 4)

    dev.off()
    cat("DAPC cluster discovery complete\\n")
    """
}
