/*
========================================================================================
    PREPARE WGS NODE DATA (WGS mode)
========================================================================================
    Merge OutFLANK outlier SNP genotypes + SV genotypes (0/1/2) into unified matrix
    for assignment classifiers at each hierarchy node.
*/

process PREPARE_WGS_NODE_DATA {
    tag "${node_name}"

    input:
    tuple val(node_name), path(outlier_dir)
    path sv_matrix

    output:
    tuple val(node_name), path("${node_name}_wgs_node"), emit: wgs_node_data

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    node_name <- "${node_name}"
    out_dir <- paste0(node_name, "_wgs_node")
    dir.create(out_dir, showWarnings = FALSE)

    # Read outlier SNP genotypes
    outlier_file <- file.path("${outlier_dir}", "outlier_genotypes.tsv")
    outlier_geno <- fread(outlier_file, header = TRUE)

    # Read SV genotype matrix
    sv_file <- "${sv_matrix}"
    sv_geno <- fread(sv_file, header = TRUE)

    cat("Node:", node_name, "\\n")
    cat("  Outlier SNP genotypes:", ncol(outlier_geno) - 1, "SNPs x", nrow(outlier_geno), "samples\\n")
    cat("  SV genotypes:", ncol(sv_geno) - 1, "SVs x", nrow(sv_geno), "samples\\n")

    # Get reference labels for this node
    labels_file <- file.path("${outlier_dir}", "reference_labels.tsv")
    ref_labels <- fread(labels_file, header = TRUE)
    ref_samples <- ref_labels[["sample_id"]]

    # Subset SV matrix to reference samples
    sv_id_col <- colnames(sv_geno)[1]
    if (sv_id_col == "sv_id") {
        # SV matrix is transposed (SVs as rows, samples as columns)
        # Transpose: samples as rows
        sv_ids <- sv_geno[["sv_id"]]
        sv_mat <- as.matrix(sv_geno[, -1])
        rownames(sv_mat) <- sv_ids
        sv_t <- as.data.table(t(sv_mat))
        sv_t[, sample_id := colnames(sv_geno)[-1]]
        setcolorder(sv_t, c("sample_id", sv_ids))
        sv_geno <- sv_t
    }

    # Filter SV data to node samples
    sv_node <- sv_geno[sv_geno[["sample_id"]] %in% ref_samples]

    if (nrow(sv_node) > 0 && ncol(sv_node) > 1) {
        cat("  SV genotypes for node samples:", ncol(sv_node) - 1, "SVs x", nrow(sv_node), "samples\\n")
    }

    # Merge outlier SNPs and SVs by sample_id
    if (ncol(outlier_geno) > 1 && ncol(sv_node) > 1) {
        # Both have data - merge
        combined <- merge(outlier_geno, sv_node, by = "sample_id", all.x = TRUE)
    } else if (ncol(outlier_geno) > 1) {
        combined <- outlier_geno
    } else if (ncol(sv_node) > 1) {
        combined <- sv_node[sv_node[["sample_id"]] %in% ref_samples]
    } else {
        # No features at all - use outlier (may be empty)
        combined <- outlier_geno
    }

    # Ensure only reference samples
    combined <- combined[combined[["sample_id"]] %in% ref_samples]

    cat("  Combined matrix:", ncol(combined) - 1, "features x", nrow(combined), "samples\\n")

    # Write output
    fwrite(combined, file.path(out_dir, "reference_genotypes.tsv"), sep = "\\t")
    file.copy(labels_file, file.path(out_dir, "reference_labels.tsv"))

    # Write unknown genotypes if they exist in the original node data
    unk_file <- file.path("${outlier_dir}", "..", basename("${outlier_dir}"))
    # Unknowns are not in outlier dir - they need to come from the original node
    # For now, create empty placeholder if no unknowns
    # The main workflow will handle unknown routing separately

    # Write node_meta
    meta <- list(
        node_name = node_name,
        n_snps = ncol(combined) - 1,
        n_reference = nrow(combined),
        n_outlier_snps = ncol(outlier_geno) - 1,
        n_sv_features = max(0, ncol(sv_node) - 1),
        sample_counts = as.list(table(ref_labels[["group"]]))
    )
    jsonlite::write_json(meta, file.path(out_dir, "node_meta.json"),
                         pretty = TRUE, auto_unbox = TRUE)
    """
}
