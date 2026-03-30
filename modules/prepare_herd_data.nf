/*
========================================================================================
    PREPARE HERD DATA
========================================================================================
    For the herd differentiation workflow (-entry HERD_DIFF):
    Takes a pre-QC'd genotype matrix and metadata with individual herd labels.
    Subsets to overlapping samples and creates a single flat node with herds as groups.
    Applies node-level SNP QC (call rate >= 80%, remove monomorphic).
*/

process PREPARE_HERD_DATA {
    tag "herd_diff"
    label 'process_medium'

    publishDir "${params.outdir}/herd_differentiation/data", mode: 'copy'

    input:
    path genotype_matrix
    path metadata

    output:
    path "herd_genotype_matrix.tsv", emit: herd_matrix
    path "herd_node",                emit: herd_node_dir

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(jsonlite)

    geno <- fread("${genotype_matrix}", header = TRUE)
    meta <- fread("${metadata}", header = TRUE)

    # Normalize column names
    if ("Sample" %in% colnames(meta) && !"sample_id" %in% colnames(meta)) setnames(meta, "Sample", "sample_id")
    if ("Herd_Cleaned" %in% colnames(meta) && !"herd" %in% colnames(meta)) setnames(meta, "Herd_Cleaned", "herd")

    id_col <- intersect(c("sample_id", "ID", "Sample_ID"), colnames(meta))[1]
    if (!is.na(id_col) && id_col != "sample_id") setnames(meta, id_col, "sample_id")

    # Subset to samples in both genotype matrix and metadata
    shared_samples <- intersect(geno[["sample_id"]], meta[["sample_id"]])
    cat("Samples in genotype matrix:", nrow(geno), "\\n")
    cat("Samples in metadata:", nrow(meta), "\\n")
    cat("Overlapping samples:", length(shared_samples), "\\n")

    if (length(shared_samples) < 2) {
        stop("Fewer than 2 overlapping samples between genotype matrix and metadata")
    }

    # Subset and remove samples without herd labels
    meta_sub <- meta[sample_id %in% shared_samples & !is.na(herd) & herd != ""]
    geno_sub <- geno[sample_id %in% meta_sub[["sample_id"]]]

    cat("\\nSamples with herd labels:", nrow(geno_sub), "\\n")
    cat("Herds found:\\n")
    herd_counts <- table(meta_sub[["herd"]])
    for (h in names(sort(herd_counts, decreasing = TRUE))) {
        cat("  ", h, ":", herd_counts[h], "\\n")
    }

    # Output 1: Full subsetted genotype matrix (for structure discovery - PCA, DAPC, ADMIXTURE)
    fwrite(geno_sub, "herd_genotype_matrix.tsv", sep = "\\t")

    # Output 2: Node directory with node-level SNP QC
    dir.create("herd_node", showWarnings = FALSE)

    ref_mat <- as.matrix(geno_sub[, -1])

    # Node SNP QC
    snp_call_rate <- colSums(!is.na(ref_mat)) / nrow(ref_mat)
    pass_cr <- snp_call_rate >= ${params.node_snp_call_rate}

    allele_freq <- colMeans(ref_mat, na.rm = TRUE) / 2
    is_mono <- allele_freq == 0 | allele_freq == 1
    is_mono[is.na(is_mono)] <- TRUE

    keep <- pass_cr & !is_mono
    kept_snps <- colnames(geno_sub)[-1][keep]

    cat("\\nSNP QC within herd node:\\n")
    cat("  Total SNPs:", ncol(ref_mat), "\\n")
    cat("  Failed call rate:", sum(!pass_cr), "\\n")
    cat("  Monomorphic:", sum(is_mono & pass_cr), "\\n")
    cat("  Retained:", sum(keep), "\\n")

    # Write reference genotypes
    ref_out <- geno_sub[, c("sample_id", kept_snps), with = FALSE]
    fwrite(ref_out, file.path("herd_node", "reference_genotypes.tsv"), sep = "\\t")

    # Write reference labels
    ref_labels <- meta_sub[, c("sample_id", "herd")]
    setnames(ref_labels, "herd", "group")
    fwrite(ref_labels, file.path("herd_node", "reference_labels.tsv"), sep = "\\t")

    # Write node metadata JSON
    sample_counts <- as.list(table(ref_labels[["group"]]))
    node_meta <- list(
        node_name = "herd_differentiation",
        level = 1L,
        column = "herd",
        groups = names(sample_counts),
        n_reference = nrow(ref_out),
        n_unknown = 0L,
        n_snps = length(kept_snps),
        min_class_size = as.integer(min(unlist(sample_counts))),
        sample_counts = sample_counts
    )
    write_json(node_meta, file.path("herd_node", "node_meta.json"),
               pretty = TRUE, auto_unbox = TRUE)

    cat("\\nHerd node prepared successfully\\n")
    """
}
