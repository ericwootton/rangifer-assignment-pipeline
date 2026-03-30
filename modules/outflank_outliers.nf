/*
========================================================================================
    OUTFLANK OUTLIER DETECTION (WGS mode, per node)
========================================================================================
    Compute per-SNP Fst, fit neutral distribution with OutFLANK,
    identify outlier SNPs (q-value < 0.05).
*/

process OUTFLANK_OUTLIERS {
    tag "${node_name}"
    label 'process_medium'

    publishDir "${params.outdir}/wgs_analyses/outflank/${node_name}", mode: 'copy'

    input:
    tuple val(node_name), path(node_dir), path(node_meta)

    output:
    tuple val(node_name), path("${node_name}_outliers"), emit: outlier_data

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(jsonlite)
    library(OutFLANK)

    meta <- fromJSON("${node_meta}")
    node_name <- meta[["node_name"]]

    # Read node data
    ref_geno <- fread(file.path("${node_dir}", "reference_genotypes.tsv"), header = TRUE)
    ref_labels <- fread(file.path("${node_dir}", "reference_labels.tsv"), header = TRUE)

    # Build genotype matrix
    geno_mat <- as.matrix(ref_geno[, -1])
    snp_names <- colnames(ref_geno)[-1]
    labels_map <- setNames(ref_labels[["group"]], ref_labels[["sample_id"]])
    pop <- labels_map[ref_geno[["sample_id"]]]

    cat("Node:", node_name, "\\n")
    cat("  Samples:", nrow(geno_mat), "SNPs:", ncol(geno_mat), "\\n")

    # Compute Fst using OutFLANK's MakeDiploidFSTMat
    fst_mat <- MakeDiploidFSTMat(
        SNPmat = geno_mat,
        locusNames = snp_names,
        popNames = pop
    )

    cat("  Fst matrix computed for", nrow(fst_mat), "SNPs\\n")

    # Run OutFLANK
    out_result <- OutFLANK(
        FstDataFrame = fst_mat,
        LeftTrimFraction = 0.05,
        RightTrimFraction = 0.05,
        Hmin = 0.1,
        NumberOfSamples = length(unique(pop)),
        qthreshold = 0.05
    )

    # Extract outliers (q-value < 0.05)
    outlier_df <- out_result[["results"]]
    outlier_df[["qvalues"]] <- ifelse(is.na(outlier_df[["qvalues"]]), 1, outlier_df[["qvalues"]])
    outliers <- outlier_df[outlier_df[["qvalues"]] < 0.05, ]

    n_outliers <- nrow(outliers)
    cat("  Outlier SNPs (q < 0.05):", n_outliers, "\\n")

    # Create output directory
    out_dir <- paste0(node_name, "_outliers")
    dir.create(out_dir, showWarnings = FALSE)

    # Write outlier SNP list
    if (n_outliers > 0) {
        outlier_snps <- outliers[["LocusName"]]
        writeLines(outlier_snps, file.path(out_dir, "outlier_snps.txt"))

        # Write filtered reference genotypes (outlier SNPs only)
        ref_out <- ref_geno[, c("sample_id", outlier_snps), with = FALSE]
        fwrite(ref_out, file.path(out_dir, "outlier_genotypes.tsv"), sep = "\\t")
    } else {
        # No outliers - write empty files
        writeLines(character(0), file.path(out_dir, "outlier_snps.txt"))
        fwrite(ref_geno[, "sample_id"], file.path(out_dir, "outlier_genotypes.tsv"), sep = "\\t")
    }

    # Write full Fst results for diagnostics
    fwrite(as.data.table(outlier_df), file.path(out_dir, "outflank_results.tsv"), sep = "\\t")

    # Copy labels
    file.copy(file.path("${node_dir}", "reference_labels.tsv"),
              file.path(out_dir, "reference_labels.tsv"))

    # Write summary
    summary_info <- list(
        node_name = node_name,
        n_snps_tested = nrow(fst_mat),
        n_outliers = n_outliers,
        outlier_fraction = round(n_outliers / max(nrow(fst_mat), 1), 4)
    )
    write_json(summary_info, file.path(out_dir, "outflank_summary.json"),
               pretty = TRUE, auto_unbox = TRUE)
    """
}
