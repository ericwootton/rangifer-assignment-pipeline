/*
========================================================================================
    HERD DIFFERENTIATION SUMMARY
========================================================================================
    Compiles all herd differentiation results into a single summary report.
*/

process HERD_DIFF_SUMMARY {
    tag "summary"

    publishDir "${params.outdir}/herd_differentiation", mode: 'copy'

    input:
    path fst_matrix
    path fst_long
    path pca_scores
    path pca_eigenvalues
    path admixture_results
    path admixture_choosek
    path dapc_results
    path assignpop_results
    path node_admixture_results
    path metadata

    output:
    path "herd_differentiation_summary.txt", emit: summary

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    cat_section <- function(title) {
        cat("\\n", paste(rep("=", 70), collapse = ""), "\\n", sep = "")
        cat(" ", title, "\\n", sep = "")
        cat(paste(rep("=", 70), collapse = ""), "\\n\\n", sep = "")
    }

    sink("herd_differentiation_summary.txt")

    cat("HERD DIFFERENTIATION ANALYSIS SUMMARY\\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")

    # ---- Sample counts ----
    cat_section("SAMPLE COUNTS BY HERD")
    meta <- fread("${metadata}", header = TRUE)
    if ("Sample" %in% colnames(meta)) setnames(meta, "Sample", "sample_id")
    if ("Herd_Cleaned" %in% colnames(meta)) setnames(meta, "Herd_Cleaned", "herd")
    id_col <- intersect(c("sample_id", "ID", "Sample_ID"), colnames(meta))[1]
    if (!is.na(id_col) && id_col != "sample_id") setnames(meta, id_col, "sample_id")

    pca <- fread("${pca_scores}", header = TRUE)
    meta_in <- meta[sample_id %in% pca[["sample_id"]]]
    herd_tab <- sort(table(meta_in[["herd"]]), decreasing = TRUE)
    for (h in names(herd_tab)) cat(sprintf("  %-25s %d\\n", h, herd_tab[h]))
    cat(sprintf("  %-25s %d\\n", "TOTAL", sum(herd_tab)))

    # ---- PCA ----
    cat_section("PCA VARIANCE EXPLAINED")
    eig <- fread("${pca_eigenvalues}", header = TRUE)
    for (i in 1:min(10, nrow(eig))) {
        cat(sprintf("  PC%-3d %6.2f%%  (cumulative: %6.2f%%)\\n",
                    i, eig[["var_percent"]][i], eig[["cumulative_var"]][i]))
    }

    # ---- Pairwise Fst ----
    cat_section("PAIRWISE FST")
    fst <- fread("${fst_long}", header = TRUE)
    fst <- fst[order(-fst)]
    cat("Top 10 most differentiated pairs:\\n")
    for (i in 1:min(10, nrow(fst))) {
        cat(sprintf("  %-20s vs %-20s  Fst = %.4f\\n",
                    fst[["pop1"]][i], fst[["pop2"]][i], fst[["fst"]][i]))
    }
    cat("\\nLeast differentiated pairs:\\n")
    fst_asc <- fst[order(fst)]
    for (i in 1:min(5, nrow(fst_asc))) {
        cat(sprintf("  %-20s vs %-20s  Fst = %.4f\\n",
                    fst_asc[["pop1"]][i], fst_asc[["pop2"]][i], fst_asc[["fst"]][i]))
    }
    cat(sprintf("\\nOverall Fst range: %.4f - %.4f\\n", min(fst[["fst"]]), max(fst[["fst"]])))
    cat(sprintf("Mean pairwise Fst: %.4f\\n", mean(fst[["fst"]])))

    # ---- ADMIXTURE ----
    cat_section("ADMIXTURE (GENOME-WIDE)")
    choosek <- readLines("${admixture_choosek}")
    for (line in choosek) cat(" ", line, "\\n")

    # ---- Node ADMIXTURE ----
    cat_section("ADMIXTURE (HERD NODE)")
    adm_dir <- "${node_admixture_results}"
    adm_choosek <- file.path(adm_dir, "choosek.txt")
    if (file.exists(adm_choosek)) {
        adm_lines <- readLines(adm_choosek)
        for (line in adm_lines) cat(" ", line, "\\n")
    }

    # ---- DAPC ----
    cat_section("DAPC CLASSIFICATION (LOOCV)")
    dapc_dir <- "${dapc_results}"
    loocv_file <- file.path(dapc_dir, "loocv_results.tsv")
    if (file.exists(loocv_file)) {
        loocv <- fread(loocv_file, header = TRUE)
        overall_acc <- mean(loocv[["correct"]])
        cat(sprintf("Overall LOOCV accuracy: %.1f%% (%d/%d)\\n\\n",
                    overall_acc * 100, sum(loocv[["correct"]]), nrow(loocv)))

        cat("Per-herd accuracy:\\n")
        for (grp in sort(unique(loocv[["true_group"]]))) {
            sub <- loocv[true_group == grp]
            acc <- mean(sub[["correct"]])
            cat(sprintf("  %-25s %.1f%% (%d/%d)\\n",
                        grp, acc * 100, sum(sub[["correct"]]), nrow(sub)))
        }

        # Confusion summary
        cat("\\nMisclassifications:\\n")
        wrong <- loocv[correct == FALSE]
        if (nrow(wrong) > 0) {
            conf <- wrong[, .N, by = c("true_group", "predicted_group")]
            conf <- conf[order(-N)]
            for (i in 1:nrow(conf)) {
                cat(sprintf("  %-20s -> %-20s (%d samples)\\n",
                            conf[["true_group"]][i],
                            conf[["predicted_group"]][i],
                            conf[["N"]][i]))
            }
        } else {
            cat("  None\\n")
        }
    }

    # ---- assignPOP ----
    cat_section("ASSIGNPOP CLASSIFICATION (LOOCV)")
    ap_dir <- "${assignpop_results}"
    summary_file <- file.path(ap_dir, "summary.tsv")
    if (file.exists(summary_file)) {
        ap_sum <- fread(summary_file, header = TRUE)
        cat("Method performance:\\n")
        for (i in 1:nrow(ap_sum)) {
            for (col in colnames(ap_sum)) {
                cat(sprintf("  %-20s %s\\n", paste0(col, ":"), ap_sum[[col]][i]))
            }
            cat("\\n")
        }
    }

    # Per-method LOOCV details
    for (method in c("svm", "bayesian", "randomforest")) {
        lf <- file.path(ap_dir, paste0("loocv_", method, ".tsv"))
        if (file.exists(lf)) {
            loocv <- fread(lf, header = TRUE)
            # Find correct column (varies)
            correct_col <- intersect(c("correct", "Correct"), colnames(loocv))[1]
            pred_col <- intersect(c("predicted_group", "predicted", "Predicted"), colnames(loocv))[1]
            true_col <- intersect(c("true_group", "true", "True"), colnames(loocv))[1]
            if (!is.na(correct_col)) {
                acc <- mean(loocv[[correct_col]])
                cat(sprintf("  %s LOOCV: %.1f%%\\n", toupper(method), acc * 100))
            }
        }
    }

    cat("\\n")
    sink()

    cat("Summary written to herd_differentiation_summary.txt\\n")
    """
}
