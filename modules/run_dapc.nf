/*
========================================================================================
    DAPC CLASSIFICATION
========================================================================================
    1. Use xvalDapc() to select optimal number of PCA axes
    2. Train DAPC with optimal axes
    3. Extract posterior probabilities and SNP loadings
*/

process RUN_DAPC {
    tag "${node_dir.name}"
    label 'process_medium'

    publishDir "${params.outdir}/classification/dapc", mode: 'copy'

    input:
    tuple path(node_dir), path(node_meta)

    output:
    tuple val("${node_dir.name}"), path("${node_dir.name}_dapc_results"), emit: results

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(adegenet)
    library(jsonlite)

    node_name <- basename("${node_dir}")
    out_dir <- paste0(node_name, "_dapc_results")
    dir.create(out_dir, showWarnings = FALSE)

    meta <- fromJSON("${node_meta}")
    cat("DAPC at node:", node_name, "\\n")
    cat("Groups:", paste(meta\$groups, collapse = ", "), "\\n")

    # Load reference data
    ref_geno <- fread(file.path("${node_dir}", "reference_genotypes.tsv"), header = TRUE)
    ref_labels <- fread(file.path("${node_dir}", "reference_labels.tsv"), header = TRUE)

    geno_mat <- as.matrix(ref_geno[, -1])
    rownames(geno_mat) <- ref_geno\$sample_id

    # Mean-impute
    for (j in 1:ncol(geno_mat)) {
        m <- mean(geno_mat[, j], na.rm = TRUE)
        geno_mat[is.na(geno_mat[, j]), j] <- round(m)
    }

    # Set up groups
    grp <- factor(ref_labels\$group[match(ref_geno\$sample_id, ref_labels\$sample_id)])
    names(grp) <- ref_geno\$sample_id

    n_groups <- length(unique(grp))
    cat("Number of groups:", n_groups, "\\n")

    # Remove columns that are all NA (would produce NaN after imputation)
    all_na <- apply(geno_mat, 2, function(x) all(is.na(x)))
    if (any(all_na)) {
        cat("Removing", sum(all_na), "all-NA columns\\n")
        geno_mat <- geno_mat[, !all_na, drop = FALSE]
    }

    if (n_groups < 2) {
        # Single group - DAPC not applicable, assign all with probability 1.0
        cat("Only 1 group at this node -- skipping DAPC, assigning all to:", levels(grp)[1], "\\n")
        single_grp <- as.character(levels(grp)[1])

        loocv_results <- data.table(
            sample_id = ref_geno\$sample_id,
            true_group = as.character(grp),
            predicted_group = as.character(grp),
            posterior_max = 1.0,
            correct = TRUE
        )
        accuracy <- 1.0

        # Handle unknowns
        unk_file <- file.path("${node_dir}", "unknown_genotypes.tsv")
        if (file.exists(unk_file)) {
            unk_geno <- fread(unk_file, header = TRUE)
            if (nrow(unk_geno) > 0) {
                unknown_assignments <- data.table(
                    sample_id = unk_geno\$sample_id,
                    predicted_group = single_grp,
                    posterior_max = 1.0
                )
                post_col <- paste0("posterior_", single_grp)
                unknown_assignments[[post_col]] <- 1.0
                fwrite(unknown_assignments, file.path(out_dir, "unknown_predictions.tsv"), sep = "\\t")
            }
        }

        snp_importance <- data.table(snp = character(), total_loading = numeric())
        fwrite(loocv_results, file.path(out_dir, "loocv_results.tsv"), sep = "\\t")
        fwrite(snp_importance, file.path(out_dir, "snp_importance.tsv"), sep = "\\t")

        summary_dt <- data.table(
            node = node_name, method = "DAPC",
            n_reference = nrow(geno_mat), n_snps = ncol(geno_mat),
            n_groups = 1, optimal_npca = NA_integer_, loocv_accuracy = 1.0
        )
        fwrite(summary_dt, file.path(out_dir, "summary.tsv"), sep = "\\t")
        cat("DAPC complete for", node_name, "(single group)\\n")
    } else {

    # Cross-validation to find optimal n.pca
    # Cap at n_samples/3 to prevent overfitting (binary posteriors)
    max_npca <- min(nrow(geno_mat) - 1, ncol(geno_mat) - 1, as.integer(nrow(geno_mat) / 3), 150)
    if (max_npca < 10) {
        n_pca_range <- seq_len(max(max_npca, 1))
    } else {
        n_pca_range <- unique(c(seq(10, max_npca, by = 10), max_npca))
    }
    n_pca_range <- n_pca_range[n_pca_range > 0 & n_pca_range <= max_npca]

    cat("Cross-validating DAPC with n.pca range:", range(n_pca_range), "\\n")

    set.seed(42)
    xval <- tryCatch({
        xvalDapc(geno_mat, grp, n.pca = n_pca_range, n.rep = 30,
                 training.set = 0.9, result = "groupMean")
    }, error = function(e) {
        cat("xvalDapc failed:", e\$message, "\\n")
        cat("Using n.pca = min(n-1, 50)\\n")
        NULL
    })

    if (!is.null(xval)) {
        msa <- xval[["Mean Successful Assignment"]]
        best_idx <- which.max(msa)
        if (length(best_idx) == 1) {
            optimal_npca <- as.integer(names(best_idx))
        } else {
            optimal_npca <- NA_integer_
        }
        # Fallback if extraction failed
        if (is.na(optimal_npca) || length(optimal_npca) != 1) {
            optimal_npca <- min(nrow(geno_mat) - 1, 50)
            cat("xvalDapc returned no clear optimum, falling back to n.pca =", optimal_npca, "\\n")
        } else {
            cat("Optimal n.pca:", optimal_npca, "\\n")
        }
    } else {
        optimal_npca <- min(nrow(geno_mat) - 1, 50)
    }

    # Train final DAPC
    n_da <- n_groups - 1
    dapc_result <- dapc(geno_mat, grp, n.pca = optimal_npca, n.da = n_da)

    # LOOCV for accuracy assessment
    cat("Running LOOCV...\\n")
    loocv_results <- data.table(
        sample_id = character(),
        true_group = character(),
        predicted_group = character(),
        posterior_max = numeric(),
        correct = logical()
    )

    for (i in 1:nrow(geno_mat)) {
        train_mat <- geno_mat[-i, , drop = FALSE]
        train_grp <- grp[-i]
        test_mat <- geno_mat[i, , drop = FALSE]

        dapc_loo <- tryCatch({
            dapc(train_mat, train_grp, n.pca = optimal_npca, n.da = n_da)
        }, error = function(e) NULL)

        if (!is.null(dapc_loo)) {
            pred <- predict(dapc_loo, newdata = test_mat)
            pred_group <- as.character(pred\$assign)
            post_max <- max(pred\$posterior)
        } else {
            pred_group <- NA
            post_max <- NA
        }

        loocv_results <- rbind(loocv_results, data.table(
            sample_id = rownames(geno_mat)[i],
            true_group = as.character(grp[i]),
            predicted_group = pred_group,
            posterior_max = post_max,
            correct = pred_group == as.character(grp[i])
        ))
    }

    accuracy <- mean(loocv_results\$correct, na.rm = TRUE)
    cat("LOOCV accuracy:", round(accuracy * 100, 1), "%\\n")

    # Predict unknowns
    unk_file <- file.path("${node_dir}", "unknown_genotypes.tsv")
    if (file.exists(unk_file)) {
        unk_geno <- fread(unk_file, header = TRUE)
        if (nrow(unk_geno) > 0) {
            unk_mat <- as.matrix(unk_geno[, -1])
            rownames(unk_mat) <- unk_geno\$sample_id
            for (j in 1:ncol(unk_mat)) {
                m <- mean(geno_mat[, j], na.rm = TRUE)
                unk_mat[is.na(unk_mat[, j]), j] <- round(m)
            }
            unk_pred <- predict(dapc_result, newdata = unk_mat)

            unknown_assignments <- data.table(
                sample_id = unk_geno\$sample_id,
                predicted_group = as.character(unk_pred\$assign),
                posterior_max = apply(unk_pred\$posterior, 1, max)
            )

            # Add all posterior probabilities
            post_dt <- as.data.table(unk_pred\$posterior)
            setnames(post_dt, paste0("posterior_", colnames(unk_pred\$posterior)))
            unknown_assignments <- cbind(unknown_assignments, post_dt)

            fwrite(unknown_assignments, file.path(out_dir, "unknown_predictions.tsv"), sep = "\\t")
        }
    }

    # SNP loadings
    loadings <- dapc_result\$var.contr
    top_snps <- head(order(-rowSums(abs(loadings))), 100)
    snp_importance <- data.table(
        snp = rownames(loadings)[top_snps],
        total_loading = rowSums(abs(loadings))[top_snps]
    )

    # Save outputs
    fwrite(loocv_results, file.path(out_dir, "loocv_results.tsv"), sep = "\\t")
    fwrite(snp_importance, file.path(out_dir, "snp_importance.tsv"), sep = "\\t")
    saveRDS(dapc_result, file.path(out_dir, "dapc_model.rds"))

    # Write summary -- force scalars to prevent 0-row data.table
    summary_dt <- data.table(
        node = as.character(node_name)[1],
        method = "DAPC",
        n_reference = as.integer(nrow(geno_mat))[1],
        n_snps = as.integer(ncol(geno_mat))[1],
        n_groups = as.integer(n_groups)[1],
        optimal_npca = as.integer(optimal_npca)[1],
        loocv_accuracy = round(as.numeric(accuracy)[1], 4)
    )
    stopifnot(nrow(summary_dt) == 1)
    fwrite(summary_dt, file.path(out_dir, "summary.tsv"), sep = "\\t")

    cat("DAPC complete for", node_name, "\\n")
    } # end else (n_groups >= 2)
    """
}
