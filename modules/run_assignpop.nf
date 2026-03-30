/*
========================================================================================
    assignPOP CLASSIFICATION
========================================================================================
    Run all three models (Bayesian, Random Forest, SVM) with LOO cross-validation.
    Select best-performing model per node.
    Extract LOD scores and variable importance.
*/

process RUN_ASSIGNPOP {
    tag "${node_dir.name}"
    label 'process_medium'

    publishDir "${params.outdir}/classification/assignpop", mode: 'copy'

    input:
    tuple path(node_dir), path(node_meta)

    output:
    tuple val("${node_dir.name}"), path("${node_dir.name}_assignpop_results"), emit: results

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(MASS)
    library(e1071)
    library(randomForest)
    library(jsonlite)

    node_name <- basename("${node_dir}")
    out_dir <- paste0(node_name, "_assignpop_results")
    dir.create(out_dir, showWarnings = FALSE)

    meta <- fromJSON("${node_meta}")
    cat("assignPOP at node:", node_name, "\\n")

    # Load reference data
    ref_geno <- fread(file.path("${node_dir}", "reference_genotypes.tsv"), header = TRUE)
    ref_labels <- fread(file.path("${node_dir}", "reference_labels.tsv"), header = TRUE)

    geno_mat <- as.matrix(ref_geno[, -1])
    rownames(geno_mat) <- ref_geno\$sample_id

    # For assignPOP: need complete data (no NAs)
    # Restrict to SNPs with no missing in reference
    complete_snps <- colSums(is.na(geno_mat)) == 0
    if (sum(complete_snps) < 100) {
        # If too few complete SNPs, mean-impute instead
        cat("Only", sum(complete_snps), "complete SNPs; using mean imputation\\n")
        for (j in 1:ncol(geno_mat)) {
            m <- mean(geno_mat[, j], na.rm = TRUE)
            geno_mat[is.na(geno_mat[, j]), j] <- round(m)
        }
    } else {
        geno_mat <- geno_mat[, complete_snps, drop = FALSE]
        cat("Using", ncol(geno_mat), "complete SNPs (no missing)\\n")
    }

    # Set up group labels
    labels <- ref_labels\$group[match(ref_geno\$sample_id, ref_labels\$sample_id)]

    # Remove all-NA columns
    all_na <- apply(geno_mat, 2, function(x) all(is.na(x)))
    if (any(all_na)) {
        cat("Removing", sum(all_na), "all-NA columns\\n")
        geno_mat <- geno_mat[, !all_na, drop = FALSE]
    }

    # In SNP mode, assignPOP receives Fst-filtered data (top 5000 per node).
    # In other modes, all SNPs are passed through.
    cat("Using", ncol(geno_mat), "SNPs for assignPOP\\n")

    n_groups <- length(unique(labels))
    cat("Number of groups:", n_groups, "\\n")

    if (n_groups < 2) {
        single_grp <- unique(labels)[1]
        cat("Only 1 group -- skipping classification, assigning all to:", single_grp, "\\n")

        loocv <- data.table(
            sample_id = ref_geno\$sample_id, true_group = labels,
            predicted_group = labels, probability = 1.0, correct = TRUE
        )
        for (mn in c("bayesian", "svm", "randomforest")) {
            fwrite(loocv, file.path(out_dir, paste0("loocv_", mn, ".tsv")), sep = "\\t")
        }

        unk_file <- file.path("${node_dir}", "unknown_genotypes.tsv")
        if (file.exists(unk_file)) {
            unk_geno <- fread(unk_file, header = TRUE)
            if (nrow(unk_geno) > 0) {
                unk_preds <- data.table(
                    sample_id = unk_geno\$sample_id,
                    predicted_group = single_grp,
                    probability = 1.0, lod = Inf
                )
                fwrite(unk_preds, file.path(out_dir, "unknown_predictions.tsv"), sep = "\\t")
                fwrite(unk_preds, file.path(out_dir, "unknown_predictions_svm.tsv"), sep = "\\t")
            }
        }

        summary_dt <- data.table(
            node = node_name, method = "assignPOP", best_model = "Bayesian",
            n_reference = nrow(geno_mat), n_snps = ncol(geno_mat), n_groups = 1,
            loocv_accuracy_lda = 1.0, loocv_accuracy_svm = 1.0,
            loocv_accuracy_rf = 1.0, best_accuracy = 1.0
        )
        fwrite(summary_dt, file.path(out_dir, "summary.tsv"), sep = "\\t")
        cat("assignPOP complete for", node_name, "(single group)\\n")
        quit(save = "no", status = 0)
    }

    # Create assignPOP-compatible data frame
    assign_df <- data.frame(pop = factor(labels), geno_mat, check.names = FALSE)

    # Run LOOCV with all three methods
    methods <- c("lda", "svm", "randomForest")
    method_names <- c("Bayesian", "SVM", "RandomForest")

    all_results <- list()
    all_accuracy <- numeric()

    for (m in seq_along(methods)) {
        cat("\\nRunning", method_names[m], "LOOCV...\\n")

        loocv <- data.table(
            sample_id = character(),
            true_group = character(),
            predicted_group = character(),
            probability = numeric(),
            correct = logical()
        )

        for (i in 1:nrow(geno_mat)) {
            train_df <- assign_df[-i, , drop = FALSE]
            test_df <- assign_df[i, -1, drop = FALSE]

            pred <- tryCatch({
                if (methods[m] == "lda") {
                    model <- MASS::lda(pop ~ ., data = train_df)
                    p <- predict(model, test_df)
                    list(class = as.character(p\$class), prob = max(p\$posterior))
                } else if (methods[m] == "svm") {
                    model <- e1071::svm(pop ~ ., data = train_df, probability = TRUE)
                    p <- predict(model, test_df, probability = TRUE)
                    probs <- attr(p, "probabilities")
                    list(class = as.character(p), prob = max(probs))
                } else {
                    # Use matrix interface to avoid formula overhead / stack overflow
                    # ntree=100 for LOOCV (sufficient for CV accuracy; final model uses 500)
                    model <- randomForest::randomForest(
                        x = as.matrix(train_df[, -1, drop = FALSE]),
                        y = train_df\$pop,
                        ntree = 100
                    )
                    p <- predict(model, as.matrix(test_df), type = "prob")
                    list(class = colnames(p)[which.max(p)], prob = max(p))
                }
            }, error = function(e) {
                list(class = NA, prob = NA)
            })

            loocv <- rbind(loocv, data.table(
                sample_id = rownames(geno_mat)[i],
                true_group = as.character(labels[i]),
                predicted_group = pred\$class,
                probability = pred\$prob,
                correct = !is.na(pred\$class) && pred\$class == as.character(labels[i])
            ))
        }

        acc <- mean(loocv\$correct, na.rm = TRUE)
        all_accuracy[m] <- acc
        all_results[[method_names[m]]] <- loocv
        cat(method_names[m], "accuracy:", round(acc * 100, 1), "%\\n")

        fwrite(loocv, file.path(out_dir, paste0("loocv_", tolower(method_names[m]), ".tsv")), sep = "\\t")
    }

    # Select best method
    best_idx <- which.max(all_accuracy)
    best_method <- method_names[best_idx]
    cat("\\nBest method:", best_method, "(", round(all_accuracy[best_idx] * 100, 1), "%)\\n")

    # Train final model with best method on all data
    if (methods[best_idx] == "lda") {
        final_model <- MASS::lda(pop ~ ., data = assign_df)
    } else if (methods[best_idx] == "svm") {
        final_model <- e1071::svm(pop ~ ., data = assign_df, probability = TRUE)
    } else {
        final_model <- randomForest::randomForest(
            x = as.matrix(assign_df[, -1, drop = FALSE]),
            y = assign_df\$pop,
            ntree = 500,
            importance = TRUE
        )
    }

    # Always save SVM model for Shiny app (ensemble always needs SVM)
    if (methods[best_idx] == "svm") {
        svm_final <- final_model
    } else {
        svm_final <- e1071::svm(pop ~ ., data = assign_df, probability = TRUE)
    }
    saveRDS(svm_final, file.path(out_dir, "svm_model.rds"))
    saveRDS(colnames(geno_mat), file.path(out_dir, "svm_snp_names.rds"))
    ref_means <- colMeans(geno_mat, na.rm = TRUE)
    saveRDS(ref_means, file.path(out_dir, "svm_ref_means.rds"))
    cat("Saved SVM model + metadata for Shiny app\\n")

    # Predict unknowns
    unk_file <- file.path("${node_dir}", "unknown_genotypes.tsv")
    if (file.exists(unk_file)) {
        unk_geno <- fread(unk_file, header = TRUE)
        if (nrow(unk_geno) > 0) {
            unk_mat <- as.matrix(unk_geno[, -1])
            # Match columns
            shared_snps <- intersect(colnames(geno_mat), colnames(unk_mat))
            unk_mat <- unk_mat[, shared_snps, drop = FALSE]
            # Impute missing
            for (j in 1:ncol(unk_mat)) {
                m <- mean(geno_mat[, shared_snps[j]], na.rm = TRUE)
                unk_mat[is.na(unk_mat[, j]), j] <- round(m)
            }

            unk_df <- data.frame(unk_mat, check.names = FALSE)

            unk_preds <- tryCatch({
                if (methods[best_idx] == "lda") {
                    p <- predict(final_model, unk_df)
                    data.table(
                        sample_id = unk_geno\$sample_id,
                        predicted_group = as.character(p\$class),
                        probability = apply(p\$posterior, 1, max)
                    )
                } else if (methods[best_idx] == "svm") {
                    p <- predict(final_model, unk_df, probability = TRUE)
                    probs <- attr(p, "probabilities")
                    data.table(
                        sample_id = unk_geno\$sample_id,
                        predicted_group = as.character(p),
                        probability = apply(probs, 1, max)
                    )
                } else {
                    p <- predict(final_model, as.matrix(unk_df), type = "prob")
                    data.table(
                        sample_id = unk_geno\$sample_id,
                        predicted_group = colnames(p)[apply(p, 1, which.max)],
                        probability = apply(p, 1, max)
                    )
                }
            }, error = function(e) {
                cat("Prediction failed:", e\$message, "\\n")
                NULL
            })

            if (!is.null(unk_preds)) {
                # Calculate LOD score: -log10((1 - Pi) / (1 - Pj_k))
                # where Pi = max prob, Pj_k = sum of remaining probs
                unk_preds[, lod := -log10((1 - probability) / probability)]
                fwrite(unk_preds, file.path(out_dir, "unknown_predictions.tsv"), sep = "\\t")
            }
        }
    }

    # Always produce SVM predictions for unknowns (ensemble always uses SVM)
    svm_pred_file <- file.path(out_dir, "unknown_predictions_svm.tsv")
    if (methods[best_idx] == "svm") {
        # Best model was SVM -- just copy the existing predictions
        if (file.exists(file.path(out_dir, "unknown_predictions.tsv"))) {
            file.copy(file.path(out_dir, "unknown_predictions.tsv"), svm_pred_file)
        }
    } else if (file.exists(unk_file)) {
        # Best model was NOT SVM -- train an SVM on all reference data and predict
        unk_geno <- fread(unk_file, header = TRUE)
        if (nrow(unk_geno) > 0) {
            svm_model <- e1071::svm(pop ~ ., data = assign_df, probability = TRUE)

            unk_mat_svm <- as.matrix(unk_geno[, -1])
            shared_snps_svm <- intersect(colnames(geno_mat), colnames(unk_mat_svm))
            unk_mat_svm <- unk_mat_svm[, shared_snps_svm, drop = FALSE]
            for (j in 1:ncol(unk_mat_svm)) {
                m_val <- mean(geno_mat[, shared_snps_svm[j]], na.rm = TRUE)
                unk_mat_svm[is.na(unk_mat_svm[, j]), j] <- round(m_val)
            }
            unk_df_svm <- data.frame(unk_mat_svm, check.names = FALSE)

            svm_preds <- tryCatch({
                p <- predict(svm_model, unk_df_svm, probability = TRUE)
                probs <- attr(p, "probabilities")
                data.table(
                    sample_id = unk_geno\$sample_id,
                    predicted_group = as.character(p),
                    probability = apply(probs, 1, max)
                )
            }, error = function(e) {
                cat("SVM prediction failed:", e\$message, "\\n")
                NULL
            })

            if (!is.null(svm_preds)) {
                svm_preds[, lod := -log10((1 - probability) / probability)]
                fwrite(svm_preds, svm_pred_file, sep = "\\t")
            }
        }
    }

    # Extract variable importance (RF only)
    if (methods[best_idx] == "randomForest") {
        imp <- randomForest::importance(final_model)
        imp_dt <- data.table(snp = rownames(imp), importance = imp[, "MeanDecreaseAccuracy"])
        imp_dt <- imp_dt[order(-importance)]
        fwrite(head(imp_dt, 100), file.path(out_dir, "snp_importance.tsv"), sep = "\\t")
    }

    # Write summary
    summary_dt <- data.table(
        node = node_name,
        method = "assignPOP",
        best_model = best_method,
        n_reference = nrow(geno_mat),
        n_snps = ncol(geno_mat),
        n_groups = length(unique(labels)),
        loocv_accuracy_lda = round(all_accuracy[1], 4),
        loocv_accuracy_svm = round(all_accuracy[2], 4),
        loocv_accuracy_rf = round(all_accuracy[3], 4),
        best_accuracy = round(all_accuracy[best_idx], 4)
    )
    fwrite(summary_dt, file.path(out_dir, "summary.tsv"), sep = "\\t")

    cat("assignPOP complete for", node_name, "\\n")
    """
}
