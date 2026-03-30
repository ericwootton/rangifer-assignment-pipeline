/*
========================================================================================
    GENERATE NODE DIAGNOSTIC PLOTS
========================================================================================
    Generates comprehensive per-node diagnostic plots as a multi-page PDF.
    Combines outputs from DAPC, assignPOP, popfinder, ADMIXTURE, and ensemble
    to produce 12 diagnostic pages per multi-group node.

    Uses base R graphics only (no ggplot2).
*/

process GENERATE_NODE_PLOTS {
    tag "${node_name}"
    label 'process_medium'
    publishDir "${params.outdir}/classification/node_plots", mode: 'copy'

    input:
    tuple val(node_name), path(node_dir), path(dapc_dir), path(assignpop_dir),
          path(popfinder_dir), path(admixture_dir), path(ensemble_file)
    path data_source

    output:
    path "${node_name}_plots.pdf", emit: plots

    script:
    def has_popfinder = popfinder_dir && popfinder_dir.name != 'NO_POPFINDER' && popfinder_dir.name != 'null'
    def has_admixture = admixture_dir && admixture_dir.name != 'NO_ADMIXTURE' && admixture_dir.name != 'null'
    """
    #!/usr/bin/env Rscript
    # v2: data source shapes (chip vs WGS)

    library(data.table)
    library(jsonlite)
    library(adegenet)

    node_name <- "${node_name}"
    has_popfinder <- ${has_popfinder ? 'TRUE' : 'FALSE'}
    has_admixture <- ${has_admixture ? 'TRUE' : 'FALSE'}

    cat("Generating diagnostic plots for node:", node_name, "\\n")

    # -------------------------------------------------------------------------
    # Helper: safe file read
    # -------------------------------------------------------------------------
    safe_fread <- function(path, ...) {
        tryCatch({
            if (file.exists(path)) fread(path, header = TRUE, ...) else NULL
        }, error = function(e) { cat("Warning: cannot read", path, ":", e[["message"]], "\\n"); NULL })
    }

    safe_readLines <- function(path) {
        tryCatch({
            if (file.exists(path)) readLines(path) else character(0)
        }, error = function(e) character(0))
    }

    # -------------------------------------------------------------------------
    # Helper: color palette for groups
    # -------------------------------------------------------------------------
    make_group_colors <- function(groups) {
        n <- length(groups)
        if (n <= 8) {
            pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                     "#FF7F00", "#A65628", "#F781BF", "#999999")
        } else {
            pal <- rainbow(n, s = 0.7, v = 0.85)
        }
        cols <- pal[seq_len(n)]
        names(cols) <- groups
        cols
    }

    # -------------------------------------------------------------------------
    # Helper: draw confusion matrix heatmap
    # -------------------------------------------------------------------------
    draw_confusion <- function(loocv_dt, method_label) {
        if (is.null(loocv_dt) || nrow(loocv_dt) == 0) {
            plot.new(); text(0.5, 0.5, paste(method_label, "- no data"), cex = 1.2)
            return(invisible(NULL))
        }
        true_col <- intersect(c("true_group", "true_pop"), colnames(loocv_dt))[1]
        pred_col <- intersect(c("predicted_group", "predicted_pop"), colnames(loocv_dt))[1]
        if (is.na(true_col) || is.na(pred_col)) {
            plot.new(); text(0.5, 0.5, paste(method_label, "- missing columns"), cex = 1.2)
            return(invisible(NULL))
        }
        true_vals <- as.character(loocv_dt[[true_col]])
        pred_vals <- as.character(loocv_dt[[pred_col]])
        all_groups <- sort(unique(c(true_vals, pred_vals)))
        cm <- table(factor(true_vals, levels = all_groups),
                    factor(pred_vals, levels = all_groups))
        n <- length(all_groups)
        acc <- sum(diag(cm)) / sum(cm)

        # Color ramp: white to steelblue
        max_val <- max(cm)
        if (max_val == 0) max_val <- 1
        col_ramp <- colorRampPalette(c("white", "steelblue"))(100)

        par(mar = c(5, 6, 4, 2))
        image(1:n, 1:n, t(as.matrix(cm))[, n:1, drop = FALSE],
              col = col_ramp, axes = FALSE,
              xlab = "Predicted", ylab = "",
              main = paste0(method_label, " (Acc: ", round(acc * 100, 1), "%)"))
        axis(1, at = 1:n, labels = all_groups, las = 2, cex.axis = 0.8)
        axis(2, at = 1:n, labels = rev(all_groups), las = 1, cex.axis = 0.8)
        mtext("True", side = 2, line = 4.5)

        # Add count text
        for (r in 1:n) {
            for (cc in 1:n) {
                val <- cm[n - cc + 1, r]
                text(r, cc, val, cex = max(0.6, min(1.2, 8 / n)),
                     col = if (val > max_val * 0.6) "white" else "black")
            }
        }
        box()
    }

    # -------------------------------------------------------------------------
    # Load node metadata
    # -------------------------------------------------------------------------
    meta <- tryCatch(
        fromJSON(file.path("${node_dir}", "node_meta.json")),
        error = function(e) list(groups = character(0), n_groups = 0, level = NA, n_snps = NA)
    )

    groups <- if (!is.null(meta[["groups"]])) meta[["groups"]] else character(0)
    n_groups <- if (!is.null(meta[["n_groups"]])) meta[["n_groups"]] else length(groups)
    node_level <- if (!is.null(meta[["level"]])) meta[["level"]] else "unknown"
    n_snps_meta <- if (!is.null(meta[["n_snps"]])) meta[["n_snps"]] else NA

    # Load reference labels for sample counts
    ref_labels <- safe_fread(file.path("${node_dir}", "reference_labels.tsv"))
    group_counts <- if (!is.null(ref_labels)) table(ref_labels[["group"]]) else NULL
    group_colors <- if (length(groups) > 0) make_group_colors(groups) else character(0)

    # Load data source mapping (chip vs WGS)
    source_dt <- safe_fread("${data_source}")
    if (!is.null(source_dt)) {
        # Collapse FinalReport/NWT into "Chip"
        source_map <- source_dt[["source"]]
        source_map[source_map != "WGS"] <- "Chip"
        names(source_map) <- source_dt[["sample_id"]]
    } else {
        source_map <- character(0)
    }
    # Shape mapping: Chip = circle (16), WGS = triangle (17)
    source_pch <- c("Chip" = 16, "WGS" = 17)

    # =========================================================================
    # Single-group node: minimal PDF
    # =========================================================================
    if (n_groups <= 1) {
        cat("Single-group node -- generating minimal PDF\\n")
        pdf("${node_name}_plots.pdf", width = 12, height = 10)
        plot.new()
        text(0.5, 0.8, paste("Node:", node_name), cex = 2.2, font = 2)
        text(0.5, 0.65, paste("Level:", node_level), cex = 1.4)
        text(0.5, 0.50, paste("Groups:", paste(groups, collapse = ", ")), cex = 1.4)
        n_total <- if (!is.null(ref_labels)) nrow(ref_labels) else 0
        text(0.5, 0.38, paste("Reference samples:", n_total), cex = 1.4)
        text(0.5, 0.22, "Classification is trivial (single group).", cex = 1.6, col = "grey40", font = 3)
        dev.off()
        cat("Done (single-group)\\n")
        quit(save = "no", status = 0)
    }

    # =========================================================================
    # Multi-group node: full 12-page PDF
    # =========================================================================

    # Pre-load all data sources
    dapc_loocv   <- safe_fread(file.path("${dapc_dir}", "loocv_results.tsv"))
    dapc_unk     <- safe_fread(file.path("${dapc_dir}", "unknown_predictions.tsv"))
    dapc_snps    <- safe_fread(file.path("${dapc_dir}", "snp_importance.tsv"))
    dapc_summary <- safe_fread(file.path("${dapc_dir}", "summary.tsv"))

    ap_loocv_svm <- safe_fread(file.path("${assignpop_dir}", "loocv_svm.tsv"))
    ap_loocv_bay <- safe_fread(file.path("${assignpop_dir}", "loocv_bayesian.tsv"))
    ap_loocv_rf  <- safe_fread(file.path("${assignpop_dir}", "loocv_randomforest.tsv"))
    ap_summary   <- safe_fread(file.path("${assignpop_dir}", "summary.tsv"))

    pf_loocv   <- if (has_popfinder) safe_fread(file.path("${popfinder_dir}", "loocv_results.tsv")) else NULL
    pf_summary <- if (has_popfinder) safe_fread(file.path("${popfinder_dir}", "summary.tsv")) else NULL

    adm_cv     <- if (has_admixture) safe_readLines(file.path("${admixture_dir}", "cv_errors.txt")) else character(0)
    adm_Q      <- if (has_admixture) safe_fread(file.path("${admixture_dir}", "optimal_Q_labeled.tsv"), header = FALSE) else NULL
    adm_choose <- if (has_admixture) safe_readLines(file.path("${admixture_dir}", "choosek.txt")) else character(0)

    ensemble <- safe_fread("${ensemble_file}")

    cat("Data loaded. Opening PDF...\\n")

    pdf("${node_name}_plots.pdf", width = 12, height = 10)

    # =========================================================================
    # PAGE 1: Title / Summary
    # =========================================================================
    tryCatch({
        plot.new()
        text(0.5, 0.88, paste("Node:", node_name), cex = 2.5, font = 2)
        text(0.5, 0.78, paste("Level:", node_level), cex = 1.5)
        text(0.5, 0.70, paste("Groups:", paste(groups, collapse = ", ")), cex = 1.3)
        if (!is.na(n_snps_meta)) text(0.5, 0.63, paste("Total SNPs:", n_snps_meta), cex = 1.3)

        if (!is.null(group_counts) && length(group_counts) > 0) {
            for (i in seq_along(group_counts)) {
                y_pos <- 0.55 - (i - 1) * 0.05
                if (y_pos > 0.1) {
                    text(0.5, y_pos,
                         paste0("  ", names(group_counts)[i], ": n = ", group_counts[i]),
                         cex = 1.1)
                }
            }
        }

        # Small bar chart of reference sample sizes in bottom portion
        if (!is.null(group_counts) && length(group_counts) > 0) {
            par(fig = c(0.15, 0.85, 0.05, 0.38), new = TRUE)
            par(mar = c(5, 4, 2, 1))
            bp <- barplot(as.numeric(group_counts), names.arg = names(group_counts),
                          col = group_colors[names(group_counts)],
                          main = "Reference Samples per Group", ylab = "Count",
                          las = 2, cex.names = 0.8, border = NA)
            text(bp, as.numeric(group_counts), labels = as.numeric(group_counts),
                 pos = 3, cex = 0.9)
            par(fig = c(0, 1, 0, 1), new = FALSE)
        }
    }, error = function(e) {
        cat("Page 1 (Title) error:", e[["message"]], "\\n")
        plot.new(); text(0.5, 0.5, paste("Title page error:", e[["message"]]), cex = 0.9)
    })

    # =========================================================================
    # PAGE 2: LOOCV Confusion Heatmaps
    # =========================================================================
    tryCatch({
        n_panels <- 2 + as.integer(!is.null(pf_loocv) && nrow(pf_loocv) > 0)
        layout_mat <- if (n_panels == 3) matrix(1:3, nrow = 1) else matrix(1:2, nrow = 1)
        layout(layout_mat)

        draw_confusion(dapc_loocv, "DAPC")
        draw_confusion(ap_loocv_svm, "SVM (assignPOP)")
        if (!is.null(pf_loocv) && nrow(pf_loocv) > 0) {
            draw_confusion(pf_loocv, "popfinder")
        }

        layout(1)
    }, error = function(e) {
        cat("Page 2 (Confusion) error:", e[["message"]], "\\n")
        plot.new(); text(0.5, 0.5, paste("Confusion heatmap error:", e[["message"]]), cex = 0.9)
        layout(1)
    })

    # =========================================================================
    # PAGE 3: LOOCV Accuracy Comparison
    # =========================================================================
    tryCatch({
        acc_names <- character(0)
        acc_vals  <- numeric(0)

        # DAPC accuracy
        if (!is.null(dapc_summary) && "loocv_accuracy" %in% colnames(dapc_summary)) {
            acc_names <- c(acc_names, "DAPC")
            acc_vals  <- c(acc_vals, as.numeric(dapc_summary[["loocv_accuracy"]][1]))
        }

        # assignPOP accuracies
        if (!is.null(ap_summary)) {
            if ("loocv_accuracy_lda" %in% colnames(ap_summary)) {
                acc_names <- c(acc_names, "LDA/Bayesian")
                acc_vals  <- c(acc_vals, as.numeric(ap_summary[["loocv_accuracy_lda"]][1]))
            }
            if ("loocv_accuracy_svm" %in% colnames(ap_summary)) {
                acc_names <- c(acc_names, "SVM")
                acc_vals  <- c(acc_vals, as.numeric(ap_summary[["loocv_accuracy_svm"]][1]))
            }
            if ("loocv_accuracy_rf" %in% colnames(ap_summary)) {
                acc_names <- c(acc_names, "Random Forest")
                acc_vals  <- c(acc_vals, as.numeric(ap_summary[["loocv_accuracy_rf"]][1]))
            }
        }

        # popfinder accuracy
        if (!is.null(pf_summary) && "loocv_accuracy" %in% colnames(pf_summary)) {
            acc_names <- c(acc_names, "popfinder")
            acc_vals  <- c(acc_vals, as.numeric(pf_summary[["loocv_accuracy"]][1]))
        }

        if (length(acc_vals) > 0) {
            bar_cols <- c("DAPC" = "#377EB8", "LDA/Bayesian" = "#984EA3",
                          "SVM" = "#E41A1C", "Random Forest" = "#4DAF4A",
                          "popfinder" = "#FF7F00")
            cols <- ifelse(acc_names %in% names(bar_cols), bar_cols[acc_names], "grey50")

            par(mar = c(7, 5, 4, 2))
            bp <- barplot(acc_vals * 100, names.arg = acc_names, col = cols,
                          main = paste0("LOOCV Accuracy Comparison - ", node_name),
                          ylab = "Accuracy (%)", ylim = c(0, 105),
                          las = 2, cex.names = 0.9, border = NA)
            abline(h = 95, lty = 2, col = "red", lwd = 1.5)
            text(max(bp) + 0.5, 96, "95%", col = "red", cex = 0.8, pos = 4)
            text(bp, acc_vals * 100, labels = paste0(round(acc_vals * 100, 1), "%"),
                 pos = 3, cex = 0.85)
        } else {
            plot.new()
            text(0.5, 0.5, "No LOOCV accuracy data available", cex = 1.3)
        }
    }, error = function(e) {
        cat("Page 3 (Accuracy) error:", e[["message"]], "\\n")
        plot.new(); text(0.5, 0.5, paste("Accuracy comparison error:", e[["message"]]), cex = 0.9)
    })

    # =========================================================================
    # PAGE 4: DAPC Scatterplot
    # =========================================================================
    tryCatch({
        dapc_model_file <- file.path("${dapc_dir}", "dapc_model.rds")
        if (file.exists(dapc_model_file)) {
            dapc_obj <- readRDS(dapc_model_file)

            # Merge in true groups
            ref_grp <- ref_labels[["group"]]
            names(ref_grp) <- ref_labels[["sample_id"]]
            sample_grps <- ref_grp[rownames(dapc_obj[["ind.coord"]])]

            n_da <- ncol(dapc_obj[["ind.coord"]])

            if (n_da >= 2) {
                # DF1 vs DF2 scatterplot
                par(mar = c(5, 5, 4, 10), xpd = NA)
                df1 <- dapc_obj[["ind.coord"]][, 1]
                df2 <- dapc_obj[["ind.coord"]][, 2]
                grp_factor <- factor(sample_grps, levels = groups)
                pt_cols <- group_colors[as.character(grp_factor)]

                # Data source shapes
                ref_ids <- rownames(dapc_obj[["ind.coord"]])
                ref_src <- ifelse(ref_ids %in% names(source_map),
                                  source_map[ref_ids], "Chip")
                pt_pch <- source_pch[ref_src]

                plot(df1, df2, col = pt_cols, pch = pt_pch, cex = 1.2,
                     xlab = "Discriminant Function 1",
                     ylab = "Discriminant Function 2",
                     main = paste0("DAPC - ", node_name))

                # Project unknowns if they exist
                unk_file <- file.path("${node_dir}", "unknown_genotypes.tsv")
                if (file.exists(unk_file)) {
                    unk_geno <- safe_fread(unk_file)
                    if (!is.null(unk_geno) && nrow(unk_geno) > 0) {
                        unk_mat <- as.matrix(unk_geno[, -1, with = FALSE])
                        rownames(unk_mat) <- unk_geno[["sample_id"]]
                        # Mean-impute using reference means
                        ref_geno <- safe_fread(file.path("${node_dir}", "reference_genotypes.tsv"))
                        if (!is.null(ref_geno)) {
                            ref_mat <- as.matrix(ref_geno[, -1, with = FALSE])
                            for (j in seq_len(ncol(unk_mat))) {
                                na_idx <- is.na(unk_mat[, j])
                                if (any(na_idx)) {
                                    col_mean <- mean(ref_mat[, j], na.rm = TRUE)
                                    unk_mat[na_idx, j] <- round(col_mean)
                                }
                            }
                        }
                        unk_pred <- tryCatch(predict(dapc_obj, newdata = unk_mat), error = function(e) NULL)
                        if (!is.null(unk_pred) && ncol(unk_pred[["ind.scores"]]) >= 2) {
                            # Unknown shapes by source
                            unk_ids <- rownames(unk_pred[["ind.scores"]])
                            unk_src <- ifelse(unk_ids %in% names(source_map),
                                              source_map[unk_ids], "Chip")
                            unk_pch <- source_pch[unk_src]
                            points(unk_pred[["ind.scores"]][, 1],
                                   unk_pred[["ind.scores"]][, 2],
                                   pch = unk_pch, cex = 1.5, lwd = 2, col = "black")
                        }
                    }
                }

                # Legend: groups (color) + data source (shape)
                legend_labels <- groups
                legend_cols <- group_colors[groups]
                legend_pch <- rep(16, length(groups))
                # Add unknown if plotted
                if (exists("unk_pred") && !is.null(unk_pred)) {
                    legend_labels <- c(legend_labels, "Unknown")
                    legend_cols <- c(legend_cols, "black")
                    legend_pch <- c(legend_pch, 16)
                }
                # Add source shape entries
                legend_labels <- c(legend_labels, "", "Chip", "WGS")
                legend_cols <- c(legend_cols, NA, "grey30", "grey30")
                legend_pch <- c(legend_pch, NA, 16, 17)
                legend("topright", inset = c(-0.15, 0), legend = legend_labels,
                       col = legend_cols, pch = legend_pch, cex = 0.8, bty = "n")

            } else {
                # Only 1 DF (2 groups): density / strip plot
                par(mar = c(5, 5, 4, 2))
                df1 <- dapc_obj[["ind.coord"]][, 1]
                grp_factor <- factor(sample_grps, levels = groups)

                # Density for each group
                plot(NULL, xlim = range(df1, na.rm = TRUE) * 1.1,
                     ylim = c(0, 1),
                     xlab = "Discriminant Function 1", ylab = "Density (scaled)",
                     main = paste0("DAPC - ", node_name, " (2 groups)"))

                for (g in groups) {
                    vals <- df1[grp_factor == g]
                    if (length(vals) > 2) {
                        d <- density(vals, na.rm = TRUE)
                        d_scaled <- d[["y"]] / max(d[["y"]])
                        lines(d[["x"]], d_scaled, col = group_colors[g], lwd = 2)
                    }
                }

                # Rug / jittered points at bottom (shapes by data source)
                ref_ids <- rownames(dapc_obj[["ind.coord"]])
                ref_src <- ifelse(ref_ids %in% names(source_map),
                                  source_map[ref_ids], "Chip")
                for (g in groups) {
                    idx <- which(grp_factor == g)
                    vals <- df1[idx]
                    pch_vals <- source_pch[ref_src[idx]]
                    points(vals, rep(0, length(vals)) + runif(length(vals), -0.03, 0.03),
                           col = group_colors[g], pch = pch_vals, cex = 0.8)
                }

                legend("topright",
                       legend = c(groups, "", "Chip", "WGS"),
                       col = c(group_colors[groups], NA, "grey30", "grey30"),
                       lwd = c(rep(2, length(groups)), NA, NA, NA),
                       pch = c(rep(16, length(groups)), NA, 16, 17),
                       cex = 0.9, bty = "n")
            }
            par(xpd = FALSE)
        } else {
            plot.new()
            text(0.5, 0.5, "DAPC model file not found", cex = 1.3)
        }
    }, error = function(e) {
        cat("Page 4 (DAPC scatter) error:", e[["message"]], "\\n")
        plot.new(); text(0.5, 0.5, paste("DAPC scatterplot error:", e[["message"]]), cex = 0.9)
    })

    # =========================================================================
    # PAGE 5: DAPC Loading Plot (top 30 SNPs)
    # =========================================================================
    tryCatch({
        if (!is.null(dapc_snps) && nrow(dapc_snps) > 0) {
            n_show <- min(30, nrow(dapc_snps))
            top <- dapc_snps[order(-dapc_snps[["total_loading"]])][seq_len(n_show)]
            # Reverse for horizontal barplot (top at top)
            top <- top[n_show:1]

            par(mar = c(5, 12, 4, 2))
            barplot(top[["total_loading"]], names.arg = top[["snp"]],
                    horiz = TRUE, las = 1, col = "steelblue", border = NA,
                    main = paste0("Top ", n_show, " SNPs by DAPC Loading - ", node_name),
                    xlab = "Total Loading",
                    cex.names = max(0.5, min(0.8, 20 / n_show)))
        } else {
            plot.new()
            text(0.5, 0.5, "No DAPC SNP importance data available", cex = 1.3)
        }
    }, error = function(e) {
        cat("Page 5 (DAPC loadings) error:", e[["message"]], "\\n")
        plot.new(); text(0.5, 0.5, paste("DAPC loading error:", e[["message"]]), cex = 0.9)
    })

    # =========================================================================
    # PAGE 6: Posterior / Probability Distributions
    # =========================================================================
    tryCatch({
        # Collect available LOOCV datasets
        classifier_data <- list()
        if (!is.null(dapc_loocv) && nrow(dapc_loocv) > 0) {
            prob_col <- if ("posterior_max" %in% colnames(dapc_loocv)) "posterior_max" else "probability"
            classifier_data[["DAPC"]] <- dapc_loocv
            attr(classifier_data[["DAPC"]], "prob_col") <- prob_col
        }
        if (!is.null(ap_loocv_svm) && nrow(ap_loocv_svm) > 0) {
            prob_col <- if ("probability" %in% colnames(ap_loocv_svm)) "probability" else "posterior_max"
            classifier_data[["SVM"]] <- ap_loocv_svm
            attr(classifier_data[["SVM"]], "prob_col") <- prob_col
        }
        if (!is.null(pf_loocv) && nrow(pf_loocv) > 0) {
            prob_col <- if ("posterior_max" %in% colnames(pf_loocv)) "posterior_max" else "probability"
            classifier_data[["popfinder"]] <- pf_loocv
            attr(classifier_data[["popfinder"]], "prob_col") <- prob_col
        }

        n_panels <- length(classifier_data)
        if (n_panels > 0) {
            par(mfrow = c(1, n_panels), mar = c(5, 5, 4, 2))

            for (clf_name in names(classifier_data)) {
                dt <- classifier_data[[clf_name]]
                pc <- attr(dt, "prob_col")
                true_col <- intersect(c("true_group", "true_pop"), colnames(dt))[1]
                correct_col <- "correct"

                if (is.na(true_col) || !pc %in% colnames(dt)) {
                    plot.new()
                    text(0.5, 0.5, paste(clf_name, "- missing cols"))
                    next
                }

                grp_vals <- as.character(dt[[true_col]])
                probs <- as.numeric(dt[[pc]])
                is_correct <- if (correct_col %in% colnames(dt)) as.logical(dt[[correct_col]]) else rep(TRUE, length(probs))
                all_grps <- sort(unique(grp_vals))
                n_g <- length(all_grps)

                plot(NULL, xlim = c(0.5, n_g + 0.5), ylim = c(0, 1),
                     xlab = "True Group", ylab = "Posterior / Probability",
                     main = clf_name, xaxt = "n")
                axis(1, at = seq_len(n_g), labels = all_grps, las = 2, cex.axis = 0.8)

                for (gi in seq_len(n_g)) {
                    g <- all_grps[gi]
                    idx <- which(grp_vals == g)
                    p_vals <- probs[idx]
                    corr <- is_correct[idx]

                    # Boxplot manually
                    if (length(p_vals) > 0) {
                        qs <- quantile(p_vals, c(0.25, 0.5, 0.75), na.rm = TRUE)
                        iqr <- qs[3] - qs[1]
                        wlo <- max(min(p_vals, na.rm = TRUE), qs[1] - 1.5 * iqr)
                        whi <- min(max(p_vals, na.rm = TRUE), qs[3] + 1.5 * iqr)
                        bw <- 0.3

                        rect(gi - bw, qs[1], gi + bw, qs[3],
                             col = adjustcolor(group_colors[g], alpha.f = 0.3),
                             border = group_colors[g])
                        segments(gi - bw, qs[2], gi + bw, qs[2],
                                 col = group_colors[g], lwd = 2)
                        segments(gi, wlo, gi, qs[1], col = group_colors[g])
                        segments(gi, qs[3], gi, whi, col = group_colors[g])
                        segments(gi - bw/2, wlo, gi + bw/2, wlo, col = group_colors[g])
                        segments(gi - bw/2, whi, gi + bw/2, whi, col = group_colors[g])

                        # Jittered points (shapes by data source)
                        jitter_x <- gi + runif(length(p_vals), -0.15, 0.15)
                        sample_ids_dt <- if ("sample_id" %in% colnames(dt)) dt[["sample_id"]][idx] else rep("", length(idx))
                        src_vals <- ifelse(sample_ids_dt %in% names(source_map),
                                           source_map[sample_ids_dt], "Chip")
                        pch_vals <- source_pch[src_vals]
                        points(jitter_x, p_vals, pch = pch_vals, cex = 0.6,
                               col = adjustcolor(group_colors[g], alpha.f = 0.7))
                    }
                }
                abline(h = 0.5, lty = 3, col = "grey50")
            }
            par(mfrow = c(1, 1))
        } else {
            plot.new()
            text(0.5, 0.5, "No posterior/probability data available", cex = 1.3)
        }
    }, error = function(e) {
        cat("Page 6 (Posterior distributions) error:", e[["message"]], "\\n")
        par(mfrow = c(1, 1))
        plot.new(); text(0.5, 0.5, paste("Posterior distributions error:", e[["message"]]), cex = 0.9)
    })

    # =========================================================================
    # PAGE 7: ADMIXTURE Structure Plot
    # =========================================================================
    tryCatch({
        if (has_admixture && !is.null(adm_Q) && nrow(adm_Q) > 0) {
            sample_ids <- as.character(adm_Q[[1]])
            q_mat <- as.matrix(adm_Q[, -1, with = FALSE])
            rownames(q_mat) <- sample_ids
            K <- ncol(q_mat)

            # Get group assignments from reference labels
            if (!is.null(ref_labels)) {
                grp_map <- ref_labels[["group"]]
                names(grp_map) <- ref_labels[["sample_id"]]
                sample_groups <- grp_map[sample_ids]
                sample_groups[is.na(sample_groups)] <- "Unknown"
            } else {
                sample_groups <- rep("Unknown", length(sample_ids))
            }

            # Sort by group, then by dominant Q within group
            dom_q <- apply(q_mat, 1, which.max)
            dom_prop <- apply(q_mat, 1, max)
            sort_order <- order(sample_groups, dom_q, -dom_prop)

            q_sorted <- q_mat[sort_order, , drop = FALSE]
            grp_sorted <- sample_groups[sort_order]

            # Color palette for K clusters
            if (K <= 8) {
                q_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                              "#FF7F00", "#A65628", "#F781BF", "#999999")[seq_len(K)]
            } else {
                q_colors <- rainbow(K, s = 0.75, v = 0.85)
            }

            par(mar = c(8, 5, 4, 2))
            barplot(t(q_sorted), col = q_colors, border = NA, space = 0,
                    main = paste0("ADMIXTURE Structure (K=", K, ") - ", node_name),
                    ylab = "Ancestry Proportion",
                    xlab = "", xaxt = "n", ylim = c(0, 1))

            # Add group dividers and labels
            unique_grps <- unique(grp_sorted)
            grp_boundaries <- cumsum(table(factor(grp_sorted, levels = unique_grps)))
            grp_midpoints <- c(0, grp_boundaries[-length(grp_boundaries)])
            grp_midpoints <- (grp_midpoints + grp_boundaries) / 2

            for (b in grp_boundaries[-length(grp_boundaries)]) {
                abline(v = b, col = "black", lwd = 1.5)
            }
            mtext(unique_grps, side = 1, at = grp_midpoints, line = 1.5,
                  cex = max(0.6, min(0.9, 10 / length(unique_grps))), las = 2)

            legend("topright", legend = paste0("K", seq_len(K)),
                   fill = q_colors, border = NA, cex = 0.7, bty = "n")
        } else {
            plot.new()
            text(0.5, 0.5, "ADMIXTURE not available for this node", cex = 1.3, col = "grey40")
        }
    }, error = function(e) {
        cat("Page 7 (ADMIXTURE structure) error:", e[["message"]], "\\n")
        plot.new(); text(0.5, 0.5, paste("ADMIXTURE structure error:", e[["message"]]), cex = 0.9)
    })

    # =========================================================================
    # PAGE 8: ADMIXTURE CV Error Plot
    # =========================================================================
    tryCatch({
        if (has_admixture && length(adm_cv) > 0) {
            # Parse "K=N  CV=X.XXXXX" format
            cv_lines <- adm_cv[grepl("K=", adm_cv)]
            if (length(cv_lines) > 0) {
                k_vals <- as.numeric(gsub(".*K=([0-9]+).*", "\\\\1", cv_lines))
                cv_vals <- as.numeric(gsub(".*CV=([0-9.]+).*", "\\\\1", cv_lines))

                valid <- !is.na(k_vals) & !is.na(cv_vals)
                k_vals <- k_vals[valid]
                cv_vals <- cv_vals[valid]

                if (length(k_vals) > 0) {
                    opt_k <- k_vals[which.min(cv_vals)]

                    par(mar = c(5, 5, 4, 2))
                    plot(k_vals, cv_vals, type = "b", pch = 19, col = "steelblue",
                         lwd = 2, cex = 1.3,
                         xlab = "K (Number of Clusters)",
                         ylab = "Cross-Validation Error",
                         main = paste0("ADMIXTURE CV Error - ", node_name))
                    points(opt_k, cv_vals[which.min(cv_vals)],
                           pch = 18, col = "red", cex = 2.5)
                    text(opt_k, cv_vals[which.min(cv_vals)],
                         labels = paste0("Optimal K=", opt_k),
                         pos = 4, col = "red", cex = 1.1, font = 2)
                    grid(col = "grey80", lty = 3)
                } else {
                    plot.new()
                    text(0.5, 0.5, "Could not parse CV error values", cex = 1.3)
                }
            } else {
                plot.new()
                text(0.5, 0.5, "No CV error lines found", cex = 1.3)
            }
        } else {
            plot.new()
            text(0.5, 0.5, "ADMIXTURE CV error not available for this node", cex = 1.3, col = "grey40")
        }
    }, error = function(e) {
        cat("Page 8 (CV error) error:", e[["message"]], "\\n")
        plot.new(); text(0.5, 0.5, paste("ADMIXTURE CV error plot error:", e[["message"]]), cex = 0.9)
    })

    # =========================================================================
    # PAGE 9: Ensemble Tier Distribution
    # =========================================================================
    tryCatch({
        if (!is.null(ensemble) && nrow(ensemble) > 0) {
            # Identify column names -- handle both old and new ensemble formats
            assign_col <- intersect(c("assignment"), colnames(ensemble))[1]
            tier_col <- intersect(c("tier"), colnames(ensemble))[1]

            if (!is.na(assign_col) && !is.na(tier_col)) {
                ens_assigned <- as.character(ensemble[[assign_col]])
                ens_tier <- as.character(ensemble[[tier_col]])
                ens_sample <- as.character(ensemble[["sample_id"]])

                valid <- !is.na(ens_assigned) & !is.na(ens_tier)
                ens_assigned <- ens_assigned[valid]
                ens_tier <- ens_tier[valid]
                ens_sample <- ens_sample[valid]

                if (length(ens_assigned) > 0) {
                    # Determine which are reference vs unknown
                    ref_ids <- if (!is.null(ref_labels)) ref_labels[["sample_id"]] else character(0)
                    is_ref <- ens_sample %in% ref_ids

                    tier_levels <- c("Strongly Supported", "Supported", "Inconclusive")
                    tier_colors <- c("Strongly Supported" = "green3",
                                     "Supported" = "steelblue",
                                     "Inconclusive" = "orange")
                    all_groups_ens <- sort(unique(ens_assigned))

                    has_both <- any(is_ref) && any(!is_ref)

                    if (has_both) {
                        par(mfrow = c(1, 2), mar = c(7, 5, 4, 2))

                        # Reference samples
                        ct_ref <- table(factor(ens_tier[is_ref], levels = tier_levels),
                                        factor(ens_assigned[is_ref], levels = all_groups_ens))
                        barplot(ct_ref, col = tier_colors[rownames(ct_ref)], beside = FALSE,
                                main = "Reference Samples", ylab = "Count",
                                las = 2, cex.names = 0.8, border = NA)

                        # Unknown samples
                        ct_unk <- table(factor(ens_tier[!is_ref], levels = tier_levels),
                                        factor(ens_assigned[!is_ref], levels = all_groups_ens))
                        barplot(ct_unk, col = tier_colors[rownames(ct_unk)], beside = FALSE,
                                main = "Unknown Samples", ylab = "Count",
                                las = 2, cex.names = 0.8, border = NA)

                        legend("topright", legend = tier_levels, fill = tier_colors[tier_levels],
                               cex = 0.8, bty = "n", border = NA)
                        par(mfrow = c(1, 1))
                    } else {
                        ct <- table(factor(ens_tier, levels = tier_levels),
                                    factor(ens_assigned, levels = all_groups_ens))
                        par(mar = c(7, 5, 4, 8), xpd = NA)
                        barplot(ct, col = tier_colors[rownames(ct)], beside = FALSE,
                                main = paste0("Ensemble Tier Distribution - ", node_name),
                                ylab = "Count", las = 2, cex.names = 0.8, border = NA)
                        legend("topright", inset = c(-0.15, 0),
                               legend = tier_levels, fill = tier_colors[tier_levels],
                               cex = 0.8, bty = "n", border = NA)
                        par(xpd = FALSE)
                    }
                } else {
                    plot.new()
                    text(0.5, 0.5, "No valid ensemble assignments", cex = 1.3)
                }
            } else {
                plot.new()
                text(0.5, 0.5, "Ensemble file missing required columns", cex = 1.3)
            }
        } else {
            plot.new()
            text(0.5, 0.5, "No ensemble data available", cex = 1.3)
        }
    }, error = function(e) {
        cat("Page 9 (Ensemble tiers) error:", e[["message"]], "\\n")
        par(mfrow = c(1, 1))
        plot.new(); text(0.5, 0.5, paste("Ensemble tier error:", e[["message"]]), cex = 0.9)
    })

    # =========================================================================
    # PAGE 10: Classifier Agreement
    # =========================================================================
    tryCatch({
        if (!is.null(ensemble) && nrow(ensemble) > 0) {
            # Detect column names for each classifier
            dapc_grp_col <- intersect(c("dapc_group"), colnames(ensemble))[1]
            svm_grp_col <- intersect(c("svm_group", "ap_group"), colnames(ensemble))[1]
            pf_grp_col  <- intersect(c("pf_group"), colnames(ensemble))[1]
            dapc_prob_col <- intersect(c("dapc_posterior"), colnames(ensemble))[1]
            svm_prob_col  <- intersect(c("svm_probability", "ap_probability"), colnames(ensemble))[1]

            # Collect available classifier columns
            clf_names <- character(0)
            clf_grp_cols <- character(0)
            if (!is.na(dapc_grp_col)) { clf_names <- c(clf_names, "DAPC"); clf_grp_cols <- c(clf_grp_cols, dapc_grp_col) }
            if (!is.na(svm_grp_col)) { clf_names <- c(clf_names, "SVM/AP"); clf_grp_cols <- c(clf_grp_cols, svm_grp_col) }
            if (!is.na(pf_grp_col))  { clf_names <- c(clf_names, "popfinder"); clf_grp_cols <- c(clf_grp_cols, pf_grp_col) }

            n_clf <- length(clf_names)

            if (n_clf >= 2) {
                layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1.5))

                # Panel 1: Pairwise agreement matrix
                agree_mat <- matrix(NA, n_clf, n_clf)
                rownames(agree_mat) <- clf_names
                colnames(agree_mat) <- clf_names

                for (a in seq_len(n_clf)) {
                    for (b in seq_len(n_clf)) {
                        g1 <- as.character(ensemble[[clf_grp_cols[a]]])
                        g2 <- as.character(ensemble[[clf_grp_cols[b]]])
                        valid_idx <- !is.na(g1) & !is.na(g2)
                        if (sum(valid_idx) > 0) {
                            agree_mat[a, b] <- round(mean(g1[valid_idx] == g2[valid_idx]) * 100, 1)
                        }
                    }
                }

                par(mar = c(5, 6, 4, 2))
                col_ramp <- colorRampPalette(c("white", "green3"))(100)
                image(1:n_clf, 1:n_clf, agree_mat, col = col_ramp,
                      axes = FALSE, xlab = "", ylab = "",
                      main = "Pairwise Agreement (%)",
                      zlim = c(0, 100))
                axis(1, at = seq_len(n_clf), labels = clf_names, las = 2, cex.axis = 0.9)
                axis(2, at = seq_len(n_clf), labels = clf_names, las = 1, cex.axis = 0.9)
                for (a in seq_len(n_clf)) {
                    for (b in seq_len(n_clf)) {
                        if (!is.na(agree_mat[a, b])) {
                            text(a, b, paste0(agree_mat[a, b], "%"), cex = 1.0,
                                 col = if (agree_mat[a, b] > 70) "white" else "black")
                        }
                    }
                }
                box()

                # Panel 2: SVM/AP probability vs DAPC posterior scatter
                if (!is.na(dapc_prob_col) && !is.na(svm_prob_col) &&
                    !is.na(dapc_grp_col) && !is.na(svm_grp_col)) {

                    dapc_p <- as.numeric(ensemble[[dapc_prob_col]])
                    svm_p  <- as.numeric(ensemble[[svm_prob_col]])
                    dapc_g <- as.character(ensemble[[dapc_grp_col]])
                    svm_g  <- as.character(ensemble[[svm_grp_col]])
                    agree  <- !is.na(dapc_g) & !is.na(svm_g) & dapc_g == svm_g
                    disagree <- !is.na(dapc_g) & !is.na(svm_g) & dapc_g != svm_g

                    valid_pts <- !is.na(dapc_p) & !is.na(svm_p)

                    par(mar = c(5, 5, 4, 2))
                    plot(dapc_p[valid_pts], svm_p[valid_pts],
                         pch = ifelse(agree[valid_pts], 19, 1),
                         col = ifelse(agree[valid_pts],
                                      adjustcolor("green3", alpha.f = 0.6),
                                      adjustcolor("red", alpha.f = 0.6)),
                         cex = 1.0,
                         xlab = "DAPC Posterior",
                         ylab = "SVM/AP Probability",
                         main = "Classifier Probability Agreement")
                    abline(0, 1, lty = 2, col = "grey50")

                    # Highlight validator_conflict samples with red X marks
                    conflict_col <- intersect(c("validator_conflict"), colnames(ensemble))[1]
                    if (!is.na(conflict_col)) {
                        conflict_flag <- as.logical(ensemble[[conflict_col]])
                        conflict_idx <- which(valid_pts & !is.na(conflict_flag) & conflict_flag)
                        if (length(conflict_idx) > 0) {
                            points(dapc_p[conflict_idx], svm_p[conflict_idx],
                                   pch = 4, col = "red", cex = 2, lwd = 2.5)
                        }
                    }

                    n_agree <- sum(agree & valid_pts, na.rm = TRUE)
                    n_disagree <- sum(disagree & valid_pts, na.rm = TRUE)
                    legend("bottomright",
                           legend = c(paste0("Agree (n=", n_agree, ")"),
                                      paste0("Disagree (n=", n_disagree, ")"),
                                      "Conflict"),
                           col = c("green3", "red", "red"),
                           pch = c(19, 1, 4), cex = 0.9, bty = "n")
                } else {
                    plot.new()
                    text(0.5, 0.5, "Insufficient data for scatter", cex = 1.1)
                }

                layout(1)
            } else {
                plot.new()
                text(0.5, 0.5, "Need at least 2 classifiers for agreement analysis", cex = 1.3)
            }
        } else {
            plot.new()
            text(0.5, 0.5, "No ensemble data for classifier agreement", cex = 1.3)
        }
    }, error = function(e) {
        cat("Page 10 (Classifier agreement) error:", e[["message"]], "\\n")
        layout(1)
        plot.new(); text(0.5, 0.5, paste("Classifier agreement error:", e[["message"]]), cex = 0.9)
    })

    # =========================================================================
    # PAGE 11: Per-Group Accuracy Breakdown
    # =========================================================================
    tryCatch({
        # Compute per-group accuracy for each method
        per_group_acc <- list()

        compute_per_group <- function(loocv_dt, method_name) {
            if (is.null(loocv_dt) || nrow(loocv_dt) == 0) return(NULL)
            true_col <- intersect(c("true_group", "true_pop"), colnames(loocv_dt))[1]
            pred_col <- intersect(c("predicted_group", "predicted_pop"), colnames(loocv_dt))[1]
            if (is.na(true_col) || is.na(pred_col)) return(NULL)
            true_vals <- as.character(loocv_dt[[true_col]])
            pred_vals <- as.character(loocv_dt[[pred_col]])
            all_grps <- sort(unique(true_vals))
            accs <- sapply(all_grps, function(g) {
                idx <- true_vals == g
                if (sum(idx) == 0) return(NA)
                mean(pred_vals[idx] == g, na.rm = TRUE)
            })
            names(accs) <- all_grps
            accs
        }

        dapc_pga <- compute_per_group(dapc_loocv, "DAPC")
        if (!is.null(dapc_pga)) per_group_acc[["DAPC"]] <- dapc_pga

        svm_pga <- compute_per_group(ap_loocv_svm, "SVM")
        if (!is.null(svm_pga)) per_group_acc[["SVM"]] <- svm_pga

        bay_pga <- compute_per_group(ap_loocv_bay, "LDA/Bayesian")
        if (!is.null(bay_pga)) per_group_acc[["LDA/Bayesian"]] <- bay_pga

        rf_pga <- compute_per_group(ap_loocv_rf, "Random Forest")
        if (!is.null(rf_pga)) per_group_acc[["Random Forest"]] <- rf_pga

        pf_pga <- compute_per_group(pf_loocv, "popfinder")
        if (!is.null(pf_pga)) per_group_acc[["popfinder"]] <- pf_pga

        if (length(per_group_acc) > 0) {
            # Unify group names across methods
            all_grps_union <- sort(unique(unlist(lapply(per_group_acc, names))))
            method_names_pg <- names(per_group_acc)

            # Build matrix: methods x groups
            acc_matrix <- matrix(NA, nrow = length(method_names_pg), ncol = length(all_grps_union))
            rownames(acc_matrix) <- method_names_pg
            colnames(acc_matrix) <- all_grps_union
            for (mn in method_names_pg) {
                matched <- match(names(per_group_acc[[mn]]), all_grps_union)
                acc_matrix[mn, matched] <- per_group_acc[[mn]]
            }

            method_colors <- c("DAPC" = "#377EB8", "SVM" = "#E41A1C",
                               "LDA/Bayesian" = "#984EA3",
                               "Random Forest" = "#4DAF4A",
                               "popfinder" = "#FF7F00")
            mc <- ifelse(method_names_pg %in% names(method_colors),
                         method_colors[method_names_pg], "grey50")

            par(mar = c(7, 5, 4, 10), xpd = NA)
            bp <- barplot(acc_matrix * 100, beside = TRUE, col = mc,
                          main = paste0("Per-Group Accuracy - ", node_name),
                          ylab = "Accuracy (%)", ylim = c(0, 110),
                          las = 2, cex.names = 0.8, border = NA)
            abline(h = 95, lty = 2, col = "red", lwd = 1)
            legend("topright", inset = c(-0.15, 0),
                   legend = method_names_pg, fill = mc,
                   cex = 0.7, bty = "n", border = NA)
            par(xpd = FALSE)
        } else {
            plot.new()
            text(0.5, 0.5, "No per-group accuracy data available", cex = 1.3)
        }
    }, error = function(e) {
        cat("Page 11 (Per-group accuracy) error:", e[["message"]], "\\n")
        plot.new(); text(0.5, 0.5, paste("Per-group accuracy error:", e[["message"]]), cex = 0.9)
    })

    # =========================================================================
    # PAGE 12: Misclassification Summary Table
    # =========================================================================
    tryCatch({
        # Collect all misclassified samples across methods
        collect_misclass <- function(loocv_dt, method_name) {
            if (is.null(loocv_dt) || nrow(loocv_dt) == 0) return(NULL)
            true_col <- intersect(c("true_group", "true_pop"), colnames(loocv_dt))[1]
            pred_col <- intersect(c("predicted_group", "predicted_pop"), colnames(loocv_dt))[1]
            if (is.na(true_col) || is.na(pred_col)) return(NULL)
            wrong <- loocv_dt[[true_col]] != loocv_dt[[pred_col]]
            wrong[is.na(wrong)] <- TRUE
            if (!any(wrong)) return(NULL)
            data.table(
                sample_id = loocv_dt[["sample_id"]][wrong],
                true_group = as.character(loocv_dt[[true_col]][wrong]),
                method = method_name,
                predicted = as.character(loocv_dt[[pred_col]][wrong])
            )
        }

        misclass_list <- list()
        mc_dapc <- collect_misclass(dapc_loocv, "DAPC")
        if (!is.null(mc_dapc)) misclass_list[["DAPC"]] <- mc_dapc
        mc_svm <- collect_misclass(ap_loocv_svm, "SVM")
        if (!is.null(mc_svm)) misclass_list[["SVM"]] <- mc_svm
        mc_bay <- collect_misclass(ap_loocv_bay, "LDA/Bayesian")
        if (!is.null(mc_bay)) misclass_list[["LDA/Bayesian"]] <- mc_bay
        mc_rf <- collect_misclass(ap_loocv_rf, "Random Forest")
        if (!is.null(mc_rf)) misclass_list[["Random Forest"]] <- mc_rf
        mc_pf <- collect_misclass(pf_loocv, "popfinder")
        if (!is.null(mc_pf)) misclass_list[["popfinder"]] <- mc_pf

        if (length(misclass_list) > 0) {
            all_misclass <- rbindlist(misclass_list)

            # Pivot to wide: one row per sample, columns for each method's prediction
            unique_samples <- unique(all_misclass[, .(sample_id, true_group)])
            method_list <- unique(all_misclass[["method"]])

            # Build a wide table
            wide <- unique_samples
            for (meth in method_list) {
                mc_sub <- all_misclass[all_misclass[["method"]] == meth]
                mc_sub_dedup <- mc_sub[, .(pred = predicted[1]), by = sample_id]
                setnames(mc_sub_dedup, "pred", meth)
                wide <- merge(wide, mc_sub_dedup, by = "sample_id", all.x = TRUE)
            }

            # Render as text on plot
            plot.new()
            par(mar = c(1, 2, 3, 1))
            title(main = paste0("Misclassification Summary - ", node_name), line = 1)

            # Table dimensions
            col_names <- c("Sample", "True Group", method_list)
            n_cols <- length(col_names)
            n_rows <- min(nrow(wide), 40)  # Cap to avoid overflow

            # Max characters per column (truncate long strings)
            max_chars <- max(8, floor(80 / n_cols))

            if (n_rows == 0) {
                text(0.5, 0.5, "No misclassifications!", cex = 1.5, col = "green4")
            } else {
                # Calculate column widths
                col_width <- 1 / n_cols
                row_height <- min(0.025, 0.85 / (n_rows + 1))
                start_y <- 0.95

                # Helper to truncate text
                trunc_text <- function(s, maxc) {
                    ifelse(nchar(s) > maxc, paste0(substr(s, 1, maxc - 2), ".."), s)
                }

                # Header
                for (ci in seq_along(col_names)) {
                    x_pos <- (ci - 1) * col_width + col_width / 2
                    text(x_pos, start_y, trunc_text(col_names[ci], max_chars), font = 2,
                         cex = max(0.5, min(0.75, 5 / n_cols)))
                }
                segments(0, start_y - row_height * 0.4, 1, start_y - row_height * 0.4, col = "grey30")

                # Data rows
                for (ri in seq_len(n_rows)) {
                    y_pos <- start_y - row_height * (ri + 0.5)
                    if (y_pos < 0.02) break

                    row_data <- c(
                        as.character(wide[["sample_id"]][ri]),
                        as.character(wide[["true_group"]][ri])
                    )
                    for (meth in method_list) {
                        val <- wide[[meth]][ri]
                        row_data <- c(row_data, if (is.na(val)) "-" else as.character(val))
                    }

                    for (ci in seq_along(row_data)) {
                        x_pos <- (ci - 1) * col_width + col_width / 2
                        text(x_pos, y_pos, trunc_text(row_data[ci], max_chars),
                             cex = max(0.45, min(0.65, 4 / n_cols)),
                             col = if (ci > 2 && row_data[ci] != "-") "red3" else "black")
                    }
                }

                if (nrow(wide) > 40) {
                    text(0.5, 0.01,
                         paste0("... and ", nrow(wide) - 40, " more misclassified samples (truncated)"),
                         cex = 0.7, col = "grey40", font = 3)
                }
            }
        } else {
            plot.new()
            text(0.5, 0.6, paste0("Misclassification Summary - ", node_name),
                 cex = 1.5, font = 2)
            text(0.5, 0.4, "No misclassifications across all methods!",
                 cex = 1.5, col = "green4", font = 3)
        }
    }, error = function(e) {
        cat("Page 12 (Misclassification table) error:", e[["message"]], "\\n")
        plot.new(); text(0.5, 0.5, paste("Misclassification table error:", e[["message"]]), cex = 0.9)
    })

    # =========================================================================
    # Close PDF
    # =========================================================================
    dev.off()
    cat("Diagnostic PDF generated:", paste0("${node_name}_plots.pdf"), "\\n")
    """
}
