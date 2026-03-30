/*
========================================================================================
    ENSEMBLE DECISION -- SVM-PRIMARY LOGIC
========================================================================================
    Assignment always follows assignPOP's SVM (loocv_svm.tsv).
    DAPC and popfinder serve as validators, not equal voters.
    Tiers:
        Strongly Supported -- SVM + both validators agree (or SVM + DAPC when popfinder absent)
        Supported          -- SVM + one validator agrees (or SVM alone when only one
                              validator exists and it disagrees but SVM prob > upper threshold)
        Provisional        -- SVM assignment but neither validator agrees
    Threshold calibration retained for supplementary confidence info.
*/

process ENSEMBLE_DECISION {
    tag "${node_name}"
    label 'process_low'

    publishDir "${params.outdir}/classification/ensemble", mode: 'copy'

    input:
    tuple val(node_name), path(dapc_dir), path(assignpop_dir), path(popfinder_dir)

    output:
    path "${node_name}_ensemble.tsv", emit: decisions

    script:
    def has_popfinder = popfinder_dir && popfinder_dir.name != 'NO_POPFINDER' && popfinder_dir.name != 'null'
    """
    #!/usr/bin/env Rscript

    library(data.table)

    node_name <- "${node_name}"
    cat("Ensemble decision (SVM-primary) for node:", node_name, "\\n")

    # -- Load SVM LOOCV (always use SVM regardless of assignPOP best_model) --
    svm_loocv_file <- file.path("${assignpop_dir}", "loocv_svm.tsv")
    if (!file.exists(svm_loocv_file)) {
        stop("loocv_svm.tsv not found in assignPOP results for node: ", node_name)
    }
    svm_loocv <- fread(svm_loocv_file, header = TRUE)
    cat("Loaded SVM LOOCV:", nrow(svm_loocv), "reference samples\\n")

    # -- Load DAPC LOOCV --
    dapc_loocv <- fread(file.path("${dapc_dir}", "loocv_results.tsv"), header = TRUE)

    # -- Load popfinder LOOCV (optional) --
    has_popfinder <- file.exists(file.path("${popfinder_dir}", "loocv_results.tsv"))
    if (has_popfinder) {
        pf_loocv <- fread(file.path("${popfinder_dir}", "loocv_results.tsv"), header = TRUE)
        cat("popfinder LOOCV loaded:", nrow(pf_loocv), "samples\\n")
    } else {
        cat("popfinder not available at this node\\n")
    }

    # -- Calibrate thresholds (supplementary confidence info) --
    calibrate_threshold <- function(loocv_dt, prob_col = "posterior_max") {
        if (!prob_col %in% colnames(loocv_dt)) {
            prob_col <- "probability"
        }
        if (!prob_col %in% colnames(loocv_dt)) {
            return(list(upper = 0.95, middle = 0.80))
        }
        correct_probs <- loocv_dt[correct == TRUE][[prob_col]]
        if (length(correct_probs) < 3) return(list(upper = 0.95, middle = 0.80))

        upper <- quantile(correct_probs, 0.05, na.rm = TRUE)  # 95% of correct above this
        middle <- quantile(correct_probs, 0.20, na.rm = TRUE)  # 80% of correct above this

        upper <- max(upper, 0.70)
        middle <- max(middle, 0.50)

        list(upper = as.numeric(upper), middle = as.numeric(middle))
    }

    svm_thresh <- calibrate_threshold(svm_loocv, "probability")
    dapc_thresh <- calibrate_threshold(dapc_loocv, "posterior_max")
    cat("SVM thresholds  -- upper:", svm_thresh[["upper"]], " middle:", svm_thresh[["middle"]], "\\n")
    cat("DAPC thresholds -- upper:", dapc_thresh[["upper"]], " middle:", dapc_thresh[["middle"]], "\\n")

    if (has_popfinder) {
        pf_thresh <- calibrate_threshold(pf_loocv, "posterior_max")
        cat("popfinder thresholds -- upper:", pf_thresh[["upper"]], " middle:", pf_thresh[["middle"]], "\\n")
    }

    # -- Load unknown predictions --
    svm_preds <- NULL
    dapc_preds <- NULL
    pf_preds   <- NULL

    # SVM unknown predictions: assignPOP writes unknown_predictions.tsv using its best model.
    # We need SVM specifically, but unknowns are predicted with best_model only.
    # The unknown_predictions.tsv from assignPOP uses the best model, so we load it
    # and note that SVM was used if best_model == SVM, otherwise we lack SVM unknowns.
    # For reference samples we use loocv_svm.tsv directly.
    ap_unk_file <- file.path("${assignpop_dir}", "unknown_predictions.tsv")
    if (file.exists(ap_unk_file)) svm_preds <- fread(ap_unk_file, header = TRUE)

    dapc_unk_file <- file.path("${dapc_dir}", "unknown_predictions.tsv")
    if (file.exists(dapc_unk_file)) dapc_preds <- fread(dapc_unk_file, header = TRUE)

    if (has_popfinder) {
        pf_unk_file <- file.path("${popfinder_dir}", "unknown_predictions.tsv")
        if (file.exists(pf_unk_file)) pf_preds <- fread(pf_unk_file, header = TRUE)
    }

    # -- Build per-sample table (reference + unknown) --
    # Reference: use LOOCV predicted groups
    ref_svm <- svm_loocv[, .(sample_id,
                              svm_group = predicted_group,
                              svm_probability = probability)]

    ref_dapc <- dapc_loocv[, .(sample_id,
                                dapc_group = predicted_group,
                                dapc_posterior = posterior_max)]

    ensemble <- merge(ref_svm, ref_dapc, by = "sample_id", all = TRUE)

    if (has_popfinder) {
        ref_pf <- pf_loocv[, .(sample_id,
                                pf_group = predicted_group,
                                pf_probability = posterior_max)]
        ensemble <- merge(ensemble, ref_pf, by = "sample_id", all = TRUE)
    }

    # Unknown samples
    if (!is.null(svm_preds) || !is.null(dapc_preds) || !is.null(pf_preds)) {
        unk_svm <- if (!is.null(svm_preds)) {
            svm_preds[, .(sample_id,
                          svm_group = predicted_group,
                          svm_probability = probability)]
        } else NULL

        unk_dapc <- if (!is.null(dapc_preds)) {
            dapc_preds[, .(sample_id,
                           dapc_group = predicted_group,
                           dapc_posterior = posterior_max)]
        } else NULL

        unk_pf <- if (!is.null(pf_preds)) {
            pf_preds[, .(sample_id,
                         pf_group = predicted_group,
                         pf_probability = probability)]
        } else NULL

        unk_ensemble <- unk_svm
        if (!is.null(unk_dapc)) {
            if (is.null(unk_ensemble)) {
                unk_ensemble <- unk_dapc
            } else {
                unk_ensemble <- merge(unk_ensemble, unk_dapc, by = "sample_id", all = TRUE)
            }
        }
        if (!is.null(unk_pf) && has_popfinder) {
            if (is.null(unk_ensemble)) {
                unk_ensemble <- unk_pf
            } else {
                unk_ensemble <- merge(unk_ensemble, unk_pf, by = "sample_id", all = TRUE)
            }
        }

        if (!is.null(unk_ensemble)) {
            # Ensure same columns as reference ensemble before rbind
            for (col in setdiff(colnames(ensemble), colnames(unk_ensemble))) {
                unk_ensemble[[col]] <- NA
            }
            for (col in setdiff(colnames(unk_ensemble), colnames(ensemble))) {
                ensemble[[col]] <- NA
            }
            ensemble <- rbind(ensemble, unk_ensemble[, colnames(ensemble), with = FALSE])
        }
    }

    # Ensure popfinder columns exist even if popfinder absent
    if (!"pf_group" %in% colnames(ensemble)) {
        ensemble[, pf_group := NA_character_]
    }
    if (!"pf_probability" %in% colnames(ensemble)) {
        ensemble[, pf_probability := NA_real_]
    }

    cat("Total samples in ensemble:", nrow(ensemble), "\\n")

    # -- Apply SVM-primary decision rules --
    ensemble[, c("assignment", "tier", "n_validators_agree",
                 "validators_agreeing", "validators_dissenting",
                 "validator_conflict", "validator_consensus_group") := list(
        NA_character_, NA_character_, NA_integer_,
        NA_character_, NA_character_,
        FALSE, NA_character_
    )]

    for (i in 1:nrow(ensemble)) {
        svm_g  <- ensemble[[i, "svm_group"]]
        svm_p  <- ensemble[[i, "svm_probability"]]
        dapc_g <- ensemble[[i, "dapc_group"]]
        dapc_p <- ensemble[[i, "dapc_posterior"]]
        pf_g   <- ensemble[[i, "pf_group"]]
        pf_p   <- ensemble[[i, "pf_probability"]]

        # If SVM assignment is missing, mark as Provisional with NA
        if (is.na(svm_g)) {
            set(ensemble, i, "assignment", NA_character_)
            set(ensemble, i, "tier", "Provisional")
            set(ensemble, i, "n_validators_agree", 0L)
            next
        }

        # Assignment always follows SVM
        set(ensemble, i, "assignment", svm_g)

        # Determine which validators agree
        dapc_agrees <- !is.na(dapc_g) && dapc_g == svm_g
        pf_agrees   <- !is.na(pf_g) && pf_g == svm_g
        pf_ran      <- !is.na(pf_g)

        agreeing   <- character()
        dissenting <- character()
        if (!is.na(dapc_g)) {
            if (dapc_agrees) agreeing <- c(agreeing, "DAPC") else dissenting <- c(dissenting, "DAPC")
        }
        if (pf_ran) {
            if (pf_agrees) agreeing <- c(agreeing, "popfinder") else dissenting <- c(dissenting, "popfinder")
        }

        n_agree <- length(agreeing)
        set(ensemble, i, "n_validators_agree", as.integer(n_agree))
        set(ensemble, i, "validators_agreeing",
            if (length(agreeing) > 0) paste(agreeing, collapse = ",") else NA_character_)
        set(ensemble, i, "validators_dissenting",
            if (length(dissenting) > 0) paste(dissenting, collapse = ",") else NA_character_)

        # Validator conflict: both DAPC and popfinder agree on a DIFFERENT group than SVM
        val_conflict <- FALSE
        val_consensus <- NA_character_
        if (!is.na(dapc_g) && pf_ran) {
            if (!dapc_agrees && !pf_agrees && dapc_g == pf_g) {
                val_conflict <- TRUE
                val_consensus <- dapc_g
            }
        }
        set(ensemble, i, "validator_conflict", val_conflict)
        set(ensemble, i, "validator_consensus_group", val_consensus)

        # -- Tier assignment --
        if (pf_ran) {
            # Both validators available
            if (dapc_agrees && pf_agrees) {
                tier_val <- "Strongly Supported"
            } else if (dapc_agrees || pf_agrees) {
                tier_val <- "Supported"
            } else {
                tier_val <- "Provisional"
            }
        } else {
            # Only DAPC as validator (popfinder did not run)
            if (dapc_agrees) {
                tier_val <- "Strongly Supported"
            } else if (!is.na(dapc_g)) {
                # DAPC ran but disagrees; promote to Supported if
                # SVM probability exceeds its calibrated upper threshold
                if (!is.na(svm_p) && svm_p >= svm_thresh[["upper"]]) {
                    tier_val <- "Supported"
                } else {
                    tier_val <- "Provisional"
                }
            } else {
                # DAPC also missing -- SVM alone
                if (!is.na(svm_p) && svm_p >= svm_thresh[["upper"]]) {
                    tier_val <- "Supported"
                } else {
                    tier_val <- "Provisional"
                }
            }
        }

        set(ensemble, i, "tier", tier_val)
    }

    # -- Add node column --
    ensemble[, node := node_name]

    # -- Reorder columns to specified output --
    out_cols <- c("sample_id", "svm_group", "svm_probability",
                  "dapc_group", "dapc_posterior",
                  "pf_group", "pf_probability",
                  "assignment", "tier", "n_validators_agree",
                  "validators_agreeing", "validators_dissenting",
                  "validator_conflict", "validator_consensus_group", "node")
    setcolorder(ensemble, intersect(out_cols, colnames(ensemble)))

    # -- Summary --
    cat("\\nEnsemble results (SVM-primary):\\n")
    cat("  Strongly Supported:", sum(ensemble[["tier"]] == "Strongly Supported", na.rm = TRUE), "\\n")
    cat("  Supported:         ", sum(ensemble[["tier"]] == "Supported", na.rm = TRUE), "\\n")
    cat("  Provisional:       ", sum(ensemble[["tier"]] == "Provisional", na.rm = TRUE), "\\n")
    cat("  Validator conflicts:", sum(ensemble[["validator_conflict"]] == TRUE, na.rm = TRUE), "\\n")

    fwrite(ensemble, "${node_name}_ensemble.tsv", sep = "\\t")
    cat("Ensemble complete for", node_name, "\\n")
    """
}
