/*
========================================================================================
    COMPILE FINAL ASSIGNMENTS
========================================================================================
    Cascade through hierarchy levels, applying gate checks, and generate the final
    structured output with all fields.
*/

process COMPILE_ASSIGNMENTS {
    tag "compile"
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path ensemble_files   // collected from all nodes
    path metadata
    path sample_stats

    output:
    path "final_assignments.tsv",     emit: final_assignments
    path "assignment_summary.txt",    emit: summary

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    sink("assignment_summary.txt")
    cat("=== Final Assignment Summary ===\\n\\n")

    meta <- fread("${metadata}", header = TRUE)
    stats <- fread("${sample_stats}", header = TRUE)

    # Normalize metadata column names
    if ("Sample" %in% colnames(meta) && !"sample_id" %in% colnames(meta)) setnames(meta, "Sample", "sample_id")
    if ("Subspecies" %in% colnames(meta) && !"subspecies" %in% colnames(meta)) setnames(meta, "Subspecies", "subspecies")
    if ("Ecotype_Cleaned" %in% colnames(meta) && !"ecotype" %in% colnames(meta)) setnames(meta, "Ecotype_Cleaned", "ecotype")
    if ("Herd_Cleaned" %in% colnames(meta) && !"herd" %in% colnames(meta)) setnames(meta, "Herd_Cleaned", "herd")

    # Auto-detect ID column (fallback for other formats)
    id_col <- intersect(c("sample_id", "ID", "Sample_ID"), colnames(meta))[1]
    if (!is.na(id_col) && id_col != "sample_id") setnames(meta, id_col, "sample_id")

    # Load all ensemble results
    ens_files <- list.files(".", pattern = "_ensemble\\\\.tsv\$", full.names = TRUE)
    all_ens <- rbindlist(lapply(ens_files, fread, header = TRUE), fill = TRUE)

    cat("Loaded", length(ens_files), "node results\\n")
    cat("Total sample-node entries:", nrow(all_ens), "\\n\\n")

    # Separate by level
    l1 <- all_ens[grepl("^root\$|^L1_|^level1", node, ignore.case = TRUE)]
    l2 <- all_ens[grepl("^L2_|^level2", node, ignore.case = TRUE)]
    l3 <- all_ens[grepl("^L3_|^level3", node, ignore.case = TRUE)]

    # Get unique samples
    all_samples <- unique(all_ens\$sample_id)

    # Build final output table
    output <- data.table(sample_id = all_samples)

    # Merge sample stats
    output <- merge(output, stats[, .(sample_id, call_rate, het_rate)], by = "sample_id", all.x = TRUE)

    # Level 1
    if (nrow(l1) > 0) {
        l1_out <- l1[, .(
            level1_assignment = assignment,
            level1_tier = tier,
            level1_dapc_group = dapc_group,
            level1_dapc_posterior = dapc_posterior,
            level1_svm_group = svm_group,
            level1_svm_probability = svm_probability
        ), by = sample_id]
        if ("pf_group" %in% colnames(l1)) {
            pf_cols <- l1[, .(level1_pf_group = pf_group, level1_pf_probability = pf_probability), by = sample_id]
            l1_out <- merge(l1_out, pf_cols, by = "sample_id", all.x = TRUE)
        }
        output <- merge(output, l1_out, by = "sample_id", all.x = TRUE)
    }

    # Level 2: route each sample to the L2 node matching its L1 assignment
    if (nrow(l2) > 0 && "level1_assignment" %in% colnames(output)) {
        # Build expected L2 node name from L1 assignment
        output[, expected_l2_node := paste0("L2_", gsub(" ", "_", level1_assignment))]

        # Filter L2 results to only matching nodes
        l2_routed <- merge(l2, output[, .(sample_id, expected_l2_node)], by = "sample_id")
        l2_routed <- l2_routed[node == expected_l2_node]

        l2_out <- l2_routed[, .(
            level2_assignment = assignment[1],
            level2_tier = tier[1],
            level2_dapc_group = dapc_group[1],
            level2_svm_group = svm_group[1]
        ), by = sample_id]

        output[, expected_l2_node := NULL]
        output <- merge(output, l2_out, by = "sample_id", all.x = TRUE)

        # Gate check: clear L2 if L1 was Inconclusive or Provisional
        output[is.na(level1_tier) | level1_tier %in% c("Inconclusive", "Provisional"), c("level2_assignment", "level2_tier") := list(NA, NA)]
    }

    # Level 3: route each sample to the L3 node matching its L1+L2 assignments
    if (nrow(l3) > 0 && "level2_assignment" %in% colnames(output)) {
        output[, expected_l3_node := paste0("L3_", gsub(" ", "_", level1_assignment), "_", gsub(" ", "_", level2_assignment))]

        l3_routed <- merge(l3, output[, .(sample_id, expected_l3_node)], by = "sample_id")
        l3_routed <- l3_routed[node == expected_l3_node]

        l3_out <- l3_routed[, .(
            level3_assignment = assignment[1],
            level3_tier = tier[1],
            level3_dapc_group = dapc_group[1],
            level3_svm_group = svm_group[1]
        ), by = sample_id]

        output[, expected_l3_node := NULL]
        output <- merge(output, l3_out, by = "sample_id", all.x = TRUE)

        # Gate check: clear L3 if L2 was Inconclusive, Provisional, or NA
        output[is.na(level2_tier) | level2_tier %in% c("Inconclusive", "Provisional"), c("level3_assignment", "level3_tier") := list(NA, NA)]
    }

    # Determine deepest confident assignment (vectorized to handle multi-row groups)
    has_l3 <- "level3_tier" %in% colnames(output)
    has_l2 <- "level2_tier" %in% colnames(output)
    has_l1 <- "level1_tier" %in% colnames(output)
    confident_tiers <- c("Strongly Supported", "Supported")

    output[, deepest_confident_assignment := NA_character_]
    if (has_l1) output[!is.na(level1_tier) & level1_tier %in% confident_tiers, deepest_confident_assignment := level1_assignment]
    if (has_l2) output[!is.na(level2_tier) & level2_tier %in% confident_tiers, deepest_confident_assignment := level2_assignment]
    if (has_l3) output[!is.na(level3_tier) & level3_tier %in% confident_tiers, deepest_confident_assignment := level3_assignment]

    # Quality flags - vectorized
    output[, quality_flags := NA_character_]
    output[!is.na(call_rate) & call_rate < 0.80, quality_flags := "low_call_rate"]
    output[!is.na(het_rate) & het_rate > 0.42, quality_flags := fifelse(
        is.na(quality_flags), "possible_mixed", paste0(quality_flags, ";possible_mixed")
    )]

    # Summary statistics
    cat("\\nTotal samples:", nrow(output), "\\n")
    if ("level1_tier" %in% colnames(output)) {
        cat("\\nLevel 1 tiers:\\n")
        print(table(output\$level1_tier, useNA = "ifany"))
    }
    if ("level2_tier" %in% colnames(output)) {
        cat("\\nLevel 2 tiers:\\n")
        print(table(output\$level2_tier, useNA = "ifany"))
    }
    if ("level3_tier" %in% colnames(output)) {
        cat("\\nLevel 3 tiers:\\n")
        print(table(output\$level3_tier, useNA = "ifany"))
    }
    cat("\\nDeepest confident assignments:\\n")
    print(table(output\$deepest_confident_assignment, useNA = "ifany"))

    sink()

    fwrite(output, "final_assignments.tsv", sep = "\\t")
    cat("Final assignments compiled:", nrow(output), "samples\\n")
    """
}
