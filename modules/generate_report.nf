/*
========================================================================================
    GENERATE FINAL REPORT
========================================================================================
    Creates a comprehensive summary report of the entire pipeline run.
*/

process GENERATE_REPORT {
    tag "report"
    label 'process_medium'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path final_assignments
    path importance_summary
    path pca_results
    path hierarchy_json

    output:
    path "pipeline_report.pdf",  emit: report
    path "pipeline_summary.txt", emit: summary

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(jsonlite)

    sink("pipeline_summary.txt")
    cat("================================================================\\n")
    cat("  RANGIFER ASSIGNMENT PIPELINE - SUMMARY REPORT\\n")
    cat("================================================================\\n\\n")

    # Load assignments
    assignments <- fread("${final_assignments}", header = TRUE)
    hierarchy <- fromJSON("${hierarchy_json}")

    cat("Date:", as.character(Sys.Date()), "\\n\\n")
    cat("Total samples processed:", nrow(assignments), "\\n\\n")

    # Hierarchy summary
    cat("--- HIERARCHY STRUCTURE ---\\n\\n")
    if (!is.null(hierarchy\$level1)) {
        cat("Level 1:", paste(hierarchy\$level1\$groups, collapse = ", "), "\\n")
        cat("  Reference samples:", sum(unlist(hierarchy\$level1\$sample_counts)), "\\n\\n")
    }
    if (!is.null(hierarchy\$level2)) {
        for (key in names(hierarchy\$level2)) {
            node <- hierarchy\$level2[[key]]
            cat("Level 2 (", key, "):", paste(node\$groups, collapse = ", "), "\\n")
            cat("  Reference samples:", sum(unlist(node\$sample_counts)), "\\n")
        }
        cat("\\n")
    }

    # Assignment results
    cat("--- ASSIGNMENT RESULTS ---\\n\\n")

    if ("level1_tier" %in% colnames(assignments)) {
        cat("Level 1 Assignments:\\n")
        cat("  Strongly Supported:", sum(assignments\$level1_tier == "Strongly Supported", na.rm = TRUE), "\\n")
        cat("  Supported:", sum(assignments\$level1_tier == "Supported", na.rm = TRUE), "\\n")
        cat("  Provisional:", sum(assignments\$level1_tier == "Provisional", na.rm = TRUE), "\\n")
        cat("  Inconclusive:", sum(assignments\$level1_tier == "Inconclusive", na.rm = TRUE), "\\n")
        if ("level1_assignment" %in% colnames(assignments)) {
            cat("\\n  By group:\\n")
            print(table(assignments\$level1_assignment, assignments\$level1_tier, useNA = "ifany"))
        }
        cat("\\n")
    }

    if ("level2_tier" %in% colnames(assignments)) {
        cat("Level 2 Assignments:\\n")
        cat("  Strongly Supported:", sum(assignments\$level2_tier == "Strongly Supported", na.rm = TRUE), "\\n")
        cat("  Supported:", sum(assignments\$level2_tier == "Supported", na.rm = TRUE), "\\n")
        cat("  Provisional:", sum(assignments\$level2_tier == "Provisional", na.rm = TRUE), "\\n")
        cat("  Inconclusive:", sum(assignments\$level2_tier == "Inconclusive", na.rm = TRUE), "\\n")
        cat("\\n")
    }

    if ("level3_tier" %in% colnames(assignments)) {
        cat("Level 3 Assignments:\\n")
        cat("  Strongly Supported:", sum(assignments\$level3_tier == "Strongly Supported", na.rm = TRUE), "\\n")
        cat("  Supported:", sum(assignments\$level3_tier == "Supported", na.rm = TRUE), "\\n")
        cat("  Provisional:", sum(assignments\$level3_tier == "Provisional", na.rm = TRUE), "\\n")
        cat("  Inconclusive:", sum(assignments\$level3_tier == "Inconclusive", na.rm = TRUE), "\\n")
        cat("\\n")
    }

    # Deepest confident
    cat("Deepest Confident Assignments:\\n")
    dca <- table(assignments\$deepest_confident_assignment, useNA = "ifany")
    print(dca)

    # Quality flags
    cat("\\n--- QUALITY FLAGS ---\\n")
    if ("quality_flags" %in% colnames(assignments)) {
        flags <- assignments[!is.na(quality_flags)]
        cat("Samples with quality flags:", nrow(flags), "\\n")
        if (nrow(flags) > 0) {
            print(table(flags\$quality_flags))
        }
    }

    sink()

    # Generate PDF report
    pdf("pipeline_report.pdf", width = 12, height = 10)

    # Title page
    plot.new()
    text(0.5, 0.7, "Rangifer Assignment Pipeline", cex = 2, font = 2)
    text(0.5, 0.5, paste("Date:", Sys.Date()), cex = 1.2)
    text(0.5, 0.4, paste("Samples:", nrow(assignments)), cex = 1.2)

    # Assignment tier distribution
    if ("level1_tier" %in% colnames(assignments)) {
        tier_counts <- table(assignments\$level1_tier, useNA = "ifany")
        barplot(tier_counts, main = "Level 1 Assignment Tiers",
                col = c("Inconclusive" = "red", "Strongly Supported" = "green3",
                         "Supported" = "steelblue"),
                ylab = "Number of Samples")
    }

    # Deepest confident assignment distribution
    if ("deepest_confident_assignment" %in% colnames(assignments)) {
        dca_counts <- sort(table(assignments\$deepest_confident_assignment), decreasing = TRUE)
        if (length(dca_counts) > 0) {
            par(mar = c(5, 10, 3, 1))
            barplot(rev(dca_counts), horiz = TRUE, las = 1,
                    main = "Deepest Confident Assignments",
                    xlab = "Number of Samples", col = "steelblue")
        }
    }

    dev.off()
    cat("Report generated\\n")
    """
}
