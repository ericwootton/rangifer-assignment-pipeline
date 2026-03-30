/*
========================================================================================
    FEATURE IMPORTANCE ANALYSIS
========================================================================================
    Extract and compare informative SNPs across methods and nodes.
*/

process FEATURE_IMPORTANCE {
    tag "feature_importance"
    label 'process_medium'

    publishDir "${params.outdir}/feature_importance", mode: 'copy'

    input:
    path dapc_results
    path assignpop_results

    output:
    path "importance_summary.tsv",    emit: importance_summary
    path "cross_method_comparison.tsv", emit: comparison
    path "importance_plots.pdf",      emit: plots

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    # Collect all SNP importance files
    # list.files pattern matches filenames only, so find all snp_importance.tsv
    # then determine method from parent directory name
    all_imp_files <- list.files(".", pattern = "^snp_importance\\\\.tsv\$",
                                recursive = TRUE, full.names = TRUE)
    cat("Found", length(all_imp_files), "importance files\\n")

    all_imp <- data.table()

    for (f in all_imp_files) {
        parent_dir <- basename(dirname(f))
        if (grepl("_dapc_results", parent_dir)) {
            method <- "DAPC"
            node <- gsub("_dapc_results.*", "", parent_dir)
        } else if (grepl("_assignpop_results", parent_dir)) {
            method <- "assignPOP"
            node <- gsub("_assignpop_results.*", "", parent_dir)
        } else {
            next
        }
        dt <- fread(f, header = TRUE)
        if (nrow(dt) == 0) next
        dt[, c("method", "node") := list(method, node)]
        all_imp <- rbind(all_imp, dt, fill = TRUE)
    }

    if (nrow(all_imp) == 0) {
        cat("No importance data found\\n")
        fwrite(data.table(message = "No importance data"), "importance_summary.tsv", sep = "\\t")
        fwrite(data.table(message = "No data"), "cross_method_comparison.tsv", sep = "\\t")
        pdf("importance_plots.pdf"); plot.new(); text(0.5, 0.5, "No data"); dev.off()
        quit(save = "no")
    }

    fwrite(all_imp, "importance_summary.tsv", sep = "\\t")

    # Cross-method comparison: find SNPs ranked highly by both methods
    snp_col <- intersect(c("snp", "snp_id", "SNP"), colnames(all_imp))[1]
    if (is.na(snp_col)) snp_col <- colnames(all_imp)[1]

    comparison <- all_imp[, .(
        n_methods = uniqueN(method),
        methods = paste(unique(method), collapse = ","),
        nodes = paste(unique(node), collapse = ",")
    ), by = snp_col]
    setorder(comparison, -n_methods)

    fwrite(comparison, "cross_method_comparison.tsv", sep = "\\t")

    # Plots
    pdf("importance_plots.pdf", width = 12, height = 8)
    for (n in unique(all_imp\$node)) {
        node_imp <- all_imp[node == n]
        if (nrow(node_imp) > 0) {
            for (m in unique(node_imp\$method)) {
                method_imp <- node_imp[method == m]
                imp_col <- intersect(c("total_loading", "importance"), colnames(method_imp))[1]
                if (!is.na(imp_col) && nrow(method_imp) >= 5) {
                    top <- head(method_imp[order(-get(imp_col))], 20)
                    par(mar = c(5, 12, 3, 1))
                    barplot(rev(top[[imp_col]]), names.arg = rev(top[[snp_col]]),
                            horiz = TRUE, las = 1, cex.names = 0.6,
                            main = paste(n, "-", m, "top 20 SNPs"),
                            xlab = imp_col, col = "steelblue")
                }
            }
        }
    }
    dev.off()

    cat("Feature importance analysis complete\\n")
    """
}
