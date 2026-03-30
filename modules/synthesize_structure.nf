/*
========================================================================================
    SYNTHESIZE STRUCTURE DISCOVERY
========================================================================================
    Combines PCA, DAPC clusters, and fastSTRUCTURE results to define the
    classification hierarchy. Outputs a JSON hierarchy definition.
*/

process SYNTHESIZE_STRUCTURE {
    tag "synthesize"
    label 'process_low'

    publishDir "${params.outdir}/structure_discovery", mode: 'copy'

    input:
    path pca_results
    path cluster_results
    path structure_results
    path metadata

    output:
    path "hierarchy_definition.json", emit: hierarchy_definition
    path "synthesis_report.txt",      emit: report

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(jsonlite)

    sink("synthesis_report.txt")
    cat("=== Structure Synthesis Report ===\\n\\n")

    meta <- fread("${metadata}", header = TRUE)

    # Normalize metadata column names
    if ("Sample" %in% colnames(meta) && !"sample_id" %in% colnames(meta)) setnames(meta, "Sample", "sample_id")
    if ("Subspecies" %in% colnames(meta) && !"subspecies" %in% colnames(meta)) setnames(meta, "Subspecies", "subspecies")
    if ("Ecotype_Cleaned" %in% colnames(meta) && !"ecotype" %in% colnames(meta)) setnames(meta, "Ecotype_Cleaned", "ecotype")
    if ("Herd_Cleaned" %in% colnames(meta) && !"herd" %in% colnames(meta)) setnames(meta, "Herd_Cleaned", "herd")

    # Auto-detect columns
    id_col <- intersect(c("sample_id", "ID", "Sample_ID"), colnames(meta))[1]
    if (!is.na(id_col) && id_col != "sample_id") setnames(meta, id_col, "sample_id")

    # Find hierarchy columns in metadata
    level_cols <- list()
    for (col_name in c("lineage", "subspecies", "major_lineage", "level1")) {
        if (col_name %in% colnames(meta)) { level_cols[["level1"]] <- col_name; break }
    }
    for (col_name in c("ecotype", "subspecies_ecotype", "level2")) {
        if (col_name %in% colnames(meta)) { level_cols[["level2"]] <- col_name; break }
    }
    for (col_name in c("herd", "population", "Source_ID", "level3")) {
        if (col_name %in% colnames(meta)) { level_cols[["level3"]] <- col_name; break }
    }

    cat("Detected hierarchy columns:\\n")
    for (lvl in names(level_cols)) cat("  ", lvl, ":", level_cols[[lvl]], "\\n")

    # Build hierarchy definition
    # For each level, list the groups and sample counts
    hierarchy <- list()

    # Reference samples = those with metadata at each level
    ref_samples <- meta[!is.na(meta[[level_cols[["level1"]]]])]
    cat("\\nReference samples with Level 1 labels:", nrow(ref_samples), "\\n")

    # Level 1 node (root)
    if (!is.null(level_cols[["level1"]])) {
        l1_groups <- ref_samples[, .N, by = c(level_cols[["level1"]])]
        setnames(l1_groups, 1, "group")
        cat("\\nLevel 1 groups:\\n")
        print(l1_groups)

        hierarchy[["level1"]] <- list(
            node_name = "root",  # no spaces to sanitize
            level = 1,
            column = level_cols[["level1"]],
            groups = l1_groups\$group,
            sample_counts = setNames(l1_groups\$N, l1_groups\$group)
        )
    }

    # Level 2 nodes (one per Level 1 group)
    if (!is.null(level_cols[["level2"]])) {
        hierarchy[["level2"]] <- list()
        for (l1_grp in hierarchy[["level1"]]\$groups) {
            l1_col <- level_cols[["level1"]]
            l2_col <- level_cols[["level2"]]
            subset <- ref_samples[get(l1_col) == l1_grp & !is.na(get(l2_col))]
            if (nrow(subset) > 0) {
                l2_groups <- subset[, .N, by = c(l2_col)]
                setnames(l2_groups, 1, "group")
                cat("\\nLevel 2 within", l1_grp, ":\\n")
                print(l2_groups)

                hierarchy[["level2"]][[l1_grp]] <- list(
                    node_name = gsub(" ", "_", paste0("L2_", l1_grp)),
                    level = 2,
                    parent_group = l1_grp,
                    column = l2_col,
                    groups = l2_groups\$group,
                    sample_counts = setNames(l2_groups\$N, l2_groups\$group)
                )
            }
        }
    }

    # Level 3 nodes (one per Level 2 group)
    if (!is.null(level_cols[["level3"]])) {
        hierarchy[["level3"]] <- list()
        for (l1_grp in names(hierarchy[["level2"]])) {
            l2_node <- hierarchy[["level2"]][[l1_grp]]
            for (l2_grp in l2_node\$groups) {
                l1_col <- level_cols[["level1"]]
                l2_col <- level_cols[["level2"]]
                l3_col <- level_cols[["level3"]]
                subset <- ref_samples[get(l1_col) == l1_grp & get(l2_col) == l2_grp & !is.na(get(l3_col))]
                if (nrow(subset) > 0) {
                    l3_groups <- subset[, .N, by = c(l3_col)]
                    setnames(l3_groups, 1, "group")
                    cat("\\nLevel 3 within", l1_grp, "/", l2_grp, ":\\n")
                    print(l3_groups)

                    node_key <- gsub(" ", "_", paste0(l1_grp, "_", l2_grp))
                    hierarchy[["level3"]][[node_key]] <- list(
                        node_name = gsub(" ", "_", paste0("L3_", l1_grp, "_", l2_grp)),
                        level = 3,
                        parent_group = l2_grp,
                        grandparent_group = l1_grp,
                        column = l3_col,
                        groups = l3_groups\$group,
                        sample_counts = setNames(l3_groups\$N, l3_groups\$group)
                    )
                }
            }
        }
    }

    sink()

    # Write hierarchy definition
    write_json(hierarchy, "hierarchy_definition.json", pretty = TRUE, auto_unbox = TRUE)
    cat("Hierarchy definition written\\n")
    """
}
