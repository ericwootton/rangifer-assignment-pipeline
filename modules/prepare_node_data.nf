/*
========================================================================================
    PREPARE NODE DATA
========================================================================================
    For each node in the hierarchy, prepare a subset of the genotype matrix containing
    only the reference samples for that node, with node-specific SNP filtering:
    1. Subset to relevant reference samples
    2. Recalculate call rates within subset, remove SNPs < 80%
    3. Remove monomorphic SNPs (MAF=0 only, NO higher MAF cutoff)
*/

process PREPARE_NODE_DATA {
    tag "prepare_nodes"
    label 'process_medium'

    publishDir "${params.outdir}/classification/node_data", mode: 'copy'

    input:
    path genotype_matrix
    path metadata
    path hierarchy_json

    output:
    path "nodes/*", emit: node_data

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(jsonlite)

    geno <- fread("${genotype_matrix}", header = TRUE)
    meta <- fread("${metadata}", header = TRUE)
    hierarchy <- fromJSON("${hierarchy_json}")

    # Normalize metadata column names
    if ("Sample" %in% colnames(meta) && !"sample_id" %in% colnames(meta)) setnames(meta, "Sample", "sample_id")
    if ("Subspecies" %in% colnames(meta) && !"subspecies" %in% colnames(meta)) setnames(meta, "Subspecies", "subspecies")
    if ("Ecotype_Cleaned" %in% colnames(meta) && !"ecotype" %in% colnames(meta)) setnames(meta, "Ecotype_Cleaned", "ecotype")
    if ("Herd_Cleaned" %in% colnames(meta) && !"herd" %in% colnames(meta)) setnames(meta, "Herd_Cleaned", "herd")

    # Auto-detect ID column (fallback for other formats)
    id_col <- intersect(c("sample_id", "ID", "Sample_ID"), colnames(meta))[1]
    if (!is.na(id_col) && id_col != "sample_id") setnames(meta, id_col, "sample_id")

    dir.create("nodes", showWarnings = FALSE)

    # Process each node in the hierarchy
    process_node <- function(node_info, node_key) {
        level <- node_info\$level
        column <- node_info\$column
        groups <- node_info\$groups
        node_name <- node_info\$node_name

        cat("\\nProcessing node:", node_name, "(Level", level, ")\\n")
        cat("  Groups:", paste(groups, collapse = ", "), "\\n")

        # Get reference samples for this node
        ref_ids <- meta[get(column) %in% groups]\$sample_id
        ref_ids <- intersect(ref_ids, geno\$sample_id)

        # Get unknown samples (those without labels at this level)
        all_ids <- geno\$sample_id
        unknown_ids <- setdiff(all_ids, meta[!is.na(get(column))]\$sample_id)

        cat("  Reference samples:", length(ref_ids), "\\n")
        cat("  Unknown samples:", length(unknown_ids), "\\n")

        if (length(ref_ids) < 2) {
            cat("  SKIPPING: insufficient reference samples\\n")
            return(invisible(NULL))
        }

        # Subset genotype matrix
        ref_geno <- geno[sample_id %in% ref_ids]
        ref_mat <- as.matrix(ref_geno[, -1])

        # Node-specific SNP QC
        # 1. Call rate within node >= 80%
        snp_call_rate <- colSums(!is.na(ref_mat)) / nrow(ref_mat)
        pass_cr <- snp_call_rate >= ${params.node_snp_call_rate}

        # 2. Remove monomorphic (MAF = 0 only, NOT a conventional MAF filter)
        allele_freq <- colMeans(ref_mat, na.rm = TRUE) / 2
        is_mono <- allele_freq == 0 | allele_freq == 1
        is_mono[is.na(is_mono)] <- TRUE

        keep <- pass_cr & !is_mono
        cat("  SNPs after node QC:", sum(keep), "(removed", sum(!keep),
            "- call rate:", sum(!pass_cr), ", monomorphic:", sum(is_mono & pass_cr), ")\\n")

        kept_snps <- colnames(geno)[-1][keep]

        # Prepare node directory
        node_dir <- file.path("nodes", node_name)
        dir.create(node_dir, showWarnings = FALSE, recursive = TRUE)

        # Write reference genotype matrix
        ref_out <- geno[sample_id %in% ref_ids, c("sample_id", kept_snps), with = FALSE]
        fwrite(ref_out, file.path(node_dir, "reference_genotypes.tsv"), sep = "\\t")

        # Write unknown genotype matrix
        if (length(unknown_ids) > 0) {
            unk_out <- geno[sample_id %in% unknown_ids, c("sample_id", kept_snps), with = FALSE]
            fwrite(unk_out, file.path(node_dir, "unknown_genotypes.tsv"), sep = "\\t")
        }

        # Write reference labels
        ref_labels <- meta[sample_id %in% ref_ids, c("sample_id", column), with = FALSE]
        setnames(ref_labels, column, "group")
        fwrite(ref_labels, file.path(node_dir, "reference_labels.tsv"), sep = "\\t")

        # Write node metadata JSON
        min_class_size <- min(table(ref_labels\$group))
        node_meta <- list(
            node_name = node_name,
            level = level,
            column = column,
            groups = groups,
            n_reference = length(ref_ids),
            n_unknown = length(unknown_ids),
            n_snps = length(kept_snps),
            min_class_size = as.integer(min_class_size),
            sample_counts = as.list(table(ref_labels\$group))
        )
        write_json(node_meta, file.path(node_dir, "node_meta.json"), pretty = TRUE, auto_unbox = TRUE)

        cat("  Node data written to:", node_dir, "\\n")
    }

    # Process Level 1
    if (!is.null(hierarchy\$level1)) {
        process_node(hierarchy\$level1, "level1")
    }

    # Process Level 2 nodes
    if (!is.null(hierarchy\$level2)) {
        for (key in names(hierarchy\$level2)) {
            process_node(hierarchy\$level2[[key]], paste0("level2_", key))
        }
    }

    # Process Level 3 nodes
    if (!is.null(hierarchy\$level3)) {
        for (key in names(hierarchy\$level3)) {
            process_node(hierarchy\$level3[[key]], paste0("level3_", key))
        }
    }

    cat("\\nAll nodes prepared\\n")
    """
}
