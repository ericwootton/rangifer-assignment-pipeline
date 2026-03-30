/*
========================================================================================
    CONVERT NWT GENOTYPE MATRIX TO ALLELE MATRIX
========================================================================================
    The NWT data uses two-letter allele codes (CC, AG, TT, NN=missing).
    These are already on some strand - likely the same Illumina TOP strand as
    the FinalReport since they come from the same chip platform.

    Column names are mapped to CHROM:POS coordinates using the chip_snp_map.
    Allele pairs are passed through as-is for strand resolution in the merge step.
*/

process CONVERT_NWT_GENOTYPES {
    tag "convert_nwt"
    label 'process_low'

    publishDir "${params.outdir}/data_integration/nwt", mode: 'copy'

    input:
    path nwt_csv
    path snp_map

    output:
    path "nwt_alleles.tsv", emit: allele_matrix
    path "nwt_stats.txt",   emit: stats

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    cat("=== NWT Genotype Conversion ===\\n\\n")

    # Load SNP name -> CHROM:POS map
    snp_map <- fread("${snp_map}", header = FALSE, col.names = c("snp_name", "chrom", "pos", "snp_id"))
    var_lookup <- unique(snp_map[grepl("^Var", snp_name)], by = "snp_name")
    setkey(var_lookup, snp_name)
    cat("Var-name map entries:", nrow(var_lookup), "\\n")

    nwt <- fread("${nwt_csv}", header = TRUE)
    cat("NWT raw:", nrow(nwt), "rows x", ncol(nwt), "columns\\n")

    setnames(nwt, 1, "sample_id")
    nwt <- nwt[!duplicated(sample_id)]
    cat("After dedup:", nrow(nwt), "samples\\n")

    snp_cols <- colnames(nwt)[-1]

    # Map column names to CHROM:POS
    col_mapping <- var_lookup[snp_cols, on = "snp_name", nomatch = NA]
    mapped_mask <- !is.na(col_mapping\$snp_id)
    cat("Mapped:", sum(mapped_mask), "Unmapped:", sum(!mapped_mask), "\\n")

    mapped_snps <- snp_cols[mapped_mask]
    new_names <- col_mapping\$snp_id[mapped_mask]

    # Subset and rename
    result <- nwt[, c("sample_id", mapped_snps), with = FALSE]
    setnames(result, c("sample_id", mapped_snps), c("sample_id", new_names))

    # Convert NN to NA, keep all other allele pairs as-is
    for (col in new_names) {
        vals <- result[[col]]
        vals[vals == "NN" | vals == "" | nchar(vals) < 2] <- NA_character_
        set(result, j = col, value = vals)
    }

    cat("Output:", nrow(result), "samples x", ncol(result) - 1, "SNPs\\n")

    fwrite(result, "nwt_alleles.tsv", sep = "\\t", na = "NA", quote = FALSE)

    writeLines(c(
        "NWT Conversion Stats",
        paste("Samples:", nrow(result)),
        paste("SNPs mapped to coordinates:", ncol(result) - 1),
        "Alleles on chip reporting strand (likely Illumina TOP)"
    ), "nwt_stats.txt")
    """
}
