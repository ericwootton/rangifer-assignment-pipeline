/*
========================================================================================
    CONVERT ILLUMINA FINALREPORT TO ALLELE MATRIX
========================================================================================
    Parses Illumina FinalReport CSV and outputs actual A/C/G/T allele pairs per sample
    per SNP. Column names are converted to CHROM:POS using the chip_snp_map.

    Does NOT code to 012 here - that happens in the merge step where we can match
    alleles to VCF REF/ALT accounting for strand differences.

    Output format: TSV with sample_id as rows, CHROM:POS as columns.
    Cell values are two-character allele pairs: "AG", "CC", "TT", or NA for missing.
*/

process CONVERT_FINALREPORT {
    tag "convert_finalreport"
    label 'process_medium'

    publishDir "${params.outdir}/data_integration/finalreport", mode: 'copy'

    input:
    path finalreport
    path samplesheet
    path snp_map

    output:
    path "finalreport_alleles.tsv", emit: allele_matrix
    path "finalreport_stats.txt",   emit: stats

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    cat("=== FinalReport Conversion ===\\n\\n")

    # Load SNP name -> CHROM:POS map
    snp_map <- fread("${snp_map}", header = FALSE, col.names = c("snp_name", "chrom", "pos", "snp_id"))
    snp_lookup <- unique(snp_map[, .(snp_name, snp_id)], by = "snp_name")
    cat("SNP map entries:", nrow(snp_lookup), "\\n")

    # Read FinalReport, skip header lines until [Data]
    lines <- readLines("${finalreport}", n = 20)
    skip_n <- grep("^\\\\[Data\\\\]", lines)
    if (length(skip_n) == 0) stop("Cannot find [Data] section in FinalReport")

    fr <- fread("${finalreport}", skip = skip_n, header = TRUE)
    cat("FinalReport rows:", nrow(fr), "\\n")

    # Rename columns with spaces/special chars to avoid backtick issues
    if ("Sample ID" %in% colnames(fr)) setnames(fr, "Sample ID", "SampleID")
    if ("SNP Name" %in% colnames(fr)) setnames(fr, "SNP Name", "SNPName")
    if ("Allele1 - Top" %in% colnames(fr)) setnames(fr, "Allele1 - Top", "Allele1Top")
    if ("Allele2 - Top" %in% colnames(fr)) setnames(fr, "Allele2 - Top", "Allele2Top")

    # Extract clean sample ID (second field of compound ID)
    fr[, clean_id := sapply(strsplit(as.character(SampleID), "#"), function(x) {
        if (length(x) >= 2) x[2] else x[1]
    })]

    # Map SNP names to CHROM:POS
    fr <- merge(fr, snp_lookup, by.x = "SNPName", by.y = "snp_name", all.x = TRUE)
    cat("Mapped:", sum(!is.na(fr\$snp_id)), "Unmapped:", sum(is.na(fr\$snp_id)), "\\n")
    fr <- fr[!is.na(snp_id)]

    # Encode allele pair as two-character string; "-" means missing
    fr[, allele_pair := ifelse(
        Allele1Top == "-" | Allele2Top == "-",
        NA_character_,
        paste0(Allele1Top, Allele2Top)
    )]

    # Deduplicate: multiple probe designs (e.g. _modl variants) can map to same position
    fr <- unique(fr, by = c("clean_id", "snp_id"))

    # Pivot to wide: samples as rows, CHROM:POS as columns, allele pairs as values
    geno_wide <- dcast(fr, clean_id ~ snp_id, value.var = "allele_pair")
    setnames(geno_wide, "clean_id", "sample_id")

    cat("Output:", nrow(geno_wide), "samples x", ncol(geno_wide) - 1, "SNPs\\n")

    fwrite(geno_wide, "finalreport_alleles.tsv", sep = "\\t", na = "NA", quote = FALSE)

    writeLines(c(
        "FinalReport Conversion Stats",
        paste("Samples:", nrow(geno_wide)),
        paste("SNPs mapped to coordinates:", ncol(geno_wide) - 1),
        "Alleles reported on Illumina TOP strand"
    ), "finalreport_stats.txt")
    """
}
