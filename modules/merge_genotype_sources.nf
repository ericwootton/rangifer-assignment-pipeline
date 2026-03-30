/*
========================================================================================
    MERGE GENOTYPE SOURCES
========================================================================================
    Combines the extracted WGS genotypes (TSV) with chip allele matrices into a
    unified 012-coded genotype matrix.

    All sources use CHROM:POS column names (converters mapped chip names to coords).
    Chip converters output actual allele pairs (e.g., "AG", "CC", "TT") on Illumina
    TOP strand, which may differ from the VCF reference forward strand.

    For each shared SNP, chip alleles are matched to VCF REF/ALT by:
      1. Direct match: chip alleles ⊂ {REF, ALT}
      2. Complement match: complement(chip alleles) ⊂ {REF, ALT}
         (A↔T, C↔G for Illumina TOP → forward strand conversion)
      3. If neither matches, the SNP is flagged and excluded

    012 coding: 0 = homozygous REF, 1 = heterozygous, 2 = homozygous ALT
*/

process MERGE_GENOTYPE_SOURCES {
    tag "merge_sources"
    label 'process_high'

    publishDir "${params.outdir}/data_integration", mode: 'copy'

    input:
    path wgs_genotypes
    path wgs_snp_info
    path chip_matrices

    output:
    path "merged_genotype_matrix.tsv", emit: genotype_matrix
    path "snp_info.tsv",              emit: snp_info
    path "data_source.tsv",           emit: data_source
    path "merge_stats.txt",           emit: stats

    script:
    def chip_files = chip_matrices instanceof List ? chip_matrices.join(',') : chip_matrices
    """
    #!/usr/bin/env Rscript

    library(data.table)

    sink("merge_stats.txt")
    cat("=== Genotype Source Merge Report ===\\n\\n")

    # -----------------------------------------------
    # 1. Load WGS genotypes (already extracted from VCF)
    # -----------------------------------------------
    gt_raw <- fread("${wgs_genotypes}", header = TRUE)
    n_snps <- nrow(gt_raw)

    snp_info <- gt_raw[, .(chrom = CHROM, pos = POS, ref = REF, alt = ALT)]
    snp_info[, snp_id := paste0(chrom, ":", pos)]

    sample_cols <- colnames(gt_raw)[5:ncol(gt_raw)]
    cat("WGS samples:", length(sample_cols), "\\n")
    cat("WGS SNPs:", n_snps, "\\n")

    # Convert GT strings to 012
    wgs_geno <- matrix(NA_integer_, nrow = n_snps, ncol = length(sample_cols))
    for (j in seq_along(sample_cols)) {
        gt <- gt_raw[[sample_cols[j]]]
        wgs_geno[, j] <- ifelse(gt == "0/0" | gt == "0|0", 0L,
            ifelse(gt == "0/1" | gt == "1/0" | gt == "0|1" | gt == "1|0", 1L,
                ifelse(gt == "1/1" | gt == "1|1", 2L, NA_integer_)))
    }

    # Transpose: samples as rows, SNPs as columns
    wgs_dt <- data.table(sample_id = sample_cols, t(wgs_geno))
    setnames(wgs_dt, c("sample_id", snp_info\$snp_id))
    cat("WGS matrix:", nrow(wgs_dt), "samples x", ncol(wgs_dt) - 1, "SNPs\\n")

    # Build REF/ALT lookup keyed by snp_id
    ref_lookup <- setNames(snp_info\$ref, snp_info\$snp_id)
    alt_lookup <- setNames(snp_info\$alt, snp_info\$snp_id)

    # -----------------------------------------------
    # 2. Load chip allele matrices (if any)
    # -----------------------------------------------
    chip_file_list <- strsplit("${chip_files}", ",")[[1]]
    chip_file_list <- chip_file_list[chip_file_list != "" & !grepl("NO_CHIP_DATA", chip_file_list)]

    # Complement function for Illumina TOP → forward strand
    complement <- function(base) {
        comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
        comp[base]
    }

    if (length(chip_file_list) > 0) {
        chip_dts <- list()
        for (f in chip_file_list) {
            if (file.exists(f) && file.size(f) > 0) {
                dt <- fread(f, header = TRUE)
                cat("Chip source", basename(f), ":", nrow(dt), "samples x", ncol(dt) - 1, "SNPs\\n")
                chip_dts[[f]] <- dt
            }
        }

        if (length(chip_dts) > 0) {
            # Find shared SNPs between chip sources and VCF
            vcf_snps <- snp_info\$snp_id

            chip_combined_list <- list()
            for (f in names(chip_dts)) {
                dt <- chip_dts[[f]]
                chip_snps <- colnames(dt)[-1]
                shared <- intersect(chip_snps, vcf_snps)
                cat("  Shared with VCF:", length(shared), "of", length(chip_snps), "\\n")
                chip_combined_list[[f]] <- dt[, c("sample_id", shared), with = FALSE]
            }

            chip_combined <- rbindlist(chip_combined_list, fill = TRUE)

            # Ensure all SNP columns are character and strip any embedded quotes
            snp_cols_chip <- colnames(chip_combined)[-1]
            for (col in snp_cols_chip) {
                v <- as.character(chip_combined[[col]])
                v <- gsub('"', '', v, fixed = TRUE)
                set(chip_combined, j = col, value = v)
            }

            # Remove chip samples that duplicate WGS sample IDs
            overlap <- intersect(chip_combined\$sample_id, wgs_dt\$sample_id)
            if (length(overlap) > 0) {
                cat("Overlapping sample IDs (chip removed):", length(overlap), "\\n")
                chip_combined <- chip_combined[!sample_id %in% overlap]
            }
            cat("Chip samples after dedup:", nrow(chip_combined), "\\n")

            # -----------------------------------------------
            # 3. Convert chip allele pairs to 012 using
            #    complement-aware base matching against VCF REF/ALT
            # -----------------------------------------------
            shared_snps <- intersect(colnames(chip_combined)[-1], vcf_snps)
            cat("\\nConverting", length(shared_snps), "shared SNPs from allele pairs to 012...\\n")

            n_direct <- 0L
            n_complement <- 0L
            n_ambiguous <- 0L
            n_failed <- 0L
            n_all_na <- 0L
            failed_snps <- character(0)
            ambiguous_snps <- character(0)

            for (snp in shared_snps) {
                ref <- ref_lookup[snp]
                alt <- alt_lookup[snp]
                allele_pairs <- chip_combined[[snp]]

                # Extract individual alleles from two-character pairs
                a1 <- substr(allele_pairs, 1, 1)
                a2 <- substr(allele_pairs, 2, 2)

                # Skip SNPs with no data (set to integer NA for type safety)
                non_na <- which(!is.na(allele_pairs))
                if (length(non_na) == 0) {
                    n_all_na <- n_all_na + 1L
                    set(chip_combined, j = snp, value = NA_integer_)
                    next
                }

                # Check if this is an ambiguous A/T or C/G SNP
                # For these, complement maps to the same allele set,
                # so strand cannot be determined from base identity alone
                is_ambiguous_snp <- (ref == "A" & alt == "T") |
                                    (ref == "T" & alt == "A") |
                                    (ref == "C" & alt == "G") |
                                    (ref == "G" & alt == "C")

                # Gather observed alleles from multiple samples for robust detection
                sample_check <- non_na[1:min(20, length(non_na))]
                all_alleles <- unique(c(a1[sample_check], a2[sample_check]))
                all_alleles <- all_alleles[!is.na(all_alleles)]

                use_complement <- FALSE
                if (all(all_alleles %in% c(ref, alt))) {
                    # Direct match
                    n_direct <- n_direct + 1L
                    if (is_ambiguous_snp) {
                        n_ambiguous <- n_ambiguous + 1L
                        ambiguous_snps <- c(ambiguous_snps, snp)
                    }
                } else if (all(complement(all_alleles) %in% c(ref, alt))) {
                    # Complement match
                    use_complement <- TRUE
                    n_complement <- n_complement + 1L
                } else {
                    n_failed <- n_failed + 1L
                    failed_snps <- c(failed_snps, snp)
                    set(chip_combined, j = snp, value = NA_integer_)
                    next
                }

                # Apply complement if needed
                if (use_complement) {
                    a1 <- complement(a1)
                    a2 <- complement(a2)
                }

                # Code to 012: count ALT alleles
                alt_count <- ifelse(is.na(a1) | is.na(a2), NA_integer_,
                    (a1 == alt) + (a2 == alt))

                # Validate: alleles should be ref or alt
                invalid <- !is.na(a1) & !is.na(a2) & !(a1 %in% c(ref, alt) & a2 %in% c(ref, alt))
                alt_count[invalid] <- NA_integer_

                set(chip_combined, j = snp, value = as.integer(alt_count))
            }

            cat("Direct allele match:", n_direct, "SNPs\\n")
            cat("Complement match:", n_complement, "SNPs\\n")
            cat("A/T or C/G ambiguous SNPs (strand indeterminate):", n_ambiguous, "\\n")
            cat("All-NA SNPs:", n_all_na, "\\n")
            cat("Failed to match:", n_failed, "SNPs\\n")
            if (n_failed > 0 && n_failed <= 20) {
                cat("  Failed SNPs:", paste(failed_snps, collapse = ", "), "\\n")
            }

            # Remove failed SNPs from chip data
            if (n_failed > 0) {
                keep_snps <- setdiff(shared_snps, failed_snps)
                chip_combined <- chip_combined[, c("sample_id", keep_snps), with = FALSE]
            }

            # -----------------------------------------------
            # 4. Merge WGS + chip
            # -----------------------------------------------
            wgs_snp_cols <- colnames(wgs_dt)[-1]
            chip_aligned <- chip_combined[, c("sample_id", intersect(colnames(chip_combined)[-1], wgs_snp_cols)), with = FALSE]

            merged <- rbindlist(list(wgs_dt, chip_aligned), fill = TRUE)
        } else {
            merged <- wgs_dt
        }
    } else {
        merged <- wgs_dt
    }

    cat("\\n=== Final Merged Matrix ===\\n")
    cat("Samples:", nrow(merged), "\\n")
    cat("SNPs:", ncol(merged) - 1, "\\n")

    sink()

    # Build data_source table mapping each sample to its genotyping platform
    wgs_ids <- sample_cols
    ds_list <- list(data.table(sample_id = wgs_ids, source = "WGS"))

    if (exists("chip_combined") && nrow(chip_combined) > 0) {
        # Determine which chip source each sample came from
        for (f in names(chip_dts)) {
            chip_ids <- intersect(chip_dts[[f]]\$sample_id, merged\$sample_id)
            src_label <- ifelse(grepl("nwt", basename(f), ignore.case = TRUE), "NWT", "FinalReport")
            ds_list[[length(ds_list) + 1]] <- data.table(sample_id = chip_ids, source = src_label)
        }
    }
    data_source <- rbindlist(ds_list)
    # Deduplicate (WGS takes priority if overlap)
    data_source <- data_source[!duplicated(sample_id)]
    data_source <- data_source[sample_id %in% merged\$sample_id]

    fwrite(merged, "merged_genotype_matrix.tsv", sep = "\\t")
    fwrite(snp_info, "snp_info.tsv", sep = "\\t")
    fwrite(data_source, "data_source.tsv", sep = "\\t")
    """
}
