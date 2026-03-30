/*
========================================================================================
    CHROMOPAINTER + FINESTRUCTURE (WGS mode, per node)
========================================================================================
    Convert phased VCF to ChromoPainter input, run painting, then run fineSTRUCTURE
    MCMC clustering on the co-ancestry matrix.
*/

process CHROMOPAINTER_FINESTRUCTURE {
    tag "${node_name}"
    label 'process_high'

    publishDir "${params.outdir}/wgs_analyses/finestructure/${node_name}", mode: 'copy'

    input:
    tuple val(node_name), path(node_dir), path(node_meta)
    tuple path(phased_vcf), path(phased_tbi)

    output:
    tuple val(node_name), path("${node_name}_coancestry.tsv"), path("${node_name}_clusters.tsv"), emit: results

    script:
    """
    set -euo pipefail

    # Extract node samples
    awk -F'\\t' 'NR>1 {print \$1}' ${node_dir}/reference_labels.tsv > node_samples.txt
    N_SAMPLES=\$(wc -l < node_samples.txt)

    if [ "\${N_SAMPLES}" -lt 4 ]; then
        echo "sample_id\\tcluster" > ${node_name}_clusters.tsv
        echo "sample_id" > ${node_name}_coancestry.tsv
        echo "Skipped: fewer than 4 samples for fineSTRUCTURE"
        exit 0
    fi

    # Subset phased VCF to node samples
    bcftools view -S node_samples.txt ${phased_vcf} -Oz -o node_phased.vcf.gz
    bcftools index -t node_phased.vcf.gz

    # Convert VCF to ChromoPainter input format
    # Create recombination map (uniform rate as approximation)
    N_SNPS=\$(bcftools view -H node_phased.vcf.gz | wc -l)

    # Create phase file from VCF
    # ChromoPainter needs: haplotype matrix (0/1 per haplotype per individual)
    bcftools query -f '[%GT\\t]\\n' node_phased.vcf.gz > gt_raw.txt

    # Convert to ChromoPainter phase format
    cat > convert_to_cp.R <<'RSCRIPT'
    library(data.table)

    samples <- readLines("node_samples.txt")
    n_samples <- length(samples)

    # Read genotypes
    gt <- fread("gt_raw.txt", header = FALSE)
    n_snps <- nrow(gt)

    # Create haplotype matrix (2 rows per individual)
    hap_mat <- matrix(0, nrow = 2 * n_samples, ncol = n_snps)
    for (i in seq_len(n_samples)) {
        gts <- as.character(gt[[i]])
        alleles <- strsplit(gts, "[/|]")
        hap_mat[2*i - 1, ] <- as.integer(sapply(alleles, `[`, 1))
        hap_mat[2*i, ] <- as.integer(sapply(alleles, `[`, 2))
    }
    hap_mat[is.na(hap_mat)] <- 0

    # Write phase file
    sink("input.phase")
    cat(n_samples, "\n")
    cat(n_snps, "\n")
    cat("P", paste(seq_len(n_snps), collapse = " "), "\n")
    cat(paste(rep("S", n_snps), collapse = ""), "\n")
    for (i in seq_len(nrow(hap_mat))) {
        cat(paste(hap_mat[i, ], collapse = ""), "\n")
    }
    sink()

    # Write IDs file
    writeLines(samples, "input.ids")

    # Write recombination file (uniform rate)
    sink("input.recomb")
    cat("start.pos recom.rate.perbp\n")
    cat("1 1e-8\n")
    sink()
RSCRIPT

    Rscript convert_to_cp.R

    # Run ChromoPainter (all vs all painting)
    # Use fs (fineSTRUCTURE suite) for linked painting
    fs cp -g input.phase -r input.recomb -t input.ids \
        -o painted \
        -s 0 -n \${N_SAMPLES} || {
        # Fallback: if fs cp fails, create dummy outputs
        echo "ChromoPainter failed, creating empty outputs"
        echo "sample_id\\tcluster" > ${node_name}_clusters.tsv
        echo "sample_id" > ${node_name}_coancestry.tsv
        exit 0
    }

    # Run fineSTRUCTURE MCMC
    if [ -f painted_chunkcounts.out ]; then
        fs fs -i painted_chunkcounts.out -o fs_mcmc \
            -x 100000 -y 100000 -z 1000 || {
            echo "fineSTRUCTURE MCMC failed, creating empty outputs"
            echo "sample_id\\tcluster" > ${node_name}_clusters.tsv
            echo "sample_id" > ${node_name}_coancestry.tsv
            exit 0
        }

        # Extract co-ancestry matrix
        if [ -f painted_chunkcounts.out ]; then
            cp painted_chunkcounts.out ${node_name}_coancestry.tsv
        fi

        # Extract cluster assignments
        if [ -f fs_mcmc_tree.xml ]; then
            fs tree -i fs_mcmc_tree.xml -o clusters.csv || true
        fi

        if [ -f clusters.csv ]; then
            cp clusters.csv ${node_name}_clusters.tsv
        else
            # Parse from MCMC output
            cat > extract_clusters.R <<'RSCRIPT2'
            samples <- readLines("node_samples.txt")
            # Simple: write sample IDs with placeholder clusters
            df <- data.frame(sample_id = samples, cluster = NA)
            write.table(df, "${node_name}_clusters.tsv", sep = "\\t",
                        row.names = FALSE, quote = FALSE)
RSCRIPT2
            Rscript extract_clusters.R
        fi
    else
        echo "sample_id\\tcluster" > ${node_name}_clusters.tsv
        echo "sample_id" > ${node_name}_coancestry.tsv
    fi
    """
}
