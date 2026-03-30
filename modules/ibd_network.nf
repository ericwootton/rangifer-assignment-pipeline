/*
========================================================================================
    IBD NETWORK (WGS mode, per node)
========================================================================================
    Uses IBDseq to detect IBD segments from phased VCF, then builds an IBD-sharing
    network adjacency matrix.
*/

process IBD_NETWORK {
    tag "${node_name}"
    label 'process_medium'

    publishDir "${params.outdir}/wgs_analyses/ibd/${node_name}", mode: 'copy'

    input:
    tuple val(node_name), path(node_dir), path(node_meta)
    tuple path(phased_vcf), path(phased_tbi)

    output:
    tuple val(node_name), path("${node_name}_ibd_matrix.tsv"), path("${node_name}_ibd_segments.tsv"), emit: results

    script:
    """
    set -euo pipefail

    # Extract node samples
    awk -F'\\t' 'NR>1 {print \$1}' ${node_dir}/reference_labels.tsv > node_samples.txt
    N_SAMPLES=\$(wc -l < node_samples.txt)

    if [ "\${N_SAMPLES}" -lt 2 ]; then
        # Not enough samples for IBD
        echo "sample_i\\tsample_j\\ttotal_ibd_length" > ${node_name}_ibd_matrix.tsv
        echo "sample1\\tsample2\\tchrom\\tstart\\tend\\tLOD" > ${node_name}_ibd_segments.tsv
        echo "Skipped: fewer than 2 samples"
        exit 0
    fi

    # Subset phased VCF to node samples
    bcftools view -S node_samples.txt ${phased_vcf} -Oz -o node_phased.vcf.gz
    bcftools index -t node_phased.vcf.gz

    # Run IBDseq
    java -jar \$(which ibdseq.jar 2>/dev/null || echo /opt/ibdseq/ibdseq.jar) \\
        LOD=3.0 \\
        ibdtrim=0 \\
        in=node_phased.vcf.gz \\
        out=ibd_result || true

    # Build IBD adjacency matrix
    if [ -f ibd_result.ibd ]; then
        cp ibd_result.ibd ${node_name}_ibd_segments.tsv
    else
        echo "sample1\\tsample2\\tchrom\\tstart\\tend\\tLOD" > ${node_name}_ibd_segments.tsv
    fi

    # Compute pairwise total IBD length
    cat > build_matrix.R <<'RSCRIPT'
    library(data.table)

    segments_file <- "${node_name}_ibd_segments.tsv"
    samples_file <- "node_samples.txt"

    samples <- readLines(samples_file)
    n <- length(samples)

    if (file.info(segments_file)\$size > 100) {
        segs <- tryCatch(
            fread(segments_file, header = FALSE),
            error = function(e) data.table()
        )
        if (nrow(segs) > 0 && ncol(segs) >= 5) {
            colnames(segs)[1:5] <- c("sample1", "sample2", "chrom", "start", "end")
            segs[, length := end - start]

            # Sum IBD length per pair
            pair_ibd <- segs[, .(total_ibd = sum(length)), by = .(sample1, sample2)]

            # Build matrix
            mat <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))
            for (i in seq_len(nrow(pair_ibd))) {
                s1 <- pair_ibd[["sample1"]][i]
                s2 <- pair_ibd[["sample2"]][i]
                if (s1 %in% samples && s2 %in% samples) {
                    mat[s1, s2] <- pair_ibd[["total_ibd"]][i]
                    mat[s2, s1] <- pair_ibd[["total_ibd"]][i]
                }
            }
        } else {
            mat <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))
        }
    } else {
        mat <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))
    }

    out <- data.table(sample_id = rownames(mat), as.data.table(mat))
    fwrite(out, "${node_name}_ibd_matrix.tsv", sep = "\\t")
RSCRIPT

    Rscript build_matrix.R
    """
}
