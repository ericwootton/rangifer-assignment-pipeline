/*
========================================================================================
    NODE-LEVEL ADMIXTURE
========================================================================================
    Runs ADMIXTURE at each hierarchy node to estimate ancestry proportions
    and select optimal K via cross-validation error.

    Outputs Q matrices for each K and the optimal K selection.
*/

process RUN_NODE_ADMIXTURE {
    tag "${node_name}"
    label 'process_medium'

    publishDir "${params.outdir}/classification/admixture_nodes", mode: 'copy'

    input:
    tuple val(node_name), path(bed), path(bim), path(fam)

    output:
    tuple val(node_name), path("${node_name}_admixture"), emit: results

    script:
    """
    #!/usr/bin/env bash
    set -euo pipefail

    mkdir -p ${node_name}_admixture

    N_SAMP=\$(wc -l < ${fam})
    N_SNPS=\$(wc -l < ${bim})
    echo "Node: ${node_name}"
    echo "Samples: \$N_SAMP, SNPs: \$N_SNPS"

    # Determine K range
    # Min K=2, max K = min(n_groups * 2, n_samples - 1, max_k_structure)
    MAX_K=${params.max_k_structure}
    if [ \$MAX_K -ge \$N_SAMP ]; then
        MAX_K=\$((N_SAMP - 1))
    fi
    # Need at least K=2
    if [ \$MAX_K -lt 2 ]; then
        echo "Too few samples for ADMIXTURE (n=\$N_SAMP). Skipping."
        echo "optimal_k=1" > ${node_name}_admixture/choosek.txt
        echo "Skipped: too few samples" > ${node_name}_admixture/summary.txt
        exit 0
    fi

    BEST_K=2
    BEST_CV=999999

    for K in \$(seq 2 \$MAX_K); do
        echo "Running K=\$K"
        admixture --cv -j\${SLURM_CPUS_PER_TASK:-4} --seed=42 ${bed} \$K \
            > ${node_name}_admixture/admixture_K\${K}.log 2>&1 || {
            echo "K=\$K failed, skipping"
            continue
        }

        # Move output files
        mv ${bed.baseName}.\${K}.Q ${node_name}_admixture/ 2>/dev/null || true
        mv ${bed.baseName}.\${K}.P ${node_name}_admixture/ 2>/dev/null || true

        # Parse CV error
        CV=\$(grep "CV error" ${node_name}_admixture/admixture_K\${K}.log | awk '{print \$NF}')
        if [ -n "\$CV" ]; then
            echo "K=\$K  CV=\$CV" >> ${node_name}_admixture/cv_errors.txt
            LESS=\$(awk "BEGIN {print (\$CV < \$BEST_CV)}")
            if [ "\$LESS" = "1" ]; then
                BEST_CV=\$CV
                BEST_K=\$K
            fi
        fi
    done

    echo "optimal_k=\$BEST_K" > ${node_name}_admixture/choosek.txt
    echo "cv_error=\$BEST_CV" >> ${node_name}_admixture/choosek.txt

    # Copy optimal Q matrix with a clear name
    if [ -f "${node_name}_admixture/${bed.baseName}.\${BEST_K}.Q" ]; then
        cp ${node_name}_admixture/${bed.baseName}.\${BEST_K}.Q \
           ${node_name}_admixture/optimal_Q.txt
    fi

    # Attach sample IDs to optimal Q matrix
    if [ -f "${node_name}_admixture/optimal_Q.txt" ]; then
        awk 'NR==FNR {ids[NR]=\$1; next} {print ids[FNR] "\\t" \$0}' \
            ${fam} ${node_name}_admixture/optimal_Q.txt \
            > ${node_name}_admixture/optimal_Q_labeled.tsv
    fi

    echo "Optimal K=\$BEST_K (CV error=\$BEST_CV)" | tee ${node_name}_admixture/summary.txt
    """
}
