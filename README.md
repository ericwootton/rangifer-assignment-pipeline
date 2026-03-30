# Rangifer Assignment Pipeline

A Nextflow pipeline for hierarchical probabilistic assignment of *Rangifer tarandus* (caribou/reindeer) samples to subspecies, ecotype, and herd using multi-method ensemble classification.

## Overview

This pipeline performs hierarchical population assignment through four stages:

1. **Data Integration** - Merges genotype data from multiple sources (WGS VCF, Illumina FinalReport, NWT chip format) with strand-aware allele resolution
2. **Quality Control** - SNP and sample filtering, duplicate detection, and cross-platform batch correction
3. **Structure Discovery** - PCA, DAPC cluster analysis, and ADMIXTURE to define the classification hierarchy from metadata
4. **Hierarchical Classification** - Per-node ensemble classification using DAPC, SVM/LDA/Random Forest (assignPOP-style), neural network (popfinder), and ADMIXTURE, with tiered confidence assignments

## Modes

The pipeline supports four execution modes via `--mode`:

| Mode | Description |
|------|-------------|
| `snp` (default) | SNP chip data with multi-platform integration and Fst-based top-5000 SNP selection per node |
| `mito` | Mitochondrial variants with haploid-appropriate QC (no strand filter, relaxed het thresholds) |
| `wgs` | Whole-genome data with SV calling (Delly), phasing (Eagle), Fst outlier detection (OutFLANK), windowed Fst/Dxy (pixy), IBD networks (IBDseq), and fineSTRUCTURE |
| `herd_diff` | Post-QC herd-level differentiation analysis (pairwise Fst, PCA, DAPC, ADMIXTURE, classifier accuracy) |

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.04.0
- [Singularity/Apptainer](https://apptainer.org/) for container execution
- SLURM scheduler (for HPC execution via the `rorqual` profile)

### Containers

The pipeline expects Singularity `.sif` containers in the directory specified by `--container_dir`. Required containers:

| Container | Used by |
|-----------|---------|
| `bcftools.sif` | VCF extraction, SV merging |
| `plink.sif` | BED format conversion |
| `plink2.sif` | (optional) |
| `r-popgen-dt.sif` | R-based analyses (QC, PCA, DAPC, assignPOP, ensemble, reporting) |
| `r-popgen-outflank.sif` | OutFLANK Fst outlier detection (WGS mode) |
| `popfinder.sif` | Neural network classifier |
| `admixture.sif` | ADMIXTURE ancestry estimation |
| `python-scientific.sif` | General Python utilities |
| `delly.sif` | Structural variant calling (WGS mode) |
| `eagle.sif` | Reference-free phasing (WGS mode) |
| `pixy.sif` | Windowed Fst/Dxy (WGS mode) |
| `ibdseq.sif` | IBD segment detection (WGS mode) |
| `finestructure.sif` | ChromoPainter/fineSTRUCTURE (WGS mode) |

## Quick Start

### SNP mode (multi-platform chip data)

```bash
nextflow run main.nf \
    -profile singularity \
    --mode snp \
    --wgs_vcf /path/to/chip_variants.vcf.gz \
    --finalreport /path/to/FinalReport.csv \
    --nwt_genotypes /path/to/NWT_genotypes.csv \
    --chip_manifest /path/to/chip_manifest.csv \
    --fasta /path/to/reference.fna \
    --metadata /path/to/metadata.csv \
    --container_dir /path/to/containers \
    --outdir results
```

### Mito mode

```bash
nextflow run main.nf \
    -profile singularity \
    --mode mito \
    --wgs_vcf /path/to/mito_variants.vcf.gz \
    --metadata /path/to/metadata.csv \
    --container_dir /path/to/containers \
    --outdir results \
    --min_het 0.0 --max_het 0.05 \
    --skip_strand_filter true
```

### WGS mode

```bash
nextflow run main.nf \
    -profile singularity \
    --mode wgs \
    --wgs_vcf /path/to/wgs_variants.vcf.gz \
    --bam_dir /path/to/bams \
    --fasta /path/to/reference.fna \
    --metadata /path/to/metadata.csv \
    --container_dir /path/to/containers \
    --outdir results \
    --skip_strand_filter true
```

### Herd differentiation mode

```bash
nextflow run main.nf \
    -profile singularity \
    --mode herd_diff \
    --genotype_matrix /path/to/batch_corrected_matrix.tsv \
    --metadata /path/to/herd_metadata.csv \
    --container_dir /path/to/containers \
    --outdir results
```

## Input Files

### Metadata CSV

A CSV file with the following columns:

| Column | Description |
|--------|-------------|
| `Sample` | Sample identifier (must match VCF sample IDs) |
| `Subspecies` | Level 1 classification label |
| `Ecotype_Cleaned` | Level 2 classification label |
| `Herd_Cleaned` | Level 3 classification label |

The pipeline auto-detects alternative column names (`sample_id`, `ID`, `Source_ID`, etc.).

### VCF

A bgzipped, tabix-indexed VCF file (`.vcf.gz` + `.tbi`). For SNP mode, this contains chip-derived variants; for mito mode, mitochondrial variants; for WGS mode, genome-wide variants.

## Parameters

### QC thresholds

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_call_rate` | 0.70 | Minimum sample call rate (platform-adjusted) |
| `--max_het` | 0.50 | Maximum sample heterozygosity |
| `--min_het` | 0.15 | Minimum sample heterozygosity |
| `--snp_call_rate` | 0.70 | Minimum SNP call rate |
| `--max_snp_het` | 0.49 | Maximum SNP heterozygosity (paralogy filter) |
| `--ibs_threshold` | 0.95 | IBS threshold for duplicate detection |
| `--skip_strand_filter` | false | Skip strand-ambiguous SNP removal (use for single-platform runs) |

### Classification

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_samples_dapc` | 5 | Minimum samples per class for DAPC |
| `--min_samples_assignpop` | 5 | Minimum samples per class for assignPOP |
| `--min_samples_popfinder` | 5 | Minimum samples per class for popfinder |
| `--max_k_structure` | 15 | Maximum K for ADMIXTURE |
| `--max_k_dapc` | 20 | Maximum K for DAPC cluster discovery |
| `--batch_alpha` | 0.01 | Significance threshold for batch correction |
| `--node_snp_call_rate` | 0.80 | Per-node SNP call rate threshold |

## Output

Key output files in `--outdir`:

```
results/
  data_integration/         # Merged genotype matrix, SNP map
  qc/                       # QC reports, filtered matrices
  structure_discovery/       # PCA, DAPC clusters, ADMIXTURE, hierarchy definition
  classification/
    node_data/              # Per-node genotype subsets
    dapc/                   # DAPC LOOCV results per node
    assignpop/              # SVM/LDA/RF LOOCV results per node
    popfinder/              # Neural network results per node
    ensemble/               # Tiered ensemble assignments per node
    node_plots/             # 12-page diagnostic PDFs per node
  wgs_analyses/             # (WGS mode) SV, phasing, Fst/Dxy, IBD, fineSTRUCTURE
  herd_differentiation/     # (herd_diff mode) Pairwise Fst, summary
  final_assignments.tsv     # Final hierarchical assignments with confidence tiers
  pipeline_report.pdf       # Summary report
  pipeline_info/            # Nextflow execution reports
```

### Confidence Tiers

Each sample assignment receives a confidence tier:

| Tier | Criteria |
|------|----------|
| **Strongly Supported** | SVM + both validators (DAPC and popfinder) agree, or SVM + DAPC agree when popfinder is absent |
| **Supported** | SVM + one validator agrees |
| **Provisional** | SVM assignment with no validator agreement |

## Profiles

| Profile | Description |
|---------|-------------|
| `singularity` | Enable Singularity container execution |
| `rorqual` | SLURM execution on Alliance Canada HPC clusters |
| `test` | Reduced resource limits for testing |
| `debug` | Print hostname before each process |

## License

GPL-3.0
