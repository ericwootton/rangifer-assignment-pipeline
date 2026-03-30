/*
========================================================================================
    BUILD SNP COORDINATE MAP
========================================================================================
    Builds the mapping from chip SNP names to genomic coordinates (CHROM:POS).

    Chip SNP names encode scaffold_position:
      - FinalReport: "1_1002319" = scaffold 1, position 1002319
      - NWT CSV: "Var100_1027221" = scaffold 100, position 1027221

    FASTA headers map scaffold numbers to accessions:
      ">JAHWTM010000001.1 Rangifer tarandus ... scaffold1, ..."

    Output: TSV with columns [snp_name, chrom, pos, snp_id]
    where snp_id = "CHROM:POS" matching the VCF format.
*/

process BUILD_SNP_MAP {
    tag "build_snp_map"
    label 'process_low'

    publishDir "${params.outdir}/data_integration", mode: 'copy'

    input:
    path manifest
    path fasta

    output:
    path "chip_snp_map.tsv",    emit: snp_map
    path "scaffold_map.tsv",    emit: scaffold_map

    script:
    """
    # Step 1: Build scaffold number -> accession mapping from FASTA headers
    # Header format: >JAHWTM010000001.1 Rangifer tarandus caribou isolate CA-31 scaffold1, ...
    grep '^>' ${fasta} | tr -d '\\r' | awk '{
        accession = substr(\$1, 2)
        for (i = 2; i <= NF; i++) {
            if (\$i ~ /^scaffold[0-9]/) {
                gsub(/,/, "", \$i)
                gsub(/^scaffold/, "", \$i)
                print \$i "\\t" accession
                break
            }
        }
    }' > scaffold_map.tsv

    echo "Scaffold map entries: \$(wc -l < scaffold_map.tsv)"

    # Step 2: Parse manifest to build SNP name -> CHROM:POS mapping
    # Manifest Name column (field 2): scaffoldNum_position or scaffoldNum_position_modl...
    # Also generate Var-prefixed names for NWT format
    tr -d '\\r' < ${manifest} | awk -F',' '
        NR == FNR {
            # scaffold_map.tsv is tab-separated; -F is comma for manifest
            split(\$0, a, "\\t")
            map[a[1]] = a[2]
            next
        }
        /^\\[Assay\\]/ { in_assay = 1; next }
        /^\\[/ { in_assay = 0; next }
        in_assay && \$2 != "" && \$2 != "Name" {
            name = \$2
            n = split(name, parts, "_")
            if (parts[1] !~ /^[0-9]+\$/) next
            scaffold_num = parts[1]
            position = parts[2]
            if (scaffold_num in map) {
                chrom = map[scaffold_num]
                snp_id = chrom ":" position
                # Chip name (FinalReport format): full name as-is
                print name "\\t" chrom "\\t" position "\\t" snp_id
                # Var name short (NWT format): Var{scaffold}_{position}
                var_name = "Var" scaffold_num "_" position
                print var_name "\\t" chrom "\\t" position "\\t" snp_id
                # Var name full (NWT format with suffixes): Var + full chip name
                if (n > 2) {
                    var_full = "Var" name
                    print var_full "\\t" chrom "\\t" position "\\t" snp_id
                }
            }
        }
    ' scaffold_map.tsv - | sort -u > chip_snp_map.tsv

    total=\$(wc -l < chip_snp_map.tsv)
    echo "Total SNP name mappings (chip + Var): \${total}"
    """
}
