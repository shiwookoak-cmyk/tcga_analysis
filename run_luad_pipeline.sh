#!/usr/bin/env bash

# ==============================================================================
#  Title: Run Full LUAD TCGA Analysis Pipeline
#  2026-07-03 | Shiwoo Koak | Seoul National University
#
#  Description:
#    Run the LUAD example workflow: download TCGA data if needed, perform
#    gene-wise survival analysis, run pathway analysis, map significant KEGG metabolic
#    pathways, run ssGSEA immune signature scoring, run example expression
#    analyses, create an oncoprint, and copy generated figures for GitHub.
#
#  Usage:
#    bash run_luad_pipeline.sh
#
#  Output:
#    data/TCGA_LUAD/
#    data/LUAD_coxph_<method>.tsv
#    data/LUAD_coxph_<method>_gsea_*.tsv
#    data/LUAD_coxph_<method>_gsea_kegg_map.svg
#    data/LUAD_ssgsea_scores.tsv
#    data/LUAD_ssgsea_*_heatmap.svg
#    data/LUAD_*_box.svg
#    data/LUAD_*_scatter*.svg
#    data/LUAD_*_gene_group_heatmap.svg
#    data/LUAD_oncoprint_fisher.tsv
#    data/LUAD_oncoprint_fisher_heatmap.svg
#    data/LUAD_oncoprint.svg
#    figures/*.svg
#    data/logs/LUAD_pipeline_<timestamp>.log
# ==============================================================================

set -euo pipefail

CANCER="LUAD"
SURVIVAL_METHOD="continuous"
GMT_FILE="immune_signature.gmt"
KEGG_PADJ_CUTOFF="0.05"
DRIVER_GENES=("TP53" "KRAS" "STK11" "KEAP1" "EGFR" "NF1")
FIGURE_DIR="figures"

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PROJECT_DIR}"

mkdir -p data data/logs "${FIGURE_DIR}"

LOG_FILE="data/logs/${CANCER}_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "Starting full TCGA pipeline for ${CANCER}"
echo "  - Survival method: ${SURVIVAL_METHOD}"
echo "  - Driver genes: ${DRIVER_GENES[*]}"
echo "  - GMT file: ${GMT_FILE}"
echo "  - Project directory: ${PROJECT_DIR}"
echo "  - Log file: ${LOG_FILE}"

REQUIRED_DATA=(
  "data/TCGA_${CANCER}/${CANCER}_rnaseq_counts.tsv"
  "data/TCGA_${CANCER}/${CANCER}_rnaseq_log2tpm.tsv"
  "data/TCGA_${CANCER}/${CANCER}_clinical.tsv"
  "data/TCGA_${CANCER}/${CANCER}_snv.tsv"
  "data/TCGA_${CANCER}/${CANCER}_cnv.tsv"
)

DATA_READY="yes"
for data_file in "${REQUIRED_DATA[@]}"; do
  if [[ ! -s "${data_file}" ]]; then
    DATA_READY="no"
  fi
done

if [[ "${DATA_READY}" == "yes" ]]; then
  echo "Step 1/8: TCGA input files already exist. Skipping download."
else
  echo "Step 1/8: Downloading TCGA data..."
  Rscript src/01_download_tcga.R "${CANCER}"
fi

echo "Step 2/8: Running gene-wise survival analysis..."
Rscript src/02_gene_survival_analysis.R "${CANCER}" "${SURVIVAL_METHOD}"

SURVIVAL_FILE="data/${CANCER}_coxph_${SURVIVAL_METHOD}.tsv"
KEGG_GSEA_FILE="data/${CANCER}_coxph_${SURVIVAL_METHOD}_gsea_kegg.tsv"

echo "Step 3/8: Running pathway analysis..."
Rscript src/03_pathway_analysis.R "${SURVIVAL_FILE}"

echo "Step 4/8: Mapping significant KEGG metabolic pathways..."
if Rscript src/04_map_kegg_pathways.R "${KEGG_GSEA_FILE}" "${KEGG_PADJ_CUTOFF}"; then
  echo "KEGG metabolic map complete."
else
  echo "KEGG metabolic map was not created. Continuing with remaining analyses."
fi

echo "Step 5/8: Running ssGSEA immune signature scoring..."
Rscript src/05_run_ssgsea.R "${CANCER}" "${GMT_FILE}"

echo "Step 6/8: Running example mRNA expression analysis..."
python src/06_analyze_expression.py

echo "Step 7/8: Creating oncoprint..."
python src/07_oncoprint.py "${CANCER}" "${DRIVER_GENES[@]}"

echo "Step 8/8: Copying figures for GitHub upload..."
FIGURE_FILES=(
  "data/${CANCER}_coxph_${SURVIVAL_METHOD}_gsea_kegg_map.svg"
  "data/${CANCER}_ssgsea_signature_correlation_heatmap.svg"
  "data/${CANCER}_ssgsea_sample_score_heatmap.svg"
  "data/${CANCER}_PLA2G4A_paired_tumor_normal_box.svg"
  "data/${CANCER}_PLA2G4A_STK11_mutation_group_box.svg"
  "data/${CANCER}_PLA2G4A_PTGS2_scatter_sample_type.svg"
  "data/${CANCER}_PLA2G4A_PTGS2_STK11_scatter.svg"
  "data/${CANCER}_STK11_gene_group_heatmap.svg"
  "data/${CANCER}_oncoprint_fisher_heatmap.svg"
  "data/${CANCER}_oncoprint.svg"
)

for figure_file in "${FIGURE_FILES[@]}"; do
  if [[ -s "${figure_file}" ]]; then
    cp "${figure_file}" "${FIGURE_DIR}/"
    echo "  - Copied ${figure_file} to ${FIGURE_DIR}/"
  else
    echo "  - Figure not found, skipping: ${figure_file}"
  fi
done

echo "Pipeline complete."
