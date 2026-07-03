#!/usr/bin/env bash

# ==============================================================================
#  Title: Run Full LUAD TCGA Analysis Pipeline
#  2026-07-02 | Shiwoo Koak | Seoul National University
#
#  Description:
#    Run the LUAD example workflow: download TCGA data if needed, perform
#    gene-wise survival analysis, run GSEA, map significant KEGG metabolic
#    pathways, run ssGSEA immune signature scoring, and create an oncoprint.
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
#    data/LUAD_oncoprint_fisher.tsv
#    data/LUAD_oncoprint_fisher_heatmap.svg
#    data/LUAD_oncoprint.svg
#    data/logs/LUAD_pipeline_<timestamp>.log
# ==============================================================================

set -euo pipefail

CANCER="LUAD"
SURVIVAL_METHOD="median"
GMT_FILE="immune_signature.gmt"
KEGG_PADJ_CUTOFF="0.05"
DRIVER_GENES=("TP53" "KRAS" "STK11" "KEAP1" "EGFR" "NF1")

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${PROJECT_DIR}"

mkdir -p data data/logs

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
  echo "Step 1/6: TCGA input files already exist. Skipping download."
else
  echo "Step 1/6: Downloading TCGA data..."
  Rscript src/01_download_tcga.R "${CANCER}"
fi

echo "Step 2/6: Running gene-wise survival analysis..."
Rscript src/02_gene_survival_analysis.R "${CANCER}" "${SURVIVAL_METHOD}"

SURVIVAL_FILE="data/${CANCER}_coxph_${SURVIVAL_METHOD}.tsv"
KEGG_GSEA_FILE="data/${CANCER}_coxph_${SURVIVAL_METHOD}_gsea_kegg.tsv"

echo "Step 3/6: Running GSEA..."
Rscript src/03_run_gsea.R "${SURVIVAL_FILE}"

echo "Step 4/6: Mapping significant KEGG metabolic pathways..."
if Rscript src/04_map_kegg_pathways.R "${KEGG_GSEA_FILE}" "${KEGG_PADJ_CUTOFF}"; then
  echo "KEGG metabolic map complete."
else
  echo "KEGG metabolic map was not created. Continuing with remaining analyses."
fi

echo "Step 5/6: Running ssGSEA immune signature scoring..."
Rscript src/05_run_ssgsea.R "${CANCER}" "${GMT_FILE}"

echo "Step 6/6: Creating oncoprint..."
python src/07_oncoprint.py "${CANCER}" "${DRIVER_GENES[@]}"

echo "Pipeline complete."
