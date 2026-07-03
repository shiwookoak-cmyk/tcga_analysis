# ==============================================================================
#  Title: Download TCGA Data
#  2026-06-16 | Shiwoo Koak | Seoul National University
#
#  Description:
#    Download TCGA RNA-seq, clinical, somatic mutation, and GISTIC gene-level
#    copy-number data for downstream quality control and analysis.
#
#  Usage:
#    Rscript src/01_download_tcga.R <cancer_type>
#
#  Output:
#    data/TCGA_<cancer_type>/<cancer_type>_rnaseq_counts.tsv
#    data/TCGA_<cancer_type>/<cancer_type>_rnaseq_log2tpm.tsv
#    data/TCGA_<cancer_type>/<cancer_type>_clinical.tsv
#    data/TCGA_<cancer_type>/<cancer_type>_snv.tsv
#    data/TCGA_<cancer_type>/<cancer_type>_cnv.tsv
# ==============================================================================

################################################################################
# Step 0: Load Packages --------------------------------------------------------
################################################################################

library(dplyr)
library(readr)
library(TCGAbiolinks)
library(SummarizedExperiment)

################################################################################
# Step 1: Utilities ------------------------------------------------------------
################################################################################

normalize_barcode <- function(x) {
  substr(x, 1, 16)
}

################################################################################
# Step 2: Download and Process RNA-seq -----------------------------------------
################################################################################

get_rnaseq <- function(cancer_type) {
  
  query_rna <- GDCquery(
    project       = paste0("TCGA-", cancer_type),
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  GDCdownload(query_rna)
  rna_seq <- GDCprepare(query_rna)
  
  # Build gene metadata (retaining raw versioned IDs)
  gene_metadata <- data.frame(
    gene_name = rowData(rna_seq)$gene_name,
    gene_id   = rowData(rna_seq)$gene_id,
    gene_type = rowData(rna_seq)$gene_type,
    stringsAsFactors = FALSE
  )
  
  # Extract matrices
  count_mat <- assays(rna_seq)[["unstranded"]]
  tpm_mat   <- assays(rna_seq)[["tpm_unstrand"]]
  
  # 1. Sort alphabetically by full barcode (puts 11R before 21R)
  count_mat <- count_mat[, order(colnames(count_mat)), drop = FALSE]
  tpm_mat   <- tpm_mat[, order(colnames(tpm_mat)), drop = FALSE]
  
  # 2. Truncate to 16 characters (TCGA-XX-YYYY-01A)
  colnames(count_mat) <- normalize_barcode(colnames(count_mat))
  colnames(tpm_mat)   <- normalize_barcode(colnames(tpm_mat))
  
  # 3. Drop duplicate aliquots after barcode truncation (keeps first occurrence)
  count_mat <- count_mat[, !duplicated(colnames(count_mat)), drop = FALSE]
  tpm_mat   <- tpm_mat[, !duplicated(colnames(tpm_mat)), drop = FALSE]
  
  # Build tidy data frames
  count_df <- as.data.frame(count_mat, check.names = FALSE) %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(gene_metadata, by = "gene_id") %>%
    relocate(gene_name, gene_id, gene_type)
  
  tpm_df <- as.data.frame(tpm_mat, check.names = FALSE) %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(gene_metadata, by = "gene_id") %>%
    relocate(gene_name, gene_id, gene_type)
  
  # Apply log2(TPM + 1) transformation
  log2tpm_df <- tpm_df %>%
    mutate(across(where(is.numeric), ~ log2(. + 1)))
  
  list(
    count   = count_df,
    log2tpm = log2tpm_df
  )
}

################################################################################
# Step 3: Download Clinical Data -----------------------------------------------
################################################################################

get_clinical <- function(cancer_type) {
  
  clinical_data <- GDCquery_clinic(
    project = paste0("TCGA-", cancer_type),
    type    = "clinical"
  )
  
  clinical_out <- clinical_data %>%
    mutate(
      age = as.numeric(age_at_index),
      sex = sex_at_birth,
      OS.time = coalesce(
        as.numeric(days_to_death),
        as.numeric(days_to_last_follow_up)
      ),
      OS = as.integer(toupper(vital_status) == "DEAD")
    ) %>%
    dplyr::select(
      sample = submitter_id,
      age,
      sex,
      OS.time,
      OS,
      primary_diagnosis
    ) %>%
    distinct()
  
  clinical_out
}

################################################################################
# Step 4: Download Somatic Mutation Data ---------------------------------------
################################################################################

get_snv <- function(cancer_type) {
  
  query_snv <- GDCquery(
    project       = paste0("TCGA-", cancer_type),
    data.category = "Simple Nucleotide Variation",
    data.type     = "Masked Somatic Mutation"
  )
  
  GDCdownload(query_snv)
  snv_data <- GDCprepare(query_snv)
  
  protein_altering_variants <- c(
    "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site",
    "Translation_Start_Site", "Nonstop_Mutation", "In_Frame_Del",
    "In_Frame_Ins", "Missense_Mutation"
  )
  
  snv_out <- snv_data %>%
    mutate(sample = normalize_barcode(Tumor_Sample_Barcode)) %>%
    transmute(
      sample,
      gene_name = Hugo_Symbol,
      gene_id   = Gene,
      variant_classification = Variant_Classification,
      all_effects
    ) %>%
    filter(variant_classification %in% protein_altering_variants)
  
  snv_out
}

################################################################################
# Step 5: Download GISTIC Gene-Level CNV ---------------------------------------
################################################################################

get_cnv <- function(cancer_type) {
  
  cnv_data <- getGistic(disease = cancer_type, type = "thresholded")
  
  colnames(cnv_data)[1]  <- "gene_name"
  colnames(cnv_data)[-1] <- normalize_barcode(colnames(cnv_data)[-1])
  
  cnv_data
}

################################################################################
# Step 6: Main Execution -------------------------------------------------------
################################################################################

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript src/01_download_tcga.R <cancer_type>")
}

cancer_type <- toupper(args[1])
output_base <- "data"

# Create output directory relative to the provided base
output_dir <- file.path(output_base, paste0("TCGA_", cancer_type))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message(paste("Starting pipeline for", cancer_type))
message(paste("Output will be saved to", output_dir))

# 1) RNA-seq data
message("Downloading and processing RNA-seq data...")
rnaseq_data <- get_rnaseq(cancer_type = cancer_type)

write_tsv(
  rnaseq_data$count,
  file.path(output_dir, paste0(cancer_type, "_rnaseq_counts.tsv"))
)

write_tsv(
  rnaseq_data$log2tpm,
  file.path(output_dir, paste0(cancer_type, "_rnaseq_log2tpm.tsv"))
)

# 2) Clinical data
message("Downloading clinical data...")
clinical_df <- get_clinical(cancer_type)
write_tsv(
  clinical_df,
  file.path(output_dir, paste0(cancer_type, "_clinical.tsv"))
)

# 3) Somatic mutation data
message("Downloading somatic mutation data...")
snv_df <- get_snv(cancer_type)
write_tsv(
  snv_df,
  file.path(output_dir, paste0(cancer_type, "_snv.tsv"))
)

# 4) GISTIC gene-level CNV data
message("Downloading CNV data...")
cnv_df <- get_cnv(cancer_type)
write_tsv(
  cnv_df,
  file.path(output_dir, paste0(cancer_type, "_cnv.tsv"))
)

message("Pipeline complete!")
