# ==============================================================================
#  Title: Run ssGSEA Immune Signature Scoring
#  2026-07-02 | Shiwoo Koak | Seoul National University
#
#  Description:
#    Calculate single-sample GSEA (ssGSEA) scores from TCGA RNA-seq expression
#    data using a GMT immune signature file. The script saves the ssGSEA score
#    matrix and two heatmaps: one showing Pearson correlations between
#    signatures, and one showing signature scores across TCGA samples. Samples
#    are restricted to primary tumors and low-library-size samples are removed.
#    Lowly expressed genes are filtered using raw counts without restricting to
#    protein-coding annotations. Genes are retained when they have at least 15
#    raw counts in at least 20% of retained samples.
#    ssGSEA is run without global score normalization, then scores are z-scored
#    separately for each gene set across samples.
#
#    The example immune_signature.gmt file used with this script was generated
#    from the immune signature collection reported in:
#    https://www.cell.com/cancer-cell/fulltext/S1535-6108(21)00222-1
#
#  Usage:
#    Rscript src/05_run_ssgsea.R <cancer_type> <gmt_file>
#
#  Output:
#    data/<cancer_type>_ssgsea_scores.tsv
#    data/<cancer_type>_ssgsea_signature_correlation_heatmap.svg
#    data/<cancer_type>_ssgsea_sample_score_heatmap.svg
#
#  Notes:
#    The GMT file is expected to contain signature name, signature category,
#    and gene symbols in columns 1, 2, and 3+, respectively.
#    The saved score table contains per-signature z-scored ssGSEA scores.
# ==============================================================================

################################################################################
# Step 0: Load Packages --------------------------------------------------------
################################################################################

library(dplyr)
library(readr)
library(GSVA)
library(pheatmap)
library(svglite)

################################################################################
# Step 1: Utilities ------------------------------------------------------------
################################################################################

calculate_mad <- function(x) {
  median(abs(x - median(x, na.rm = TRUE)), na.rm = TRUE)
}

read_gmt <- function(gmt_file) {
  
  gmt_lines <- readLines(gmt_file)
  gmt_list  <- strsplit(gmt_lines, "\t")
  
  signature_names <- vapply(gmt_list, `[`, character(1), 1)
  gene_sets <- lapply(gmt_list, function(x) {
    genes <- unique(x[-c(1, 2)])
    genes[genes != ""]
  })
  names(gene_sets) <- signature_names
  
  gene_sets
}

################################################################################
# Step 2: Prepare and Filter Input Data ----------------------------------------
################################################################################

filter_samples <- function(expression_data, counts_data) {
  
  sample_cols <- grep("^TCGA-", colnames(expression_data), value = TRUE)
  metadata_cols <- setdiff(colnames(expression_data), sample_cols)
  
  message("  - Samples in expression/counts data: ", length(sample_cols))
  
  sample_cols <- sample_cols[substr(sample_cols, 14, 16) == "01A"]
  message("  - Primary tumor samples (01A): ", length(sample_cols))
  
  sample_info <- data.frame(
    sample = sample_cols,
    stringsAsFactors = FALSE
  ) %>%
    arrange(sample)
  
  if (nrow(sample_info) == 0) {
    stop("No primary tumor samples were found.", call. = FALSE)
  }
  
  counts_matrix <- as.matrix(counts_data[, sample_info$sample])
  storage.mode(counts_matrix) <- "numeric"
  
  library_sizes <- colSums(counts_matrix, na.rm = TRUE)
  log2_library_sizes <- log2(library_sizes)
  
  library_median <- median(log2_library_sizes, na.rm = TRUE)
  library_mad    <- calculate_mad(log2_library_sizes)
  library_cutoff <- library_median - 3 * library_mad
  
  n_before_library_filter <- nrow(sample_info)
  
  sample_info <- sample_info %>%
    mutate(
      library_size      = library_sizes[sample],
      log2_library_size = log2_library_sizes[sample]
    ) %>%
    filter(
      is.finite(log2_library_size),
      log2_library_size >= library_cutoff
    )
  
  message(
    "  - Samples after library-size filtering: ",
    nrow(sample_info),
    " (removed ",
    n_before_library_filter - nrow(sample_info),
    ")"
  )
  
  if (nrow(sample_info) == 0) {
    stop(
      "No samples remained after library-size filtering.",
      call. = FALSE
    )
  }
  
  sample_cols <- sample_info$sample
  
  list(
    expression_data = expression_data %>%
      dplyr::select(all_of(metadata_cols), all_of(sample_cols)),
    counts_data = counts_data %>%
      dplyr::select(all_of(metadata_cols), all_of(sample_cols))
  )
}

filter_genes <- function(expression_data, counts_data) {
  
  message("  - Genes in expression/counts data: ", nrow(expression_data))
  
  sample_cols <- grep("^TCGA-", colnames(counts_data), value = TRUE)
  
  counts_matrix <- as.matrix(counts_data[, sample_cols])
  
  storage.mode(counts_matrix) <- "numeric"
  
  min_samples <- ceiling(0.20 * length(sample_cols))
  count_filter <- rowSums(counts_matrix >= 15, na.rm = TRUE) >= min_samples
  
  message(
    "  - Genes with >=15 raw counts in at least 20% of samples: ",
    sum(count_filter, na.rm = TRUE)
  )
  
  if (sum(count_filter, na.rm = TRUE) == 0) {
    stop("No genes remained after count-based filtering.", call. = FALSE)
  }
  
  list(
    expression_data = expression_data[count_filter, ],
    counts_data = counts_data[count_filter, ]
  )
}

prepare_analysis_data <- function(expression_data, counts_data) {
  
  sample_filtered <- filter_samples(
    expression_data = expression_data,
    counts_data     = counts_data
  )
  
  gene_filtered <- filter_genes(
    expression_data = sample_filtered$expression_data,
    counts_data     = sample_filtered$counts_data
  )
  
  expression_data <- gene_filtered$expression_data
  sample_cols <- grep("^TCGA-", colnames(expression_data), value = TRUE)
  
  expression_matrix <- expression_data %>%
    filter(!is.na(gene_name), gene_name != "") %>%
    group_by(gene_name) %>%
    summarize(
      across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  message("  - Genes with gene symbols for ssGSEA: ", nrow(expression_matrix))
  
  expression_mat <- as.matrix(expression_matrix[, sample_cols])
  rownames(expression_mat) <- expression_matrix$gene_name
  
  expression_mat
}

filter_gene_sets <- function(gene_sets, expression_genes, min_size = 3) {
  
  filtered_gene_sets <- lapply(gene_sets, function(genes) {
    intersect(genes, expression_genes)
  })
  
  missing_gene_sets <- lapply(gene_sets, function(genes) {
    setdiff(genes, expression_genes)
  })
  
  matched_sizes <- lengths(filtered_gene_sets)
  message("Gene set matching summary:")
  for (signature in names(filtered_gene_sets)) {
    message(
      paste0(
        "  - ", signature, ": ",
        matched_sizes[[signature]], " genes matched"
      )
    )
    
    missing_genes <- missing_gene_sets[[signature]]
    if (length(missing_genes) > 0) {
      message(
        paste0(
          "    Missing genes: ",
          paste(missing_genes, collapse = ", ")
        )
      )
    }
  }
  
  filtered_gene_sets <- filtered_gene_sets[matched_sizes >= min_size]
  
  if (length(filtered_gene_sets) == 0) {
    stop("No gene sets had enough matched genes for ssGSEA.", call. = FALSE)
  }
  
  filtered_gene_sets
}

run_ssgsea <- function(expression_mat, gene_sets) {
  
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    ssgsea_param <- GSVA::ssgseaParam(
      exprData  = expression_mat,
      geneSets  = gene_sets,
      normalize = FALSE
    )
    
    ssgsea_scores <- GSVA::gsva(
      param   = ssgsea_param,
      verbose = FALSE
    )
    
  } else {
    ssgsea_scores <- GSVA::gsva(
      expr          = expression_mat,
      gset.idx.list = gene_sets,
      method        = "ssgsea",
      ssgsea.norm   = FALSE,
      verbose       = FALSE
    )
  }
  
  as.matrix(ssgsea_scores)
}

zscore_gene_sets <- function(ssgsea_scores) {
  
  score_means <- rowMeans(ssgsea_scores, na.rm = TRUE)
  score_sds   <- apply(ssgsea_scores, 1, sd, na.rm = TRUE)
  
  constant_scores <- is.na(score_sds) | score_sds == 0
  if (any(constant_scores)) {
    message(
      "  - Gene sets with constant ssGSEA scores set to z-score 0: ",
      paste(rownames(ssgsea_scores)[constant_scores], collapse = ", ")
    )
    score_sds[constant_scores] <- 1
  }
  
  zscore_mat <- sweep(ssgsea_scores, 1, score_means, "-")
  zscore_mat <- sweep(zscore_mat, 1, score_sds, "/")
  zscore_mat[constant_scores, ] <- 0
  
  zscore_mat
}

################################################################################
# Step 3: Plot Heatmaps --------------------------------------------------------
################################################################################

plot_signature_correlation <- function(ssgsea_scores, output_path) {
  
  signature_cor <- cor(
    t(ssgsea_scores),
    method = "pearson",
    use    = "pairwise.complete.obs"
  )
  
  correlation_colors <- colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(101)
  
  svglite(output_path, width = 9, height = 8)
  on.exit(dev.off(), add = TRUE)
  
  pheatmap(
    signature_cor,
    color        = correlation_colors,
    breaks       = seq(-1, 1, length.out = 102),
    fontsize     = 8,
    angle_col    = "45",
    border_color = NA
  )
}

plot_sample_scores <- function(ssgsea_scores, output_path) {
  
  plot_scores <- ssgsea_scores
  color_limits <- quantile(
    plot_scores,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  plot_scores[plot_scores < color_limits[[1]]] <- color_limits[[1]]
  plot_scores[plot_scores > color_limits[[2]]] <- color_limits[[2]]
  
  score_colors <- colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(101)
  
  svglite(output_path, width = 10, height = 7)
  on.exit(dev.off(), add = TRUE)
  
  pheatmap(
    plot_scores,
    color         = score_colors,
    show_colnames = FALSE,
    fontsize      = 8,
    border_color  = NA
  )
}

################################################################################
# Step 4: Main Execution -------------------------------------------------------
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript src/05_run_ssgsea.R <cancer_type> <gmt_file>")
}

cancer_type <- toupper(args[1])
gmt_file <- args[2]

output_dir <- "data"
expression_file <- file.path(
  output_dir,
  paste0("TCGA_", cancer_type),
  paste0(cancer_type, "_rnaseq_log2tpm.tsv")
)
counts_file <- file.path(
  output_dir,
  paste0("TCGA_", cancer_type),
  paste0(cancer_type, "_rnaseq_counts.tsv")
)

message("Starting ssGSEA analysis")
message(paste("  - Cancer type:", cancer_type))
message(paste("  - Expression file:", expression_file))
message(paste("  - Counts file:", counts_file))
message(paste("  - GMT file:", gmt_file))

expression_data <- read_tsv(expression_file, show_col_types = FALSE)
counts_data <- read_tsv(counts_file, show_col_types = FALSE)
gmt_data <- read_gmt(gmt_file)

message("Filtering samples and genes...")
expression_mat <- prepare_analysis_data(
  expression_data = expression_data,
  counts_data     = counts_data
)

gene_sets <- filter_gene_sets(
  gene_sets        = gmt_data,
  expression_genes = rownames(expression_mat),
  min_size         = 3
)

message(paste("Running unnormalized ssGSEA for", length(gene_sets), "signatures..."))
ssgsea_raw_scores <- run_ssgsea(
  expression_mat = expression_mat,
  gene_sets      = gene_sets
)

message("Z-scoring ssGSEA scores within each gene set...")
ssgsea_scores <- zscore_gene_sets(ssgsea_raw_scores)

score_file <- file.path(output_dir, paste0(cancer_type, "_ssgsea_scores.tsv"))
write_tsv(
  as.data.frame(ssgsea_scores) %>%
    tibble::rownames_to_column("signature"),
  score_file
)
message(paste("Saved ssGSEA scores to:", score_file))

correlation_heatmap_file <- file.path(
  output_dir,
  paste0(cancer_type, "_ssgsea_signature_correlation_heatmap.svg")
)
plot_signature_correlation(
  ssgsea_scores = ssgsea_scores,
  output_path   = correlation_heatmap_file
)
message(paste("Saved signature correlation heatmap to:", correlation_heatmap_file))

sample_heatmap_file <- file.path(
  output_dir,
  paste0(cancer_type, "_ssgsea_sample_score_heatmap.svg")
)
plot_sample_scores(
  ssgsea_scores = ssgsea_scores,
  output_path   = sample_heatmap_file
)
message(paste("Saved sample score heatmap to:", sample_heatmap_file))

message("Done.")
