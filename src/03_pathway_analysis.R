# ==============================================================================
#  Title: Run Pathway Analysis
#  2026-06-23 | Shiwoo Koak | Seoul National University
#
#  Description:
#    Run pathway analysis from gene-level results. This script runs Gene Set
#    Enrichment Analysis (GSEA) for KEGG, Reactome, Hallmark, and Gene Ontology
#    (GO BP, GO MF, and GO CC). It also runs Over-Representation Analysis (ORA)
#    for KEGG, Reactome, Hallmark, and GO. ORA significant genes are selected
#    using FDR < 0.05, then split into HR > 1 and HR < 1 gene sets.
#
#  Usage:
#    Rscript src/03_pathway_analysis.R <analysis_file>
#
#  Output:
#    data/<analysis_file_basename>_gsea_kegg.tsv
#    data/<analysis_file_basename>_gsea_reactome.tsv
#    data/<analysis_file_basename>_gsea_hallmark.tsv
#    data/<analysis_file_basename>_gsea_gobp.tsv
#    data/<analysis_file_basename>_gsea_gomf.tsv
#    data/<analysis_file_basename>_gsea_gocc.tsv
#    data/<analysis_file_basename>_ora_<direction>_<database>.tsv
# ==============================================================================

################################################################################
# Step 0: Load Packages --------------------------------------------------------
################################################################################

library(dplyr)
library(readr)
library(tools)
library(clusterProfiler)
library(ReactomePA)
library(msigdbr)
library(org.Hs.eg.db)
library(GO.db)

set.seed(1)

################################################################################
# Step 1: Utilities ------------------------------------------------------------
################################################################################

create_rank_metric <- function(analysis_data) {
  
  message("Using hazard ratios and p-values to create the GSEA ranking metric.")
  
  rank_data <- analysis_data %>%
    mutate(
      ensembl_id = sub("\\..*$", "", gene_id),
      score      = -log10(pvalue) * sign(log2(HR))
    )
  
  rank_data <- rank_data %>%
    filter(!is.na(ensembl_id), ensembl_id != "", is.finite(score))
  
  if (nrow(rank_data) == 0) {
    stop("No valid genes remained after creating the GSEA ranking metric.")
  }
  
  rank_data <- rank_data %>%
    group_by(ensembl_id) %>%
    slice_max(order_by = abs(score), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(score)) %>%
    dplyr::select(ensembl_id, score)
  
  gene_list <- rank_data$score
  names(gene_list) <- rank_data$ensembl_id
  
  gene_list
}

ensembl_ids_to_entrez <- function(ensembl_ids) {
  
  ensembl_ids <- unique(ensembl_ids[!is.na(ensembl_ids) & ensembl_ids != ""])
  
  if (length(ensembl_ids) == 0) {
    return(character())
  }
  
  entrez_ids <- unname(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = ensembl_ids,
      column    = "ENTREZID",
      keytype   = "ENSEMBL",
      multiVals = "first"
    )
  )
  
  unique(entrez_ids[!is.na(entrez_ids) & entrez_ids != ""])
}

ensembl_to_entrez <- function(gene_list) {
  
  if (is.null(names(gene_list)) || any(is.na(names(gene_list)) | names(gene_list) == "")) {
    stop("The ranked gene list must be named with Ensembl gene IDs.")
  }
  
  ensembl_ids <- names(gene_list)
  
  entrez_ids <- unname(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = ensembl_ids,
      column    = "ENTREZID",
      keytype   = "ENSEMBL",
      multiVals = "first"
    )
  )
  
  rank_data <- data.frame(
    entrez_id = entrez_ids,
    score     = unname(gene_list),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(entrez_id), entrez_id != "", is.finite(score)) %>%
    group_by(entrez_id) %>%
    slice_max(order_by = abs(score), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(score)) %>%
    dplyr::select(entrez_id, score)
  
  if (nrow(rank_data) == 0) {
    stop("No genes could be mapped from Ensembl IDs to Entrez IDs.")
  }
  
  entrez_gene_list <- rank_data$score
  names(entrez_gene_list) <- rank_data$entrez_id
  
  entrez_gene_list
}

convert_core_enrichment <- function(gsea_df, keytype = "ENTREZID") {
  
  if (nrow(gsea_df) == 0 || !"core_enrichment" %in% colnames(gsea_df)) {
    return(gsea_df)
  }
  
  core_gene_list <- strsplit(
    ifelse(is.na(gsea_df$core_enrichment), "", gsea_df$core_enrichment),
    split = "/",
    fixed = TRUE
  )
  
  gene_ids <- unique(unlist(core_gene_list, use.names = FALSE))
  gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  
  if (length(gene_ids) == 0) {
    return(gsea_df)
  }
  
  gene_symbols <- unname(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = gene_ids,
      column    = "SYMBOL",
      keytype   = keytype,
      multiVals = "first"
    )
  )
  names(gene_symbols) <- gene_ids
  
  convert_gene_ids <- function(ids) {
    
    ids <- ids[!is.na(ids) & ids != ""]
    
    if (length(ids) == 0) {
      return(NA_character_)
    }
    
    symbols <- gene_symbols[ids]
    symbols[is.na(symbols)] <- ids[is.na(symbols)]
    
    paste(unname(symbols), collapse = "/")
  }
  
  gsea_df$core_enrichment <- vapply(
    core_gene_list,
    convert_gene_ids,
    character(1)
  )
  
  gsea_df
}

convert_gene_id_column <- function(enrichment_df, column_name, keytype = "ENTREZID") {
  
  if (nrow(enrichment_df) == 0 || !column_name %in% colnames(enrichment_df)) {
    return(enrichment_df)
  }
  
  gene_list <- strsplit(
    ifelse(is.na(enrichment_df[[column_name]]), "", enrichment_df[[column_name]]),
    split = "/",
    fixed = TRUE
  )
  
  gene_ids <- unique(unlist(gene_list, use.names = FALSE))
  gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  
  if (length(gene_ids) == 0) {
    return(enrichment_df)
  }
  
  gene_symbols <- unname(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = gene_ids,
      column    = "SYMBOL",
      keytype   = keytype,
      multiVals = "first"
    )
  )
  names(gene_symbols) <- gene_ids
  
  convert_gene_ids <- function(ids) {
    
    ids <- ids[!is.na(ids) & ids != ""]
    
    if (length(ids) == 0) {
      return(NA_character_)
    }
    
    symbols <- gene_symbols[ids]
    symbols[is.na(symbols)] <- ids[is.na(symbols)]
    
    paste(unname(symbols), collapse = "/")
  }
  
  enrichment_df[[column_name]] <- vapply(
    gene_list,
    convert_gene_ids,
    character(1)
  )
  
  enrichment_df
}

clean_hallmark_name <- function(term_id) {
  
  term_name <- sub("^HALLMARK_", "", term_id)
  gsub("_", " ", term_name)
}

get_hallmark_gene_sets <- function() {
  
  hallmark_data <- msigdbr::msigdbr(
    species    = "Homo sapiens",
    collection = "H"
  ) %>%
    transmute(
      ID          = gs_name,
      Description = clean_hallmark_name(gs_name),
      ensembl     = ensembl_gene
    ) %>%
    filter(!is.na(ensembl), ensembl != "") %>%
    distinct()
  
  if (nrow(hallmark_data) == 0) {
    stop("No Hallmark Ensembl genes were found from msigdbr.")
  }
  
  term2gene <- hallmark_data %>%
    distinct(ID, ensembl)
  
  term2name <- hallmark_data %>%
    distinct(ID, Description)
  
  term_names <- term2name$Description
  names(term_names) <- term2name$ID
  
  list(
    term2gene = term2gene,
    term2name = term2name,
    term_names = term_names
  )
}

get_go_descriptions <- function(ids) {
  
  go_terms <- GO.db::GOTERM[ids]
  descriptions <- unname(AnnotationDbi::Term(go_terms))
  names(descriptions) <- ids
  
  descriptions
}

get_kegg_descriptions <- function(ids) {
  
  kegg_data <- clusterProfiler::download_KEGG(
    species  = "hsa",
    keggType = "KEGG",
    keyType  = "ncbi-geneid"
  )
  
  descriptions <- kegg_data$KEGGPATHID2NAME$to
  names(descriptions) <- kegg_data$KEGGPATHID2NAME$from
  
  descriptions[ids]
}

get_reactome_descriptions <- function(ids) {
  
  reactome_data <- ReactomePA::gson_Reactome("human")
  descriptions <- reactome_data@gsid2name$name
  names(descriptions) <- reactome_data@gsid2name$gsid
  
  descriptions[ids]
}

repair_gsea_metadata <- function(gsea_fit, gsea_df, ontology, description_fun) {
  
  if (nrow(gsea_df) == 0) {
    return(gsea_df)
  }
  
  if (!"ONTOLOGY" %in% colnames(gsea_df)) {
    gsea_df <- gsea_df %>%
      mutate(ONTOLOGY = ontology, .before = ID)
  }
  
  missing_rows <- which(is.na(gsea_df$ID))
  
  if (length(missing_rows) == 0) {
    return(gsea_df)
  }
  
  gene_set_sizes <- vapply(gsea_fit@geneSets, length, integer(1))
  min_gs_size <- gsea_fit@params$minGSSize
  max_gs_size <- gsea_fit@params$maxGSSize
  
  result_ids <- gsea_df$ID[!is.na(gsea_df$ID)]
  missing_ids <- setdiff(names(gsea_fit@geneSets), result_ids)
  missing_ids <- missing_ids[
    gene_set_sizes[missing_ids] >= min_gs_size &
      gene_set_sizes[missing_ids] <= max_gs_size
  ]
  
  n_fill <- min(length(missing_rows), length(missing_ids))
  
  if (n_fill == 0) {
    return(gsea_df)
  }
  
  fill_rows <- missing_rows[seq_len(n_fill)]
  fill_ids  <- missing_ids[seq_len(n_fill)]
  descriptions <- description_fun(fill_ids)
  
  gsea_df$ONTOLOGY[fill_rows]    <- ontology
  gsea_df$ID[fill_rows]          <- fill_ids
  gsea_df$Description[fill_rows] <- unname(descriptions[fill_ids])
  gsea_df$setSize[fill_rows]     <- gene_set_sizes[fill_ids]
  
  message("  - Recovered metadata for ", n_fill, " ", ontology, " terms with NA statistics")
  
  gsea_df
}

format_gsea_result <- function(gsea_fit, ontology, description_fun, core_keytype = "ENTREZID") {
  
  gsea_df <- as.data.frame(gsea_fit)
  
  if (nrow(gsea_df) == 0) {
    return(gsea_df)
  }
  
  gsea_df <- repair_gsea_metadata(
    gsea_fit        = gsea_fit,
    gsea_df         = gsea_df,
    ontology        = ontology,
    description_fun = description_fun
  )
  
  gsea_df %>%
    arrange(p.adjust, desc(NES)) %>%
    convert_core_enrichment(keytype = core_keytype) %>%
    dplyr::select(ONTOLOGY, ID, Description, setSize, everything())
}

save_result <- function(result_df, output_dir, analysis_name, suffix) {
  
  output_path <- file.path(output_dir, paste0(analysis_name, "_", suffix, ".tsv"))
  write_tsv(result_df, output_path)
  message("Saved ", suffix, " results to: ", output_path)
}

################################################################################
# Step 2: Run GSEA -------------------------------------------------------------
################################################################################

run_gsea_hallmark <- function(gene_list) {
  
  message("Running clusterProfiler GSEA: Hallmark")
  message("  - Ranked Ensembl genes: ", length(gene_list))
  
  hallmark_sets <- get_hallmark_gene_sets()
  
  message("  - Hallmark gene sets: ", n_distinct(hallmark_sets$term2gene$ID))
  
  gsea_fit <- clusterProfiler::GSEA(
    geneList     = gene_list,
    TERM2GENE    = hallmark_sets$term2gene,
    TERM2NAME    = hallmark_sets$term2name,
    pvalueCutoff = 1,
    verbose      = FALSE
  )
  
  format_gsea_result(
    gsea_fit        = gsea_fit,
    ontology        = "Hallmark",
    description_fun = function(ids) hallmark_sets$term_names[ids],
    core_keytype    = "ENSEMBL"
  )
}

run_gsea_reactome <- function(gene_list) {
  
  message("Running ReactomePA gsePathway")
  entrez_gene_list <- ensembl_to_entrez(gene_list)
  message("  - Ranked Entrez genes: ", length(entrez_gene_list))
  
  gsea_fit <- ReactomePA::gsePathway(
    geneList     = entrez_gene_list,
    organism     = "human",
    pvalueCutoff = 1,
    verbose      = FALSE
  )
  
  format_gsea_result(
    gsea_fit        = gsea_fit,
    ontology        = "Reactome",
    description_fun = get_reactome_descriptions,
    core_keytype    = "ENTREZID"
  )
}

run_gsea_kegg <- function(gene_list) {
  
  message("Running clusterProfiler gseKEGG")
  entrez_gene_list <- ensembl_to_entrez(gene_list)
  message("  - Ranked Entrez genes: ", length(entrez_gene_list))
  
  gsea_fit <- gseKEGG(
    geneList     = entrez_gene_list,
    organism     = "hsa",
    keyType      = "ncbi-geneid",
    pvalueCutoff = 1,
    verbose      = FALSE
  )
  
  format_gsea_result(
    gsea_fit        = gsea_fit,
    ontology        = "KEGG",
    description_fun = get_kegg_descriptions,
    core_keytype    = "ENTREZID"
  )
}

run_gsea_go <- function(gene_list, ontology) {
  
  message("Running clusterProfiler gseGO: ", ontology)
  message("  - Ranked Ensembl genes: ", length(gene_list))
  
  gsea_fit <- gseGO(
    geneList     = gene_list,
    OrgDb        = org.Hs.eg.db,
    keyType      = "ENSEMBL",
    ont          = ontology,
    pvalueCutoff = 1,
    verbose      = FALSE
  )
  
  format_gsea_result(
    gsea_fit        = gsea_fit,
    ontology        = ontology,
    description_fun = get_go_descriptions,
    core_keytype    = "ENSEMBL"
  )
}

################################################################################
# Step 3: Run ORA --------------------------------------------------------------
################################################################################

create_ora_gene_sets <- function(analysis_data) {
  
  ora_data <- analysis_data %>%
    mutate(ensembl_id = sub("\\..*$", "", gene_id)) %>%
    filter(
      !is.na(ensembl_id),
      ensembl_id != "",
      !is.na(HR),
      !is.na(FDR),
      is.finite(HR),
      HR > 0
    ) %>%
    group_by(ensembl_id) %>%
    slice_max(order_by = abs(log2(HR)), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  if (nrow(ora_data) == 0) {
    stop("No valid genes remained for ORA.")
  }
  
  up_ensembl <- ora_data %>%
    filter(HR > 1, FDR < 0.05) %>%
    pull(ensembl_id)
  
  down_ensembl <- ora_data %>%
    filter(HR < 1, FDR < 0.05) %>%
    pull(ensembl_id)
  
  universe_ensembl <- unique(ora_data$ensembl_id)
  universe_entrez  <- ensembl_ids_to_entrez(universe_ensembl)
  up_entrez        <- ensembl_ids_to_entrez(up_ensembl)
  down_entrez      <- ensembl_ids_to_entrez(down_ensembl)
  
  if (length(universe_entrez) == 0) {
    stop("No ORA universe genes could be mapped from Ensembl IDs to Entrez IDs.")
  }
  
  message("ORA universe Ensembl genes: ", length(universe_ensembl))
  message("ORA universe Entrez genes: ", length(universe_entrez))
  message("ORA up genes (HR > 1) with Ensembl IDs: ", length(up_ensembl))
  message("ORA down genes (HR < 1) with Ensembl IDs: ", length(down_ensembl))
  message("ORA up genes (HR > 1) with Entrez IDs: ", length(up_entrez))
  message("ORA down genes (HR < 1) with Entrez IDs: ", length(down_entrez))
  
  list(
    universe_ensembl = universe_ensembl,
    universe_entrez  = universe_entrez,
    up = list(
      ensembl = up_ensembl,
      entrez  = up_entrez
    ),
    down = list(
      ensembl = down_ensembl,
      entrez  = down_entrez
    )
  )
}

format_ora_result <- function(ora_fit, ontology, gene_keytype = "ENTREZID") {
  
  ora_df <- as.data.frame(ora_fit)
  
  if (nrow(ora_df) == 0) {
    return(ora_df)
  }
  
  if (!"ONTOLOGY" %in% colnames(ora_df)) {
    ora_df <- ora_df %>%
      mutate(ONTOLOGY = ontology, .before = ID)
  }
  
  ora_df %>%
    arrange(p.adjust, pvalue) %>%
    convert_gene_id_column(column_name = "geneID", keytype = gene_keytype) %>%
    dplyr::select(ONTOLOGY, ID, Description, everything())
}

run_ora_go <- function(entrez_genes, universe_entrez, ontology) {
  
  message("Running clusterProfiler enrichGO: ", ontology)
  message("  - Significant Entrez genes: ", length(entrez_genes))
  
  if (length(entrez_genes) == 0) {
    return(data.frame())
  }
  
  ora_fit <- enrichGO(
    gene          = entrez_genes,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = ontology,
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = FALSE
  )
  
  format_ora_result(
    ora_fit      = ora_fit,
    ontology     = ontology,
    gene_keytype = "ENTREZID"
  )
}

run_ora_kegg <- function(entrez_genes, universe_entrez) {
  
  message("Running clusterProfiler enrichKEGG")
  message("  - Significant Entrez genes: ", length(entrez_genes))
  
  if (length(entrez_genes) == 0) {
    return(data.frame())
  }
  
  ora_fit <- enrichKEGG(
    gene          = entrez_genes,
    universe      = universe_entrez,
    organism      = "hsa",
    keyType       = "ncbi-geneid",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
  
  format_ora_result(
    ora_fit      = ora_fit,
    ontology     = "KEGG",
    gene_keytype = "ENTREZID"
  )
}

run_ora_reactome <- function(entrez_genes, universe_entrez) {
  
  message("Running ReactomePA enrichPathway")
  message("  - Significant Entrez genes: ", length(entrez_genes))
  
  if (length(entrez_genes) == 0) {
    return(data.frame())
  }
  
  ora_fit <- ReactomePA::enrichPathway(
    gene          = entrez_genes,
    universe      = universe_entrez,
    organism      = "human",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = FALSE
  )
  
  format_ora_result(
    ora_fit      = ora_fit,
    ontology     = "Reactome",
    gene_keytype = "ENTREZID"
  )
}

run_ora_hallmark <- function(ensembl_genes, universe_ensembl) {
  
  message("Running clusterProfiler enricher: Hallmark")
  message("  - Significant Ensembl genes: ", length(ensembl_genes))
  
  if (length(ensembl_genes) == 0) {
    return(data.frame())
  }
  
  hallmark_sets <- get_hallmark_gene_sets()
  
  ora_fit <- clusterProfiler::enricher(
    gene          = ensembl_genes,
    universe      = universe_ensembl,
    TERM2GENE     = hallmark_sets$term2gene,
    TERM2NAME     = hallmark_sets$term2name,
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
  
  format_ora_result(
    ora_fit      = ora_fit,
    ontology     = "Hallmark",
    gene_keytype = "ENSEMBL"
  )
}

################################################################################
# Step 4: Main Execution -------------------------------------------------------
################################################################################

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript src/03_pathway_analysis.R <analysis_file>")
}

analysis_file <- args[1]
output_dir    <- "data"

message("Starting pathway analysis")
message("Analysis file: ", analysis_file)
message("Output directory: ", output_dir)

message("Loading gene-level results...")
analysis_data <- read_tsv(analysis_file, show_col_types = FALSE)

message("Creating ranked gene list...")
gene_list <- create_rank_metric(analysis_data)
message("Ranked Ensembl genes: ", length(gene_list))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

analysis_name <- file_path_sans_ext(basename(analysis_file))

hallmark_df <- run_gsea_hallmark(gene_list)
save_result(hallmark_df, output_dir, analysis_name, "gsea_hallmark")

reactome_df <- run_gsea_reactome(gene_list)
save_result(reactome_df, output_dir, analysis_name, "gsea_reactome")

kegg_df <- run_gsea_kegg(gene_list)
save_result(kegg_df, output_dir, analysis_name, "gsea_kegg")

gobp_df <- run_gsea_go(gene_list, "BP")
save_result(gobp_df, output_dir, analysis_name, "gsea_gobp")

gomf_df <- run_gsea_go(gene_list, "MF")
save_result(gomf_df, output_dir, analysis_name, "gsea_gomf")

gocc_df <- run_gsea_go(gene_list, "CC")
save_result(gocc_df, output_dir, analysis_name, "gsea_gocc")


message("Creating ORA gene sets...")
ora_genes <- create_ora_gene_sets(analysis_data)

message("Running ORA...")

ora_up_kegg_df <- run_ora_kegg(
  entrez_genes    = ora_genes$up$entrez,
  universe_entrez = ora_genes$universe_entrez
)
save_result(ora_up_kegg_df, output_dir, analysis_name, "ora_up_kegg")

ora_down_kegg_df <- run_ora_kegg(
  entrez_genes    = ora_genes$down$entrez,
  universe_entrez = ora_genes$universe_entrez
)
save_result(ora_down_kegg_df, output_dir, analysis_name, "ora_down_kegg")

ora_up_reactome_df <- run_ora_reactome(
  entrez_genes    = ora_genes$up$entrez,
  universe_entrez = ora_genes$universe_entrez
)
save_result(ora_up_reactome_df, output_dir, analysis_name, "ora_up_reactome")

ora_down_reactome_df <- run_ora_reactome(
  entrez_genes    = ora_genes$down$entrez,
  universe_entrez = ora_genes$universe_entrez
)
save_result(ora_down_reactome_df, output_dir, analysis_name, "ora_down_reactome")

ora_up_hallmark_df <- run_ora_hallmark(
  ensembl_genes    = ora_genes$up$ensembl,
  universe_ensembl = ora_genes$universe_ensembl
)
save_result(ora_up_hallmark_df, output_dir, analysis_name, "ora_up_hallmark")

ora_down_hallmark_df <- run_ora_hallmark(
  ensembl_genes    = ora_genes$down$ensembl,
  universe_ensembl = ora_genes$universe_ensembl
)
save_result(ora_down_hallmark_df, output_dir, analysis_name, "ora_down_hallmark")

ora_up_gobp_df <- run_ora_go(
  entrez_genes    = ora_genes$up$entrez,
  universe_entrez = ora_genes$universe_entrez,
  ontology        = "BP"
)
save_result(ora_up_gobp_df, output_dir, analysis_name, "ora_up_gobp")

ora_down_gobp_df <- run_ora_go(
  entrez_genes    = ora_genes$down$entrez,
  universe_entrez = ora_genes$universe_entrez,
  ontology        = "BP"
)
save_result(ora_down_gobp_df, output_dir, analysis_name, "ora_down_gobp")

ora_up_gomf_df <- run_ora_go(
  entrez_genes    = ora_genes$up$entrez,
  universe_entrez = ora_genes$universe_entrez,
  ontology        = "MF"
)
save_result(ora_up_gomf_df, output_dir, analysis_name, "ora_up_gomf")

ora_down_gomf_df <- run_ora_go(
  entrez_genes    = ora_genes$down$entrez,
  universe_entrez = ora_genes$universe_entrez,
  ontology        = "MF"
)
save_result(ora_down_gomf_df, output_dir, analysis_name, "ora_down_gomf")

ora_up_gocc_df <- run_ora_go(
  entrez_genes    = ora_genes$up$entrez,
  universe_entrez = ora_genes$universe_entrez,
  ontology        = "CC"
)
save_result(ora_up_gocc_df, output_dir, analysis_name, "ora_up_gocc")

ora_down_gocc_df <- run_ora_go(
  entrez_genes    = ora_genes$down$entrez,
  universe_entrez = ora_genes$universe_entrez,
  ontology        = "CC"
)
save_result(ora_down_gocc_df, output_dir, analysis_name, "ora_down_gocc")

message("Pipeline complete!")
