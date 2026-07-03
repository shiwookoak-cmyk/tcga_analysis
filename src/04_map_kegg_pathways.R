# ==============================================================================
#  Title: Map KEGG GSEA Pathways on the Global Metabolic Map
#  2026-06-25 | Shiwoo Koak | Seoul National University
#
#  Description:
#    Visualize significant metabolic KEGG GSEA pathways from 03_run_gsea.R on the
#    human global metabolic map. Significant KEGG pathways are restricted to the
#    metabolism class, converted to their annotated KEGG modules, and highlighted
#    on hsa01100. Pathway names are used as labels. Positive NES pathways are
#    colored red, and negative NES pathways are colored blue.
#
#  Usage:
#    Rscript src/04_map_kegg_pathways.R <gsea_file> [padj_cutoff]
#
#  Output:
#    data/<gsea_file_basename>_map.svg
# ==============================================================================

################################################################################
# Step 0: Load Packages --------------------------------------------------------
################################################################################

library(dplyr)
library(readr)
library(tools)
library(ggplot2)
library(ggrepel)
library(ggraph)
library(ggfx)
library(KEGGREST)
library(tidygraph)
library(ggkegg)
library(svglite)

################################################################################
# Step 1: Utilities ------------------------------------------------------------
################################################################################

collapse_ids <- function(ids, max_ids = 30) {
  
  ids <- unique(ids[!is.na(ids) & ids != ""])
  
  if (length(ids) == 0) {
    return("none")
  }
  
  if (length(ids) > max_ids) {
    return(
      paste0(
        paste(ids[seq_len(max_ids)], collapse = ", "),
        ", ... (",
        length(ids),
        " total)"
      )
    )
  }
  
  paste(ids, collapse = ", ")
}

wrap_label <- function(label, width = 22) {
  paste(strwrap(label, width = width), collapse = "\n")
}

clean_module_name <- function(module_name) {
  sub("\\s+\\[PATH:[^]]+\\]$", "", module_name)
}

standardize_module_id <- function(module_id) {
  sub("^hsa_", "", module_id)
}

has_module_highlight <- function(data, module_ids) {
  
  module_ids <- intersect(module_ids, colnames(data))
  
  if (length(module_ids) == 0) {
    return(rep(FALSE, nrow(data)))
  }
  
  rowSums(data[, module_ids, drop = FALSE], na.rm = TRUE) > 0
}

is_artificial_origin <- function(x, y) {
  is.finite(x) & is.finite(y) & x == 0 & y == 0
}

################################################################################
# Step 2: Prepare GSEA and KEGG Module Data ------------------------------------
################################################################################

filter_significant_gsea <- function(gsea_data, padj_cutoff) {
  
  if (is.na(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff > 1) {
    stop("padj_cutoff must be greater than 0 and less than or equal to 1.")
  }
  
  if (!"qvalue" %in% colnames(gsea_data)) {
    gsea_data$qvalue <- NA_real_
  }
  
  significant_pathways <- gsea_data %>%
    filter(
      !is.na(ID),
      ID != "",
      !is.na(NES),
      !is.na(p.adjust),
      p.adjust <= padj_cutoff,
      NES != 0
    ) %>%
    mutate(
      direction = ifelse(NES > 0, "positive", "negative")
    ) %>%
    arrange(p.adjust, desc(abs(NES)))
  
  if (nrow(significant_pathways) == 0) {
    stop("No significant pathways were found.")
  }
  
  significant_pathways
}

get_pathway_class <- function(pathway_id) {
  
  pathway_data <- tryCatch(
    KEGGREST::keggGet(pathway_id)[[1]],
    error = function(e) {
      message("  - Skipping pathway ", pathway_id, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(pathway_data) || is.null(pathway_data$CLASS)) {
    return(NA_character_)
  }
  
  paste(pathway_data$CLASS, collapse = "; ")
}

filter_metabolic_pathways <- function(significant_pathways) {
  
  pathway_classes <- vapply(
    significant_pathways$ID,
    get_pathway_class,
    character(1)
  )
  
  metabolic_pathways <- significant_pathways %>%
    mutate(kegg_class = pathway_classes) %>%
    filter(!is.na(kegg_class), grepl("^Metabolism", kegg_class))
  
  if (nrow(metabolic_pathways) == 0) {
    stop("No significant metabolic KEGG pathways were found.")
  }
  
  metabolic_pathways
}

get_pathway_modules <- function(metabolic_pathways) {
  
  module_data <- lapply(seq_len(nrow(metabolic_pathways)), function(i) {
    
    pathway_id <- metabolic_pathways$ID[i]
    
    pathway_data <- tryCatch(
      KEGGREST::keggGet(pathway_id)[[1]],
      error = function(e) {
        message("  - Skipping pathway ", pathway_id, ": ", conditionMessage(e))
        NULL
      }
    )
    
    if (is.null(pathway_data) || is.null(pathway_data$MODULE)) {
      message("  - No KEGG modules found for pathway ", pathway_id)
      return(NULL)
    }
    
    module_ids <- names(pathway_data$MODULE)
    module_names <- unname(pathway_data$MODULE)
    
    data.frame(
      direction       = metabolic_pathways$direction[i],
      ID              = pathway_id,
      Description     = metabolic_pathways$Description[i],
      NES             = metabolic_pathways$NES[i],
      pvalue          = metabolic_pathways$pvalue[i],
      p.adjust        = metabolic_pathways$p.adjust[i],
      qvalue          = metabolic_pathways$qvalue[i],
      kegg_class      = metabolic_pathways$kegg_class[i],
      module_id       = standardize_module_id(module_ids),
      human_module_id = module_ids,
      module_name     = clean_module_name(module_names),
      row.names       = NULL
    )
  })
  
  module_table <- bind_rows(module_data) %>%
    distinct() %>%
    arrange(match(direction, c("positive", "negative")), p.adjust, ID, module_id)
  
  if (nrow(module_table) == 0) {
    stop("No significant metabolic KEGG pathways had KEGG module annotations.")
  }
  
  module_table
}

################################################################################
# Step 3: Build Metabolic Map --------------------------------------------------
################################################################################

add_module_highlights <- function(graph, module_table) {
  
  highlighted_modules <- character()
  module_ids <- unique(module_table$module_id)
  
  for (module_id in module_ids) {
    
    graph <- tryCatch(
      graph %>%
        highlight_module(module(module_id)),
      error = function(e) {
        message("  - Could not highlight module ", module_id, ": ", conditionMessage(e))
        graph
      }
    )
    
    node_data <- graph %>%
      activate(nodes) %>%
      as_tibble()
    
    edge_data <- graph %>%
      activate(edges) %>%
      as_tibble()
    
    module_visible <- (
      module_id %in% colnames(node_data) &&
        any(node_data[[module_id]], na.rm = TRUE)
    ) || (
      module_id %in% colnames(edge_data) &&
        any(edge_data[[module_id]], na.rm = TRUE)
    )
    
    if (module_visible) {
      highlighted_modules <- c(highlighted_modules, module_id)
    }
  }
  
  attr(graph, "highlighted_modules") <- unique(highlighted_modules)
  
  graph
}

combine_module_flags <- function(graph, module_ids, output_column) {
  
  module_flags <- graph %>%
    as_tibble()
  
  module_ids <- intersect(module_ids, colnames(module_flags))
  
  if (length(module_ids) == 0) {
    graph <- graph %>%
      mutate("{output_column}" := FALSE)
    
    return(graph)
  }
  
  module_flags <- module_flags %>%
    dplyr::select(all_of(module_ids))
  
  graph %>%
    mutate("{output_column}" := rowSums(module_flags, na.rm = TRUE) > 0)
}

get_edge_coordinates <- function(edge_data, node_data) {
  
  node_lookup <- node_data %>%
    mutate(node_id = row_number()) %>%
    dplyr::select(node_id, x, y)
  
  edge_data %>%
    dplyr::select(from, to) %>%
    left_join(
      node_lookup %>% rename(from_x = x, from_y = y),
      by = c("from" = "node_id")
    ) %>%
    left_join(
      node_lookup %>% rename(to_x = x, to_y = y),
      by = c("to" = "node_id")
    ) %>%
    filter(
      is.finite(from_x),
      is.finite(from_y),
      is.finite(to_x),
      is.finite(to_y)
    )
}

get_visible_modules <- function(module_ids, node_data, edge_data) {
  
  visible_modules <- vapply(module_ids, function(module_id) {
    
    finite_nodes <- node_data %>%
      filter(is.finite(x), is.finite(y))
    
    node_visible <- any(
      has_module_highlight(finite_nodes, module_id),
      na.rm = TRUE
    )
    
    matched_edges <- edge_data[
      has_module_highlight(edge_data, module_id),
      ,
      drop = FALSE
    ]
    
    if ("drawable_edge" %in% colnames(matched_edges)) {
      matched_edges <- matched_edges[
        matched_edges$drawable_edge %in% TRUE,
        ,
        drop = FALSE
      ]
    }
    
    edge_visible <- matched_edges %>%
      get_edge_coordinates(node_data) %>%
      nrow() > 0
    
    node_visible | edge_visible
  }, logical(1))
  
  module_ids[visible_modules]
}

create_pathway_summary <- function(module_table, visible_modules) {
  
  module_table %>%
    mutate(module_visible = module_id %in% visible_modules) %>%
    group_by(
      direction,
      ID,
      Description,
      NES,
      pvalue,
      p.adjust,
      qvalue,
      kegg_class
    ) %>%
    summarise(
      n_modules          = n_distinct(module_id),
      n_visible_modules  = n_distinct(module_id[module_visible]),
      module_ids         = paste(sort(unique(module_id)), collapse = "; "),
      visible_module_ids = paste(
        sort(unique(module_id[module_visible])),
        collapse = "; "
      ),
      module_names       = paste(sort(unique(module_name)), collapse = "; "),
      used_for_plot      = n_visible_modules > 0,
      .groups            = "drop"
    ) %>%
    arrange(match(direction, c("positive", "negative")), p.adjust, desc(abs(NES)))
}

prepare_metabolic_map <- function(module_table) {
  
  positive_pathways <- module_table %>%
    filter(direction == "positive") %>%
    pull(ID) %>%
    unique()
  
  negative_pathways <- module_table %>%
    filter(direction == "negative") %>%
    pull(ID) %>%
    unique()
  
  positive_modules <- module_table %>%
    filter(direction == "positive") %>%
    pull(module_id) %>%
    unique()
  
  negative_modules <- module_table %>%
    filter(direction == "negative") %>%
    pull(module_id) %>%
    unique()
  
  message("  - Positive NES pathways: ", collapse_ids(positive_pathways))
  message("  - Negative NES pathways: ", collapse_ids(negative_pathways))
  message("  - Positive NES modules: ", collapse_ids(positive_modules))
  message("  - Negative NES modules: ", collapse_ids(negative_modules))
  
  metabolic_map <- pathway("hsa01100") %>%
    process_line()
  
  message("Highlighting KEGG modules on hsa01100...")
  metabolic_map <- add_module_highlights(metabolic_map, module_table)
  highlighted_modules <- attr(metabolic_map, "highlighted_modules")
  
  positive_highlight_modules <- intersect(positive_modules, highlighted_modules)
  negative_highlight_modules <- intersect(negative_modules, highlighted_modules)
  
  if (length(c(positive_highlight_modules, negative_highlight_modules)) == 0) {
    stop("No KEGG modules could be highlighted on hsa01100.")
  }
  
  metabolic_map <- metabolic_map %>%
    activate(nodes) %>%
    combine_module_flags(positive_highlight_modules, "positive_highlight") %>%
    combine_module_flags(negative_highlight_modules, "negative_highlight") %>%
    mutate(
      compound      = convert_id("compound"),
      drawable_node = is.finite(x) & is.finite(y)
    ) %>%
    activate(edges) %>%
    combine_module_flags(positive_highlight_modules, "positive_highlight") %>%
    combine_module_flags(negative_highlight_modules, "negative_highlight")
  
  node_data <- metabolic_map %>%
    activate(nodes) %>%
    as_tibble()
  
  edge_data <- metabolic_map %>%
    activate(edges) %>%
    as_tibble()
  
  from_x <- node_data$x[edge_data$from]
  from_y <- node_data$y[edge_data$from]
  to_x <- node_data$x[edge_data$to]
  to_y <- node_data$y[edge_data$to]
  
  edge_dx <- abs(from_x - to_x)
  edge_dy <- abs(from_y - to_y)
  edge_length <- sqrt(edge_dx^2 + edge_dy^2)
  
  artificial_origin_edge <- is_artificial_origin(from_x, from_y) |
    is_artificial_origin(to_x, to_y)
  
  long_diagonal_edge <- edge_data$type == "line" &
    edge_length > 100 &
    edge_dx > 25 &
    edge_dy > 25
  
  drawable_edge <- is.finite(from_x) &
    is.finite(from_y) &
    is.finite(to_x) &
    is.finite(to_y) &
    !artificial_origin_edge &
    !long_diagonal_edge
  
  drawable_edge[is.na(drawable_edge)] <- FALSE
  
  metabolic_map <- metabolic_map %>%
    activate(edges) %>%
    mutate(drawable_edge = drawable_edge)
  
  edge_data <- metabolic_map %>%
    activate(edges) %>%
    as_tibble()
  
  visible_modules <- get_visible_modules(
    module_ids = unique(c(positive_highlight_modules, negative_highlight_modules)),
    node_data  = node_data,
    edge_data  = edge_data
  )
  
  pathway_summary <- create_pathway_summary(
    module_table     = module_table,
    visible_modules  = visible_modules
  )
  
  if (sum(pathway_summary$used_for_plot) == 0) {
    stop("No significant pathway modules had visible coordinates on hsa01100.")
  }
  
  message(
    "  - Pathways with visible module coordinates: ",
    n_distinct(pathway_summary$ID[pathway_summary$used_for_plot]),
    " / ",
    n_distinct(pathway_summary$ID)
  )
  message("  - Modules highlighted on hsa01100: ", length(visible_modules))
  message(
    "  - Positive NES highlighted nodes/edges: ",
    sum(node_data$positive_highlight, na.rm = TRUE),
    " / ",
    sum(edge_data$positive_highlight, na.rm = TRUE)
  )
  message(
    "  - Negative NES highlighted nodes/edges: ",
    sum(node_data$negative_highlight, na.rm = TRUE),
    " / ",
    sum(edge_data$negative_highlight, na.rm = TRUE)
  )
  
  list(
    metabolic_map = metabolic_map,
    pathway_table = pathway_summary,
    module_table  = module_table,
    node_data     = node_data,
    edge_data     = edge_data
  )
}

################################################################################
# Step 4: Label and Plot Metabolic Map -----------------------------------------
################################################################################

create_pathway_labels <- function(pathway_table, node_data) {
  
  plotted_pathways <- pathway_table %>%
    filter(used_for_plot) %>%
    arrange(match(direction, c("positive", "negative")), p.adjust, desc(abs(NES)))
  
  if (nrow(plotted_pathways) == 0) {
    return(plotted_pathways)
  }
  
  label_data <- lapply(unique(plotted_pathways$ID), function(pathway_id) {
    
    pathway_data <- plotted_pathways %>%
      filter(.data$ID == .env$pathway_id)
    
    pathway_node <- node_data[
      is.finite(node_data$x) &
        is.finite(node_data$y) &
        node_data$type == "map" &
        node_data$name == paste0("path:", pathway_id),
      ,
      drop = FALSE
    ] %>%
      dplyr::select(x, y)
    
    if (nrow(pathway_node) == 0) {
      return(NULL)
    }
    
    data.frame(
      direction   = pathway_data$direction[1],
      ID          = pathway_data$ID[1],
      Description = pathway_data$Description[1],
      NES         = pathway_data$NES[1],
      p.adjust    = pathway_data$p.adjust[1],
      x           = median(pathway_node$x, na.rm = TRUE),
      y           = median(pathway_node$y, na.rm = TRUE),
      label       = wrap_label(pathway_data$Description[1]),
      row.names   = NULL
    )
  })
  
  label_data <- bind_rows(label_data)
  
  if (nrow(label_data) == 0) {
    return(label_data)
  }
  
  label_data %>%
    arrange(match(direction, c("positive", "negative")), p.adjust, desc(abs(NES)))
}

plot_metabolic_map <- function(metabolic_data, output_path) {
  
  metabolic_map <- metabolic_data$metabolic_map
  pathway_table <- metabolic_data$pathway_table
  node_data     <- metabolic_data$node_data
  
  pathway_label_data <- create_pathway_labels(
    pathway_table = pathway_table,
    node_data     = node_data
  ) %>%
    mutate(label_color = ifelse(direction == "positive", "red", "blue"))
  
  node_coordinates <- node_data %>%
    filter(is.finite(x), is.finite(y), !is_artificial_origin(x, y))
  
  plot_x_limits <- range(node_coordinates$x)
  plot_y_limits <- range(node_coordinates$y)
  plot_x_margin <- diff(plot_x_limits) * 0.08
  plot_y_margin <- diff(plot_y_limits) * 0.07
  label_x_limits <- plot_x_limits + c(-plot_x_margin, plot_x_margin)
  label_y_limits <- plot_y_limits + c(-plot_y_margin, plot_y_margin)
  
  message("  - Highlighted pathway labels: ", nrow(pathway_label_data))
  
  plot_obj <- metabolic_map %>%
    ggraph(x = x, y = y) +
    geom_node_point(
      size = 1,
      aes(
        color = I(fgcolor),
        filter = drawable_node & fgcolor != "none" & type != "line"
      )
    ) +
    geom_edge_link(
      width = 0.1,
      aes(
        color = I(fgcolor),
        filter = drawable_edge & type == "line" & fgcolor != "none"
      )
    ) +
    ggfx::with_outer_glow(
      geom_edge_link(
        width = 1,
        aes(
          color = I(fgcolor),
          filter = drawable_edge & fgcolor != "none" & positive_highlight
        )
      ),
      colour = "red",
      expand = 3
    ) +
    ggfx::with_outer_glow(
      geom_node_point(
        size = 2,
        aes(
          color = I(fgcolor),
          filter = drawable_node & fgcolor != "none" & positive_highlight
        )
      ),
      colour = "red",
      expand = 3
    ) +
    ggfx::with_outer_glow(
      geom_edge_link(
        width = 1,
        aes(
          color = I(fgcolor),
          filter = drawable_edge & fgcolor != "none" & negative_highlight
        )
      ),
      colour = "blue",
      expand = 3
    ) +
    ggfx::with_outer_glow(
      geom_node_point(
        size = 2,
        aes(
          color = I(fgcolor),
          filter = drawable_node & fgcolor != "none" & negative_highlight
        )
      ),
      colour = "blue",
      expand = 3
    ) +
    ggrepel::geom_label_repel(
      data = pathway_label_data,
      aes(
        x = x,
        y = y,
        label = label,
        color = I(label_color)
      ),
      inherit.aes        = FALSE,
      size               = 3.2,
      lineheight         = 0.88,
      fontface           = "bold",
      fill               = "#FFFFFFDD",
      label.size         = 0.15,
      label.padding      = grid::unit(0.08, "lines"),
      box.padding        = 0.18,
      point.padding      = 0.05,
      min.segment.length = Inf,
      segment.color      = NA,
      force              = 4,
      force_pull         = 2,
      max.iter           = 100000,
      max.time           = 10,
      max.overlaps       = Inf,
      xlim               = label_x_limits,
      ylim               = label_y_limits,
      seed               = 1
    ) +
    scale_x_continuous(limits = label_x_limits, expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = label_y_limits, expand = expansion(mult = 0.01)) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(10, 10, 10, 10)
    )
  
  ggsave(
    filename = output_path,
    plot     = plot_obj,
    width    = 15,
    height   = 8.5,
    device   = svglite::svglite,
    bg       = "white"
  )
}

################################################################################
# Step 5: Main Execution -------------------------------------------------------
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1 || length(args) > 2) {
  stop(
    "Usage: Rscript src/04_map_kegg_pathways.R ",
    "<gsea_file> [padj_cutoff]"
  )
}

gsea_file <- args[1]
output_dir <- "data"
padj_cutoff <- if (length(args) == 2) as.numeric(args[2]) else 0.05

message("Starting KEGG GSEA metabolic map visualization")
message("GSEA file: ", gsea_file)
message("Output directory: ", output_dir)
message("Adjusted p-value cutoff: ", padj_cutoff)

if (!file.exists(gsea_file)) {
  stop("GSEA file does not exist: ", gsea_file)
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

gsea_name <- file_path_sans_ext(basename(gsea_file))
plot_path <- file.path(
  output_dir,
  paste0(gsea_name, "_map.svg")
)

message("Loading GSEA results...")
gsea_data <- read_tsv(gsea_file, show_col_types = FALSE)

message("Filtering significant pathways...")
significant_pathways <- filter_significant_gsea(
  gsea_data   = gsea_data,
  padj_cutoff = padj_cutoff
)

message("  - Significant pathways: ", nrow(significant_pathways))
message("  - Positive NES: ", sum(significant_pathways$direction == "positive"))
message("  - Negative NES: ", sum(significant_pathways$direction == "negative"))

message("Restricting to KEGG metabolism pathways...")
metabolic_pathways <- filter_metabolic_pathways(significant_pathways)

message("  - Significant metabolic pathways: ", nrow(metabolic_pathways))
message("  - Positive NES: ", sum(metabolic_pathways$direction == "positive"))
message("  - Negative NES: ", sum(metabolic_pathways$direction == "negative"))

message("Identifying KEGG modules under significant metabolic pathways...")
module_table <- get_pathway_modules(metabolic_pathways)

message(
  "  - Metabolic pathways with KEGG modules: ",
  n_distinct(module_table$ID),
  " / ",
  nrow(metabolic_pathways)
)
message("  - KEGG modules identified: ", n_distinct(module_table$module_id))

message("Building and highlighting the human metabolic map...")
metabolic_data <- prepare_metabolic_map(module_table)

message("Plotting metabolic map...")
plot_metabolic_map(
  metabolic_data = metabolic_data,
  output_path    = plot_path
)

message("Saved metabolic map to: ", plot_path)
message("Pipeline complete!")
