## -----------------------------------------------------------------------------
## Script Name: 00_recovery_liver_packages_functions.R
## Project: ANDY recovery liver
## -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------
##  Package Installation Helper 
# ------------------------------------------------------------------------------

# required_packages <- c(
#   "tidyverse", "biomaRt",
#   "SummarizedExperiment", "DESeq2", "apeglm",
#   "ComplexHeatmap",
#   "dendsort", "dendextend",
#   "GeneOverlap", "clusterProfiler",
#   "org.Mm.eg.db",
#   "circlize", "grid",
#   "sva"
# )
#
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
#
# for (pkg in required_packages) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     tryCatch({
#       BiocManager::install(pkg, ask = FALSE, update = TRUE)
#     }, error = function(e) {
#       message("Falling back to install.packages() for ", pkg)
#       install.packages(pkg)
#     })
#   }
# }
#
# rm(pkg, required_packages)



# ------------------------------------------------------------------------------
##  Load Required Libraries
# ------------------------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(SummarizedExperiment)
library(biomaRt)
library(ComplexHeatmap)
library(dendsort)
library(dendextend)
library(GeneOverlap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(circlize)
library(grid)
library(sva)





# ------------------------------------------------------------------------------
##  Function: summarize_DEGs
##  Purpose: Summarize DESeq2 differential expression results by counting
##           significant DEGs (padj ≤ 0.05), and reporting the number that
##           are upregulated and downregulated.
## ------------------------------------------------------------------------------

summarize_DEGs <- function(df, contrast_name = "Contrast") {
  
  # Convert to tibble and carry gene symbols as a column
  df <- df %>% 
    dplyr::as_tibble(rownames = "Symbol") %>%
    dplyr::filter(!is.na(padj))  
  
  # DEG summaries
  total <- df %>% 
    dplyr::filter(padj <= 0.05) %>% 
    nrow()
  
  up <- df %>% 
    dplyr::filter(padj <= 0.05, log2FoldChange > 0) %>% 
    nrow()
  
  down <- df %>% 
    dplyr::filter(padj <= 0.05, log2FoldChange < 0) %>% 
    nrow()
  
  # Return summary table
  tibble::tibble(
    Contrast = contrast_name,
    Total_DEGs = total,
    Upregulated = up,
    Downregulated = down
  )
}




# ------------------------------------------------------------------------------
##  Helper Functions for Multi-Contrast DEG Comparisons
# ------------------------------------------------------------------------------

# Prep DEG data frame: ensure tibble format, Symbol column, and uniqueness
prep_deg_df <- function(x, contrast_name) {
  
  # Convert to tibble with Symbol as rownames if needed
  out <- if (tibble::is_tibble(x)) {
    x
  } else if (inherits(x, "DataFrame")) {
    tibble::as_tibble(as.data.frame(x), rownames = "Symbol")
  } else {
    tibble::as_tibble(x, rownames = "Symbol")
  }
  
  # Ensure Symbol is unique
  if ("Symbol" %in% names(out)) {
    out <- dplyr::distinct(out, Symbol, .keep_all = TRUE)
  }
  
  # Add contrast identifier column
  out <- dplyr::mutate(out, Contrast = contrast_name)
  
  return(out)
}



# ------------------------------------------------------------------------------
##  Compare DEGs Across Multiple Contrasts
##  dfs: list of DEG data frames (must include Symbol and log2FoldChange)
##  contrast_names: vector of contrast names matching dfs order
##  direction:
##     "all" = all DEGs shared (no filtering on direction)
##     "up" = consistently upregulated across contrasts
##     "down" = consistently downregulated across contrasts
##     "opposite" = genes with mixed directions (discordant LFC signs)
# ------------------------------------------------------------------------------

compare_contrasts_multi <- function(dfs, contrast_names,
                                    direction = c("all", "up", "down", "opposite")) {
  
  # Input validation
  stopifnot(length(dfs) >= 2, length(dfs) == length(contrast_names))
  direction <- match.arg(direction)
  
  # Prefix columns for each contrast (except Symbol)
  dfs_named <- Map(
    function(df, nm) {
      df <- tibble::as_tibble(df)
      if (!"Symbol" %in% names(df))
        stop("Each df must contain a 'Symbol' column.")
      dplyr::rename_with(df, .fn = ~ paste0(nm, "_", .x), .cols = -Symbol)
    },
    dfs, contrast_names
  )
  
  # Inner join → retains only shared DEGs
  joined <- Reduce(function(x, y) dplyr::inner_join(x, y, by = "Symbol"), dfs_named)
  
  # Identify relevant log2FC columns
  logfc_cols <- paste0(contrast_names, "_log2FoldChange")
  
  if (!all(logfc_cols %in% names(joined))) {
    missing_cols <- setdiff(logfc_cols, names(joined))
    stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Determine sign consistency
  lfc_mat <- as.matrix(joined[, logfc_cols, drop = FALSE])
  all_up <- apply(lfc_mat, 1, function(v) all(v > 0))
  all_down <- apply(lfc_mat, 1, function(v) all(v < 0))
  opposite <- !(all_up | all_down)
  
  # Apply direction filter
  keep <- switch(
    direction,
    all = rep(TRUE, nrow(joined)),
    up = all_up,
    down = all_down,
    opposite = opposite
  )
  
  return(joined[keep, , drop = FALSE])
}


# ------------------------------------------------------------------------------
##  Identify DEGs Unique to Each Contrast (No Overlap with Others)
##  dfs: A *named list* of DEG data frames or tibbles. Each element must contain a "Symbol" column.
##  Purpose: For each contrast, return the genes that are present only in that contrast
##     and absent from all other contrasts in the list.
##  Output: A combined tibble containing:
##        - unique DEGs
##        - contrast label (Contrast)
##        - all DEG metadata from original tables
# ------------------------------------------------------------------------------

get_unique_DEGs_multi <- function(dfs) {
  
  # Require a named list with ≥ 2 tibbles
  stopifnot(is.list(dfs), length(dfs) >= 2L)
  
  if (is.null(names(dfs)) || any(names(dfs) == "")) {
    stop("Please provide a *named* list of tibbles (names = contrast labels).")
  }
  
  # Get sets of Symbols from each tibble
  symbol_sets <- lapply(dfs, function(df) df$Symbol)
  
  # For each contrast, keep genes not present in any other contrast
  unique_lists <- lapply(names(dfs), function(nm) {
    
    this_df <- dfs[[nm]]
    other_syms <- unlist(symbol_sets[names(symbol_sets) != nm], use.names = FALSE)
    uniq_syms  <- setdiff(this_df$Symbol, unique(other_syms))
    
    # Filter to unique symbols and add Contrast column
    this_df %>%
      dplyr::filter(Symbol %in% uniq_syms) %>%
      dplyr::mutate(Contrast = nm, .before = 1)
  })
  
  dplyr::bind_rows(unique_lists)
}



# ------------------------------------------------------------------------------
##  Summarize Unique DEGs Per Contrast
##  unique_df: Output from get_unique_DEGs_multi().
##  Purpose: Generate a tidy summary table reporting:
##        - Total number of unique DEGs
##        - Upregulated unique DEGs
##        - Downregulated unique DEGs
##  Output: A tibble with one row per contrast, sorted by total unique DEGs.
# ------------------------------------------------------------------------------

summarize_unique_DEGs_multi <- function(unique_df) {
  
  unique_df %>%
    dplyr::group_by(Contrast) %>%
    dplyr::summarise(
      Total         = dplyr::n(),
      Upregulated   = sum(log2FoldChange > 0, na.rm = TRUE),
      Downregulated = sum(log2FoldChange < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(Total))
}



# ------------------------------------------------------------------------------
##  Dotplot Function for Single-Contrast clusterProfiler ORA Results (GO)
##
##  Purpose:
##     Generate a regulation-aware dotplot for a single contrast, using
##     clusterProfiler::compareCluster output (e.g., GO MF results).
##
##  Arguments:
##     cp_res           compareCluster result object (class compareClusterResult)
##     show_path_number Number of enriched pathways to display (default = 30)
##     plot_title       Optional plot title
##
##  Notes:
##     - Removes species suffix from pathway names.
##     - Orders pathways within Up/Down groups by adjusted p-value, then GeneRatio.
##     - Produces dotplots consistent with manuscript figures.
# ------------------------------------------------------------------------------

clusterProfiler_ora_cc_single_dotplot <- function(cp_res,
                                                  show_path_number = 30,
                                                  plot_title = NULL) {

  # Remove species suffix from pathway descriptions
  cp_res@compareClusterResult$Description <- stringr::str_replace(
    cp_res@compareClusterResult$Description,
    " - Mus musculus \\(house mouse\\)", ""
  )
  
  # Reorder pathways by significance and GeneRatio within each Regulation group
  cp_res@compareClusterResult <- cp_res@compareClusterResult %>%
    dplyr::group_by(Regulation) %>%
    dplyr::arrange(p.adjust, dplyr::desc(GeneRatio), .by_group = TRUE) %>%
    dplyr::ungroup()
  
  # Generate dotplot
  clusterProfiler::dotplot(
    cp_res,
    x = "Regulation",
    showCategory = show_path_number,
    label_format = function(x) stringr::str_wrap(x, width = 25)
  ) +
    scale_x_discrete(labels = c("Upregulated" = "Up", "Downregulated" = "Down")) +
    scale_fill_continuous(high = "blue", low = "red") +
    guides(
      size = guide_legend(order = 1, title = "Gene Ratio"),
      fill = guide_colorbar(order = 2, title = "Adj. p-value")
    ) +
    ggtitle(plot_title) +
    theme(
      plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(face = "bold", size = 12),
      axis.text.y  = element_text(face = "bold", size = 9),
      strip.text   = element_text(size = 12)
    )
}



# ------------------------------------------------------------------------------
##  Add Gene Symbols to clusterProfiler compareCluster Results
##
##  Arguments:
##     cp_res   A compareClusterResult object from clusterProfiler
##     OrgDb    Annotation database (default = org.Mm.eg.db)
##
##  Output:
##     A tibble replicating cp_res@compareClusterResult but with an
##     added `geneSymbol` column for downstream pathway interpretation.
# ------------------------------------------------------------------------------

add_symbols_to_compareCluster <- function(cp_res, OrgDb = org.Mm.eg.db) {
  
  # Extract enrichment results
  df <- cp_res@compareClusterResult
  
  # Expand geneID column
  df_expanded <- df %>%
    tidyr::separate_rows(geneID, sep = "/")
  
  # Map Entrez to Symbol
  entrez2symbol <- clusterProfiler::bitr(
    df_expanded$geneID,
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = OrgDb
  ) %>%
    dplyr::distinct()
  
  # Add SYMBOL column
  df_expanded <- df_expanded %>%
    dplyr::left_join(entrez2symbol, by = c("geneID" = "ENTREZID"))
  
  # Collapse back to pathway-level
  df_with_symbols <- df_expanded %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(
      geneID = paste(geneID, collapse = "/"),
      geneSymbol = paste(na.omit(SYMBOL), collapse = "/"),
      .groups = "drop"
    ) %>%
    dplyr::right_join(df, by = "ID") %>%
    
    # Keep only the correct geneID column
    dplyr::rename(geneID = geneID.x) %>%
    dplyr::select(-geneID.y) %>%
    
    # Move geneID and geneSymbol to the END of the table
    dplyr::relocate(geneID, geneSymbol, .after = dplyr::last_col())
  
  return(df_with_symbols)
}


# ------------------------------------------------------------------------------
## Summarize the clusterProfiler results
# ------------------------------------------------------------------------------
summarize_ora_results <- function(df) {
  df %>%
    dplyr::distinct(ID, .keep_all = TRUE) %>%   # ensure each pathway counted once
    dplyr::summarise(
      total_pathways = dplyr::n(),
      upregulated = sum(Regulation == "Upregulated", na.rm = TRUE),
      downregulated  = sum(Regulation == "Downregulated", na.rm = TRUE)
    )
}




# ------------------------------------------------------------------------------
##  Dotplot for Shared Enriched Pathways Across Multiple Contrasts
##
##  Arguments:
##     cp_df            A merged ORA results data frame (from multiple contrasts)
##     show_path_number Number of pathways to show per facet (default = 30)
##
##  Output:
##     A faceted multi-contrast dotplot.
# ------------------------------------------------------------------------------

dotplot_from_df_multi_contrasts <- function(cp_df, show_path_number = 30, names_list = NULL, plot_title = NULL) {
  # Clean up pathway names (optional: remove species tag)
  cp_df$Description <- stringr::str_replace(
    cp_df$Description,
    " - Mus musculus \\(house mouse\\)", ""
  )

  # Reorder by Adj P-value and then by GeneRatio within each Contrast & Regulation ----
  cp_df <- cp_df %>%
    dplyr::group_by(Contrast, Regulation) %>%
    dplyr::arrange(p.adjust, dplyr::desc(GeneRatio), .by_group = TRUE) %>%
    dplyr::ungroup()

  # Convert names_list into a labeller using ggplot2's as_labeller
  labeller_fn <- if (!is.null(names_list)) ggplot2::as_labeller(names_list) else "label_value"

  # Create dotplot using compareClusterResult class
  clusterProfiler::dotplot(
    new("compareClusterResult", compareClusterResult = cp_df),
    x = "Regulation",
    showCategory = show_path_number,
    label_format = function(x) stringr::str_wrap(x, width = 25)
  ) +
    ggplot2::facet_grid(~Contrast, labeller = labeller_fn) +
    ggplot2::scale_x_discrete(labels = c("Upregulated" = "Up", "Downregulated" = "Down")) +
    ggplot2::scale_fill_continuous(high = "blue", low = "red") +
    ggplot2::guides(
      size = ggplot2::guide_legend(order = 1, title = "Gene Ratio"),
      fill = ggplot2::guide_colorbar(order = 2, title = "Adj. p-value")
    ) +
    ggtitle(plot_title) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(face = "bold", size = 10),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      strip.text = ggplot2::element_text(size = 12)
    )
}





# ------------------------------------------------------------------------------
##  Summarize Shared Pathways Across Two Contrasts
##  Arguments:
##     df1, df2          Data frames of ORA results containing:
##                       ID, Regulation, and other clusterProfiler fields
##     contrast1_name    Name of the first contrast (string)
##     contrast2_name    Name of the second contrast (string)
##
##  Output:
##     A tibble summarizing:
##        - Shared Upregulated pathways
##        - Shared Downregulated pathways
##        - Total Shared pathways
##        - Shared pathways with Opposite regulation
# ------------------------------------------------------------------------------

summarize_shared_pathways <- function(df1, df2,
                                      contrast1_name = "Contrast1",
                                      contrast2_name = "Contrast2") {
  
  # Select pathway ID and regulation
  df1_ids <- df1 %>% dplyr::select(ID, Regulation)
  df2_ids <- df2 %>% dplyr::select(ID, Regulation)
  
  # Shared pathways with SAME regulation direction only
  shared <- df1_ids %>%
    dplyr::inner_join(df2_ids, by = c("ID", "Regulation"))
  
  # Count shared Up / Down
  shared_summary <- shared %>%
    dplyr::count(Regulation, name = "SharedCount") %>%
    dplyr::mutate(
      Contrast1 = contrast1_name,
      Contrast2 = contrast2_name
    )
  
  # Total shared (same ID + same regulation)
  total_shared <- tibble::tibble(
    Regulation  = "Total Shared",
    SharedCount = nrow(shared),
    Contrast1 = contrast1_name,
    Contrast2 = contrast2_name
  )
  
  # Combine summaries
  dplyr::bind_rows(shared_summary, total_shared)
}



# ------------------------------------------------------------------------------
##  Identify Unique clusterProfiler Pathways Across Multiple Contrasts
##
##  Arguments:
##     dfs        A *named* list of clusterProfiler ORA result tibbles.
##                Names must correspond to contrasts (e.g., NFDvsCD, RecvsCD).
##     key_cols   Columns that define pathway identity. Default: c("ID", "Regulation")
##
##  Output:
##     A tibble containing all original annotation columns for pathways that
##     are unique to exactly one contrast. Includes a Contrast column.
# ------------------------------------------------------------------------------

unique_cp_pathways <- function(dfs, key_cols = c("ID", "Regulation")) {
  
  # Input validation
  stopifnot(is.list(dfs), length(dfs) >= 2L)
  if (is.null(names(dfs)) || any(names(dfs) == "")) {
    stop("Provide a *named* list: e.g., list(NFDvsCD=..., RecvsCD=...).")
  }
  
  # Validate each df has key columns
  invisible(lapply(names(dfs), function(nm) {
    df <- dfs[[nm]]
    if (!tibble::is_tibble(df)) {
      stop(sprintf("Input '%s' must be a tibble.", nm))
    }
    missing <- setdiff(key_cols, names(df))
    if (length(missing) > 0) {
      stop(sprintf("Input '%s' is missing columns: %s",
                   nm, paste(missing, collapse = ", ")))
    }
    TRUE
  }))
  
  # Bind all inputs and tag with Contrast label
  bound_full <- purrr::imap_dfr(
    dfs,
    ~ dplyr::mutate(.x, Contrast = .y, .before = 1)
  )
  
  # Each contrast × pathway identity combination (unique)
  key_per_contrast <- bound_full %>%
    dplyr::distinct(Contrast, !!!rlang::syms(key_cols))
  
  # Count in how many contrasts each pathway identity appears
  key_counts <- key_per_contrast %>%
    dplyr::count(!!!rlang::syms(key_cols), name = "n_contrasts")
  
  # Keep pathways appearing in exactly one contrast
  unique_keys <- key_counts %>%
    dplyr::filter(n_contrasts == 1L) %>%
    dplyr::select(!!!rlang::syms(key_cols))
  
  # Return full annotation rows for unique pathways
  bound_full %>%
    dplyr::semi_join(unique_keys, by = key_cols)
}


# ------------------------------------------------------------------------------
##  Summarize Unique clusterProfiler Pathways per Contrast
##
##  Arguments:
##     unique_df   Output from unique_cp_pathways(); must contain
##                 columns: Contrast, ID, Regulation
##
##  Output:
##     A tidy tibble with rows for:
##        Contrast × (Upregulated / Downregulated / Total)
# ------------------------------------------------------------------------------

summarize_unique_cp_pathways <- function(unique_df) {
  
  if (!("Contrast" %in% names(unique_df))) stop("Expected a 'Contrast' column.")
  if (!("Regulation" %in% names(unique_df))) stop("Expected a 'Regulation' column.")
  if (!("ID" %in% names(unique_df))) stop("Expected an 'ID' column.")
  
  purrr::map_dfr(unique(unique_df$Contrast), function(nm) {
    
    df <- unique_df %>% dplyr::filter(Contrast == nm)
    
    # Count Up/Down per contrast
    reg_counts <- df %>%
      dplyr::distinct(ID, Regulation) %>%
      dplyr::count(Regulation, name = "Count") %>%
      dplyr::mutate(Contrast = nm)
    
    # Total unique pathways per contrast
    total_count <- df %>%
      dplyr::distinct(ID, Regulation) %>%
      dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(Regulation = "Total", Contrast = nm)
    
    # Combine and format
    dplyr::bind_rows(reg_counts, total_count) %>%
      dplyr::rename(Pathways = Regulation) %>%
      dplyr::select(Contrast, Pathways, Count)
  })
}



# ------------------------------------------------------------------------------
##  Dotplot for Unique Enriched Pathways (Single Contrast)
##
##  Arguments:
##     cp_df            Tibble of unique pathway results; must contain:
##                      ID, Description, Regulation, p.adjust, GeneRatio, etc.
##     show_path_number Number of pathways to display (default = 30)
##     plot_title       Optional plot title
##
##  Output:
##     A ggplot2 dotplot showing Up/Down unique pathways for a contrast.
# ------------------------------------------------------------------------------

clusterProfiler_ora_cc_unique_dotplot <- function(cp_df,
                                                  show_path_number = 30,
                                                  plot_title = NULL) {

  # Ensure tibble becomes plain data.frame
  cp_df <- as.data.frame(cp_df)

  # Clean pathway Description names
  cp_df$Description <- stringr::str_replace(
    cp_df$Description,
    " - Mus musculus \\(house mouse\\)",
    ""
  )

  # Reorder by significance and GeneRatio within each regulation group
  cp_df <- cp_df %>%
    dplyr::group_by(Regulation) %>%
    dplyr::arrange(p.adjust, dplyr::desc(GeneRatio), .by_group = TRUE) %>%
    dplyr::ungroup()

  # Convert into compareClusterResult object
  cp_res <- methods::new(
    "compareClusterResult",
    compareClusterResult = cp_df
  )

  # Dotplot
  clusterProfiler::dotplot(
    cp_res,
    x = "Regulation",
    showCategory = show_path_number,
    label_format = function(x) stringr::str_wrap(x, width = 25)
  ) +
    scale_x_discrete(
      labels = c("Upregulated" = "Up",
                                "Downregulated" = "Down")) +
    scale_fill_continuous(
      high = "blue", low = "red",
      labels = scales::label_number(accuracy = 0.001)) +

    scale_size_continuous(
      labels = scales::label_number(accuracy = 0.001)
    ) +

    guides(
      size = guide_legend(order = 1, title = "Gene Ratio"),
      fill = guide_colorbar(order = 2, title = "Adj. p-value")
    ) +
    ggtitle(plot_title) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(face = "bold", size = 12),
      axis.text.y  = element_text(face = "bold", size = 10),
      strip.text   = element_text(size = 12)
    )
}