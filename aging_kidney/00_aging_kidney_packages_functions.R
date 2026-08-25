## -----------------------------------------------------------------------------
## Script Name: 00_aging_kidney_packages_functions.R
## Project: ANDY aging kidney
## -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------
##  Package Installation Helper
# ------------------------------------------------------------------------------

# required_packages <- c(
#   "tidyverse", "biomaRt",
#   "SummarizedExperiment", "DESeq2", "apeglm",
#   "ComplexHeatmap",
#   "dendsort", "dendextend",
#   "GeneOverlap",
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
