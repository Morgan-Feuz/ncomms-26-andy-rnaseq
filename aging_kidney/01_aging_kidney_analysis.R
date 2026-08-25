## -----------------------------------------------------------------------------
## Script Name: 01_aging_kidney_analysis.R
## Project: Bulk RNA-seq Analysis of ANDY CD/NFD and young/old WT mice
## Study: NFD vs CD and Old vs Young DEG-level comparison
## Analysis Type: Set-up, differential expression analysis, heatmap
##                visualization
## Tissues: Kidney
## Dependencies: Requires 00_aging_kidney_packages_functions.R script.
##               The NFD and CD samples used in this script come from the
##               recovery cohort (a different sequencing batch from the
##               Old/Young samples analyzed here) -- run
##               recovery_kidney/01_recovery_kidney_analysis.R first, in the
##               same R session, so its DESeq2 results and normalized counts
##               are available before continuing.
## Last Updated: 2026-08-25
## -----------------------------------------------------------------------------


# ------------------------------------------------------------------------------
##  Data Preprocessing
# ------------------------------------------------------------------------------

# Read in raw counts matrix from featureCounts output
counts_file <- "path/to/counts_matrix.csv"  # GEO-deposited counts matrix
counts_data_old_young_kidney <- read.csv(counts_file)

# Contains 57,132 genes across 12 samples
dim(counts_data_old_young_kidney)



# ------------------------------------------------------------------------------
##  Map Ensembl Gene IDs to Gene Symbols (biomaRt)
# ------------------------------------------------------------------------------

# Extract Ensembl gene IDs (with version)
ensembl.ids <- counts_data_old_young_kidney$Geneid
head(ensembl.ids)

# Load Ensembl BioMart (version 112)
ensembl <- useEnsembl(
  biomart = "genes",
  version = "112",
  dataset = "mmusculus_gene_ensembl"
)

# Retrieve mapping to external gene names (symbols)
symbol <- getBM(
  attributes = c("ensembl_gene_id_version", "external_gene_name"),
  filters = "ensembl_gene_id_version",
  values = ensembl.ids,
  mart = ensembl
)

# Merge gene symbols into count matrix
counts_data_old_young_kidney <- counts_data_old_young_kidney %>%
  dplyr::rename(ensembl_gene_id_version = Geneid) %>%
  dplyr::left_join(symbol, by = "ensembl_gene_id_version") %>%
  dplyr::filter(!(is.na(external_gene_name) | external_gene_name == ""))

# 56,648 genes retained after filtering out rows without symbols
dim(counts_data_old_young_kidney)



# ------------------------------------------------------------------------------
##  Handle Duplicated Gene Symbols
# ------------------------------------------------------------------------------

# Identify duplicated gene symbols
duplicated_symbols <- counts_data_old_young_kidney$external_gene_name[
  duplicated(counts_data_old_young_kidney$external_gene_name)
] %>% sort()

length(duplicated_symbols)  # 235 duplicated symbols

# Inspect duplicated symbols
dup_symbols_df <- counts_data_old_young_kidney %>%
  dplyr::filter(external_gene_name %in% duplicated_symbols)

# Remove duplicates (keep first occurrence)
counts_data_old_young_kidney <- counts_data_old_young_kidney %>%
  dplyr::distinct(external_gene_name, .keep_all = TRUE)

# 56,413 genes x 14 columns after deduplication
dim(counts_data_old_young_kidney)



# ------------------------------------------------------------------------------
##  Prepare Final Count Matrix
# ------------------------------------------------------------------------------

# Set gene symbols as rownames
rownames(counts_data_old_young_kidney) <- counts_data_old_young_kidney$external_gene_name

# Remove Ensembl ID and gene symbol columns
counts_data_old_young_kidney <- counts_data_old_young_kidney[, -c(1, 14)]

# Result should be 56,413 genes x 12 samples
dim(counts_data_old_young_kidney)



# ------------------------------------------------------------------------------
##  Load and Align Sample Metadata
# ------------------------------------------------------------------------------

# Read metadata (sample info table)
sampleinfo_file <- "path/to/sample_info.csv"  # sample metadata, deposited with GEO submission
colData_old_young_kidney <- read.csv(
  sampleinfo_file,
  header = TRUE
)

head(colData_old_young_kidney)

# Move sample names into rownames
rownames(colData_old_young_kidney) <- colData_old_young_kidney[, 1]

# Ensure sample names match between metadata and count matrix
# Required for DESeq2
all(colnames(counts_data_old_young_kidney) %in% rownames(colData_old_young_kidney)) # TRUE
all(colnames(counts_data_old_young_kidney) == rownames(colData_old_young_kidney)) # TRUE 

# Distribution across conditions
table(colData_old_young_kidney$Treatment)
# Old: 6 Young: 6



# ------------------------------------------------------------------------------
##  Convert Counts to Matrix for DESeq2
# ------------------------------------------------------------------------------

counts_data_matrix_old_young_kidney <- data.matrix(counts_data_old_young_kidney)

head(counts_data_matrix_old_young_kidney)


# ------------------------------------------------------------------------------
##  Create a SummarizedExperiment Object
# ------------------------------------------------------------------------------

# Construct SummarizedExperiment for the NFD recovery liver samples
se_old_young_kidney <- SummarizedExperiment(
  assays = list(counts = counts_data_matrix_old_young_kidney),
  colData = colData_old_young_kidney
)

se_old_young_kidney


# ------------------------------------------------------------------------------
##  Prefiltering
# ------------------------------------------------------------------------------

# Create a DESeqDataSet for the NFD recovery study
dds_old_young_kidney <- DESeqDataSet(
  se_old_young_kidney, 
  design = ~ Treatment)

# Confirm initial gene count (expected 56,413)
nrow(dds_old_young_kidney)



# ------------------------------------------------------------------------------
##  Prefilter Lowly Expressed Genes
# ------------------------------------------------------------------------------

# Determine smallest group size (used for filtering threshold)
smallestGroupSize_old_young_kidney <- min(table(colData(dds_old_young_kidney)$Treatment))
print(smallestGroupSize_old_young_kidney) # 6

# Keep genes with count ≥ 10 in at least `smallestGroupSize` samples
keep_old_young_kidney <- rowSums(counts(dds_old_young_kidney) >= 10) >= smallestGroupSize_old_young_kidney
dds_old_young_kidney <- dds_old_young_kidney[keep_old_young_kidney, ]

# 16,558 genes should remain
nrow(dds_old_young_kidney)



# ------------------------------------------------------------------------------
##  Estimate Size Factors (Normalization)
# ------------------------------------------------------------------------------

# Normalize for sequencing depth differences
dds_old_young_kidney <- estimateSizeFactors(dds_old_young_kidney)



# ------------------------------------------------------------------------------
##  Differential Expression Analysis with DESeq2
# ------------------------------------------------------------------------------

# Set reference level ("Young") and run DESeq2 pipeline
# ------------------------------------------------------------------------------
dds_old_young_kidney$Treatment <- relevel(
  dds_old_young_kidney$Treatment, 
  ref = "Young")

dds_old_young_kidney <- DESeq(dds_old_young_kidney)

# Review size factors 
sizeFactors(dds_old_young_kidney)


# Inspect model coefficients / contrast names
resultsNames(dds_old_young_kidney)
# Treatment_Old_vs_Young



# ------------------------------------------------------------------------------
##  Contrast: Old vs Young (kidney)
# ------------------------------------------------------------------------------

# Extract DESeq2 results
Old_Young_kidney <- results(
  dds_old_young_kidney, 
  contrast = c("Treatment", "Old", "Young"))

# Apply LFC shrinkage (apeglm)
Old_Young_shrunk_kidney <- lfcShrink(
  dds_old_young_kidney,
  coef = "Treatment_Old_vs_Young",
  res  = Old_Young_kidney,
  type = "apeglm"
) %>% 
  as.data.frame()

# Identify significant DEGs (FDR ≤ 0.05)
Old_Young_adjp_kidney <- subset(Old_Young_shrunk_kidney, padj <= 0.05)

# Summary table using summarize_DEGs function
summary_Old_Young_kidney <- summarize_DEGs(
  Old_Young_shrunk_kidney, 
  contrast_name = "Old vs Young")

summary_Old_Young_kidney
# Total_DEGs = 8260
# Upregulated = 4233
# Downregulated = 4027




# ------------------------------------------------------------------------------
## ComplexHeatmap Visualization of Top DEGs
# ------------------------------------------------------------------------------

# Extract top DEGs (Old vs Young, sorted by adjusted p-value)
# ------------------------------------------------------------------------------
# Sort by padj and select top 1000 most significant DEGs (Old vs Young)
Old_Young_adjp_kidney <- Old_Young_adjp_kidney[order(Old_Young_adjp_kidney$padj), ]
Old_Young_adjp_kidney_top1000 <- Old_Young_adjp_kidney[1:1000, ]

dim(Old_Young_adjp_kidney_top1000)  # Should be 1000 × 5


# Subset normalized expression matrix for those genes
# ------------------------------------------------------------------------------
# Retrieve the normalized counts
normalized_counts_old_young_kidney <- counts(dds_old_young_kidney, normalized = TRUE)
head(normalized_counts_old_young_kidney)

top_genes_old_young <- rownames(Old_Young_adjp_kidney_top1000)

norm_cts_sub_old_young <- normalized_counts_old_young_kidney[rownames(
  normalized_counts_old_young_kidney) %in% top_genes_old_young, ]

dim(norm_cts_sub_old_young) # 1000 x 12

# Scale expression matrix (row-wise z-scores for heatmap)
heat_old_young <- t(scale(t(norm_cts_sub_old_young)))


# ------------------------------------------------------------------------------
# Define annotations and colors
# ------------------------------------------------------------------------------

# Rename treatment labels for clarity
condition_labels_old_young <- colData_old_young_kidney$Treatment
condition_labels_old_young <- factor(
  condition_labels_old_young, 
  levels = c("Young", "Old"),
  labels = c("Young","Old")
  )

# Column annotations
column_ha_old_young <- HeatmapAnnotation(
  Condition = condition_labels_old_young,
  annotation_name_side = "right",
  show_annotation_name = FALSE,
  col = list(Condition = c(
    "Old" = '#808080',
    "Young" = '#ff6200'
  ))
)

# ------------------------------------------------------------------------------
# Rotate column dendrogram for visual grouping
# ------------------------------------------------------------------------------

# Build dendrogram
dend_old_young <- hclust(dist(t(heat_old_young))) %>% as.dendrogram()

# Plot the heatmap
# ------------------------------------------------------------------------------
hmap_old_young <- Heatmap(
  heat_old_young,
  name = 'Gene\nZ-Score',
  top_annotation = column_ha_old_young,
  cluster_rows = TRUE,
  cluster_columns = dend_old_young,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  row_dend_reorder = TRUE,
  column_dend_reorder = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 12)
)

draw(
  hmap_old_young,
  heatmap_legend_side = 'right',
  annotation_legend_side = 'right',
  merge_legend = TRUE
)





# ------------------------------------------------------------------------------
# Heatmap using top 1000 DEGs from NFD vs CD (recovery study) contrast
# ------------------------------------------------------------------------------
# Note: The NFD and CD samples analyzed in this section come from the recovery
# cohort -- a separate sequencing batch from the Old/Young samples analyzed
# above. Rather than re-reading exported files, this section reuses the
# DESeq2 results and normalized counts already computed by
# recovery_kidney/01_recovery_kidney_analysis.R (NFD_CD_adjp_kidney,
# NFD_CD_shrunk_kidney, colData_rec_main_kidney, normalized_counts_main_rec_kidney).
# Run that script first, in the same R session, before continuing here.

length(NFD_CD_adjp_kidney$Symbol) # 2957

# Build sample metadata with an explicit Sample_ID column
# ------------------------------------------------------------------------------
nfd_cd_colData_kidney <- colData_rec_main_kidney %>%
  tibble::rownames_to_column("Sample_ID")

head(nfd_cd_colData_kidney)

# Normalized counts matrix (genes x samples) from the recovery kidney analysis
# ------------------------------------------------------------------------------
norm_cts_nfd_cd_mat_kidney <- normalized_counts_main_rec_kidney

# Check result
norm_cts_nfd_cd_mat_kidney[1:5, 1:5]

# Keep only CD and NFD samples
nfd_cd_filtered <- nfd_cd_colData_kidney %>%
  dplyr::filter(Condition %in% c("CD", "NFD"))

# Match sample IDs in the same order as the matrix columns
keep_samples_nfd_cd <- colnames(norm_cts_nfd_cd_mat_kidney)[
  colnames(norm_cts_nfd_cd_mat_kidney) %in% nfd_cd_filtered$Sample_ID
]

# Subset counts matrix
norm_cts_nfd_cd_mat_kidney <- norm_cts_nfd_cd_mat_kidney[, keep_samples_nfd_cd]

# Scale after subsetting
heat_norm_cts_nfd_cd_mat <- t(scale(t(norm_cts_nfd_cd_mat_kidney)))

# Subset colData in the same order
nfd_cd_filtered <- nfd_cd_filtered %>%
  dplyr::filter(Sample_ID %in% keep_samples_nfd_cd) %>%
  dplyr::arrange(match(Sample_ID, keep_samples_nfd_cd))

# Check dimensions
dim(nfd_cd_filtered) # 10 x 4
dim(heat_norm_cts_nfd_cd_mat) # 16305 x 10


# ------------------------------------------------------------------------------
# Extract top DEGs (NFD vs CD, sorted by adjusted p-value)
# ------------------------------------------------------------------------------

# Sort the adjusted p-values
NFD_CD_adjp_kidney <- NFD_CD_adjp_kidney[order(NFD_CD_adjp_kidney$padj), ]
NFD_CD_adjp_top1000 <- NFD_CD_adjp_kidney[1:1000, ]

dim(NFD_CD_adjp_top1000)  # Should be 1000 × 6

# Subset the heatmap matrix using the top 1000 DEGs from the NFD vs CD contrast
top_genes_nfd_cd <- NFD_CD_adjp_top1000$Symbol

heat_norm_cts_nfd_cd_mat_sub <- heat_norm_cts_nfd_cd_mat[
  rownames(heat_norm_cts_nfd_cd_mat) %in% top_genes_nfd_cd, ]

dim(heat_norm_cts_nfd_cd_mat_sub) # e.g., 1000 x 10


# Build condition labels
condition_labels_nfd_cd <- factor(
  nfd_cd_filtered$Condition,
  levels = c("CD", "NFD")
)

column_ha_nfd_cd <- HeatmapAnnotation(
  Condition = condition_labels_nfd_cd,
  annotation_name_side = "right",
  show_annotation_name = FALSE,
  col = list(Condition = c(
    "CD"  = '#000080',  
    "NFD" = '#a00000'  
  ))
)


# Cluster columns (samples)
dend_nfd_cd <- hclust(dist(t(heat_norm_cts_nfd_cd_mat_sub))) %>% as.dendrogram()


# Build NFD vs CD heatmap
hmap_nfd_cd <- Heatmap(
  heat_norm_cts_nfd_cd_mat_sub,
  name = 'Gene\nZ-Score',
  top_annotation = column_ha_nfd_cd,
  cluster_rows = TRUE,
  cluster_columns = dend_nfd_cd,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
  row_dend_reorder = TRUE,
  column_dend_reorder = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 12)
)

# Draw the heatmap
draw(
  hmap_nfd_cd,
  heatmap_legend_side = 'right',
  annotation_legend_side = 'right',
  merge_legend = TRUE
)



# ------------------------------------------------------------------------------
# Top 1000 DEGs from Old vs Young plotted in NFD vs CD heatmap
# ------------------------------------------------------------------------------

# Draw Old vs Young heatmap to initialize clustering
# (This ensures row_order() reflects the actual displayed clustering)
hmap_old_young <- draw(hmap_old_young)

# Extract the clustered gene order from Old vs Young heatmap
row_order_old_young <- unlist(ComplexHeatmap::row_order(hmap_old_young))
gene_order_old_young <- rownames(heat_old_young)[row_order_old_young]

# Subset NFD vs CD normalized counts using that same gene order
norm_cts_sub_nfd_cd_from_old_young <- norm_cts_nfd_cd_mat_kidney[
  gene_order_old_young[gene_order_old_young %in% rownames(norm_cts_nfd_cd_mat_kidney)],
]

dim(norm_cts_sub_nfd_cd_from_old_young)  # e.g., 987 x 10


# Match sample IDs in the same order as metadata
keep_samples <- colnames(norm_cts_sub_nfd_cd_from_old_young)[
  colnames(norm_cts_sub_nfd_cd_from_old_young) %in% nfd_cd_filtered$Sample_ID
]

# Subset counts matrix
norm_cts_sub_nfd_cd_from_old_young <- norm_cts_sub_nfd_cd_from_old_young[, keep_samples]

# Scale expression values (Z-score per gene)
heat_nfd_cd_from_old_young <- t(scale(t(norm_cts_sub_nfd_cd_from_old_young)))

# Subset and order metadata to match matrix columns
nfd_cd_filtered <- nfd_cd_filtered %>%
  dplyr::filter(Sample_ID %in% keep_samples) %>%
  dplyr::arrange(match(Sample_ID, keep_samples))

# Build condition labels and annotation
condition_labels_nfd_cd <- factor(
  nfd_cd_filtered$Condition,
  levels = c("CD", "NFD")
)

# Set heatmap column annotations
column_ha_nfd_cd <- HeatmapAnnotation(
  Condition = condition_labels_nfd_cd,
  annotation_name_side = "right",
  show_annotation_name = FALSE,
  col = list(
    Condition = c(
      "CD"  = "#000080",  
      "NFD" = "#a00000" 
    )
  )
)

# Cluster columns (samples) but keep gene rows fixed
dend_nfd_cd <- hclust(dist(t(heat_nfd_cd_from_old_young))) %>% as.dendrogram()
dend_nfd_cd <- dendextend::rotate(dend_nfd_cd, order = c(5:10, 1:4))  

# Ensure gene order only includes genes present in this matrix
gene_order_old_young <- gene_order_old_young[gene_order_old_young %in% rownames(heat_nfd_cd_from_old_young)]


# ------------------------------------------------------------------------------
# Plot heatmap — same gene order as Old vs Young
# ------------------------------------------------------------------------------

# Define heatmap params
hmap_nfd_cd_from_old_young <- Heatmap(
  heat_nfd_cd_from_old_young,
  name = "Gene\nZ-Score",
  top_annotation = column_ha_nfd_cd,
  cluster_rows = FALSE,                  # don't recluster — use given order
  row_order = gene_order_old_young,      # enforce Old vs Young gene order
  cluster_columns = dend_nfd_cd,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_dend_reorder = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 12)
)

# Draw the heatmap
draw(
  hmap_nfd_cd_from_old_young,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)



# Sanity checks ---
# How many genes are in each heatmap
nrow(heat_old_young) # 1000
nrow(heat_nfd_cd_from_old_young) # 987


# Check overlap of gene sets
sum(rownames(norm_cts_sub_nfd_cd_from_old_young) %in% gene_order_old_young) # 987
# Should equal nrow(norm_cts_sub_nfd_cd_from_old_young)

# Check sample matching
all(colnames(norm_cts_sub_nfd_cd_from_old_young) == nfd_cd_filtered$Sample_ID) # TRUE
# Should return TRUE




# ------------------------------------------------------------------------------
# Top 1000 DEGs from NFD vs CD  plotted in Old vs Young Subset heatmap
# ------------------------------------------------------------------------------

# Draw once to initialize clustering
hmap_nfd_cd <- draw(hmap_nfd_cd)

# Extract clustered gene order from NFD vs CD heatmap
row_order_nfd_cd <- unlist(ComplexHeatmap::row_order(hmap_nfd_cd))
gene_order_nfd_cd <- rownames(heat_norm_cts_nfd_cd_mat_sub)[row_order_nfd_cd]



# Use NFD vs CD heatmap gene order to plot the Old vs Young dataset
# ------------------------------------------------------------------------------

# Create Old vs Young metadata from DESeq2 object
old_young_filtered <- as.data.frame(SummarizedExperiment::colData(dds_old_young_kidney)) %>%
  dplyr::filter(Treatment %in% c("Old", "Young")) %>%
  dplyr::arrange(Treatment)


# Subset normalized counts to same genes, preserving NFD vs CD gene order
norm_cts_sub_old_young_from_nfd_cd <- normalized_counts_old_young_kidney[
  gene_order_nfd_cd[gene_order_nfd_cd %in% rownames(normalized_counts_old_young_kidney)],
]

# Determine shared samples between count matrix and metadata
keep_samples_old_young <- intersect(
  colnames(norm_cts_sub_old_young_from_nfd_cd),
  old_young_filtered$Sample_ID
)

# Subset count matrix to those samples
norm_cts_sub_old_young_from_nfd_cd <- norm_cts_sub_old_young_from_nfd_cd[, keep_samples_old_young]

# Filter and arrange metadata in same order
old_young_filtered <- old_young_filtered %>%
  dplyr::filter(Sample_ID %in% keep_samples_old_young) %>%
  dplyr::arrange(match(Sample_ID, keep_samples_old_young))

# Scale expression per gene
heat_old_young_from_nfd_cd <- t(scale(t(norm_cts_sub_old_young_from_nfd_cd)))

# Build annotation for Old vs Young
condition_labels_old_young <- factor(
  old_young_filtered$Treatment,
  levels = c("Young", "Old")
)

# Define heatmap column annotations
column_ha_old_young <- HeatmapAnnotation(
  Condition = condition_labels_old_young,
  annotation_name_side = "right",
  show_annotation_name = FALSE,
  col = list(
    Condition = c(
      "Young" = "#ff6200",  # orange
      "Old"   = "#808080"   # gray
    )
  )
)

# Cluster columns (samples) but keep same gene order as NFD vs CD
dend_old_young <- hclust(dist(t(heat_old_young_from_nfd_cd))) %>% as.dendrogram()

# Ensure gene order only includes genes present in this matrix
gene_order_nfd_cd <- gene_order_nfd_cd[gene_order_nfd_cd %in% rownames(heat_old_young_from_nfd_cd)]

# Build Old vs Young heatmap using same row order
hmap_old_young_from_nfd_cd <- Heatmap(
  heat_old_young_from_nfd_cd,
  name = "Gene\nZ-Score",
  top_annotation = column_ha_old_young,
  cluster_rows = FALSE,
  row_order = gene_order_nfd_cd,
  cluster_columns = dend_old_young,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_dend_reorder = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 12)
)

# Draw the heatmap
draw(
  hmap_old_young_from_nfd_cd,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)


# ------------------------------------------------------------------------------
# Prepare contrast-specific DEG data for contrast comparisons
# ------------------------------------------------------------------------------

# Use prep_deg_df function
NFD_CD_kidney <- prep_deg_df(NFD_CD_adjp_kidney, "NFDvsCD")
Old_Young_kidney <- prep_deg_df(Old_Young_adjp_kidney, "OldvsYoung")


# ------------------------------------------------------------------------------
##  Contrast Comparisons: NFD vs CD and Old vs Young
# ------------------------------------------------------------------------------
# Using: compare_contrasts_multi()
# ------------------------------------------------------------------------------
##  Shared DEGs (any direction; no filtering by sign)
# ------------------------------------------------------------------------------
shared_all_NFDvsCD_OldvsYoung_kidney <- compare_contrasts_multi(
  dfs = list(NFD_CD_kidney, Old_Young_kidney),
  contrast_names = c("NFDvsCD", "OldvsYoung"),
  direction = "all"
)

length(shared_all_NFDvsCD_OldvsYoung_kidney$Symbol)
# 2148 shared DEGs




# ------------------------------------------------------------------------------
##  Shared Upregulated DEGs (consistent direction)
# ------------------------------------------------------------------------------
shared_up_NFDvsCD_OldvsYoung_kidney <- compare_contrasts_multi(
  dfs = list(NFD_CD_kidney, Old_Young_kidney),
  contrast_names = c("NFDvsCD", "OldvsYoung"),
  direction = "up"
)

length(shared_up_NFDvsCD_OldvsYoung_kidney$Symbol)
# 857 shared upregulated DEGs




# ------------------------------------------------------------------------------
##  Shared Downregulated DEGs (consistent direction)
# ------------------------------------------------------------------------------
shared_down_NFDvsCD_OldvsYoung_kidney <- compare_contrasts_multi(
  dfs = list(NFD_CD_kidney, Old_Young_kidney),
  contrast_names = c("NFDvsCD", "OldvsYoung"),
  direction = "down"
)

length(shared_down_NFDvsCD_OldvsYoung_kidney$Symbol)
# 1226 shared downregulated DEGs




# ------------------------------------------------------------------------------
##  Shared DEGs with Opposite Directions (discordant log2FC)
# ------------------------------------------------------------------------------
shared_opposite_NFDvsCD_OldvsYoung_kidney <- compare_contrasts_multi(
  dfs = list(NFD_CD_kidney, Old_Young_kidney),
  contrast_names = c("NFDvsCD", "OldvsYoung"),
  direction = "opposite"
)

length(shared_opposite_NFDvsCD_OldvsYoung_kidney$Symbol)
# 65 genes with opposite directions



# ------------------------------------------------------------------------------
## Tidy summary of shared DEGs across contrasts
# ------------------------------------------------------------------------------

shared_deg_summary_old_young_kidney <- tibble::tibble(
  Comparison = "NFDvsCD + OldvsYoung",
  Shared_Total = length(shared_all_NFDvsCD_OldvsYoung_kidney$Symbol),
  Shared_Up = length(shared_up_NFDvsCD_OldvsYoung_kidney$Symbol),
  Shared_Down = length(shared_down_NFDvsCD_OldvsYoung_kidney$Symbol),
  Shared_Opposite = length(shared_opposite_NFDvsCD_OldvsYoung_kidney$Symbol)
)

shared_deg_summary_old_young_kidney
# Comparison        Shared_Total Shared_Up Shared_Down Shared_Opposite
# NFDvsCD + OldvsYoung         2148       857        1226              65



# ------------------------------------------------------------------------------
##  Identify and Summarize DEGs Unique to Each Contrast
# ------------------------------------------------------------------------------

dfs_unique_old_young_kidney <- list(
  NFDvsCD = NFD_CD_kidney,
  OldvsYoung = Old_Young_kidney
)

# Use the get_unique_DEGs_multi and summarize_unique_DEGs_multi functions
unique_deg_rows_old_young_kidney <- get_unique_DEGs_multi(dfs_unique_old_young_kidney)
unique_deg_summary_old_young_kidney <- summarize_unique_DEGs_multi(unique_deg_rows_old_young_kidney)

print(unique_deg_summary_old_young_kidney)
# Contrast Total Upregulated Downregulated
# 1 OldvsYoung  6112        3345          2767
# 2 NFDvsCD      809         527           282



# ------------------------------------------------------------------------------
## Heatmap of Log2 Fold Changes for Shared DEGs 
# ------------------------------------------------------------------------------

# Prepare pivoted matrix of log2FC values
# ------------------------------------------------------------------------------

# Use the output from compare_contrasts()
# Add contrast labels explicitly for clarity
logFC_heat_df_kidney <- shared_all_NFDvsCD_OldvsYoung_kidney %>%
  dplyr::select(Symbol, 
                OldvsYoung_log2FoldChange, 
                NFDvsCD_log2FoldChange)

dim(logFC_heat_df_kidney) # 2148


# Convert to matrix with gene symbols as row names
# Ensure rownames = Symbol column
logFC_heat_matrix <- logFC_heat_df_kidney %>%
  tibble::column_to_rownames(var = "Symbol") %>%
  as.matrix()


# Plot heatmap of log2FC values
# Set custom column labels for display
col_hm_labels <- c("Old vs Young", "NFD vs CD")

# Define heatmap 
hmap_lfc <- Heatmap(
  logFC_heat_matrix,
  name = "Gene\nLog2FC",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  row_dend_reorder = TRUE,
  show_column_names = TRUE,
  column_labels = col_hm_labels,
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    title_position = "topcenter"
  )
)

# Draw heatmap
draw(
  hmap_lfc,
  heatmap_legend_side = "right",
  padding = unit(c(2, 2, 2, 12), "mm") 
)



# ------------------------------------------------------------------------------
## Hypergeometric Test of DEG Overlap Using GeneOverlap
# ------------------------------------------------------------------------------

# Define shared gene universe (background set)
# NFD_CD_shrunk_kidney is reused from recovery_kidney/01_recovery_kidney_analysis.R
dim(NFD_CD_shrunk_kidney) # 16305

# Extract gene symbols from results tables (gene symbols are stored as rownames)
universe_nfd_cd_kidney  <- rownames(NFD_CD_shrunk_kidney)
universe_old_young_kidney <- rownames(Old_Young_shrunk_kidney)

# Check universe lengths
length(universe_nfd_cd_kidney) # 16305
length(universe_old_young_kidney) # 16558

# Intersect to get common background
background_genes_shared_kidney <- intersect(
  universe_nfd_cd_kidney, 
  universe_old_young_kidney
  )

universe_size_shared_kidney <- length(background_genes_shared_kidney)


# Report size
cat("Total background genes (common to both):", length(background_genes_shared_kidney), "\n") 
# Total background genes (common to both): 15944


# Define gene lists
# ------------------------------------------------------------------------------
# Ensure Symbol column exists
NFD_genes_kidney <- NFD_CD_kidney$Symbol
Old_genes_kidney <- Old_Young_kidney$Symbol

# Check DEG lengths
cat("NFDvsCD DEGs:", length(NFD_genes_kidney), "\n") # 2957
cat("OldvsYoung DEGs:", length(Old_genes_kidney), "\n") # 8260


# Perform hypergeometric test with GeneOverlap
# ------------------------------------------------------------------------------
go.obj_shared_kidney <- newGeneOverlap(
  listA = NFD_genes_kidney,
  listB = Old_genes_kidney,
  genome.size = universe_size_shared_kidney
)

go.obj_shared_kidney <- testGeneOverlap(go.obj_shared_kidney)

# View result summary
print(go.obj_shared_kidney)

# Detailed information about this GeneOverlap object:
# listA size=2957, e.g. Cyp4a14 Azgp1 Ttc38
# listB size=8260, e.g. Sox17 Mrpl15 Lypla1
# Intersection size=2148, e.g. Cyp4a14 Azgp1 Ugt1a2
# Union size=9069, e.g. Cyp4a14 Azgp1 Ttc38
# Genome size=15944
# # Contingency Table:
# notA  inA
# notB 6875  809
# inB  6112 2148
# Overlapping p-value=4.4e-144
# Odds ratio=3.0
# Overlap tested using Fisher's exact test (alternative=greater)
# Jaccard Index=0.2



