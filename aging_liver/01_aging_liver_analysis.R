## -----------------------------------------------------------------------------
## Script Name: 01_aging_liver_analysis.R
## Project: Bulk RNA-seq Analysis of ANDY CD/NFD and young/old WT mice
## Study: NFD vs CD and Old vs Young DEG-level comparison
## Analysis Type: Set-up, differential expression analysis, heatmap
##                visualization
## Tissues: Liver
## Dependencies: Requires 00_aging_liver_packages_functions.R script
## Last Updated: 2026-08-25
## -----------------------------------------------------------------------------

# ==============================================================================
## NFD vs CD and Old vs Young (DEG level analysis)
# ==============================================================================
# ------------------------------------------------------------------------------
##  Data Preprocessing
# ------------------------------------------------------------------------------

# Read in raw counts matrix from featureCounts output
counts_file <- "path/to/counts_matrix.csv"  # GEO-deposited counts matrix
counts_data_orig_liver <- read.csv(counts_file)
dim(counts_data_orig_liver)
# Contains counts for 57,132 genes across 43 samples


# ------------------------------------------------------------------------------
##  Map Ensembl gene IDs to gene symbols using biomaRt 
# ------------------------------------------------------------------------------
# Extract Ensembl gene IDs for annotation
ensembl.ids <- counts_data_orig_liver$Geneid
head(counts_data_orig_liver)


# Load mouse dataset from Ensembl BioMart (version 112)
# Use listEnsembl() to fine 'mmusculus_gene_ensembl' dataset
ensembl <- useEnsembl(biomart = "genes", 
                      version = "112", 
                      dataset = 'mmusculus_gene_ensembl')

# Retrieve mapping from Ensembl gene IDs (with version) to external gene names (symbols)
symbol <- getBM(attributes = c('ensembl_gene_id_version', 'external_gene_name'), 
                filters = "ensembl_gene_id_version", 
                values = ensembl.ids,
                mart = ensembl)

# Merge gene symbols into counts matrix
counts_data_orig_liver <- counts_data_orig_liver %>% 
  dplyr::rename("ensembl_gene_id_version" = "Geneid") %>% 
  dplyr::left_join(symbol, by = "ensembl_gene_id_version") %>% 
  dplyr::filter(!(is.na(external_gene_name) | external_gene_name == "")) # remove rows without gene symbols

dim(counts_data_orig_liver)
# Result: 56,648 genes retained after filtering out missing symbols


# ------------------------------------------------------------------------------
##  Handle duplicated gene symbols 
# ------------------------------------------------------------------------------

# Check how many gene symbols are duplicated
duplicated_symbols <- counts_data_orig_liver$external_gene_name[duplicated(counts_data_orig_liver$external_gene_name)] %>% 
  sort()

length(duplicated_symbols) # 235 symbols appear more than once


# Optional: view or export duplicated entries for manual inspection
dup_symbols_df <- counts_data_orig_liver %>% 
  dplyr::filter(external_gene_name %in% duplicated_symbols)


# Remove duplicates by keeping the first occurrence per gene symbol
# Rationale: downstream tools require unique rownames (gene identifiers)
counts_data_orig_liver <- counts_data_orig_liver %>% 
  dplyr::distinct(external_gene_name, .keep_all = TRUE)

dim(counts_data_orig_liver)
# 56413 x 44



# ------------------------------------------------------------------------------
##  Prepare final count matrix 
# ------------------------------------------------------------------------------

# Set gene symbols as row names and drop unneeded columns
rownames(counts_data_orig_liver) <- counts_data_orig_liver$external_gene_name
counts_data_orig_liver <- counts_data_orig_liver[, -c(1, 44)]  # remove Ensembl ID and symbol columns

# Confirm matrix dimensions: 56,413 genes x 42 samples
dim(counts_data_orig_liver)


# ------------------------------------------------------------------------------
##  Load and align sample metadata 
# ------------------------------------------------------------------------------

# Read in sample information (e.g., group, batch, treatment)
sampleinfo_file <- "path/to/sample_info.csv"  # sample metadata, deposited with GEO submission
colData_orig_liver <- read.csv(sampleinfo_file, header = T)
head(colData_orig_liver)

# Create matrix from data frame
rownames(colData_orig_liver) <- colData_orig_liver[, 1]
colData_orig_liver <- colData_orig_liver[, -1]  # remove redundant column

head(colData_orig_liver)


# Check sample name consistency between metadata and counts matrix
all(colnames(counts_data_orig_liver) %in% rownames(colData_orig_liver)) # should be TRUE
all(colnames(counts_data_orig_liver) == rownames(colData_orig_liver)) # should be TRUE and in order


# ------------------------------------------------------------------------------
##  Final conversion to matrix 
# ------------------------------------------------------------------------------

# Convert to numeric matrix format for downstream analysis (e.g., DESeq2)
counts_data_matrix_orig_liver <- data.matrix(counts_data_orig_liver)


# ------------------------------------------------------------------------------
## Correct for Batch Effects with ComBat-seq and Create a SummarizedExperiment Object
# ------------------------------------------------------------------------------

# Subset: 'Young' and 'CD2W' samples for targeted batch correction
# ------------------------------------------------------------------------------
# These samples were sequenced in different batches, so we batch-correct only this subset
# to avoid confounding true biological differences in other contrasts (e.g., NFD vs CD1D)

batch_corr_sub <- counts_data_orig_liver  %>% 
  dplyr::select(M15L, M20L, M25L, M40L, M45L, M5L, # CD2W
                YL1, YL2, YL3, YL4, YL5, YL6) # Young

# Convert to matrix dtype
batch_corr_sub_matrix <- data.matrix(batch_corr_sub)


# Extract metadata (batch assignments) for batch correction
batch_corr_conditions <- c("CD2W", "Young")
batch_corr_sub_colData <- colData_orig_liver%>%
  dplyr::filter(Treatment %in% batch_corr_conditions)



# Remaining samples (main contrasts: CD1D, Old, NFD etc.) remain uncorrected
# These were all sequenced together and do not require batch correction
counts_data_main_grps <- counts_data_orig_liver %>% 
  dplyr::select(M28L, M33L, M38L, M3L, M43L, M8L, 
                M7512L, M7513L, M7514L, M7515L, M7516L, M7517L, 
                M16L, M1L, M21L, M41L, M46L, M6L,
                M101L, M104L, M106L, M107L, M108L, M109L,
                M17L, M27L, M2L, M37L, M42L, M7L)

# Convert to matrix dtype
counts_data_main_grps_matrix <- data.matrix(counts_data_main_grps)


# Note: This dataset includes additional treatment conditions (ND1DS, ND2D)
# that are not reported in the manuscript. They are retained here (rather than
# dropped) because DESeq2's dispersion estimation uses all samples in the
# design; removing them would change the DESeq2 results for the contrasts
# that are used (NFD vs CD, Old vs Young).
main_sample_conditions <- c("CD1D", "Old", "ND1DL", "ND1DS", "ND2D")
main_grps_colData <- colData_orig_liver%>%
  dplyr::filter(Treatment %in% main_sample_conditions)


# ------------------------------------------------------------------------------
## Apply ComBat-seq for batch correction
# ------------------------------------------------------------------------------

# ComBat_seq models batch effects in raw count data using a negative binomial framework
# It preserves biological variation while adjusting for unwanted batch effects
adjusted_counts <- ComBat_seq(batch_corr_sub_matrix, 
                              batch = batch_corr_sub_colData$Batch)


# Sanity checks: make sure sample metadata aligns with corrected matrix
stopifnot(all(colnames(adjusted_counts) == rownames(batch_corr_sub_colData)))


# ------------------------------------------------------------------------------
## Merge adjusted batch-corrected subset with uncorrected main groups
# ------------------------------------------------------------------------------
# This creates the full expression matrix for downstream DESeq2 analysis
final_counts_matrix_orig_liver <- cbind(adjusted_counts, counts_data_main_grps_matrix)
dim(final_counts_matrix_orig_liver)
# 56413 x 42

final_colData_orig_liver <- rbind(batch_corr_sub_colData, main_grps_colData)


# Confirm column and row alignment
stopifnot(all(colnames(final_counts_matrix_orig_liver) == rownames(final_colData_orig_liver)))


# -----------------------------------------------------------------------------
## Dynamically assign 'Condition' from Treatment column
# ------------------------------------------------------------------------------
# Directly copy clean Treatment
final_colData_orig_liver$Condition <- as.factor(final_colData_orig_liver$Treatment)

# Inspect group counts
table(final_colData_orig_liver$Condition)
# CD1D  CD2W ND1DL ND1DS  ND2D   Old Young 
# 6     6     6     6     6     6     6 


# ------------------------------------------------------------------------------
## Create SummarizedExperiment object
# ------------------------------------------------------------------------------
# Using batch-corrected samples
se_orig_liver <- SummarizedExperiment(
  assays = list(counts = final_counts_matrix_orig_liver), 
  colData = final_colData_orig_liver)



# ------------------------------------------------------------------------------
## Prefiltering for DESeq2
# ------------------------------------------------------------------------------

# Create a DESeqDataSet from the SummarizedExperiment
# ------------------------------------------------------------------------------
# The design formula determines how DESeq2 models gene expression.
dds_orig_liver <- DESeqDataSet(se_orig_liver, design = ~ Treatment)

# Confirm initial gene count
nrow(dds_orig_liver)  # Should be 56,413 genes from the raw count matrix



# Prefilter lowly expressed genes
# ------------------------------------------------------------------------------
# DESeq2 recommends filtering out genes with very low counts across all samples
# These genes contribute little statistical power and increase the multiple testing burden

# Dynamically determine the minimum number of samples in the smallest group
smallestGroupSize_orig_liver <- min(table(colData(dds_orig_liver)$Treatment)) 
print(smallestGroupSize_orig_liver) # 6

# Keep genes with count ≥10 in at least `smallestGroupSize` samples
keep_orig_liver <- rowSums(counts(dds_orig_liver) >= 10) >= smallestGroupSize_orig_liver
dds_orig_liver <- dds_orig_liver[keep_orig_liver, ]
nrow(dds_orig_liver) 
# 15300 genes



# Estimate size factors for normalization
# ------------------------------------------------------------------------------
# Size factors account for differences in sequencing depth across samples
dds_orig_liver <- estimateSizeFactors(dds_orig_liver)



# ------------------------------------------------------------------------------
##  Differential Expression Analysis with DESeq2
# ------------------------------------------------------------------------------

# Set reference level ("CD1D") and run DESeq2 pipeline
# -----------------------------------------------------------------------------
# Use batch-corrected, prefiltered DESeqDataSet object `dds`
dds_orig_liver$Treatment <- relevel(dds_orig_liver$Treatment, ref = "CD1D")

dds_orig_liver <- DESeq(dds_orig_liver)

# View size factors and dispersion plot
sizeFactors(dds_orig_liver)

plotDispEsts(dds_orig_liver)  # Check that dispersions follow expected curve

# View available contrast names (i.e., model coefficients)
resultsNames(dds_orig_liver)


# Differential Expression – ND1DL vs CD1D (NFD vs CD)
# ------------------------------------------------------------------------------

# Contrast: ND1DL (niacin-deficient) vs CD1D (control-fed)
ND1DL_CD1D_liver <- results(dds_orig_liver, contrast = c('Treatment', 'ND1DL', 'CD1D'))

# Apply shrinkage to log2 fold changes for more stable estimates
ND1DL_CD1D_shrunk_liver <- lfcShrink(
  dds_orig_liver,
  coef = 'Treatment_ND1DL_vs_CD1D',
  res = ND1DL_CD1D_liver,
  type = "apeglm"
  ) %>% 
  as.data.frame()

# Subset significantly differentially expressed genes (FDR ≤ 0.05)
ND1DL_CD1D_adjp_liver <- subset(ND1DL_CD1D_shrunk_liver, padj <= 0.05)

# Summary table using summarize_DEGs function
summary_ND1DL_CD1D_liver <- summarize_DEGs(ND1DL_CD1D_shrunk_liver, contrast_name = "NFD vs CD")
summary_ND1DL_CD1D_liver
# Contrast  Total_DEGs Upregulated Downregulated
# NFD vs CD       1165         658           507


# ------------------------------------------------------------------------------
# Relevel for Old vs Young analysis
# ------------------------------------------------------------------------------
# Change reference level so Young is baseline
dds_orig_liver$Treatment <- relevel(dds_orig_liver$Treatment, ref = "Young")

# Rerun DESeq with updated factor levels
dds_orig_liver <- DESeq(dds_orig_liver)
resultsNames(dds_orig_liver)


# Differential Expression – Old vs Young
# ------------------------------------------------------------------------------
Old_Young_liver <- results(dds_orig_liver, contrast = c('Treatment', 'Old', 'Young'))

Old_Young_shrunk_liver <- lfcShrink(
  dds = dds_orig_liver,
  coef = "Treatment_Old_vs_Young",
  res = Old_Young_liver, 
  type = "apeglm"
  ) %>% 
  as.data.frame()

# Subset DEGs (padj ≤ 0.05)
Old_Young_adjp_liver <- subset(Old_Young_shrunk_liver, padj <= 0.05)

# Summary table using summarize_DEGs function
summary_Old_Young_liver <- summarize_DEGs(Old_Young_shrunk_liver, contrast_name = "Old vs Young")
summary_Old_Young_liver
# Contrast  Total_DEGs Upregulated Downregulated
# Old vs Young       7114        4061          3053




# ------------------------------------------------------------------------------
##  ComplexHeatmap Visualization of Top DEGs (Old vs Young–ordered)
# ------------------------------------------------------------------------------

# Extract top DEGs (ND1DL vs CD1D, sorted by adjusted p-value)
# ------------------------------------------------------------------------------
# Sort by padj and select top 1000 most significant DEGs (ND1DLvsCD1D)
ND1DL_CD1D_adjp_liver <- ND1DL_CD1D_adjp_liver[order(ND1DL_CD1D_adjp_liver$padj), ]
ND1DL_CD1D_adjp_sub_liver <- ND1DL_CD1D_adjp_liver[1:1000, ]

dim(ND1DL_CD1D_adjp_sub_liver)  # Should be 1000 × 5

# Sort by padj and select top 1000 most significant DEGs (OldvsYoung)
Old_Young_adjp_liver <- Old_Young_adjp_liver[order(Old_Young_adjp_liver$padj), ]
Old_Young_adjp_sub_liver <- Old_Young_adjp_liver[1:1000, ]

dim(Old_Young_adjp_sub_liver)  # Should be 1000 × 5


# Subset normalized expression matrix for those genes
# ------------------------------------------------------------------------------

# Retrieve the normalized counts
normalized_counts_orig_liver <- counts(dds_orig_liver, normalized = TRUE)
head(normalized_counts_orig_liver)

# Note: the reciprocal heatmap (ordered by the NFD vs CD DEG set) is built
# further below, using the same normalized counts and sample subset.
top_genes_orig_liver <- rownames(Old_Young_adjp_sub_liver)
norm_cts_sub_orig_liver <- normalized_counts_orig_liver[rownames(
  normalized_counts_orig_liver) %in% top_genes_orig_liver, ]

dim(norm_cts_sub_orig_liver) # 1000 x 42

# Scale expression matrix (row-wise z-scores for heatmap)
heat_orig_liver <- t(scale(t(norm_cts_sub_orig_liver)))


# Subset to conditions of interest (ND1DL, CD1D, Young, Old)
# ------------------------------------------------------------------------------

# Define target conditions
treatment_sub_orig_liver <- c("Old", "ND1DL", "CD1D", "Young")

# Subset metadata for those treatments
norm_col_orig_liver <- final_colData_orig_liver[final_colData_orig_liver$Treatment %in% treatment_sub_orig_liver, ]
sample_names_sub_orig_liver <- rownames(norm_col_orig_liver)

# Subset expression matrix for selected samples
mat_sub_orig_liver <- norm_cts_sub_orig_liver[, sample_names_sub_orig_liver]
heat_sub_orig_liver <- t(scale(t(mat_sub_orig_liver)))  # z-score rows again on subset


# Define annotations and colors
# ------------------------------------------------------------------------------

# Rename treatment labels for clarity
condition_labels_orig_liver <- norm_col_orig_liver$Treatment
condition_labels_orig_liver <- factor(
  condition_labels_orig_liver, 
  levels = c("Young", "CD1D", "Old", "ND1DL"),
  labels = c("Young", "CD", "Old", "NFD"))

# Column annotations
column_ha_sub_orig_liver <- HeatmapAnnotation(
  Condition = condition_labels_orig_liver,
  annotation_name_side = "right",
  show_annotation_name = FALSE,
  col = list(Condition = c(
    "NFD" = '#a00000',
    "Old" = '#808080',
    "CD" = '#000080',
    "Young" = '#ff6200'
  ))
)


# Rotate column dendrogram for visual grouping
# ------------------------------------------------------------------------------
# Build and rotate dendrogram
# (the NFD vs CD ordering used for the reciprocal heatmap below is applied
# separately to dend_nfdcd_liver, since it's built from a different gene set)
dend_orig_liver <- hclust(dist(t(heat_sub_orig_liver))) %>% as.dendrogram()

dend_orig_liver <- dendextend::rotate(
  dend_orig_liver, 
  order = c(3:8, 16:21, 22:24, 14, 15, 10:13, 9, 1, 2))  # Old vs Young (1000)

# Plot the heatmap
# ------------------------------------------------------------------------------
hmap_sub_orig_liver <- Heatmap(
  heat_sub_orig_liver,
  name = 'Gene\nZ-Score',
  top_annotation = column_ha_sub_orig_liver,
  cluster_rows = TRUE,
  cluster_columns = dend_orig_liver,
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
  hmap_sub_orig_liver,
  heatmap_legend_side = 'right',
  annotation_legend_side = 'right',
  merge_legend = TRUE
)



# ------------------------------------------------------------------------------
##  ComplexHeatmap Visualization of Top DEGs (NFD vs CD–ordered, reciprocal heatmap)
# ------------------------------------------------------------------------------

# Subset normalized expression matrix for the NFD vs CD top DEGs
# ------------------------------------------------------------------------------
top_genes_nfdcd_liver <- rownames(ND1DL_CD1D_adjp_sub_liver)
norm_cts_sub_nfdcd_liver <- normalized_counts_orig_liver[rownames(
  normalized_counts_orig_liver) %in% top_genes_nfdcd_liver, ]

dim(norm_cts_sub_nfdcd_liver) # 1000 x 42

# Scale expression matrix (row-wise z-scores for heatmap)
heat_nfdcd_liver <- t(scale(t(norm_cts_sub_nfdcd_liver)))


# Subset to the same sample/condition subset used above (ND1DL, CD1D, Young, Old)
# ------------------------------------------------------------------------------
mat_sub_nfdcd_liver <- norm_cts_sub_nfdcd_liver[, sample_names_sub_orig_liver]
heat_sub_nfdcd_liver <- t(scale(t(mat_sub_nfdcd_liver)))  # z-score rows again on subset


# Rotate column dendrogram for visual grouping
# ------------------------------------------------------------------------------
# Build and rotate dendrogram
dend_nfdcd_liver <- hclust(dist(t(heat_sub_nfdcd_liver))) %>% as.dendrogram()

dend_nfdcd_liver <- dendextend::rotate(
  dend_nfdcd_liver,
  order = c(4:13, 2:3, 23:24, 14:17, 18:22, 1))  # NFD vs CD (1000)

# Plot the heatmap
# ------------------------------------------------------------------------------
# Reuses column_ha_sub_orig_liver, since the displayed samples/conditions are
# identical to the Old vs Young–ordered heatmap above — only the gene set and
# column clustering order differ.
hmap_sub_nfdcd_liver <- Heatmap(
  heat_sub_nfdcd_liver,
  name = 'Gene\nZ-Score',
  top_annotation = column_ha_sub_orig_liver,
  cluster_rows = TRUE,
  cluster_columns = dend_nfdcd_liver,
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
  hmap_sub_nfdcd_liver,
  heatmap_legend_side = 'right',
  annotation_legend_side = 'right',
  merge_legend = TRUE
)




# ------------------------------------------------------------------------------
# Prepare contrast-specific DEG data for contrast comparisons
# ------------------------------------------------------------------------------

# Use prep_deg_df function
ND1DL_CD1D_liver <- prep_deg_df(ND1DL_CD1D_adjp_liver, "NFDvsCD")
Old_Young_liver <- prep_deg_df(Old_Young_adjp_liver, "OldvsYoung")


# ------------------------------------------------------------------------------
##  Contrast Comparisons: NFD vs CD and Old vs Young
# ------------------------------------------------------------------------------
# Using: compare_contrasts_multi()
# Inputs: ND1DL_CD1D_liver (NFD vs CD DEGs), Old_Young_liver (Old vs Young DEGs)
# ------------------------------------------------------------------------------
##  Shared DEGs (any direction; no filtering by sign)
# ------------------------------------------------------------------------------
shared_all_ND1DLvsCD1DandOldvsYoung <- compare_contrasts_multi(
  dfs = list(ND1DL_CD1D_liver, Old_Young_liver),
  contrast_names = c("NFDvsCD", "OldvsYoung"),
  direction = "all"
)

length(shared_all_ND1DLvsCD1DandOldvsYoung$Symbol)
# 712 shared DEGs


# ------------------------------------------------------------------------------
##  Shared Upregulated DEGs (consistent direction)
# ------------------------------------------------------------------------------
shared_up_ND1DLvsCD1DandOldvsYoung <- compare_contrasts_multi(
  dfs = list(ND1DL_CD1D_liver, Old_Young_liver),
  contrast_names = c("NFDvsCD", "OldvsYoung"),
  direction = "up"
)

length(shared_up_ND1DLvsCD1DandOldvsYoung$Symbol)
# 391 shared upregulated DEGs


# ------------------------------------------------------------------------------
##  Shared Downregulated DEGs (consistent direction)
# ------------------------------------------------------------------------------
shared_down_ND1DLvsCD1DandOldvsYoung <- compare_contrasts_multi(
  dfs = list(ND1DL_CD1D_liver, Old_Young_liver),
  contrast_names = c("NFDvsCD", "OldvsYoung"),
  direction = "down"
)

length(shared_down_ND1DLvsCD1DandOldvsYoung$Symbol)
# 272 shared downregulated DEGs


# ------------------------------------------------------------------------------
##  Shared DEGs with Opposite Directions (discordant log2FC)
# ------------------------------------------------------------------------------
shared_opposite_ND1DLvsCD1DandOldvsYoung <- compare_contrasts_multi(
  dfs = list(ND1DL_CD1D_liver, Old_Young_liver),
  contrast_names = c("NFDvsCD", "OldvsYoung"),
  direction = "opposite"
)

length(shared_opposite_ND1DLvsCD1DandOldvsYoung$Symbol)
# 49 genes with opposite directions




# ------------------------------------------------------------------------------
## Tidy summary of shared DEGs across contrasts
# ------------------------------------------------------------------------------

shared_deg_summary_orig_liver <- tibble::tibble(
  Comparison = "NFDvsCD + OldvsYoung",
  Shared_Total = length(shared_all_ND1DLvsCD1DandOldvsYoung$Symbol),
  Shared_Up = length(shared_up_ND1DLvsCD1DandOldvsYoung$Symbol),
  Shared_Down = length(shared_down_ND1DLvsCD1DandOldvsYoung$Symbol),
  Shared_Opposite = length(shared_opposite_ND1DLvsCD1DandOldvsYoung$Symbol)
)

shared_deg_summary_orig_liver
# Comparison        Shared_Total Shared_Up Shared_Down Shared_Opposite
# NFDvsCD + OldvsYoung          712       391         272              49



# ------------------------------------------------------------------------------
##  Identify and Summarize DEGs Unique to Each Contrast
# ------------------------------------------------------------------------------

dfs_unique_orig_liver <- list(
  NFDvsCD = ND1DL_CD1D_liver,
  OldvsYoung = Old_Young_liver
)

# Use the get_unique_DEGs_multi and summarize_unique_DEGs_multi functions
unique_deg_rows_orig_liver <- get_unique_DEGs_multi(dfs_unique_orig_liver)
unique_deg_summary_orig_liver <- summarize_unique_DEGs_multi(unique_deg_rows_orig_liver)

print(unique_deg_summary_orig_liver)
# Contrast Total Upregulated Downregulated
# OldvsYoung  6402        3645          2757
# NFDvsCD      453         243           210



# ------------------------------------------------------------------------------
##  Hypergeometric Test of DEG Overlap Using GeneOverlap
##  Contrast tested here:
##     NFD vs CD  ∩  Old vs Young
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
##  Define Universe (Background Gene Set)
# ------------------------------------------------------------------------------
# Use rownames of filtered DESeq2 object as expressed gene universe
background_genes_orig_liver <- rownames(dds_orig_liver)
universe_size_orig_liver <- length(background_genes_orig_liver)

cat("Total background genes:", universe_size_orig_liver, "\n")
# Expected: 15,300 genes after prefiltering


# ------------------------------------------------------------------------------
##  Extract DEG Lists for Testing
# ------------------------------------------------------------------------------
# Ensure Symbol column exists for each DEG table
ND1DL_CD1D_genes <- ND1DL_CD1D_liver$Symbol
Old_Young_genes <- Old_Young_liver$Symbol

cat("NFD vs CD DEGs:", length(ND1DL_CD1D_genes), "\n") # 1165
cat("Old vs Young DEGs:", length(Old_Young_genes), "\n") # 7114


# ------------------------------------------------------------------------------
##  Perform Hypergeometric Test (GeneOverlap)
# ------------------------------------------------------------------------------

go.obj_orig_liver <- newGeneOverlap(
  listA = ND1DL_CD1D_genes,
  listB = Old_Young_genes,
  genome.size = universe_size_orig_liver
)

go.obj_orig_liver <- testGeneOverlap(go.obj_orig_liver)

# Display statistical results
print(go.obj_orig_liver)

# Detailed information about this GeneOverlap object:
# listA size=1165, e.g. Mup9 Mup-ps4 Ltf
# listB size=7114, e.g. Gm6368 Gm43655 Gm4948
# Intersection size=712, e.g. Mup9 Mup-ps4 Ltf
# Union size=7567, e.g. Mup9 Mup-ps4 Ltf
# Genome size=15300
# # Contingency Table:
# notA inA
# notB 7733 453
# inB  6402 712
# Overlapping p-value=1.5e-25
# Odds ratio=1.9
# Overlap tested using Fisher's exact test (alternative=greater)
# Jaccard Index=0.1




# ------------------------------------------------------------------------------
## Heatmap of Log2 Fold Changes for Shared DEGs 
# ------------------------------------------------------------------------------

# Prepare pivoted matrix of log2FC values
# ------------------------------------------------------------------------------
# Use the output from compare_contrasts()
# Columns: Symbol, logFC_1 (ND1DLvsCD1D), logFC_2 (OldvsYoung)

# Add contrast labels explicitly for clarity
# Order to select: Symbol, OldvsYoung, NFDvsCD
logFC_heat_df_orig_liver <- shared_all_ND1DLvsCD1DandOldvsYoung %>%
  dplyr::select(Symbol,
                OldvsYoung_log2FoldChange,
                NFDvsCD_log2FoldChange)


# Convert to matrix with gene symbols as row names
# ------------------------------------------------------------------------------
# Ensure rownames = Symbol column
logFC_heat_matrix_orig_liver <- logFC_heat_df_orig_liver %>%
  tibble::column_to_rownames(var = "Symbol") %>%
  as.matrix()


# Plot heatmap of log2FC values
# ------------------------------------------------------------------------------
# Set custom column labels for display
col_hm_labels_orig_liver <- c("Old vs Young", "NFD vs CD")

hmap_lfc_orig_liver <- Heatmap(
  logFC_heat_matrix_orig_liver,
  name = "Gene\nLog2FC",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  row_dend_reorder = TRUE,
  show_column_names = TRUE,
  column_labels = col_hm_labels_orig_liver,
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    title_position = "topcenter"
  )
)

# Draw heatmap
draw(
  hmap_lfc_orig_liver,
  heatmap_legend_side = "right"
)



sessionInfo()
# > sessionInfo()
# R version 4.5.1 (2025-06-13)
# Platform: x86_64-pc-linux-gnu
# Running under: Linux Mint 21.3
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
