## -----------------------------------------------------------------------------
## Script Name: 01_recovery_kidney_analysis.R
## Project: Bulk RNA-seq Analysis of ANDY CD/NFD and young/old WT mice
## Study: NFD recovery
## Analysis Type: Set-up, differential expression analysis, heatmap
##                visualization, pathway (GO MF) over-representation analysis
## Tissues: Kidney
## Dependencies: Requires 00_recovery_kidney_packages_functions.R script
## Last Updated: 2026-08-25
## -----------------------------------------------------------------------------


# ==============================================================================
## NFD recovery study (main sequencing analysis; NFD, CD, and Rec samples)
# ==============================================================================
# ------------------------------------------------------------------------------
##  Data Preprocessing
# ------------------------------------------------------------------------------

# Read in raw counts matrix from featureCounts output
counts_file <- "path/to/counts_matrix.csv"  # GEO-deposited counts matrix

counts_data_rec_main_kidney <- read.csv(counts_file)

# Contains 57,132 genes across 15 samples
dim(counts_data_rec_main_kidney)



# ------------------------------------------------------------------------------
##  Map Ensembl Gene IDs to Gene Symbols (biomaRt)
# ------------------------------------------------------------------------------

# Extract Ensembl gene IDs (with version)
ensembl.ids <- counts_data_rec_main_kidney$Geneid
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
counts_data_rec_main_kidney <- counts_data_rec_main_kidney %>%
  dplyr::rename(ensembl_gene_id_version = Geneid) %>%
  dplyr::left_join(symbol, by = "ensembl_gene_id_version") %>%
  dplyr::filter(!(is.na(external_gene_name) | external_gene_name == ""))

# 56,648 genes retained after filtering out rows without symbols
dim(counts_data_rec_main_kidney)



# ------------------------------------------------------------------------------
##  Handle Duplicated Gene Symbols
# ------------------------------------------------------------------------------

# Identify duplicated gene symbols
duplicated_symbols <- counts_data_rec_main_kidney$external_gene_name[
  duplicated(counts_data_rec_main_kidney$external_gene_name)
] %>% sort()

length(duplicated_symbols)  # 235 duplicated symbols

# Inspect duplicated symbols
dup_symbols_df <- counts_data_rec_main_kidney %>%
  dplyr::filter(external_gene_name %in% duplicated_symbols)

# Remove duplicates (keep first occurrence)
counts_data_rec_main_kidney <- counts_data_rec_main_kidney %>%
  dplyr::distinct(external_gene_name, .keep_all = TRUE)

# 56,413 genes x 16 columns after deduplication
dim(counts_data_rec_main_kidney)



# ------------------------------------------------------------------------------
##  Prepare Final Count Matrix
# ------------------------------------------------------------------------------

# Set gene symbols as rownames
rownames(counts_data_rec_main_kidney) <- counts_data_rec_main_kidney$external_gene_name

# Remove Ensembl ID and gene symbol columns
counts_data_rec_main_kidney <- counts_data_rec_main_kidney[, -c(1, 16)]

# Result should be 56,413 genes x 14 samples
dim(counts_data_rec_main_kidney)



# ------------------------------------------------------------------------------
##  Load and Align Sample Metadata
# ------------------------------------------------------------------------------

# Read metadata (sample info table)
sampleinfo_file <- "path/to/sample_info.csv"  # sample metadata, deposited with GEO submission

colData_rec_main_kidney <- read.csv(
  sampleinfo_file,
  header = TRUE
)

head(colData_rec_main_kidney)

# Move sample names into rownames
rownames(colData_rec_main_kidney) <- colData_rec_main_kidney[, 1]
colData_rec_main_kidney <- colData_rec_main_kidney[, -1]

# Ensure sample names match between metadata and count matrix
# Required for DESeq2
all(colnames(counts_data_rec_main_kidney) %in% rownames(colData_rec_main_kidney)) # TRUE
all(colnames(counts_data_rec_main_kidney) == rownames(colData_rec_main_kidney)) # TRUE 

# Distribution across conditions
table(colData_rec_main_kidney$Condition)
# CD: 6    NFD: 4    Recovery: 4



# ------------------------------------------------------------------------------
##  Convert Counts to Matrix for DESeq2
# ------------------------------------------------------------------------------

counts_data_matrix_rec_kidney <- data.matrix(counts_data_rec_main_kidney)

head(counts_data_matrix_rec_kidney)


# ------------------------------------------------------------------------------
##  Create a SummarizedExperiment Object
# ------------------------------------------------------------------------------

# Construct SummarizedExperiment for the NFD recovery liver samples
se_main_rec_kidney <- SummarizedExperiment(
  assays = list(counts = counts_data_matrix_rec_kidney),
  colData = colData_rec_main_kidney
)


# ------------------------------------------------------------------------------
##  Prefiltering
# ------------------------------------------------------------------------------

# Create a DESeqDataSet for the NFD recovery study
dds_main_rec_kidney <- DESeqDataSet(se_main_rec_kidney, design = ~ Condition)

# Confirm initial gene count (expected 56,413)
nrow(dds_main_rec_kidney)



# ------------------------------------------------------------------------------
##  Prefilter Lowly Expressed Genes
# ------------------------------------------------------------------------------

# Determine smallest group size (used for filtering threshold)
smallestGroupSize_rec_kidney <- min(table(colData(dds_main_rec_kidney)$Condition))
print(smallestGroupSize_rec_kidney) # 4

# Keep genes with count ≥ 10 in at least `smallestGroupSize` samples
keep_rec_kidney <- rowSums(counts(dds_main_rec_kidney) >= 10) >= smallestGroupSize_rec_kidney
dds_main_rec_kidney <- dds_main_rec_kidney[keep_rec_kidney, ]

# 16,305 genes should remain
nrow(dds_main_rec_kidney)



# ------------------------------------------------------------------------------
##  Estimate Size Factors (Normalization)
# ------------------------------------------------------------------------------

# Normalize for sequencing depth differences
dds_main_rec_kidney <- estimateSizeFactors(dds_main_rec_kidney)



# ------------------------------------------------------------------------------
##  Differential Expression Analysis with DESeq2
# ------------------------------------------------------------------------------

# Set reference level ("CD") and run DESeq2 pipeline
# ------------------------------------------------------------------------------
dds_main_rec_kidney$Condition <- relevel(dds_main_rec_kidney$Condition, ref = "CD")

dds_main_rec_kidney <- DESeq(dds_main_rec_kidney)

# Review size factors and dispersion estimates
sizeFactors(dds_main_rec_kidney)
plotDispEsts(dds_main_rec_kidney)

# Inspect model coefficients / contrast names
resultsNames(dds_main_rec_kidney)
# Condition_NFD_vs_CD
# Condition_Recovery_vs_CD


# ------------------------------------------------------------------------------
##  Contrast: NFD vs CD (recovery study)
# ------------------------------------------------------------------------------

# Extract DESeq2 results
NFD_CD_kidney <- results(dds_main_rec_kidney, contrast = c("Condition", "NFD", "CD"))

# Apply LFC shrinkage (apeglm)
NFD_CD_shrunk_kidney <- lfcShrink(
  dds_main_rec_kidney,
  coef = "Condition_NFD_vs_CD",
  res  = NFD_CD_kidney,
  type = "apeglm"
) %>% 
  as.data.frame()

# Identify significant DEGs (FDR ≤ 0.05)
NFD_CD_adjp_kidney <- subset(NFD_CD_shrunk_kidney, padj <= 0.05)

# Summary table using summarize_DEGs function
summary_NFD_CD_kidney <- summarize_DEGs(NFD_CD_shrunk_kidney, contrast_name = "NFD vs CD")
summary_NFD_CD_kidney
# Total_DEGs = 2957
# Upregulated = 1418
# Downregulated = 1539


# ------------------------------------------------------------------------------
##  Contrast: Recovery vs CD
# ------------------------------------------------------------------------------

# Extract DESeq2 results
Rec_CD_kidney <- results(dds_main_rec_kidney, contrast = c("Condition", "Recovery", "CD"))

# Apply LFC shrinkage with apeglm
Rec_CD_shrunk_kidney <- lfcShrink(
  dds_main_rec_kidney,
  coef = "Condition_Recovery_vs_CD",
  res  = Rec_CD_kidney,
  type = "apeglm"
) %>% 
  as.data.frame()

# Significant DEGs (padj ≤ 0.05)
Rec_CD_adjp_kidney <- subset(Rec_CD_shrunk_kidney, padj <= 0.05)


# Summary table using summarize_DEGs function
summary_Rec_CD_kidney <- summarize_DEGs(Rec_CD_shrunk_kidney, contrast_name = "Recovery vs CD")
summary_Rec_CD_kidney
# Total_DEGs = 124
# Upregulated = 96
# Downregulated = 28


# ------------------------------------------------------------------------------
##  ComplexHeatmap Visualization of Top NFD vs CD DEGs
# ------------------------------------------------------------------------------

# Select top 1000 DEGs ranked by adjusted p-value
# ------------------------------------------------------------------------------
NFD_CD_adjp_kidney <- NFD_CD_adjp_kidney[order(NFD_CD_adjp_kidney$padj), ]
NFD_CD_top1000_kidney <- NFD_CD_adjp_kidney[1:1000, ]
dim(NFD_CD_top1000_kidney)  # 1000 × 5 columns



# ------------------------------------------------------------------------------
##  Subset Normalized Counts for Heatmap
# ------------------------------------------------------------------------------

# Retrieve normalized counts from DESeq2 object
normalized_counts_main_rec_kidney <- counts(dds_main_rec_kidney, normalized = TRUE)

# Subset to top DEG list
top_genes_kidney <- rownames(NFD_CD_top1000_kidney)

norm_cts_sub_kidney <- normalized_counts_main_rec_kidney[
  rownames(normalized_counts_main_rec_kidney) %in% top_genes_kidney,
]

dim(norm_cts_sub_kidney)  # 1000 × 14 samples

# Row-wise z-score transformation (required for heatmap visualization)
heat_mat_kidney <- t(scale(t(norm_cts_sub_kidney)))



# ------------------------------------------------------------------------------
##  Column Annotation (Condition)
# ------------------------------------------------------------------------------

condition_labels_main_rec_kidney <- colData_rec_main_kidney$Condition_Sub
condition_labels_main_rec_kidney <- factor(
  condition_labels_main_rec_kidney,
  levels = c("CD", "NFD", "Rec"),
  labels = c("CD", "NFD", "Rec")
)

column_ha_kidney <- HeatmapAnnotation(
  Condition = condition_labels_main_rec_kidney,
  annotation_name_side = "right",
  show_annotation_name = FALSE,
  col = list(
    Condition = c(
      "CD"  = "#000080",
      "NFD" = "#a00000",
      "Rec" = "#eab77c"
    )
  )
)



# ------------------------------------------------------------------------------
##  Column Clustering (optional rotation for interpretability)
# ------------------------------------------------------------------------------

# Build hierarchical clustering dendrogram
dend_kidney <- hclust(dist(t(heat_mat_kidney))) %>% as.dendrogram()

# ------------------------------------------------------------------------------
##  Plot Heatmap
# ------------------------------------------------------------------------------

hmap_main_rec_kidney <- Heatmap(
  heat_mat_kidney,
  name = "Gene\nZ-Score",
  top_annotation = column_ha_kidney,
  cluster_rows = TRUE,
  cluster_columns = dend_kidney,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  row_dend_reorder = TRUE,
  column_dend_reorder = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 12)
)

# Draw the heatmap
draw(
  hmap_main_rec_kidney,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)



# ------------------------------------------------------------------------------
# Prepare contrast-specific DEG data for contrast comparisons
# ------------------------------------------------------------------------------

# Use prep_deg_df function
NFD_CD_kidney <- prep_deg_df(NFD_CD_adjp_kidney, "NFDvsCD")
Rec_CD_kidney <- prep_deg_df(Rec_CD_adjp_kidney, "RecvsCD")


# ------------------------------------------------------------------------------
##  Contrast Comparisons: NFD vs CD and Recovery vs CD
# ------------------------------------------------------------------------------
# Using: compare_contrasts_multi()
# Inputs: NFD_CD (NFDvsCD DEGs), Rec_CD (RecvsCD DEGs)
# ------------------------------------------------------------------------------
##  Shared DEGs (any direction; no filtering by sign)
# ------------------------------------------------------------------------------
shared_all_NFDvsCD_RecvsCD_kidney <- compare_contrasts_multi(
  dfs = list(NFD_CD_kidney, Rec_CD_kidney),
  contrast_names = c("NFDvsCD", "RecvsCD"),
  direction = "all"
)

length(shared_all_NFDvsCD_RecvsCD_kidney$Symbol)
# 35 shared DEGs




# ------------------------------------------------------------------------------
##  Shared Upregulated DEGs (consistent direction)
# ------------------------------------------------------------------------------
shared_up_NFDvsCD_RecvsCD_kidney <- compare_contrasts_multi(
  dfs = list(NFD_CD_kidney, Rec_CD_kidney),
  contrast_names = c("NFDvsCD", "RecvsCD"),
  direction = "up"
)

length(shared_up_NFDvsCD_RecvsCD_kidney$Symbol)
# 12 shared upregulated DEGs




# ------------------------------------------------------------------------------
##  Shared Downregulated DEGs (consistent direction)
# ------------------------------------------------------------------------------
shared_down_NFDvsCD_RecvsCD_kidney <- compare_contrasts_multi(
  dfs = list(NFD_CD_kidney, Rec_CD_kidney),
  contrast_names = c("NFDvsCD", "RecvsCD"),
  direction = "down"
)

length(shared_down_NFDvsCD_RecvsCD_kidney$Symbol)
# 10 shared downregulated DEGs




# ------------------------------------------------------------------------------
##  Shared DEGs with Opposite Directions (discordant log2FC)
# ------------------------------------------------------------------------------
shared_opposite_NFDvsCD_RecvsCD_kidney <- compare_contrasts_multi(
  dfs = list(NFD_CD_kidney, Rec_CD_kidney),
  contrast_names = c("NFDvsCD", "RecvsCD"),
  direction = "opposite"
)

length(shared_opposite_NFDvsCD_RecvsCD_kidney$Symbol)
# 13 genes with opposite directions



# ------------------------------------------------------------------------------
## Tidy summary of shared DEGs across contrasts
# ------------------------------------------------------------------------------

shared_deg_summary_rec_kidney <- tibble::tibble(
  Comparison = "NFDvsCD + RecvsCD",
  Shared_Total = length(shared_all_NFDvsCD_RecvsCD_kidney$Symbol),
  Shared_Up = length(shared_up_NFDvsCD_RecvsCD_kidney$Symbol),
  Shared_Down = length(shared_down_NFDvsCD_RecvsCD_kidney$Symbol),
  Shared_Opposite = length(shared_opposite_NFDvsCD_RecvsCD_kidney$Symbol)
)

shared_deg_summary_rec_kidney
# Comparison        Shared_Total Shared_Up Shared_Down Shared_Opposite
# NFDvsCD + RecvsCD          35       12         10               13


# ------------------------------------------------------------------------------
##  Identify and Summarize DEGs Unique to Each Contrast
# ------------------------------------------------------------------------------

dfs_unique_main_rec_kidney <- list(
  NFDvsCD = NFD_CD_kidney,
  RecvsCD = Rec_CD_kidney
)

# Use the get_unique_DEGs_multi and summarize_unique_DEGs_multi functions
unique_deg_rows_main_rec_kidney <- get_unique_DEGs_multi(dfs_unique_main_rec_kidney)
unique_deg_summary_main_rec_kidney <- summarize_unique_DEGs_multi(unique_deg_rows_main_rec_kidney)

print(unique_deg_summary_main_rec_kidney)
# Contrast Total Upregulated Downregulated
# NFDvsCD   2922        1406          1516
# RecvsCD     89          71            18



# ------------------------------------------------------------------------------
##  Hypergeometric Test of DEG Overlap Using GeneOverlap
##  Contrast tested here:
##     NFD vs CD  ∩  Recovery vs CD
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
##  Define Universe (Background Gene Set)
# ------------------------------------------------------------------------------
# Use rownames of filtered DESeq2 object as expressed gene universe
background_genes_rec_kidney <- rownames(dds_main_rec_kidney)
universe_size_rec_kidney <- length(background_genes_rec_kidney)

cat("Total background genes:", universe_size_rec_kidney, "\n")
# Expected: 16,305 genes after prefiltering


# ------------------------------------------------------------------------------
##  Extract DEG Lists for Testing
# ------------------------------------------------------------------------------
# Ensure Symbol column exists for each DEG table
NFD_CD_genes_kidney <- NFD_CD_kidney$Symbol
Rec_CD_genes_kidney <- Rec_CD_kidney$Symbol

cat("NFD vs CD DEGs:", length(NFD_CD_genes_kidney), "\n") # 2957
cat("Rec vs CD DEGs:", length(Rec_CD_genes_kidney), "\n") # 124


# ------------------------------------------------------------------------------
##  Perform Hypergeometric Test (GeneOverlap)
# ------------------------------------------------------------------------------

go.obj_rec_kidney <- newGeneOverlap(
  listA = NFD_CD_genes_kidney,
  listB = Rec_CD_genes_kidney,
  genome.size = universe_size_rec_kidney
)

go.obj_rec_kidney <- testGeneOverlap(go.obj_rec_kidney)

# Display statistical results
print(go.obj_rec_kidney)

# Detailed information about this GeneOverlap object:
# listA size=2957, e.g. Cyp4a14 Azgp1 Ttc38
# listB size=124, e.g. Tsc2 Itih3 Apobec3
# Intersection size=35, e.g. Ttc38 Cyp2d9 Fam210a
# Union size=3046, e.g. Cyp4a14 Azgp1 Ttc38
# Genome size=16305
# # Contingency Table:
# notA  inA
# notB 13259 2922
# inB     89   35
# Overlapping p-value=3.7e-03
# Odds ratio=1.8
# Overlap tested using Fisher's exact test (alternative=greater)
# Jaccard Index=0.0



# ------------------------------------------------------------------------------
##  Heatmap of Reversible Genes (NFD ≠ CD, Recovery ≈ CD)
##  Categories:
##     - Reversible: significant in NFD vs CD, not significant in Rec vs CD
##     - Persistent: significant in both contrasts (same direction)
##     - Recovery-specific: not significant in NFD vs CD, significant in Rec vs CD
##     - Other: all remaining genes
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
##  Build a combined DEG table across contrasts
# ------------------------------------------------------------------------------

padj_cutoff_kidney <- 0.05

deg_wide_kidney <- NFD_CD_shrunk_kidney %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Symbol") %>%
  dplyr::select(Symbol, L_NFD = log2FoldChange, padj_NFD = padj) %>%
  dplyr::full_join(
    Rec_CD_shrunk_kidney %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Symbol") %>%
      dplyr::select(Symbol, L_REC = log2FoldChange, padj_REC = padj),
    by = "Symbol"
  ) %>%
  dplyr::mutate(
    # Replace missing padj values
    padj_NFD = ifelse(is.na(padj_NFD), 1, padj_NFD),
    padj_REC = ifelse(is.na(padj_REC), 1, padj_REC),
    
    # Category assignment
    category = dplyr::case_when(
      padj_NFD <= padj_cutoff & padj_REC >  padj_cutoff ~ "Reversible",
      padj_NFD <= padj_cutoff & padj_REC <= padj_cutoff &
        sign(L_NFD) == sign(L_REC)                       ~ "Persistent",
      padj_NFD >  padj_cutoff & padj_REC <= padj_cutoff  ~ "Recovery-specific",
      TRUE                                               ~ "Other"
    )
  )


# ------------------------------------------------------------------------------
##  Summary of category sizes
# ------------------------------------------------------------------------------
category_counts_kidney <- deg_wide_kidney %>%
  dplyr::count(category)

print(category_counts_kidney)
# category     n
# Other 13272
# Persistent   22 (same log2FC sign)
# Recovery-specific    89
# Reversible  2922


# ------------------------------------------------------------------------------
##  Extract reversible gene list
# ------------------------------------------------------------------------------
reversible_genes_kidney <- deg_wide_kidney %>%
  dplyr::filter(category == "Reversible") %>%
  dplyr::pull(Symbol)

length(reversible_genes_kidney)  # 2922



# ------------------------------------------------------------------------------
##  Build log2FC matrix for heatmap
# ------------------------------------------------------------------------------

lfc_tbl_kidney <- tibble::tibble(
  Symbol = reversible_genes_kidney,
  `NFD vs CD` = NFD_CD_shrunk_kidney[reversible_genes_kidney, "log2FoldChange"],
  `Rec vs CD` = Rec_CD_shrunk_kidney[reversible_genes_kidney, "log2FoldChange"]
) %>%
  dplyr::filter(!is.na(`NFD vs CD`) & !is.na(`Rec vs CD`))

logFC_heat_matrix_kidney <- lfc_tbl_kidney %>%
  tibble::column_to_rownames("Symbol") %>%
  as.matrix()



# ------------------------------------------------------------------------------
##  Heatmap of reversible gene log2 fold-changes
# ------------------------------------------------------------------------------

col_hm_labels_kidney <- c("NFD vs CD", "Rec vs CD")

hmap_lfc_kidney <- Heatmap(
  logFC_heat_matrix_kidney,
  name = "Gene\nLog2FC",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_row_names = FALSE,
  row_dend_reorder = TRUE,
  show_column_names = TRUE,
  column_labels = col_hm_labels_kidney,
  column_names_rot = 0,
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_width = unit(4, "cm")
  )
)

draw(
  hmap_lfc_kidney,
  heatmap_legend_side = "right"
)



# ------------------------------------------------------------------------------
##  Pathway Over-representation Analysis (ORA) Using clusterProfiler
##  Notes:
##     - Universe = all genes retained after DESeq2 prefiltering.
##     - GO requires Entrez IDs.
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
##  Define Universe (Background Gene Set)
# ------------------------------------------------------------------------------

universe_gene_list_rec_kidney <- rownames(NFD_CD_shrunk_kidney)
length(universe_gene_list_rec_kidney)   # 16305

universe_gene_list_entrez_rec_kidney <- clusterProfiler::bitr(
  universe_gene_list_rec_kidney,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

length(universe_gene_list_entrez_rec_kidney)  # 15786



# ------------------------------------------------------------------------------
##  Prepare Input DEG Data (Ensure tibble format with Symbol column)
##
##  Note:
##     DEG objects in this project (NFD_CD_adjp and Rec_CD_adjp) are already
##     stored as data frames with rownames = gene symbols. The code below shows
##     how to enforce this format if needed for other datasets.
# ------------------------------------------------------------------------------
# Convert adjp outs to tibbles with Symbol column
NFD_CD_adjp_kidney <- NFD_CD_adjp_kidney %>%
  tibble::as_tibble(rownames = "Symbol")

Rec_CD_adjp_kidney <- Rec_CD_adjp_kidney %>%
  tibble::as_tibble(rownames = "Symbol")



# ------------------------------------------------------------------------------
##  Convert DEG Symbols to Entrez IDs + Annotate Regulation
# ------------------------------------------------------------------------------
### NFD vs CD ------------------------------------------------------------------
NFD_CD_entrez_kidney <- clusterProfiler::bitr(
  NFD_CD_adjp_kidney$Symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

NFD_CD_cc_df_kidney <- NFD_CD_adjp_kidney %>%
  dplyr::mutate(
    Regulation = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"),
    Contrast = "NFDvsCD"
  ) %>%
  dplyr::rename(SYMBOL = Symbol) %>%
  dplyr::inner_join(NFD_CD_entrez_kidney, by = "SYMBOL")



### Rec vs CD ------------------------------------------------------------------
Rec_CD_entrez_kidney <- clusterProfiler::bitr(
  Rec_CD_adjp_kidney$Symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

Rec_CD_cc_df_kidney <- Rec_CD_adjp_kidney %>%
  dplyr::mutate(
    Regulation = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"),
    Contrast = "RecvsCD"
  ) %>%
  dplyr::rename(SYMBOL = Symbol) %>%
  dplyr::inner_join(Rec_CD_entrez_kidney, by = "SYMBOL")



# ------------------------------------------------------------------------------
##  ORA: GO (MF) for Each Contrast
# ------------------------------------------------------------------------------

### GO MF ORA: NFD vs CD --------------------------------------------------------
NFD_CD_cc_mf_ora_res_kidney <- clusterProfiler::compareCluster(
  ENTREZID ~ Regulation,
  data = NFD_CD_cc_df_kidney,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  minGSSize = 3,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  universe = universe_gene_list_entrez_rec_kidney
)



### GO MF ORA: Rec vs CD ---------------------------------------------------------
Rec_CD_cc_mf_ora_res_kidney <- clusterProfiler::compareCluster(
  ENTREZID ~ Regulation,
  data = Rec_CD_cc_df_kidney,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  minGSSize = 3,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  universe = universe_gene_list_entrez_rec_kidney
)




# ------------------------------------------------------------------------------
## Dotplot visualization of enriched  pathways (single contrast)
# ------------------------------------------------------------------------------
# Use clusterProfiler_ora_cc_single_dotplot to plot with title
# Display top 10 up and down pathways
# NFD vs CD (GO MF)
clusterProfiler_ora_cc_single_dotplot(
  cp_res = NFD_CD_cc_mf_ora_res_kidney,
  show_path_number = 10,
  plot_title = "NFD vs CD"
)

# Rec vs CD (GO MF)
clusterProfiler_ora_cc_single_dotplot(
  cp_res = Rec_CD_cc_mf_ora_res_kidney,
  show_path_number = 10,
  plot_title = "Rec vs CD"
)





# ------------------------------------------------------------------------------
## Add gene symbols back to the compareCluster results
# ------------------------------------------------------------------------------
# NFD vs CD (GO MF)
NFD_CD_cc_mf_ora_res_df_with_symbols_kidney <-
  add_symbols_to_compareCluster(NFD_CD_cc_mf_ora_res_kidney)

# Rec vs CD (GO MF)
Rec_CD_cc_mf_ora_res_df_with_symbols_kidney <-
  add_symbols_to_compareCluster(Rec_CD_cc_mf_ora_res_kidney)


# ------------------------------------------------------------------------------
## Summarize the results
# ------------------------------------------------------------------------------

# NFD vs CD (MF ORA)
summarize_ora_results(NFD_CD_cc_mf_ora_res_df_with_symbols_kidney)
# 97 total pathways, 20 up, 77 down

# Rec vs CD (MF ORA)
summarize_ora_results(Rec_CD_cc_mf_ora_res_df_with_symbols_kidney)


#-------------------------------------------------------------------------------
## Compare and Plot GO MF Pathways using single contrasts
# ------------------------------------------------------------------------------
# Prepare and label contrast-specific results
# Label each result with its contrast name
NFD_CD_cc_mf_ora_res_df_with_symbols_kidney$Contrast <- "NFD vs CD"
Rec_CD_cc_mf_ora_res_df_with_symbols_kidney$Contrast <- "Rec vs CD"

# Combine into one data frame
combined_df_rec_kidney <- bind_rows(
  NFD_CD_cc_mf_ora_res_df_with_symbols_kidney,
  Rec_CD_cc_mf_ora_res_df_with_symbols_kidney
)

# Identify shared pathways (same ID + Regulation across contrasts)
shared_ids_rec_kidney <- intersect(
  dplyr::select(NFD_CD_cc_mf_ora_res_df_with_symbols_kidney, ID, Regulation),
  dplyr::select(Rec_CD_cc_mf_ora_res_df_with_symbols_kidney, ID, Regulation)
)

# Keep only rows matching shared IDs, retain all columns
filtered_combined_rec_kidney <- semi_join(
  combined_df_rec_kidney,
  shared_ids_rec_kidney,
  by = c("ID", "Regulation"))



# ------------------------------------------------------------------------------
# Plot the results using dotplot_from_df_multi_contrasts
#-------------------------------------------------------------------------------

dotplot_from_df_multi_contrasts(
  cp_df = filtered_combined_rec_kidney,
  show_path_number = 15,
  plot_title = "Non-reversible pathways"
)



# Summarize the shared pathway results using the summarize_shared_pathways fx
summary_mf_overlap_rec_kidney <- summarize_shared_pathways(
  df1 = NFD_CD_cc_mf_ora_res_df_with_symbols_kidney,
  df2 = Rec_CD_cc_mf_ora_res_df_with_symbols_kidney,
  contrast1_name = "NFDvsCD",
  contrast2_name = "RecvsCD"
)

print(summary_mf_overlap_rec_kidney)



# ------------------------------------------------------------------------------
# Identify unique pathways between two contrasts
# ------------------------------------------------------------------------------

# Named list of GO MF ORA result tables
dfs_mf_rec_kidney <- list(
  NFDvsCD = NFD_CD_cc_mf_ora_res_df_with_symbols_kidney,
  RecvsCD = Rec_CD_cc_mf_ora_res_df_with_symbols_kidney
)

# Unique pathways per contrast using unique_cp_pathways fx
mf_unique_rec_kidney <- unique_cp_pathways(
  dfs_mf_rec_kidney,
  key_cols = c("ID", "Regulation")
)

# Summary table using summarize_unique_cp_pathways fx
summary_mf_unique_rec_kidney <- summarize_unique_cp_pathways(mf_unique_rec_kidney)

print(summary_mf_unique_rec_kidney)



# Plot the results -------------------------------------------------------------
# Filter unique pathways for one contrast
mf_unique_filtered_nfd_cd_kidney <- dplyr::filter(
  mf_unique_rec_kidney,
  Contrast == "NFDvsCD")

# Plot using clusterProfiler_ora_cc_unique_dotplot
clusterProfiler_ora_cc_unique_dotplot(
  mf_unique_filtered_nfd_cd_kidney,
  show_path_number = 10,
  plot_title = "Reversible pathways"
)



mf_unique_filtered_rec_cd_kidney <- dplyr::filter(
  mf_unique_rec_kidney,
  Contrast == "RecvsCD")

# Plot using clusterProfiler_ora_cc_unique_dotplot
clusterProfiler_ora_cc_unique_dotplot(
  mf_unique_filtered_rec_cd_kidney,
  show_path_number = 10,
  plot_title = "Recovery-specific pathways"
)



# sessionInfo()
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
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C 