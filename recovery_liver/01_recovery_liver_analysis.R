## -----------------------------------------------------------------------------
## Script Name: 01_recovery_liver_analysis.R
## Project: Bulk RNA-seq Analysis of ANDY CD/NFD and young/old WT mice
## Study: NFD recovery
## Analysis Type: Set-up, differential expression analysis, heatmap
##                visualization, pathway (GO MF) over-representation analysis
## Tissues: Liver
## Dependencies: Requires 00_recovery_liver_packages_functions.R script
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

counts_data_rec_main_liver <- read.csv(counts_file)

# Contains 57,132 genes across 16 samples
dim(counts_data_rec_main_liver)



# ------------------------------------------------------------------------------
##  Map Ensembl Gene IDs to Gene Symbols (biomaRt)
# ------------------------------------------------------------------------------

# Extract Ensembl gene IDs (with version)
ensembl.ids <- counts_data_rec_main_liver$Geneid
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
counts_data_rec_main_liver <- counts_data_rec_main_liver %>%
  dplyr::rename(ensembl_gene_id_version = Geneid) %>%
  dplyr::left_join(symbol, by = "ensembl_gene_id_version") %>%
  dplyr::filter(!(is.na(external_gene_name) | external_gene_name == ""))

# 56,648 genes retained after filtering out rows without symbols
dim(counts_data_rec_main_liver)



# ------------------------------------------------------------------------------
##  Handle Duplicated Gene Symbols
# ------------------------------------------------------------------------------

# Identify duplicated gene symbols
duplicated_symbols <- counts_data_rec_main_liver$external_gene_name[
  duplicated(counts_data_rec_main_liver$external_gene_name)
] %>% sort()

length(duplicated_symbols)  # 235 duplicated symbols

# Inspect duplicated symbols
dup_symbols_df <- counts_data_rec_main_liver %>%
  dplyr::filter(external_gene_name %in% duplicated_symbols)

# Remove duplicates (keep first occurrence)
counts_data_rec_main_liver <- counts_data_rec_main_liver %>%
  dplyr::distinct(external_gene_name, .keep_all = TRUE)

# 56,413 genes x 17 columns after deduplication
dim(counts_data_rec_main_liver)



# ------------------------------------------------------------------------------
##  Prepare Final Count Matrix
# ------------------------------------------------------------------------------

# Set gene symbols as rownames
rownames(counts_data_rec_main_liver) <- counts_data_rec_main_liver$external_gene_name

# Remove Ensembl ID and gene symbol columns
counts_data_rec_main_liver <- counts_data_rec_main_liver[, -c(1, 17)]

# Result should be 56,413 genes x 15 samples
dim(counts_data_rec_main_liver)



# ------------------------------------------------------------------------------
##  Load and Align Sample Metadata
# ------------------------------------------------------------------------------

# Read metadata (sample info table)
sampleinfo_file <- "path/to/sample_info.csv"  # sample metadata, deposited with GEO submission

colData_rec_main_liver <- read.csv(
  sampleinfo_file,
  header = TRUE
)

head(colData_rec_main_liver)

# Move sample names into rownames
rownames(colData_rec_main_liver) <- colData_rec_main_liver[, 1]
colData_rec_main_liver <- colData_rec_main_liver[, -1]

# Ensure sample names match between metadata and count matrix
# Required for DESeq2
all(colnames(counts_data_rec_main_liver) %in% rownames(colData_rec_main_liver)) # TRUE
all(colnames(counts_data_rec_main_liver) == rownames(colData_rec_main_liver)) # TRUE 

# Distribution across conditions
table(colData_rec_main_liver$Condition)
# CD: 6    NFD: 4    Recovery: 5



# ------------------------------------------------------------------------------
##  Convert Counts to Matrix for DESeq2
# ------------------------------------------------------------------------------

counts_data_matrix_rec_liver <- data.matrix(counts_data_rec_main_liver)

head(counts_data_matrix_rec_liver)


# ------------------------------------------------------------------------------
##  Create a SummarizedExperiment Object
# ------------------------------------------------------------------------------

# Construct SummarizedExperiment for the NFD recovery liver samples
se_main_rec_liver <- SummarizedExperiment(
  assays = list(counts = counts_data_matrix_rec_liver),
  colData = colData_rec_main_liver
)


# ------------------------------------------------------------------------------
##  Prefiltering
# ------------------------------------------------------------------------------

# Create a DESeqDataSet for the NFD recovery study
dds_main_rec_liver <- DESeqDataSet(se_main_rec_liver, design = ~ Condition)

# Confirm initial gene count (expected 56,413)
nrow(dds_main_rec_liver)



# ------------------------------------------------------------------------------
##  Prefilter Lowly Expressed Genes
# ------------------------------------------------------------------------------

# Determine smallest group size (used for filtering threshold)
smallestGroupSize_rec_liver <- min(table(colData(dds_main_rec_liver)$Condition))
print(smallestGroupSize_rec_liver) # 4

# Keep genes with count ≥ 10 in at least `smallestGroupSize` samples
keep_rec_liver <- rowSums(counts(dds_main_rec_liver) >= 10) >= smallestGroupSize_rec_liver
dds_main_rec_liver <- dds_main_rec_liver[keep_rec_liver, ]

# 14,500 genes should remain
nrow(dds_main_rec_liver)



# ------------------------------------------------------------------------------
##  Estimate Size Factors (Normalization)
# ------------------------------------------------------------------------------

# Normalize for sequencing depth differences
dds_main_rec_liver <- estimateSizeFactors(dds_main_rec_liver)



# ------------------------------------------------------------------------------
##  Differential Expression Analysis with DESeq2
# ------------------------------------------------------------------------------

# Set reference level ("CD") and run DESeq2 pipeline
# ------------------------------------------------------------------------------
dds_main_rec_liver$Condition <- relevel(dds_main_rec_liver$Condition, ref = "CD")

dds_main_rec_liver <- DESeq(dds_main_rec_liver)

# Review size factors and dispersion estimates
sizeFactors(dds_main_rec_liver)
plotDispEsts(dds_main_rec_liver)

# Inspect model coefficients / contrast names
resultsNames(dds_main_rec_liver)
# Condition_NFD_vs_CD
# Condition_Recovery_vs_CD


# ------------------------------------------------------------------------------
##  Contrast: NFD vs CD (recovery study)
# ------------------------------------------------------------------------------

# Extract DESeq2 results
NFD_CD <- results(dds_main_rec_liver, contrast = c("Condition", "NFD", "CD"))

# Apply LFC shrinkage (apeglm)
NFD_CD_shrunk <- lfcShrink(
  dds_main_rec_liver,
  coef = "Condition_NFD_vs_CD",
  res  = NFD_CD,
  type = "apeglm"
) %>% 
  as.data.frame()

# Identify significant DEGs (FDR ≤ 0.05)
NFD_CD_adjp <- subset(NFD_CD_shrunk, padj <= 0.05)

# Summary table using summarize_DEGs function
summary_NFD_CD <- summarize_DEGs(NFD_CD_shrunk, contrast_name = "NFD vs CD")
summary_NFD_CD
# Total_DEGs = 3695
# Upregulated = 1809
# Downregulated = 1886


# ------------------------------------------------------------------------------
##  Contrast: Recovery vs CD
# ------------------------------------------------------------------------------

# Extract DESeq2 results
Rec_CD <- results(dds_main_rec_liver, contrast = c("Condition", "Recovery", "CD"))

# Apply LFC shrinkage with apeglm
Rec_CD_shrunk <- lfcShrink(
  dds_main_rec_liver,
  coef = "Condition_Recovery_vs_CD",
  res  = Rec_CD,
  type = "apeglm"
) %>% 
  as.data.frame()

# Significant DEGs (padj ≤ 0.05)
Rec_CD_adjp <- subset(Rec_CD_shrunk, padj <= 0.05)


# Summary table using summarize_DEGs function
summary_Rec_CD <- summarize_DEGs(Rec_CD_shrunk, contrast_name = "Recovery vs CD")
summary_Rec_CD
# Total_DEGs = 343
# Upregulated = 144
# Downregulated = 199


# ------------------------------------------------------------------------------
##  ComplexHeatmap Visualization of Top NFD vs CD DEGs
# ------------------------------------------------------------------------------

# Select top 1000 DEGs ranked by adjusted p-value
# ------------------------------------------------------------------------------
NFD_CD_adjp <- NFD_CD_adjp[order(NFD_CD_adjp$padj), ]
NFD_CD_top1000 <- NFD_CD_adjp[1:1000, ]
dim(NFD_CD_top1000)  # 1000 × 5 columns



# ------------------------------------------------------------------------------
##  Subset Normalized Counts for Heatmap
# ------------------------------------------------------------------------------

# Retrieve normalized counts from DESeq2 object
normalized_counts_main_rec_liver <- counts(dds_main_rec_liver, normalized = TRUE)

# Subset to top DEG list
top_genes <- rownames(NFD_CD_top1000)

norm_cts_sub <- normalized_counts_main_rec_liver[
  rownames(normalized_counts_main_rec_liver) %in% top_genes,
]

dim(norm_cts_sub)  # 1000 × 15 samples

# Row-wise z-score transformation (required for heatmap visualization)
heat_mat <- t(scale(t(norm_cts_sub)))



# ------------------------------------------------------------------------------
##  Column Annotation (Condition)
# ------------------------------------------------------------------------------

condition_labels_main_rec_liver <- colData_rec_main_liver$Condition_Sub
condition_labels_main_rec_liver <- factor(
  condition_labels_main_rec_liver,
  levels = c("CD", "NFD", "Rec"),
  labels = c("CD", "NFD", "Rec")
)

column_ha <- HeatmapAnnotation(
  Condition = condition_labels_main_rec_liver,
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
dend <- hclust(dist(t(heat_mat))) %>% as.dendrogram()

# Rotate dendrogram (dataset-specific)
dend <- dendextend::rotate(dend, order = c(5:15, 1:4))

# ------------------------------------------------------------------------------
##  Plot Heatmap
# ------------------------------------------------------------------------------

hmap_main_rec_liver <- Heatmap(
  heat_mat,
  name = "Gene\nZ-Score",
  top_annotation = column_ha,
  cluster_rows = TRUE,
  cluster_columns = dend,
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
  hmap_main_rec_liver,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)



# ------------------------------------------------------------------------------
# Prepare contrast-specific DEG data for contrast comparisons
# ------------------------------------------------------------------------------

# Use prep_deg_df function
NFD_CD <- prep_deg_df(NFD_CD_adjp, "NFDvsCD")
Rec_CD <- prep_deg_df(Rec_CD_adjp, "RecvsCD")


# ------------------------------------------------------------------------------
##  Contrast Comparisons: NFD vs CD and Recovery vs CD
# ------------------------------------------------------------------------------
# Using: compare_contrasts_multi()
# Inputs: NFD_CD (NFDvsCD DEGs), Rec_CD (RecvsCD DEGs)
# ------------------------------------------------------------------------------
##  Shared DEGs (any direction; no filtering by sign)
# ------------------------------------------------------------------------------
shared_all_NFDvsCD_RecvsCD <- compare_contrasts_multi(
  dfs = list(NFD_CD, Rec_CD),
  contrast_names = c("NFDvsCD", "RecvsCD"),
  direction = "all"
)

length(shared_all_NFDvsCD_RecvsCD$Symbol)
# 286 shared DEGs



# ------------------------------------------------------------------------------
##  Shared Upregulated DEGs (consistent direction)
# ------------------------------------------------------------------------------
shared_up_NFDvsCD_RecvsCD <- compare_contrasts_multi(
  dfs = list(NFD_CD, Rec_CD),
  contrast_names = c("NFDvsCD", "RecvsCD"),
  direction = "up"
)

length(shared_up_NFDvsCD_RecvsCD$Symbol)
# 119 shared upregulated DEGs



# ------------------------------------------------------------------------------
##  Shared Downregulated DEGs (consistent direction)
# ------------------------------------------------------------------------------
shared_down_NFDvsCD_RecvsCD <- compare_contrasts_multi(
  dfs = list(NFD_CD, Rec_CD),
  contrast_names = c("NFDvsCD", "RecvsCD"),
  direction = "down"
)

length(shared_down_NFDvsCD_RecvsCD$Symbol)
# 167 shared downregulated DEGs




# ------------------------------------------------------------------------------
##  Shared DEGs with Opposite Directions (discordant log2FC)
# ------------------------------------------------------------------------------
shared_opposite_NFDvsCD_RecvsCD <- compare_contrasts_multi(
  dfs = list(NFD_CD, Rec_CD),
  contrast_names = c("NFDvsCD", "RecvsCD"),
  direction = "opposite"
)

length(shared_opposite_NFDvsCD_RecvsCD$Symbol)
# 0 genes with opposite directions




# ------------------------------------------------------------------------------
## Tidy summary of shared DEGs across contrasts
# ------------------------------------------------------------------------------

shared_deg_summary_rec_liver <- tibble::tibble(
  Comparison = "NFDvsCD + RecvsCD",
  Shared_Total = length(shared_all_NFDvsCD_RecvsCD$Symbol),
  Shared_Up = length(shared_up_NFDvsCD_RecvsCD$Symbol),
  Shared_Down = length(shared_down_NFDvsCD_RecvsCD$Symbol),
  Shared_Opposite = length(shared_opposite_NFDvsCD_RecvsCD$Symbol)
)

shared_deg_summary_rec_liver
# Comparison        Shared_Total Shared_Up Shared_Down Shared_Opposite
# NFDvsCD + RecvsCD          286       119         167               0


# ------------------------------------------------------------------------------
##  Identify and Summarize DEGs Unique to Each Contrast
# ------------------------------------------------------------------------------

dfs_unique_main_rec_liver <- list(
  NFDvsCD = NFD_CD,
  RecvsCD = Rec_CD
)

# Use the get_unique_DEGs_multi and summarize_unique_DEGs_multi functions
unique_deg_rows_main_rec_liver <- get_unique_DEGs_multi(dfs_unique_main_rec_liver)
unique_deg_summary_main_rec_liver <- summarize_unique_DEGs_multi(unique_deg_rows_main_rec_liver)

print(unique_deg_summary_main_rec_liver)
# Contrast Total Upregulated Downregulated
# NFDvsCD   3409        1690          1719
# RecvsCD     57          25            32



# ------------------------------------------------------------------------------
##  Hypergeometric Test of DEG Overlap Using GeneOverlap
##  Contrast tested here:
##     NFD vs CD  ∩  Recovery vs CD
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
##  Define Universe (Background Gene Set)
# ------------------------------------------------------------------------------
# Use rownames of filtered DESeq2 object as expressed gene universe
background_genes_rec_liver <- rownames(dds_main_rec_liver)
universe_size_rec_liver <- length(background_genes_rec_liver)

cat("Total background genes:", universe_size_rec_liver, "\n")
# Expected: 14,500 genes after prefiltering


# ------------------------------------------------------------------------------
##  Extract DEG Lists for Testing
# ------------------------------------------------------------------------------
# Ensure Symbol column exists for each DEG table
NFD_CD_genes <- NFD_CD$Symbol
Rec_CD_genes <- Rec_CD$Symbol

cat("NFD vs CD DEGs:", length(NFD_CD_genes), "\n") # 3695
cat("Rec vs CD DEGs:", length(Rec_CD_genes), "\n") # 343


# ------------------------------------------------------------------------------
##  Perform Hypergeometric Test (GeneOverlap)
# ------------------------------------------------------------------------------

go.obj_rec_liver <- newGeneOverlap(
  listA = NFD_CD_genes,
  listB = Rec_CD_genes,
  genome.size = universe_size_rec_liver
)

go.obj_rec_liver <- testGeneOverlap(go.obj_rec_liver)

# Display statistical results
print(go.obj_rec_liver)

# Detailed information about this GeneOverlap object:
# listA size=3695, e.g. Cyp2c40 Olig1 C3
# listB size=343, e.g. Slc22a18 Sox9 Spg7
# Intersection size=286, e.g. Cyp2c40 Olig1 Plxnb1
# Union size=3752, e.g. Cyp2c40 Olig1 C3
# Genome size=14500
# # Contingency Table:
# notA  inA
# notB 10748 3409
# inB     57  286
# Overlapping p-value=3e-115
# Odds ratio=15.8
# Overlap tested using Fisher's exact test (alternative=greater)
# Jaccard Index=0.1



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

padj_cutoff <- 0.05

deg_wide <- NFD_CD_shrunk %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Symbol") %>%
  dplyr::select(Symbol, L_NFD = log2FoldChange, padj_NFD = padj) %>%
  dplyr::full_join(
    Rec_CD_shrunk %>%
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
category_counts <- deg_wide %>%
  dplyr::count(category)

print(category_counts)
# category     n
# Other 10748
# Persistent   286
# Recovery-specific    57
# Reversible  3409


# ------------------------------------------------------------------------------
##  Extract reversible gene list
# ------------------------------------------------------------------------------
reversible_genes <- deg_wide %>%
  dplyr::filter(category == "Reversible") %>%
  dplyr::pull(Symbol)

length(reversible_genes)  # 3409



# ------------------------------------------------------------------------------
##  Build log2FC matrix for heatmap
# ------------------------------------------------------------------------------

lfc_tbl <- tibble::tibble(
  Symbol = reversible_genes,
  `NFD vs CD` = NFD_CD_shrunk[reversible_genes, "log2FoldChange"],
  `Rec vs CD` = Rec_CD_shrunk[reversible_genes, "log2FoldChange"]
) %>%
  dplyr::filter(!is.na(`NFD vs CD`) & !is.na(`Rec vs CD`))

logFC_heat_matrix <- lfc_tbl %>%
  tibble::column_to_rownames("Symbol") %>%
  as.matrix()



# ------------------------------------------------------------------------------
##  Heatmap of reversible gene log2 fold-changes
# ------------------------------------------------------------------------------

col_hm_labels <- c("NFD vs CD", "Rec vs CD")

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
  column_names_rot = 0,
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_width = unit(4, "cm")
  )
)

draw(
  hmap_lfc,
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

universe_gene_list_rec_liver <- rownames(NFD_CD_shrunk)
length(universe_gene_list_rec_liver)   # 14,500

universe_gene_list_entrez_rec_liver <- clusterProfiler::bitr(
  universe_gene_list_rec_liver,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
) %>%
  dplyr::pull(ENTREZID) %>%
  unique()

length(universe_gene_list_entrez_rec_liver)  # 14,133



# ------------------------------------------------------------------------------
##  Prepare Input DEG Data (Ensure tibble format with Symbol column)
##
##  Note:
##     DEG objects in this project (NFD_CD_adjp and Rec_CD_adjp) are already
##     stored as data frames with rownames = gene symbols. The code below shows
##     how to enforce this format if needed for other datasets.
# ------------------------------------------------------------------------------
# Convert adjp outs to tibbles with Symbol column
NFD_CD_adjp <- NFD_CD_adjp %>%
  tibble::as_tibble(rownames = "Symbol")

Rec_CD_adjp <- Rec_CD_adjp %>%
  tibble::as_tibble(rownames = "Symbol")



# ------------------------------------------------------------------------------
##  Convert DEG Symbols to Entrez IDs + Annotate Regulation
# ------------------------------------------------------------------------------
### NFD vs CD ------------------------------------------------------------------
NFD_CD_entrez <- clusterProfiler::bitr(
  NFD_CD_adjp$Symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

NFD_CD_cc_df <- NFD_CD_adjp %>%
  dplyr::mutate(
    Regulation = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"),
    Contrast = "NFDvsCD"
  ) %>%
  dplyr::rename(SYMBOL = Symbol) %>%
  dplyr::inner_join(NFD_CD_entrez, by = "SYMBOL")



### Rec vs CD ------------------------------------------------------------------
Rec_CD_entrez <- clusterProfiler::bitr(
  Rec_CD_adjp$Symbol,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

Rec_CD_cc_df <- Rec_CD_adjp %>%
  dplyr::mutate(
    Regulation = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated"),
    Contrast   = "RecvsCD"
  ) %>%
  dplyr::rename(SYMBOL = Symbol) %>%
  dplyr::inner_join(Rec_CD_entrez, by = "SYMBOL")



# ------------------------------------------------------------------------------
##  ORA: GO (MF) for Each Contrast
# ------------------------------------------------------------------------------

### GO MF ORA: NFD vs CD --------------------------------------------------------
NFD_CD_cc_mf_ora_res <- clusterProfiler::compareCluster(
  ENTREZID ~ Regulation,
  data = NFD_CD_cc_df,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  minGSSize = 3,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  universe = universe_gene_list_entrez_rec_liver
)



### GO MF ORA: Rec vs CD --------------------------------------------------------
Rec_CD_cc_mf_ora_res <- clusterProfiler::compareCluster(
  ENTREZID ~ Regulation,
  data = Rec_CD_cc_df,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  minGSSize = 3,
  maxGSSize = 800,
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  universe = universe_gene_list_entrez_rec_liver
)




# ------------------------------------------------------------------------------
## Dotplot visualization of enriched pathways (single contrast)
# ------------------------------------------------------------------------------
# Use clusterProfiler_ora_cc_single_dotplot to plot with title
# Display top 10 up and down pathways
# NFD vs CD (GO MF)
clusterProfiler_ora_cc_single_dotplot(
  cp_res = NFD_CD_cc_mf_ora_res,
  show_path_number = 10,
  plot_title = "NFD vs CD"
)

# Rec vs CD (GO MF)
clusterProfiler_ora_cc_single_dotplot(
  cp_res = Rec_CD_cc_mf_ora_res,
  show_path_number = 10,
  plot_title = "Rec vs CD"
)



# ------------------------------------------------------------------------------
## Add gene symbols back to the compareCluster results
# ------------------------------------------------------------------------------
# NFD vs CD (GO MF)
NFD_CD_cc_mf_ora_res_df_with_symbols <-
  add_symbols_to_compareCluster(NFD_CD_cc_mf_ora_res)

colnames(NFD_CD_cc_mf_ora_res_df_with_symbols)

# Rec vs CD (GO MF)
Rec_CD_cc_mf_ora_res_df_with_symbols <-
  add_symbols_to_compareCluster(Rec_CD_cc_mf_ora_res)



# ------------------------------------------------------------------------------
## Summarize the results
# ------------------------------------------------------------------------------

# NFD vs CD (GO MF ORA)
summarize_ora_results(NFD_CD_cc_mf_ora_res_df_with_symbols)
# 158 total pathways, 28 up, 130 down

# Rec vs CD (GO MF ORA)
summarize_ora_results(Rec_CD_cc_mf_ora_res_df_with_symbols)



#-------------------------------------------------------------------------------
## Compare and Plot GO MF Pathways using single contrasts
# ------------------------------------------------------------------------------
# Prepare and label contrast-specific results
# Label each result with its contrast name
NFD_CD_cc_mf_ora_res_df_with_symbols$Contrast <- "NFD vs CD"
Rec_CD_cc_mf_ora_res_df_with_symbols$Contrast <- "Rec vs CD"

# Combine into one data frame
combined_df_rec_liver <- bind_rows(
  NFD_CD_cc_mf_ora_res_df_with_symbols,
  Rec_CD_cc_mf_ora_res_df_with_symbols
)

# Identify shared pathways (same ID + Regulation across contrasts)
shared_ids_rec_liver <- intersect(
  dplyr::select(NFD_CD_cc_mf_ora_res_df_with_symbols, ID, Regulation),
  dplyr::select(Rec_CD_cc_mf_ora_res_df_with_symbols, ID, Regulation)
)

# Keep only rows matching shared IDs, retain all columns
filtered_combined_rec_liver <- semi_join(
  combined_df_rec_liver,
  shared_ids_rec_liver,
  by = c("ID", "Regulation"))


# ------------------------------------------------------------------------------
# Plot the results using dotplot_from_df_multi_contrasts
# ------------------------------------------------------------------------------

dotplot_from_df_multi_contrasts(
  cp_df = filtered_combined_rec_liver,
  show_path_number = 10,
  plot_title = "Non-reversible pathways"
)


# Summarize the shared pathway results using the summarize_shared_pathways fx
summary_mf_overlap_rec_liver <- summarize_shared_pathways(
  df1 = NFD_CD_cc_mf_ora_res_df_with_symbols,
  df2 = Rec_CD_cc_mf_ora_res_df_with_symbols,
  contrast1_name = "NFDvsCD",
  contrast2_name = "RecvsCD"
)




# ------------------------------------------------------------------------------
# Identify unique pathways between two contrasts
# ------------------------------------------------------------------------------

# Named list of GO MF ORA result tables
dfs_mf_rec_liver <- list(
  NFDvsCD = NFD_CD_cc_mf_ora_res_df_with_symbols,
  RecvsCD = Rec_CD_cc_mf_ora_res_df_with_symbols
)

# Unique pathways per contrast using unique_cp_pathways fx
mf_unique_rec_liver <- unique_cp_pathways(
  dfs_mf_rec_liver,
  key_cols = c("ID", "Regulation")
)

# Summary table using summarize_unique_cp_pathways fx
summary_mf_unique_rec_liver <- summarize_unique_cp_pathways(mf_unique_rec_liver)



# Plot the results
# Filter unique pathways for one contrast
mf_unique_filtered_nfd_cd <- dplyr::filter(
  mf_unique_rec_liver,
  Contrast == "NFDvsCD")

# Plot using clusterProfiler_ora_cc_unique_dotplot
clusterProfiler_ora_cc_unique_dotplot(
  mf_unique_filtered_nfd_cd,
  show_path_number = 10,
  plot_title = "Reversible pathways"
)




# sessionInfo() platform details (package versions are listed in the repository README)
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