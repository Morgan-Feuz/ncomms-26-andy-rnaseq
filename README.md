# ANDY NAD-Deficiency Mouse Model: RNA-seq Analysis Code

Custom analysis code accompanying:

> Benitez-Rosendo, A. & Feuz, M. B. et al. Inducible chronic NAD deficiency in mice reveals multi-systemic phenotypical, metabolic, and transcriptional changes. *Nature Communications* (in revision; manuscript NCOMMS-26-053722-T). Full citation and DOI will be added here once published.

This repository contains the downstream RNA-seq differential expression, visualization, and
pathway-enrichment code used to generate the transcriptomic results in the paper. It does
**not** contain raw or processed sequencing data. See **Data availability** below for where
that lives.

## Repository structure

| Folder | Cohort | Tissue |
|---|---|---|
| `recovery_liver/` | Recovery study (CD / NFD / Rec) | Liver |
| `recovery_kidney/` | Recovery study (CD / NFD / Rec) | Kidney |
| `aging_liver/` | Aging comparison (CD / NFD / young WT / old WT) | Liver |
| `aging_kidney/` | Aging comparison (CD / NFD / young WT / old WT) | Kidney |
| `venn_diagram/` | Cross-cohort feature overlap | n/a |

Each `00_..._packages_functions.R` / `01_..._analysis.R` pair covers the full analysis for
that cohort/tissue (differential expression, visualization, and enrichment together), rather
than being tied to a single figure, so this table is organized by dataset, not by figure.

Each of the four cohort/tissue folders contains two scripts, run in order:

| Folder | Setup script | Analysis script |
|---|---|---|
| `recovery_liver/` | `00_recovery_liver_packages_functions.R` | `01_recovery_liver_analysis.R` |
| `recovery_kidney/` | `00_recovery_kidney_packages_functions.R` | `01_recovery_kidney_analysis.R` |
| `aging_liver/` | `00_aging_liver_packages_functions.R` | `01_aging_liver_analysis.R` |
| `aging_kidney/` | `00_aging_kidney_packages_functions.R` | `01_aging_kidney_analysis.R` |

- The `00_..._packages_functions.R` script loads required packages and defines any custom
  functions used by the analysis script in the same folder.
- The `01_..._analysis.R` script is the differential expression / visualization / enrichment
  analysis itself.

**`aging_kidney/01_aging_kidney_analysis.R` is the one exception to folders being
independent**: the NFD and CD samples used in that script come from the recovery cohort (a
different sequencing batch from the Old/Young samples analyzed there), so it reuses DESeq2
results and normalized counts computed by `recovery_kidney/01_recovery_kidney_analysis.R`
rather than recomputing them. Run `recovery_kidney/01_recovery_kidney_analysis.R` first, in
the same R session, before running `aging_kidney/01_aging_kidney_analysis.R`.

`venn_diagram/` contains a single standalone Python script (see **Requirements** below, since
this is the only folder that isn't R).

## Requirements

### R (all folders except `venn_diagram/`)

Analysis was run in **R v4.5.1** with the following packages. Versions below are taken
directly from the paper's Methods section; if your installed versions differ, note that in
this README rather than silently letting it drift out of sync with what's published.

| Package | Version | Used for |
|---|---|---|
| DESeq2 | 1.48.1 | Differential expression |
| biomaRt | 2.64.0 | Ensembl ID → gene symbol conversion |
| sva (ComBat-Seq) | 3.56.0 | Batch correction (liver, young WT libraries) |
| apeglm | 1.30.0 | log2FC shrinkage (DESeq2 `lfcShrink`) |
| ComplexHeatmap | 2.24.1 | Heatmap visualization |
| circlize | 0.4.16 | Heatmap color mapping |
| dendextend | 1.19.1 | Dendrogram rendering |
| dendsort | 0.3.4 | Dendrogram leaf reordering |
| clusterProfiler | 4.16.0 | Pathway over-representation analysis |
| GeneOverlap | 1.44.0 | Feature-overlap significance (Fisher's exact test) |
| ggplot2 | 3.5.2 | Enrichment dot plots |
| tidyverse | 2.0.0 | General data manipulation/visualization |

### Python (`venn_diagram/` only)

- Python 3.12.2
- matplotlib-venn 1.1.1

A `requirements.txt` in that folder listing `matplotlib-venn` is recommended so this doesn't
require cross-referencing the README to run.

## Running the code

1. Update the placeholder paths at the top of each `01_..._analysis.R` (or the Python script)
   to point to your local copy of the relevant input files (counts matrix, sample metadata,
   etc.), obtained from the repository listed under **Data availability**.
2. Within a folder, source the `00_..._packages_functions.R` script before the
   `01_..._analysis.R` script (see the table above for exact filenames).
3. `recovery_liver/`, `recovery_kidney/`, and `aging_liver/` are independent of one another
   and can be run in any order. `aging_kidney/` is the exception: run
   `recovery_kidney/01_recovery_kidney_analysis.R` first, in the same R session, before
   running `aging_kidney/01_aging_kidney_analysis.R` (see note above).

## Data availability

Raw and processed sequencing data are **not** included in this repository. They are
deposited in:

- Gene Expression Omnibus (GEO), accession **[pending]**: raw and processed RNA-seq data
  (counts matrix, sample metadata).

See the manuscript's Data Availability statement for the current accession code(s).

## License

MIT License. See `LICENSE` for the full text.
