# ==============================================================================
# 01_data_ingestion_JM.R
# ==============================================================================
# Load STARsolo count matrices, attach patch-clamp metadata, perform ambient
# RNA decontamination (decontX), doublet detection (Scrublet), assign cell
# types via label transfer, and produce the annotated Seurat object for the
# JM_patchSeq study.
#
# Corresponds to Methods: "Single-cell RNA sequencing" and
# "Integrating Patch-seq Datasets"
#   "For each study separately, ambient RNA contamination was removed using
#    decontX (celda package)... Putative doublets were identified and removed
#    using Scrublet (expected doublet rate = 0.06; cells with predicted_doublets
#    = TRUE excluded), with cells exceeding 50% ambient contamination also
#    removed."
#
# Input:
#   ./star_counts_human/                         - STARsolo output per cell
#   ./metadata_files/*_metadata.csv              - Patch-clamp metadata per plate
#   pclamp_patched_all                           - Reference patch-seq object
#     (from scripts/02_integrate_patchseq_datasets.R — required for label transfer)
#
# Output:
#   20241017_JM_islets_seurat_object_ensg_decon.RData
#     - JM_seurobj_ensg: annotated Seurat object (JM study cells only)
#       with decontX-corrected counts, Scrublet doublet scores, cell type
#       labels, QC metrics, and circadian module score
#     - JM_seurobj_islets: JM primary islet cells only (subset)
#   20241017_human_meta.RData
#     - human_meta: combined metadata data frame with cell type assignments
#
# Dependencies:
#   tidyverse, Seurat, EnsDb.Hsapiens.v86, ensembldb, celda, reticulate
#
#
# Note on cell type annotation:
#   Cell types are assigned via a two-step approach:
#   (1) Seurat label transfer from pclamp_patched_all (reference)
#   (2) Manual correction for low-confidence predictions using marker gene
#       expression thresholds (INS, GCG, SST, PPY, GHRL, CPA1, KRT19)
#   This creates a dependency on pclamp_patched_all; if reproducing from
#   scratch, scripts/02_integrate_patchseq_datasets.R must be run first.
# ==============================================================================

library("tidyverse")
library("EnsDb.Hsapiens.v86")
library("ensembldb")
library("Seurat")
library("celda")
library("reticulate")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 1: Load STARsolo count matrices
# ═══════════════════════════════════════════════════════════════════════════════
# STARsolo outputs one tab-delimited count file per cell, organized by plate.
# Files contain unstranded, forward, and reverse strand counts; we use
# unstranded counts throughout.

cellIDs <- list.files("./star_counts_human/", pattern = "*.tab", recursive = TRUE)
cells   <- unique(str_split_fixed(cellIDs, "/", 2)[, 2])

# Helper: read a single STARsolo count file
# skip=4 skips the summary lines at the top of each file
read_count_output_tab <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- read.table(
    paste0(dir, "/", name), skip = 4, row.names = NULL,
    col.names = c("gene", "unstranded", "forward", "reverse")
  )
  return(m)
}

# Build count matrix by iterating over all cell files
# Gene IDs have version suffixes (e.g. ENSG00000254647.8) — strip to base ID
JM_mat <- read_count_output_tab("./star_counts_human/", name = cellIDs[1])
JM_mat$gene <- str_split_fixed(JM_mat$gene, "\\.", 2)[, 1]
name <- str_split_fixed(cells[1], "[._]Read", 2)[1]
colnames(JM_mat)[2:4] <- paste(name, colnames(JM_mat)[2:4], sep = "-")

for (c in cellIDs[-1]) {
  tmp_mat <- read_count_output_tab("./star_counts_human/", c)
  tmp_mat$gene <- str_split_fixed(tmp_mat$gene, "\\.", 2)[, 1]
  name <- str_split_fixed(c, "/", 2)[, 2]
  name <- str_split_fixed(name, "[._]Read", 2)[1]
  colnames(tmp_mat)[2:4] <- paste(name, colnames(tmp_mat)[2:4], sep = "-")
  JM_mat <- full_join(JM_mat, tmp_mat, by = "gene")
}

genes <- JM_mat$gene

# Keep unstranded counts only; set Ensembl IDs as rownames
JM_mat <- dplyr::select(JM_mat, contains("unstranded"))
colnames(JM_mat) <- gsub("-unstranded", "", colnames(JM_mat))
rownames(JM_mat) <- genes

cat("Count matrix dimensions:", nrow(JM_mat), "genes ×", ncol(JM_mat), "cells\n")
cat("Library size summary:\n")
print(summary(colSums(JM_mat)))

save(JM_mat, file = "20240715_JM_GFP_mat_human.RData")
#load("20240715_JM_GFP_mat_human.RData")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 2: Load and combine patch-clamp metadata
# ═══════════════════════════════════════════════════════════════════════════════
# Each recording plate has its own metadata CSV. These are loaded, cleaned,
# and combined into a single data frame aligned to the count matrix.

Nov2021_Plate1_meta <- read_csv("metadata_files/2021_November_metadata.csv") %>%
  dplyr::arrange(Plate, Well) %>%
  dplyr::select(-`Series resistance`) %>%
  dplyr::filter(!is.na(NAME_FROM_MBSU))

Dec2021_Plate2_meta <- read_csv("metadata_files/December_3_2021_Metadata.csv") %>%
  dplyr::arrange(Plate, Well) %>%
  dplyr::select(-`Series resistance`) %>%
  dplyr::filter(!is.na(NAME_FROM_MBSU))

Nov2022_meta <- read_csv("metadata_files/2022_November_metadata.csv") %>%
  dplyr::filter(Species == "Human", is.na(Note)) %>%
  dplyr::arrange(Plate, Well) %>%
  dplyr::select(-`Series resistance`) %>%
  dplyr::filter(!is.na(NAME_FROM_MBSU))

Jan2023_meta <- read_csv("metadata_files/2023_January_metadata.csv") %>%
  dplyr::arrange(Plate, Well) %>%
  dplyr::select(-`Series resistance`, -Plating_Condition, -Treatment) %>%
  dplyr::filter(!is.na(NAME_FROM_MBSU))

Mar2023_meta <- read_csv("metadata_files/2023_March_metadata.csv") %>%
  dplyr::filter(Species == "Human") %>%
  dplyr::arrange(Plate, Well) %>%
  dplyr::select(-`Series resistance`, -Plating_Condition, -Treatment) %>%
  dplyr::filter(!is.na(NAME_FROM_MBSU))

Sep2023_meta <- read_csv("metadata_files/2023_September_metadata.csv") %>%
  dplyr::arrange(Plate, Well) %>%
  dplyr::select(-`Series resistance`, -`Concentration after pre-amp(ng/ul)`,
                -Plating_Condition, -Treatment) %>%
  dplyr::filter(!is.na(NAME_FROM_MBSU))

# Convert to data frames and set rownames for Seurat compatibility
for (obj_name in c("Nov2021_Plate1_meta","Dec2021_Plate2_meta","Nov2022_meta",
                   "Jan2023_meta","Mar2023_meta","Sep2023_meta")) {
  obj <- get(obj_name) %>% as.data.frame()
  rownames(obj) <- obj$NAME_FROM_MBSU
  assign(obj_name, obj)
}

# Plate-specific fixes
Sep2023_meta$Plate <- paste0("Sep2023_Plate", Sep2023_meta$Plate)
Sep2023_meta$`Flourescence (if using reporter)` <-
  as.character(Sep2023_meta$`Flourescence (if using reporter)`)

# Standardize time-from-plating column in Mar2023 (mixed text/numeric entries)
Mar2023_meta$`Time from plating (h)` <- case_when(
  Mar2023_meta$`Time from plating (h)` == "DAY 1" ~ 24,
  Mar2023_meta$`Time from plating (h)` == "DAY 2" ~ 48,
  Mar2023_meta$`Time from plating (h)` == "DAY 3" ~ 72,
  TRUE ~ as.numeric(Mar2023_meta$`Time from plating (h)`)
)

# Standardize date format in Mar2023
Mar2023_meta$Date <- as.numeric(case_when(
  Mar2023_meta$Date == "Feb 23, 2023"   ~ "230223",
  Mar2023_meta$Date == "Feb 24, 2023"   ~ "230224",
  Mar2023_meta$Date == "Feb 28, 2023"   ~ "230228",
  Mar2023_meta$Date == "March 1, 2023"  ~ "230301",
  Mar2023_meta$Date == "March 3, 2023"  ~ "230303",
  Mar2023_meta$Date == "March 7, 2023"  ~ "230307",
  Mar2023_meta$Date == "March 8, 2023"  ~ "230308",
  Mar2023_meta$Date == "March 16, 2023" ~ "230316",
  Mar2023_meta$Date == "March 17, 2023" ~ "230317",
  TRUE ~ Mar2023_meta$Date
))

human_meta <- bind_rows(
  Nov2021_Plate1_meta, Dec2021_Plate2_meta, Nov2022_meta,
  Jan2023_meta, Mar2023_meta, Sep2023_meta
)

# Sanity checks: metadata-count matrix alignment
cat("Cells in metadata but not counts:",
    length(setdiff(rownames(human_meta), colnames(JM_mat))), "\n")  # expect 0
cat("Cells with counts but no metadata:",
    length(setdiff(colnames(JM_mat), rownames(human_meta))), "\n")

human_meta <- human_meta %>% dplyr::filter(NAME_FROM_MBSU %in% colnames(JM_mat))

# ── Create Seurat object ───────────────────────────────────────────────────────
human_seurobj <- Seurat::CreateSeuratObject(
  counts = JM_mat, assay = "RNA", meta.data = human_meta
)
# Filter to human cells with valid CellIDs (excludes mouse cells and blanks)
human_seurobj <- subset(human_seurobj, Species == "Human" & CellID != "BLANK")

cat("Cells after species/blank filter:", ncol(human_seurobj), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 3: Ambient RNA decontamination (decontX)
# ═══════════════════════════════════════════════════════════════════════════════
# decontX estimates and removes ambient RNA contamination using a Bayesian
# mixture model. Cells with >50% estimated contamination are removed.
# Raw counts are replaced with decontaminated counts for all downstream use.

human_seurobj <- NormalizeData(human_seurobj)

sce <- as.SingleCellExperiment(human_seurobj)
sce <- decontX(sce, assayName = "counts")

human_seurobj <- as.Seurat(sce)
human_seurobj[["RNA"]] <- CreateAssayObject(counts = round(decontXcounts(sce)))

# Remove cells with >50% ambient contamination
cells_keep <- which(human_seurobj$decontX_contamination < 0.5 &
                    human_seurobj$nCount_RNA > 0)
human_seurobj <- subset(human_seurobj, cells = cells_keep)

cat("Cells after decontX filter:", ncol(human_seurobj), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 4: Doublet detection (Scrublet)
# ═══════════════════════════════════════════════════════════════════════════════
# Scrublet simulates synthetic doublets from observed cells and scores each
# real cell by its similarity to the simulated doublets.
use_condaenv("scrublet_env", required = TRUE)
scrublet <- import("scrublet")

counts_matrix <- human_seurobj@assays$RNA@counts
scrub <- scrublet$Scrublet(t(counts_matrix), expected_doublet_rate = 0.06)
doublet_scores <- scrub$scrub_doublets(
  min_counts             = 2,
  min_cells              = 3,
  min_gene_variability_pctl = 85,
  n_prin_comps           = as.integer(30)
)

human_seurobj$doublet_scores      <- doublet_scores[[1]]
human_seurobj$predicted_doublets  <- doublet_scores[[2]]

# Remove predicted doublets
human_seurobj <- subset(human_seurobj, predicted_doublets == FALSE)

cat("Cells after Scrublet filter:", ncol(human_seurobj), "\n")
cat("Doublet rate:", mean(doublet_scores[[2]]), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 5: QC metrics
# ═══════════════════════════════════════════════════════════════════════════════

genes_db <- ensembldb::select(EnsDb.Hsapiens.v86,
                              key     = rownames(human_seurobj),
                              columns = "SYMBOL",
                              keytype = "GENEID")

mt_genes   <- genes_db %>% dplyr::filter(str_detect(SYMBOL, "^MT-"))   %>% pull(GENEID)
ribo_genes <- genes_db %>% dplyr::filter(str_detect(SYMBOL, "^RP[LS]")) %>% pull(GENEID)

human_seurobj[["percent.mt"]]   <- PercentageFeatureSet(human_seurobj, features = mt_genes)
human_seurobj[["percent.ribo"]] <- PercentageFeatureSet(human_seurobj, features = ribo_genes)

cat("QC summary:\n")
print(summary(human_seurobj$percent.mt))
print(summary(human_seurobj$nFeature_RNA))

# ═══════════════════════════════════════════════════════════════════════════════
# Section 6: Cell type annotation via Seurat label transfer
# ═══════════════════════════════════════════════════════════════════════════════
# Cell types are assigned by projecting JM cells onto a reference PCA space
# from the integrated multi-study patch-seq object (pclamp_patched_all).
#
# ⚠ Dependency: pclamp_patched_all must exist before running this section.
# If reproducing from scratch, run scripts/02_integrate_patchseq_datasets.R
# first to generate this reference object.

load("../human_islets_scRNAseq/20241015_patchSeq_integrated_patched-only_decon.RData")

# Light QC filter on reference to remove extreme outliers before transfer
pclamp_patched_filtered <- subset(
  x = pclamp_patched_all,
  nFeature_RNA < 12500 & nFeature_RNA > 200 &
    percent.mt < summary(pclamp_patched_all$percent.mt)[5]
)
pclamp_patched_filtered <- JoinLayers(pclamp_patched_filtered)

# Find transfer anchors using PCA projection
anchors <- FindTransferAnchors(
  reference            = pclamp_patched_filtered,
  normalization.method = "LogNormalize",
  reduction            = "pcaproject",
  reference.assay      = "RNA",
  query                = human_seurobj,
  dims                 = 1:30,
  reference.reduction  = "pca",
  verbose              = TRUE
)

human_seurobj <- MapQuery(
  anchorset           = anchors,
  reference           = pclamp_patched_filtered,
  query               = human_seurobj,
  refdata             = list(celltype = "celltype"),
  reference.reduction = "ref.integrated.rpca",
  reduction.model     = "umap"
)

# ── Manual correction for low-confidence predictions ─────────────────────────
# For cells where prediction score < 0.7, override with marker gene expression.
# Thresholds were determined by visual inspection.
# Marker gene Ensembl IDs:
#   ENSG00000157005 = SST (delta)     ENSG00000091704 = CPA1 (acinar)
#   ENSG00000171345 = KRT19 (ductal)  ENSG00000157017 = GHRL (epsilon)
#   ENSG00000254647 = INS (beta)      ENSG00000115263 = GCG (alpha)
#   ENSG00000108849 = PPY (PP)

human_seurobj$celltype <- human_seurobj$predicted.celltype

correct_celltype <- function(obj, ensg, threshold, new_type, max_score = 0.7) {
  cells    <- WhichCells(obj, expression = !!sym(ensg) > threshold, slot = "data")
  cells    <- cells[cells %in% colnames(obj)]
  low_conf <- cells[obj$predicted.celltype.score[cells] < max_score]
  obj$celltype[low_conf] <- new_type
  obj
}

human_seurobj <- correct_celltype(human_seurobj, "ENSG00000157005", 7,   "delta")
human_seurobj <- correct_celltype(human_seurobj, "ENSG00000091704", 5.5, "acinar")
human_seurobj <- correct_celltype(human_seurobj, "ENSG00000171345", 4.5, "ductal")
human_seurobj <- correct_celltype(human_seurobj, "ENSG00000157017", 4,   "epsilon")
human_seurobj <- correct_celltype(human_seurobj, "ENSG00000254647", 7.5, "beta")
human_seurobj <- correct_celltype(human_seurobj, "ENSG00000115263", 7.5, "alpha")

# PPY: only one cell passed threshold after low-confidence filter
pp <- WhichCells(human_seurobj, expression = ENSG00000108849 > 7)
pp <- pp[pp %in% colnames(human_seurobj)]
pp <- pp[human_seurobj$predicted.celltype.score[pp] < 0.7]
if (length(pp) > 0) human_seurobj$celltype[pp[4]] <- "PP"

JM_seurobj_ensg <- subset(human_seurobj, User == "JM")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 7: Cell group harmonization
# ═══════════════════════════════════════════════════════════════════════════════
# Collapse the many differentiation day labels into broad stage bins, and
# collapse plate identifiers into run-level batches.
# These labels are used as covariates and grouping variables throughout.

# Remove cells with invalid Group labels and iPSCs
JM_seurobj_ensg <- subset(JM_seurobj_ensg, Group != "N/A" & !is.na(Group))
JM_seurobj_ensg$GroupBin <- case_when(
  JM_seurobj_ensg$Group == "human islets"             ~ "Human Primary",
  JM_seurobj_ensg$Group %in% c("SCB-30","SCB-35","SCB-D23","SCB-D24",
                                 "SCB-D26","SCB-D30","SCB-D31","SCB-D32",
                                 "SCB-D33","SCB-D28",
                                 "Stage 6 - Day 29 SCβ",
                                 "Stage 6 - Day 31 SCβ")  ~ "SCB D23–35",
  JM_seurobj_ensg$Group %in% c("SCB-D45","SCB-D47","SCB-D48","SCB-D49",
                                 "SCB-D51","SCB-D52",
                                 "Stage 6 - Day 54 SCβ")  ~ "SCB D45–54",
  JM_seurobj_ensg$Group == "IPS-Islet cells"           ~ "IPSC",
  TRUE ~ JM_seurobj_ensg$Group
)

JM_seurobj_ensg <- subset(JM_seurobj_ensg, GroupBin != "IPSC")

# Collapse plate labels to run-level batches for SCTransform regression
# (done in scripts/02_integrate_patchseq_datasets.R)
JM_seurobj_ensg$Run <- case_when(
  JM_seurobj_ensg$Plate == "Dec_2021_Plate2"              ~ "Dec2021",
  JM_seurobj_ensg$Plate == "Nov2021-Plate1"               ~ "Nov2021",
  JM_seurobj_ensg$Plate %in% c("2022Nov-Plate1",
                                 "2022Nov-Plate2")         ~ "Nov2022",
  JM_seurobj_ensg$Plate %in% c("Jan-2023-Plate4",
                                 "Feb-2023-Plate3")        ~ "Jan2023",
  JM_seurobj_ensg$Plate == "April-2023-Plate1"            ~ "Apr2023",
  JM_seurobj_ensg$Plate %in% c("Sep2023_Plate3",
                                 "Sep2023_Plate4")         ~ "Sep2023",
  TRUE ~ JM_seurobj_ensg$Plate
)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 8: Circadian gene module score
# ═══════════════════════════════════════════════════════════════════════════════
# Circadian gene expression varies with recording time and introduces
# plate-level batch effects. A module score is computed here and used as a
# regression covariate in SCTransform in scripts/02_integrate_patchseq_datasets.R

circadian_genes <- c("CLOCK","BMAL1","PER1","PER2","CRY1","CRY2","NR1D1","RORA")
circadian_ensg  <- ensembldb::select(
  EnsDb.Hsapiens.v86,
  key     = circadian_genes,
  columns = "GENEID",
  keytype = "SYMBOL"
)$GENEID

JM_seurobj_ensg <- NormalizeData(JM_seurobj_ensg)
JM_seurobj_ensg <- AddModuleScore(
  JM_seurobj_ensg,
  features = list(circadian_ensg),
  nbin     = 12,
  ctrl     = 100,
  name     = "circadian_features",
  seed     = 42
)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 9: Save output
# ═══════════════════════════════════════════════════════════════════════════════
ct_full <- human_seurobj$celltype
human_meta$celltype <- "Unknown"
human_meta$celltype[human_meta$NAME_FROM_MBSU %in% names(ct_full)] <-
  as.character(ct_full)

save(JM_seurobj_ensg,
     file = "20241017_JM_islets_seurat_object_ensg_decon.RData")
save(human_meta, file = "20241017_human_meta.RData")

cat("\nDone.\n")
cat("Cells in JM_seurobj_ensg:", ncol(JM_seurobj_ensg), "\n")
cat("Cell type distribution:\n")
print(table(JM_seurobj_ensg$celltype))
print(table(JM_seurobj_ensg$GroupBin))
