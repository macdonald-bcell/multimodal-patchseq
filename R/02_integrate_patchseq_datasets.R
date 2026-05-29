# ==============================================================================
# 02_integrate_patchseq_datasets.R
# ==============================================================================
# Merge the JM_patchSeq study with published primary islet patch-seq datasets,
# apply final QC filters, run SCTransform normalization, and export the
# combined object for downstream analysis and pySCENIC.
#
# Corresponds to Methods: "Integrating Patch-seq Datasets"
#   "Patch-seq scRNA-seq data from the JM_patchSeq study and five previously
#    published primary islet patch-seq datasets were integrated into a single
#    Seurat object... All downstream analyses used decontaminated counts as
#    the RNA assay input."
#
# Input:
#   20241017_JM_islets_seurat_object_ensg_decon.RData
#     - JM_seurobj_ensg: from scripts/01_data_ingestion_JM.R
#   pclamp_patched_all
#     - Integrated multi-study patch-seq object (external studies)
#     - Source: ../human_islets_scRNAseq/20241015_patchSeq_integrated_patched-only_decon.RData
#
# Output:
#   20241209_cells_combined.RData
#     - cells.combined: merged Seurat object, Ensembl ID rownames
#     - 2,865 cells × 47,081 genes after QC filtering
#   20241209_cells_combined.h5ad
#     - Same object in AnnData format (Ensembl IDs)
#   20250129_cells_combined_symbs_GFP.h5ad
#     - AnnData with HGNC symbol rownames (+ eGFP retained)
#     - Input to python/01_h5ad_to_loom.py for pySCENIC
#
# Dependencies:
#   tidyverse, Seurat, EnsDb.Hsapiens.v86, ensembldb, reticulate, sceasy
#
# Note on SCTransform:
#   SCTransform is NOT run in this script. It is run in the context of each
#   specific downstream analysis (MOFA2, DIABLO, correlations) because each
#   requires different subsets and data representations. The combined object
#   saved here contains raw decontaminated counts only.
# ==============================================================================

library("tidyverse")
library("Seurat")
library("EnsDb.Hsapiens.v86")
library("ensembldb")
library("reticulate")
library("sceasy")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 1: Load input objects
# ═══════════════════════════════════════════════════════════════════════════════

load("20241017_JM_islets_seurat_object_ensg_decon.RData")
# JM_seurobj_ensg loaded

load("../human_islets_scRNAseq/20241015_patchSeq_integrated_patched-only_decon.RData")
# pclamp_patched_all loaded

# ═══════════════════════════════════════════════════════════════════════════════
# Section 2: Prepare external patch-seq object for merging
# ═══════════════════════════════════════════════════════════════════════════════
# The external object contains cells from XQ, HPAP, and CryoPatchSeq studies.
# We strip SCT and scaled data to reduce memory before merging, keeping only
# raw RNA counts.

patched.merge <- pclamp_patched_all
DefaultAssay(patched.merge) <- "RNA"
patched.merge <- JoinLayers(patched.merge)
patched.merge[["SCT"]] <- NULL
patched.merge[["RNA"]]$data <- NULL

# Standardize metadata fields needed for downstream analysis
patched.merge$GroupBin <- "Human Primary"
patched.merge$Differentiation <- patched.merge$Donor
patched.merge$DiabetesStatus[is.na(patched.merge$DiabetesStatus)] <- "Non-Diabetic"

# Exclude diabetic donors and cells patched by Austin (QC exclusion)
# Diabetic donors are excluded to avoid confounding in the SCβ vs β comparison
patched.ND <- subset(patched.merge, DiabetesStatus == "Non-Diabetic")
patched.ND <- subset(patched.ND, Patcher != "Austin")

cat("External patch-seq cells after filtering:", ncol(patched.ND), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 3: Prepare JM object for merging
# ═══════════════════════════════════════════════════════════════════════════════

# Extract Differentiation batch labels from free-text ID field using regex
# This is needed to identify donor batches for the mixed model random effect
diff <- str_extract(
  JM_seurobj_ensg$ID.From.Fclynne,
  "(?i)(\\d*\\s*diff(\\s*\\d+(/\\d+)?|\\s*\\d{6}|\\s*start\\s*date:\\d{8})|differentiation\\s*\\d+|\\d+\\s*differentiation)"
)
JM_seurobj_ensg$Differentiation <- diff

# For primary islet cells with missing Differentiation, use donor ID as proxy
# (primary cells have one cell per donor so donor = differentiation batch)
donors <- JM_seurobj_ensg$ID.From.Fclynne[
  is.na(JM_seurobj_ensg$Donor) & JM_seurobj_ensg$Group == "Human Primary"
]
names(donors) <- colnames(JM_seurobj_ensg)[
  is.na(JM_seurobj_ensg$Donor) & JM_seurobj_ensg$Group == "Human Primary"
]
JM_seurobj_ensg$Differentiation[names(donors)] <- donors
JM_seurobj_ensg$Donor[names(donors)] <- donors

donors2 <- JM_seurobj_ensg$Donor[
  is.na(JM_seurobj_ensg$Differentiation) & JM_seurobj_ensg$Group == "Human Primary"
]
names(donors2) <- colnames(JM_seurobj_ensg)[
  is.na(JM_seurobj_ensg$Differentiation) & JM_seurobj_ensg$Group == "Human Primary"
]
JM_seurobj_ensg$Differentiation[names(donors2)] <- donors2

# Sanity check: any JM cells still missing Differentiation?
n_missing <- sum(is.na(JM_seurobj_ensg$Differentiation) &
                   JM_seurobj_ensg$User == "JM" &
                   JM_seurobj_ensg$GroupBin != "IPSC")
cat("JM cells with missing Differentiation:", n_missing, "\n")  # expect 0

# Remove all-NA metadata columns to clean up before merge
JM_meta_clean <- JM_seurobj_ensg@meta.data %>%
  dplyr::select(where(~ any(!is.na(.))))
JM_seurobj_ensg@meta.data <- JM_meta_clean

# Fix a label with a colon that causes downstream issues in some file formats
JM_seurobj_ensg$Differentiation[
  JM_seurobj_ensg$Differentiation == "35 Diff start date:20210905"
] <- "35 Diff start date 20210905"

JM_seurobj_ensg$Study <- "JM_patchSeq"

# Apply final QC filters to JM cells before merging
# nFeature_RNA >= 200: remove low-complexity cells
# percent.mt <= 80th percentile: remove high-mitochondrial cells
# iPSCs excluded: not relevant to the SCβ vs β comparison
JM_seurobj_ensg_sub <- subset(JM_seurobj_ensg, GroupBin != "IPSC")
JM_seurobj_ensg_sub <- subset(
  JM_seurobj_ensg_sub,
  nFeature_RNA >= 200 &
    percent.mt <= summary(JM_seurobj_ensg_sub$percent.mt)[5]
)

cat("JM cells after QC filter:", ncol(JM_seurobj_ensg_sub), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 4: Merge datasets
# ═══════════════════════════════════════════════════════════════════════════════

cells.combined <- merge(patched.ND, JM_seurobj_ensg_sub)

# Standardize preservation status across studies
cells.combined$Fresh.or.Cryo[cells.combined$Study == "JM_patchSeq"]   <- "Fresh"
cells.combined$Fresh.or.Cryo[cells.combined$Study == "HPAP_Patchseq"] <- "Fresh"
cells.combined$Fresh.or.Cryo[cells.combined$Study == "CryoPatchSeq"]  <-
  cells.combined$fresh.cryp.preserve[cells.combined$Study == "CryoPatchSeq"]
cells.combined$Fresh.or.Cryo[cells.combined$Fresh.or.Cryo == "cryo-preserved"] <- "Cryo"
cells.combined$Fresh.or.Cryo[cells.combined$Fresh.or.Cryo == "fresh"]          <- "Fresh"

cells.combined <- JoinLayers(cells.combined)

# Remove SCT and scaled data from the merged object (will be recomputed per analysis)
cells.combined[["SCT"]]             <- NULL
cells.combined[["RNA"]]$scale.data  <- NULL

# Remove unused metadata columns
cells.combined$Note     <- NULL
cells.combined$Staining <- NULL
cells.combined$Drugs    <- NULL

# Keep only endocrine cell types relevant to the comparison
# "Unknown" = SCβ cells that weren't in the primary islet label transfer
cells.combined <- subset(
  cells.combined,
  celltype %in% c("alpha", "beta", "delta", "PP", "Unknown")
)

cat("Cells after cell type filter:", ncol(cells.combined), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 5: Gene detection filter
# ═══════════════════════════════════════════════════════════════════════════════
# Remove genes detected in fewer than 1 cell equivalent per study-worth of cells
# Threshold: >3/ncells * 100 percent of cells (approximately 1 cell per study)
# This removes very rare transcripts that add noise without signal

counts               <- LayerData(cells.combined, layer = "counts", assay = "RNA")
n_cells              <- ncol(cells.combined)
genes.percent.expression <- rowMeans(counts > 0) * 100
genes.filter         <- names(genes.percent.expression[
  genes.percent.expression > 3 / n_cells * 100
])

cat("Genes before detection filter:", nrow(cells.combined), "\n")
cat("Genes after detection filter:", length(genes.filter), "\n")
rm(counts)

cells.combined <- subset(x = cells.combined, features = genes.filter)
cells.combined <- subset(x = cells.combined, nFeature_RNA > 100)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 6: Final QC metrics and filtering
# ═══════════════════════════════════════════════════════════════════════════════

genes_db   <- ensembldb::select(EnsDb.Hsapiens.v86,
                                key     = rownames(cells.combined),
                                columns = "SYMBOL",
                                keytype = "GENEID")
ribo_genes <- genes_db %>% dplyr::filter(str_detect(SYMBOL, "^RP[LS]")) %>% pull(GENEID)

cells.combined[["percent.ribo"]] <- PercentageFeatureSet(cells.combined,
                                                          features = ribo_genes)

# Final QC filter: remove high-ribo, high-mt, low-complexity, and very-high-count cells
# These thresholds were determined by inspection of QC distributions
highMT <- quantile(cells.combined$percent.mt, probs = 0.9)

cells.combined <- subset(
  cells.combined,
  nFeature_RNA < 15000 &
    percent.mt  <= highMT &
    percent.ribo > 0 &
    percent.ribo < 15
)

cat("Final cell count:", ncol(cells.combined), "\n")
cat("Final gene count:", nrow(cells.combined), "\n")
cat("Study breakdown:\n")
print(table(cells.combined$Study))
print(table(cells.combined$celltype))

# ═══════════════════════════════════════════════════════════════════════════════
# Section 7: Circadian module score
# ═══════════════════════════════════════════════════════════════════════════════
# Recomputed on the combined object so the score reflects the full cell
# composition. Used as a regression covariate in downstream SCTransform calls.

cells.combined <- NormalizeData(cells.combined)

circadian_genes <- c("CLOCK","BMAL1","PER1","PER2","CRY1","CRY2","NR1D1","RORA")
circadian_ensg  <- ensembldb::select(
  EnsDb.Hsapiens.v86,
  key     = circadian_genes,
  columns = "GENEID",
  keytype = "SYMBOL"
)$GENEID

cells.combined <- AddModuleScore(
  cells.combined,
  features = list(circadian = circadian_ensg),
  nbin     = 24,
  ctrl     = 100,
  name     = "circadian_score",
  seed     = 42,
  assay    = "RNA"
)

# Rename "Unknown" cell type to "SCβ" for clarity in all downstream analyses
cells.combined$celltype[cells.combined$celltype == "Unknown"] <- "SCβ"

# ═══════════════════════════════════════════════════════════════════════════════
# Section 8: Save combined object
# ═══════════════════════════════════════════════════════════════════════════════

save(cells.combined, file = "20241209_cells_combined.RData")

# Export to AnnData format (Ensembl IDs) for downstream Python tools
cells.combined.sce <- as.SingleCellExperiment(cells.combined)
use_condaenv("pyscenic2", required = TRUE)
sceasy::convertFormat(
  cells.combined.sce,
  from             = "sce",
  to               = "anndata",
  main_layer       = "counts",
  drop_single_values = FALSE,
  outFile          = "20241209_cells_combined.h5ad"
)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 9: Export with HGNC symbol rownames for pySCENIC
# ═══════════════════════════════════════════════════════════════════════════════
# pySCENIC requires HGNC gene symbols rather than Ensembl IDs.
# Convert rownames, sum duplicate symbols, filter low-detection genes,
# and export to h5ad. This file feeds python/01_h5ad_to_loom.py.

counts   <- cells.combined[["RNA"]]$counts
metadata <- cells.combined@meta.data

tr2g <- tibble(
  gene      = rownames(counts),
  gene_name = mapIds(EnsDb.Hsapiens.v86,
                     keys     = rownames(counts),
                     column   = "SYMBOL",
                     keytype  = "GENEID",
                     multiVals = "first")
)

# Explicitly retain eGFP (not in EnsDb)
tr2g$gene_name[tr2g$gene == "eGFP"] <- "eGFP"

# Remove genes without a symbol mapping
keeps  <- tr2g[!is.na(tr2g$gene_name), ]
counts <- counts[keeps$gene, ]
rownames(counts) <- keeps$gene_name

# Sum counts for duplicate symbols (e.g. different Ensembl IDs mapping to same symbol)
counts2 <- rowsum(counts, rownames(counts)) %>% as("dgCMatrix")

# Final detection filter: keep genes detected in >0.1% of cells
genes.percent.expression <- rowMeans(counts2 > 0) * 100
genes.filter             <- names(genes.percent.expression[genes.percent.expression > 0.1])
counts2                  <- counts2[genes.filter, ]

cat("Genes in pySCENIC input (HGNC symbols):", nrow(counts2), "\n")

cells.combined.sce.symbs <- SingleCellExperiment::SingleCellExperiment(
  assays  = list(counts = counts2),
  colData = metadata
)

sceasy::convertFormat(
  cells.combined.sce.symbs,
  from             = "sce",
  to               = "anndata",
  main_layer       = "counts",
  drop_single_values = FALSE,
  outFile          = "20250129_cells_combined_symbs_GFP.h5ad"
)

cat("\nDone. Outputs:\n")
cat("  20241209_cells_combined.RData   - main Seurat object\n")
cat("  20241209_cells_combined.h5ad    - AnnData (Ensembl IDs)\n")
cat("  20250129_cells_combined_symbs_GFP.h5ad - AnnData (HGNC symbols, for pySCENIC)\n")
