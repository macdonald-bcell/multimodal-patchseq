# ==============================================================================
# 07_DIABLO.R
# ==============================================================================
# Supervised multi-omics integration using DIABLO (mixOmics block.splsda),
# run across 50 imputations to identify consensus cross-block features.
# LME4 mixed models annotate DIABLO-selected features with differential
# activity effect sizes for the cross-block correlation heatmap.
#
# Corresponds to Methods: "Supervised multi-omics integration (DIABLO)" and
# "Linear mixed models for DIABLO feature annotation"
#   "...block.splsda (mixOmics) was fit on pseudobulk-aggregated mRNA, ephys,
#    and SCENIC AUC data... across 50 imputation datasets... Features present
#    in ≥80% of models with ≥90% sign consistency on component 1 were retained
#    as consensus features... Cross-block DIABLO latent variable correlations
#    on component 1 were extracted using circosPlot(cutoff=1)... averaged
#    across valid imputation models..."
#
# Input:
#   20260113_cells_plus_filtered_mi_SCENIC.rds
#     - cells.plus: integrated Seurat object with SCT, AUC, and ephys assays
#   20260120_imp_bundle.rds
#     - imp_bundle$completed_list: 50 imputed electrophysiology datasets
#   bootstrapped_DE_triage_table.csv
#     - DE triage table from scripts/04_differential_expression.R
#   250404_patch_key.xlsx
#     - Electrophysiology parameter pretty-name lookup
#   diablo_50imps/models/diablo_imp*.rds (generated within this script)
#
# Output:
#   20260123_diablo_pseudo_celltype_mi.RData
#     - diablo.pseudo, tune.diablo.pseudo: single-run DIABLO fit + tuning
#   diablo_50imps/diablo_run_summary.csv
#     - Per-imputation diagnostics (correlation, classification accuracy)
#   diablo_50imps/feature_stability_comp1to4.csv
#     - Feature frequency and sign stability across 50 imputation models
#   diablo_50imps/stable_features_top_comp1to4.csv
#     - Consensus features (freq ≥ 0.8, sign consistency ≥ 0.9)
#   SCENIC_AUC_lmer_beta_vs_SCbeta_JM.csv
#     - LME4 results for all SCENIC regulons (SCβ vs β, JM_patchSeq only)
#   20260327_DIABLO_ephys_crossblock_heatmap_Unchanged.svg
#     - Final heatmap: ephys × cross-block DIABLO correlations with LME4 annotations
#   for_manuscript/20260123_mixo_supplemental.svg
#     - plotArrow + plotIndiv diagnostic panels
#
# Dependencies:
#   tidyverse, Seurat, mixOmics, caret, lme4, lmerTest, broom.mixed,
#   ComplexHeatmap, circlize, scico, EnsDb.Hsapiens.v86, ensembldb,
#   janitor, readxl, patchwork
#
# Runtime note:
#   Running block.splsda across 50 imputations is computationally intensive.
#   Each model fit takes ~2-5 minutes; full run is several hours.
#   Models are saved to diablo_50imps/models/ as they complete.
# ==============================================================================

suppressPackageStartupMessages({
  library("tidyverse")
  library("Seurat")
  library("mixOmics")
  library("caret")
  library("lme4")
  library("lmerTest")
  library("broom.mixed")
  library("ComplexHeatmap")
  library("circlize")
  library("scico")
  library("EnsDb.Hsapiens.v86")
  library("ensembldb")
  library("janitor")
  library("readxl")
  library("patchwork")
})

# ── Parameter definitions ─────────────────────────────────────────────────────
ephys_vars <- c(
  "CellSize_pF",
  "NormalizedTotalCapacitance_fF.pF",
  "NormalizedFirstDepolarizationCapacitance_fF.pF",
  "NormalizedLateDepolarizationCapacitance",
  "CalciumIntegralNormalizedtoCellSize_pC.pF",
  "NormalizedEarlyPeakCalciumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedLateCalciumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedLVA_.60mV_pA.pF",
  "NormalizedHVA_.20mV_pA.pF",
  "NormalizedPeakSodiumCurrentAmplitudeat.10mV_pA.pF",
  "HalfInactivationofSodiumCurrent_mV",
  "VoltageforSodiumPeakCurrent_mV",
  "ReversalPotentialbyramp_mV"
)

# Current amplitude features: more negative = more active
# Signs are flipped so positive = more active in the heatmap
flip_features <- c(
  "Normalized LVA Ca²⁺ Current  (pA/pF)",
  "Normalized Late Ca²⁺ Current  (pA/pF)",
  "Normalized Early  Ca²⁺ Current  (pA/pF)",
  "Normalized HVA Ca²⁺ Current (pA/pF)",
  "Normalized Na²⁺ Current (pA/pF)",
  "Normalized Ca²⁺ Influxed  (pC/pF)"
) #\u00b2\u207a

# Unannotated transcript pattern for exclusion
exclude_pattern <- "^RP11-|^RP[0-9]+-|^AC[0-9]+\\.|^AL[0-9]+\\.|^ENSG[0-9]+"

# Colour palette
cluster_cols2 <- c(
  "#ff7db1","#e2ff50","#713be8","#669b00","#c734e9","#01a457","#e400b3",
  "#8cffb0","#010f92","#ffd26e","#015fde","#005a02","#ff6de1","#01deaf",
  "#9d0098","#4f6600","#a282ff","#867100","#170047","#b06100","#01b4e8",
  "#ff5662","#02cfcb","#ff3a85","#006c55","#ac004a","#baecff","#862200",
  "#006d95","#ff7d56","#003767","#ffcaa5","#33002d","#ffd6db","#101c00",
  "#f6bbff","#001e22","#ff9da8","#360c00","#5f3800","#001024","#d0f700"
)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 1: Load data and build pseudobulk matrices
# ═══════════════════════════════════════════════════════════════════════════════

cells.plus <- readRDS("20260113_cells_plus_filtered_mi_SCENIC.rds")

imp_bundle     <- readRDS("20260120_imp_bundle.rds")
completed_list <- imp_bundle$completed_list

triage2     <- read_csv("bootstrapped_DE_triage_table.csv")
ephys_names <- read_excel("250404_patch_key.xlsx", .name_repair = "universal") %>%
  dplyr::select(Pretty.names, name_in_xq_sheets)

# Exclude S6D1 from DIABLO (no ephys data)
cells.combined.patched <- subset(cells.plus, celltype != "S6D1" &
                                   !is.na(CellSize_pF) & CellSize_pF != 0)

# ── SCTransform for mRNA pseudobulk ──────────────────────────────────────────
# Note: SCTransform is run here on the patched-cell subset for DIABLO only.
# This is separate from any SCTransform in scripts/01-02.
DefaultAssay(cells.combined.patched) <- "RNA"
cells.combined.patched <- SCTransform(
  cells.combined.patched,
  vars.to.regress       = c("percent.mt", "nFeature_RNA", "percent.ribo"),
  variable.features.n   = 500, #low n for minimal batch correction
  residual.features     = NULL,
  return.only.var.genes = FALSE
)

# ── Pseudobulk aggregation ───────────────────────────────────────────────────
# AggregateExpression produces one pseudobulk sample per Differentiation × celltype
auc_pseudobulk <- AggregateExpression(
  object = cells.combined.patched,
  assays = "AUC",
  normalization.method = "none",
  group.by = c("Differentiation", "celltype")
)

rna_pseudobulk <- AggregateExpression(
  object = cells.combined.patched,
  assays = "SCT",
  slot   = "data",
  normalization.method = "none",
  group.by = c("Differentiation", "celltype")
)

group_names <- colnames(rna_pseudobulk$SCT)
group_sizes <- cells.combined.patched@meta.data %>%
  group_by(Differentiation, celltype) %>%
  summarize(cell_count = n(), .groups = "drop") %>%
  unite("group", Differentiation, celltype, sep = "_", remove = FALSE) %>%
  mutate(group = if_else(str_detect(group, "^[0-9]"), paste0("g", group), group))
group_sizes_vec <- setNames(group_sizes$cell_count, group_sizes$group)
group_sizes_vec <- group_sizes_vec[group_names]

# Normalize pseudobulk by cell count, then scale
rna_pseudobulk$SCT <- sweep(rna_pseudobulk$SCT, 2, group_sizes_vec, "/")
mrna_mtx           <- t(rna_pseudobulk$SCT) %>% scale()

auc_pseudobulk$AUC <- sweep(auc_pseudobulk$AUC, 2, group_sizes_vec, "/")
auc_mtx            <- t(auc_pseudobulk$AUC) %>% scale()

# Clean SCENIC regulon name formatting
colnames(auc_mtx) <- str_replace(colnames(auc_mtx), "\\-\\(\\+\\)", "\\(\\+\\)")

# ── Gene symbol mapping and filtering ────────────────────────────────────────
ensg <- colnames(mrna_mtx)
symbs <- mapIds(EnsDb.Hsapiens.v86, keys = ensg, column = "SYMBOL",
                keytype = "GENEID", multiVals = "first")
symbs[is.na(symbs)]        <- names(symbs[is.na(symbs)])
symbs[duplicated(symbs)]   <- names(symbs[duplicated(symbs)])
colnames(mrna_mtx) <- as.character(symbs)

# Remove ribosomal genes and unannotated transcripts
mrna_mtx <- mrna_mtx[, !str_detect(colnames(mrna_mtx), "^RP[LS]")]

# Remove near-zero variance features
nzv      <- caret::nearZeroVar(mrna_mtx, freqCut = 95/5, uniqueCut = 10)
mrna_mtx <- mrna_mtx[, -nzv]
cat("mRNA features after NZV filter:", ncol(mrna_mtx), "\n")

# Load LOSO stability results to exclude unstable/weak genes
loso_summary <- readRDS("loso_summary.rds")
unstable_cut <- quantile(loso_summary$loso_rel_iqr, 0.90, na.rm = TRUE)
weak_cut     <- quantile(abs(loso_summary$signed_f_full), 0.50, na.rm = TRUE)
genes_drop   <- loso_summary %>%
  dplyr::filter(loso_rel_iqr >= unstable_cut,
                abs(signed_f_full) <= weak_cut) %>%
  pull(gene_name) %>% unique()
mrna_mtx <- mrna_mtx[, !colnames(mrna_mtx) %in% genes_drop, drop = FALSE]
cat("mRNA features after LOSO filter:", ncol(mrna_mtx), "\n")

# ── Ephys pseudobulk helper ───────────────────────────────────────────────────
# Builds scaled ephys pseudobulk matrix for a given imputation dataset
.build_ephys_pseudobulk_scaled <- function(k, cells_obj, ephys_by_imp,
                                           group_sizes_vec) {
  ephys_mat    <- ephys_by_imp[[k]]
  patch_assay  <- ephys_mat %>% as.data.frame() %>% t() %>%
    CreateAssayObject()
  cells_obj[["ephys"]] <- patch_assay
  ephys_pseudobulk <- AggregateExpression(
    object = cells_obj,
    assays = "ephys",
    normalization.method = "none",
    group.by = c("Differentiation", "celltype")
  )
  ephys_pseudobulk$ephys <- sweep(ephys_pseudobulk$ephys, 2,
                                   group_sizes_vec, "/")
  t(ephys_pseudobulk$ephys) %>% scale()
}

# Build ephys matrices for all 50 imputations
ephys_by_imp <- lapply(seq_along(completed_list), function(k) {
  d <- completed_list[[k]]
  rownames(d) <- d$CellID
  d[cells.combined.patched$CellID, ephys_vars, drop = FALSE]
})

ephys_pseudobulk_mtx <- .build_ephys_pseudobulk_scaled(
  1, cells.combined.patched, ephys_by_imp, group_sizes_vec
)

# ── Design matrix and outcome ─────────────────────────────────────────────────
samp   <- rownames(mrna_mtx)
Xp     <- list(mRNA = mrna_mtx, ephys = ephys_pseudobulk_mtx, scenic = auc_mtx)
Yp     <- group_sizes$celltype

design <- matrix(1, ncol = length(Xp), nrow = length(Xp),
                 dimnames = list(names(Xp), names(Xp)))
diag(design) <- 0

stopifnot(identical(rownames(Xp$mRNA), samp))
stopifnot(identical(rownames(Xp$ephys), samp))
stopifnot(identical(rownames(Xp$scenic), samp))
stopifnot(length(Yp) == length(samp))

# ═══════════════════════════════════════════════════════════════════════════════
# Section 2: Tune and fit single-run DIABLO model
# ═══════════════════════════════════════════════════════════════════════════════
# This single-run model is used to determine ncomp and keepX settings,
# which are then applied consistently across all 50 imputation models.

# Performance evaluation to select optimal ncomp
diablo.pseudo     <- block.plsda(Xp, Yp, ncomp = 5, design = design, scale = FALSE)
perf.diablo.pseudo <- perf(diablo.pseudo, validation = "Mfold",
                           folds = 10, nrepeat = 25)
ncomp_pseudo <- perf.diablo.pseudo$choice.ncomp$WeightedVote[
  "Overall.BER", "centroids.dist"
]
ncomp_pseudo <- pmax(2, ncomp_pseudo)
cat("Selected ncomp:", ncomp_pseudo, "\n")

# Tune feature selection (keepX) per block and component
tune.diablo.pseudo <- tune.block.splsda(
  Xp, Yp, ncomp = ncomp_pseudo,
  test.keepX = list(
    mRNA   = c(5, 10, 15, 20, 25, 30, 40, 50),
    ephys  = c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13),
    scenic = c(5, 10, 15, 20, 25, 30, 40, 50)
  ),
  design = design, scale = FALSE,
  validation = "Mfold", folds = 10, nrepeat = 10
)

list.keepX.pseudo <- tune.diablo.pseudo$choice.keepX
cat("keepX selection:\n"); print(list.keepX.pseudo)

diablo.pseudo <- block.splsda(
  Xp, Yp, ncomp = ncomp_pseudo,
  keepX  = list.keepX.pseudo,
  design = design, scale = FALSE
)

save(diablo.pseudo, tune.diablo.pseudo,
     file = "20260123_diablo_pseudo_celltype_mi.RData")

# ── Diagnostic plots (representative single model) ───────────────────────────
pop_colours   <- cluster_cols2[1:5]
block_colours <- scico::scico(3, palette = "buda", end = 0.8)

p_arrow  <- plotArrow(diablo.pseudo, ind.names = FALSE, legend = TRUE,
                      comp = c(1, 2), col = pop_colours,
                      title = "DIABLO comp 1 - 2")
p_indiv  <- plotIndiv(diablo.pseudo, ind.names = FALSE, legend = TRUE,
                      col = pop_colours, comp = c(1, 2),
                      title = "DIABLO comp 1 - 2")

ggsave("./for_manuscript/20260123_mixo_supplemental.svg",
       p_arrow + p_indiv$graph, width = 15, height = 7.5)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 3: Run DIABLO across 50 imputation datasets
# ═══════════════════════════════════════════════════════════════════════════════
# Uses ncomp and keepX from the single-run tuning above.
# Each model is saved immediately in case of interruption.

models_dir <- "diablo_50imps/models"
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

ncomp <- ncomp_pseudo

# Per-model correlation function (between block variates)
.block_cor <- function(model, comp = 1) {
  blocks <- names(model$variates)
  pairs  <- combn(blocks, 2, simplify = FALSE)
  sapply(pairs, function(p) {
    A <- model$variates[[p[1]]][, comp]
    B <- model$variates[[p[2]]][, comp]
    cor(A, B)
  })
}

diag_rows <- vector("list", length(completed_list))

for (k in seq_along(completed_list)) {
  cat("Imputation", k, "/", length(completed_list), "\n")
  out_file <- file.path(models_dir, sprintf("diablo_imp%02d.rds", k))
  if (file.exists(out_file)) {
    cat("  Already exists, skipping\n"); next
  }

  ephys_k <- .build_ephys_pseudobulk_scaled(
    k, cells.combined.patched, ephys_by_imp, group_sizes_vec
  )
  Xp_k <- list(mRNA = mrna_mtx, ephys = ephys_k, scenic = auc_mtx)

  fit <- tryCatch(
    mixOmics::block.splsda(
      X = Xp_k, Y = Yp,
      ncomp  = ncomp,
      keepX  = list.keepX.pseudo,
      design = design,
      scale  = FALSE
    ),
    error = function(e) { message("  Failed: ", e$message); NULL }
  )

  if (!is.null(fit)) {
    saveRDS(fit, file = out_file)

    # Extract cross-block correlation matrix using circosPlot
    corMat_P <- tryCatch(
      mixOmics::circosPlot(fit, cutoff = 1.0, comp = 1, plot = FALSE),
      error = function(e) NULL
    )

    diag_rows[[k]] <- tibble(
      imp    = k,
      status = "ok",
      block_cors = list(.block_cor(fit, comp = 1)),
      has_cormat = !is.null(corMat_P)
    )
  } else {
    diag_rows[[k]] <- tibble(imp = k, status = "failed",
                              block_cors = list(NULL), has_cormat = FALSE)
  }
}

diablo_run_summary <- bind_rows(diag_rows)
write_csv(diablo_run_summary, file.path("diablo_50imps", "diablo_run_summary.csv"))
cat("Models complete:", sum(diablo_run_summary$status == "ok"), "/ 50\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 4: Extract consensus features across 50 models
# ═══════════════════════════════════════════════════════════════════════════════
# Consensus definition: feature selected in ≥80% of models AND sign of loading
# consistent in ≥90% of models in which it was selected.

model_files   <- list.files(models_dir, pattern = "^diablo_imp\\d+\\.rds$",
                             full.names = TRUE)
diablo_models <- setNames(lapply(model_files, readRDS),
                          nm = gsub("^diablo_imp|\\.rds$", "",
                                    basename(model_files)))
n_runs <- length(diablo_models)
cat("Models loaded:", n_runs, "\n")

# Extract selected features and loadings per run
get_selected_tbl <- function(model, run_id, ncomp = 4) {
  purrr::map(seq_len(ncomp), function(comp) {
    purrr::map(names(model$loadings), function(b) {
      tibble(run     = run_id, comp = comp, block = b,
             feature = mixOmics::selectVar(model, block = b, comp = comp)[[b]]$name)
    }) %>% bind_rows()
  }) %>% bind_rows()
}

get_loadings_tbl <- function(model, run_id, ncomp = 4) {
  purrr::map(seq_len(ncomp), function(comp) {
    purrr::map(names(model$loadings), function(b) {
      v <- model$loadings[[b]][, comp]
      tibble(run = run_id, comp = comp, block = b,
             feature = names(v), loading = as.numeric(v))
    }) %>% bind_rows()
  }) %>% bind_rows()
}

selected_tbl <- imap(diablo_models, ~ get_selected_tbl(.x, .y, ncomp)) %>%
  bind_rows()
loadings_all <- imap(diablo_models, ~ get_loadings_tbl(.x, .y, ncomp)) %>%
  bind_rows()

# Restrict to selected features only
loadings_sel <- loadings_all %>%
  inner_join(selected_tbl, by = c("run", "comp", "block", "feature"))

# ── Sign alignment across runs ─────────────────────────────────────────────────
# DIABLO loadings may flip sign between runs (arbitrary global sign per component).
# Align all runs to a reference run using correlation of loadings.
align_sign_per_block_comp <- function(df, ref_run) {
  ref <- df %>%
    dplyr::filter(run == ref_run) %>%
    dplyr::select(comp, block, feature, loading_ref = loading)
  flips <- df %>%
    left_join(ref, by = c("comp", "block", "feature")) %>%
    group_by(run, comp, block) %>%
    summarise(
      flip = {
        x  <- loading; y <- loading_ref
        ok <- is.finite(x) & is.finite(y) & !is.na(y)
        if (sum(ok) < 5) FALSE else cor(x[ok], y[ok]) < 0
      },
      .groups = "drop"
    )
  df %>%
    left_join(flips, by = c("run", "comp", "block")) %>%
    mutate(loading = if_else(flip, -loading, loading)) %>%
    dplyr::select(run, comp, block, feature, loading)
}

ref_run          <- names(diablo_models)[1]
loadings_aligned <- align_sign_per_block_comp(loadings_sel, ref_run)

# ── Feature stability metrics ─────────────────────────────────────────────────
feature_stability <- loadings_aligned %>%
  mutate(sign = sign(loading)) %>%
  group_by(comp, block, feature) %>%
  summarise(
    n_models         = n(),
    freq             = n_models / n_runs,
    n_pos            = sum(sign > 0),
    n_neg            = sum(sign < 0),
    frac_consistent  = pmax(n_pos, n_neg) / n_models,
    mean_weight      = mean(loading,     na.rm = TRUE),
    mean_abs_weight  = mean(abs(loading), na.rm = TRUE),
    sd_weight        = sd(loading,       na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(comp, block, desc(freq), desc(frac_consistent), desc(mean_abs_weight))

write_csv(feature_stability,
          "diablo_50imps/feature_stability_comp1to4.csv")

# Consensus features: ≥80% frequency, ≥90% sign consistency
stable_features <- feature_stability %>%
  dplyr::filter(freq >= 0.8, frac_consistent >= 0.9)

stable_features_top <- stable_features %>%
  group_by(comp, block) %>%
  slice_max(order_by = mean_abs_weight, n = 200, with_ties = FALSE) %>%
  ungroup()

write_csv(stable_features_top,
          "diablo_50imps/stable_features_top_comp1to4.csv")

cat("Consensus features by block (comp 1):\n")
print(stable_features %>% dplyr::filter(comp == 1) %>% count(block))

# ═══════════════════════════════════════════════════════════════════════════════
# Section 5: Build consensus cross-block correlation matrix
# ═══════════════════════════════════════════════════════════════════════════════
# Extract circosPlot correlation matrix from each model, keep features present
# in ≥20% of models, average across models.

get_corMat <- function(model) {
  circosPlot(model, comp = c(1, 2), cutoff = 1, line = FALSE, plot = FALSE)
}

cor_mats  <- purrr::map(diablo_models, safely(get_corMat))
valid_cors <- purrr::keep(cor_mats, ~ is.null(.x$error)) %>%
  purrr::map("result")
cat("Valid correlation matrices:", length(valid_cors), "/ ", n_runs, "\n")

# Features present in ≥20% of valid models
feature_counts <- purrr::map(valid_cors, rownames) %>%
  unlist() %>% table() %>% as.data.frame() %>%
  janitor::clean_names() %>%
  dplyr::rename(feature = x)

threshold_feat <- 0.2 * length(cor_mats)
common_features <- feature_counts %>%
  dplyr::filter(freq >= threshold_feat,
                !str_detect(feature, exclude_pattern)) %>%
  pull(feature) %>% as.character()

# Align and average correlation matrices
aligned_cors <- purrr::map(valid_cors, function(mat) {
  common <- intersect(rownames(mat), common_features)
  mat[common, common, drop = FALSE]
})

long_cors <- purrr::map2(
  aligned_cors, seq_along(aligned_cors),
  ~ as_tibble(as.table(.x), .name_repair = "minimal") %>%
    set_names(c("from", "to", "correlation")) %>%
    dplyr::filter(from != to) %>%
    mutate(model = .y)
) %>% bind_rows()

mean_cors <- long_cors %>%
  group_by(from, to) %>%
  summarise(n_models = n(), mean_cor = mean(correlation, na.rm = TRUE),
            .groups = "drop")

# Build consensus correlation matrix
feature_set   <- sort(unique(c(mean_cors$from, mean_cors$to)))
consensus_mat <- matrix(NA_real_, nrow = length(feature_set),
                        ncol = length(feature_set),
                        dimnames = list(feature_set, feature_set))
for (i in seq_len(nrow(mean_cors))) {
  f1 <- mean_cors$from[i]; f2 <- mean_cors$to[i]
  consensus_mat[f1, f2] <- mean_cors$mean_cor[i]
  consensus_mat[f2, f1] <- mean_cors$mean_cor[i]
}

# ── Identify cross-block edges and build heatmap matrix ───────────────────────
cor_long <- mean_cors %>%
  mutate(
    type_from = case_when(
      str_detect(from, "\\(\\+\\)$") ~ "SCENIC",
      from %in% ephys_vars           ~ "ephys",
      TRUE                           ~ "mRNA"
    ),
    type_to = case_when(
      str_detect(to, "\\(\\+\\)$") ~ "SCENIC",
      to %in% ephys_vars           ~ "ephys",
      TRUE                         ~ "mRNA"
    ),
    edge_type = ifelse(type_from == type_to, "within", "cross")
  )

# Keep cross-block edges with |r| > 0.80, at least one endpoint is ephys
active <- cor_long %>%
  dplyr::filter(
    edge_type == "cross",
    abs(mean_cor) > 0.80,
    !str_detect(from, exclude_pattern),
    !str_detect(to,   exclude_pattern),
    from %in% ephys_vars | to %in% ephys_vars
  ) %>%
  { unique(c(.$from, .$to)) }

cat("Active features in heatmap:\n")
print(table(case_when(
  str_detect(active, "\\(\\+\\)$") ~ "SCENIC",
  active %in% ephys_vars           ~ "ephys",
  TRUE                             ~ "mRNA"
)))

# Build heatmap matrix: rows = ephys, columns = non-ephys features
ephys_active <- intersect(ephys_vars, active)
other_active <- setdiff(active, ephys_vars)

mat_plot <- consensus_mat[ephys_active, other_active]

# Pretty names for ephys rows
rownames(mat_plot) <- str_replace_all(rownames(mat_plot), "-", "_")
rownames(mat_plot) <- ephys_names$Pretty.names[
  match(rownames(mat_plot), ephys_names$name_in_xq_sheets)
]

# Clean SCENIC column names
colnames(mat_plot) <- str_replace(colnames(mat_plot), "-\\(\\+\\)", "(+)")

# Flip current amplitude rows (more negative = more active → flip for visualization)
to_flip          <- intersect(flip_features, rownames(mat_plot))
mat_plot[to_flip, ] <- -mat_plot[to_flip, ]
cat("Flipped rows:", to_flip, "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 6: LME4 mixed models for DIABLO feature annotation
# ═══════════════════════════════════════════════════════════════════════════════
# SCβ vs β-cell effect sizes for SCENIC AUC, mRNA, and ephys features.
# JM_patchSeq only; Differentiation as random effect; singular → fallback to lm.
# BH adjustment run across ALL ~300 SCENIC regulons (not just DIABLO features)
# to ensure correct multiple testing correction.
# Standardised β = estimate / SD(outcome) ≈ Cohen's d; comparable across modalities.

# ── Shared LME4 helper ───────────────────────────────────────────────────────
# Note: error message references 'feature_name' generically;
# the caller passes the feature name as the first argument.
run_lmer_feature <- function(feature_name, y, dat_with_meta) {
  y_sd         <- sd(y, na.rm = TRUE)
  if (is.na(y_sd) || y_sd < 1e-10) {
    return(tibble(feature = feature_name, estimate = NA_real_,
                  std_estimate = NA_real_, se = NA_real_,
                  t_stat = NA_real_, df = NA_real_, p.value = NA_real_))
  }

  donor_counts <- dat_with_meta %>% dplyr::count(Differentiation)
  use_lmer     <- nrow(donor_counts) >= 3 && any(donor_counts$n > 1)

  fit <- tryCatch({
    if (use_lmer)
      suppressMessages(suppressWarnings(
        lmer(y ~ group + (1 | Differentiation), data = dat_with_meta,
             control = lmerControl(optimizer = "bobyqa"))
      ))
    else
      lm(y ~ group, data = dat_with_meta)
  }, error = function(e) {
    message("Model failed: ", feature_name, " — ", e$message); NULL
  })

  if (is.null(fit))
    return(tibble(feature = feature_name, estimate = NA_real_,
                  std_estimate = NA_real_, se = NA_real_,
                  t_stat = NA_real_, df = NA_real_, p.value = NA_real_))

  # Singular → fallback to lm
  if (use_lmer && isSingular(fit)) {
    fit <- tryCatch(lm(y ~ group, data = dat_with_meta),
                    error = function(e) NULL)
    if (is.null(fit))
      return(tibble(feature = feature_name, estimate = NA_real_,
                    std_estimate = NA_real_, se = NA_real_,
                    t_stat = NA_real_, df = NA_real_, p.value = NA_real_))
  }

  res <- tryCatch({
    tidy(fit, effects = "fixed") %>%
      dplyr::filter(term == "groupSCβ") %>%
      dplyr::select(estimate, std.error, statistic,
                    df = any_of("df"), p.value) %>%
      mutate(feature = feature_name, std_estimate = estimate / y_sd)
  }, error = function(e) { message("tidy() failed: ", feature_name); NULL })

  if (is.null(res) || nrow(res) == 0)
    return(tibble(feature = feature_name, estimate = NA_real_,
                  std_estimate = NA_real_, se = NA_real_,
                  t_stat = NA_real_, df = NA_real_, p.value = NA_real_))

  if (!"df" %in% names(res)) res$df <- nrow(dat_with_meta) - 2
  res
}

# ── SCENIC AUC LME4 ──────────────────────────────────────────────────────────
# Run on ALL regulons for correct BH adjustment, not just DIABLO-selected ones.
auc_matrix   <- cells.plus[["AUC"]]$counts %>% as.matrix()
auc_metadata <- cells.plus@meta.data %>%
  dplyr::select(CellID, celltype, Differentiation, Study)
cells.plus$CellID <- rownames(cells.plus@meta.data)

auc_df <- auc_matrix %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("CellID") %>%
  left_join(auc_metadata, by = "CellID") %>%
  dplyr::filter(celltype %in% c("beta", "SCβ"), Study == "JM_patchSeq") %>%
  mutate(group = factor(celltype, levels = c("beta", "SCβ")))

cat("SCENIC LME4: beta =", sum(auc_df$celltype == "beta"),
    "SCβ =", sum(auc_df$celltype == "SCβ"), "\n")

# All regulons (for correct BH adjustment across full multiple testing space)
all_regulons <- setdiff(colnames(auc_df),
                        c("CellID", "celltype", "Differentiation",
                          "Study", "group"))

scenic_lmer_results_all <- lapply(all_regulons, function(reg) {
  run_lmer_feature(reg, auc_df[[reg]], auc_df)
}) %>%
  bind_rows() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

write_csv(scenic_lmer_results_all,
          "./for_manuscript/SCENIC_AUC_lmer_beta_vs_SCbeta_JM.csv")
cat("Significant SCENIC regulons (padj < 0.05):",
    sum(scenic_lmer_results_all$p_adj < 0.05, na.rm = TRUE), "\n")

# ── mRNA LME4 ────────────────────────────────────────────────────────────────
# Only for the small number of mRNA features selected by DIABLO
mrna_feats <- cor_long %>%
  dplyr::filter(type_from == "mRNA" | type_to == "mRNA") %>%
  { unique(c(.$from[.$type_from == "mRNA"], .$to[.$type_to == "mRNA"])) }

expr_mat_sct <- cells.plus[["SCT"]]$scale.data
tr2g <- tibble(
  gene      = rownames(expr_mat_sct),
  gene_name = mapIds(EnsDb.Hsapiens.v86, keys = rownames(expr_mat_sct),
                     column = "SYMBOL", keytype = "GENEID", multiVals = "first")
)
tr2g$gene_name[is.na(tr2g$gene_name)] <- tr2g$gene[is.na(tr2g$gene_name)]
rownames(expr_mat_sct) <- tr2g$gene_name
expr_mat_sct <- expr_mat_sct[rownames(expr_mat_sct) %in% mrna_feats, ]

rna_metadata <- cells.plus@meta.data %>%
  dplyr::select(CellID, celltype, Differentiation, Study)

mrna_df <- expr_mat_sct %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("CellID") %>%
  left_join(rna_metadata, by = "CellID") %>%
  dplyr::filter(celltype %in% c("beta", "SCβ"), Study == "JM_patchSeq") %>%
  mutate(group = factor(celltype, levels = c("beta", "SCβ")))

mrna_lmer_results <- lapply(mrna_feats, function(gene) {
  run_lmer_feature(gene, mrna_df[[gene]], mrna_df)
}) %>%
  bind_rows() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH"))

# ── Ephys LME4 (JM_patchSeq within-study, pooled across imputations) ─────────
run_lmer_within_study <- function(var, imp_list, study = "JM_patchSeq") {
  results <- lapply(seq_along(imp_list), function(k) {
    dat <- imp_list[[k]] %>%
      dplyr::filter(celltype %in% c("beta", "SCβ"), Study == study) %>%
      mutate(group = factor(celltype, levels = c("beta", "SCβ")),
             y     = .data[[var]]) %>%
      dplyr::filter(!is.na(y)) %>%
      mutate(y = as.numeric(y))
    if (nrow(dat) < 10) return(NULL)
    y_sd         <- sd(dat$y, na.rm = TRUE)
    donor_counts <- dat %>% dplyr::count(Donor)
    fit <- tryCatch({
      if (nrow(donor_counts) >= 3 && any(donor_counts$n > 1))
        lmer(y ~ group + (1 | Donor), data = dat,
             control = lmerControl(optimizer = "bobyqa"))
      else
        lm(y ~ group, data = dat)
    }, error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    res <- tryCatch({
      tidy(fit, effects = "fixed") %>%
        dplyr::filter(term == "groupSCβ") %>%
        dplyr::select(estimate, std.error, statistic,
                      df = any_of("df"), p.value) %>%
        mutate(std_estimate = estimate / y_sd)
    }, error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) return(NULL)
    if (!"df" %in% names(res)) res$df <- nrow(dat) - 2
    res
  }) %>% Filter(Negate(is.null), .)
  if (length(results) == 0) return(NULL)
  results <- bind_rows(results)
  m <- nrow(results)
  Q_bar  <- mean(results$estimate)
  U_bar  <- mean(results$std.error^2)
  B      <- ifelse(m > 1, var(results$estimate), 0)
  T_var  <- U_bar + (1 + 1/m) * B
  se     <- sqrt(T_var)
  t_stat <- Q_bar / se
  df_rb  <- ifelse(B > 0,
                   (m - 1) * (1 + U_bar / ((1 + 1/m) * B))^2,
                   mean(results$df, na.rm = TRUE))
  tibble(variable     = var,
         estimate     = Q_bar,
         std_estimate = mean(results$std_estimate),
         se = se, t_stat = t_stat, df = df_rb,
         p_value      = 2 * pt(-abs(t_stat), df = df_rb),
         n_imp        = m)
}

ephys_within <- lapply(ephys_vars, run_lmer_within_study,
                       imp_list = completed_list) %>%
  bind_rows() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    direction = case_when(
      std_estimate > 0 & p_adj <= 0.05 ~ "Higher in SCβ",
      std_estimate < 0 & p_adj <= 0.05 ~ "Higher in beta",
      TRUE                              ~ "Unchanged"
    )
  ) %>%
  left_join(ephys_names, by = c("variable" = "name_in_xq_sheets")) %>%
  mutate(Pretty.names = ifelse(is.na(Pretty.names), variable, Pretty.names))

# ═══════════════════════════════════════════════════════════════════════════════
# Section 7: Cross-block correlation heatmap with LME4 annotations
# ═══════════════════════════════════════════════════════════════════════════════

# ── Column annotations (SCENIC + mRNA features) ───────────────────────────────
# Unified effect size and direction from LME4, on same scale across modalities
col_direction_unified <- bind_rows(
  # SCENIC regulons
  scenic_lmer_results_all %>%
    mutate(feature = str_replace(feature, "-\\(\\+\\)", "(+)")) %>%
    dplyr::filter(feature %in% colnames(mat_plot)) %>%
    mutate(direction = case_when(
      std_estimate >  0 & p_adj <= 0.05 ~ "Higher in SCβ",
      std_estimate <  0 & p_adj <= 0.05 ~ "Higher in beta",
      TRUE                               ~ "Unchanged"
    )) %>%
    dplyr::select(feature, std_estimate, direction),
  # mRNA features (from DE triage — same standardised β scale)
  mrna_lmer_results %>%
    mutate(direction = case_when(
      std_estimate >  0 & p_adj <= 0.05 ~ "Higher in SCβ",
      std_estimate <  0 & p_adj <= 0.05 ~ "Higher in beta",
      TRUE                               ~ "Unchanged"
    )) %>%
    dplyr::select(feature, std_estimate, direction)
) %>% column_to_rownames("feature")

col_annot_full <- data.frame(
  Direction     = col_direction_unified[colnames(mat_plot), "direction"],
  `Effect size` = col_direction_unified[colnames(mat_plot), "std_estimate"],
  row.names     = colnames(mat_plot),
  check.names   = FALSE
)

# ── Row annotations (ephys features) ─────────────────────────────────────────
ephys_row_annot <- ephys_within %>%
  dplyr::filter(Pretty.names %in% rownames(mat_plot)) %>%
  dplyr::select(Pretty.names, std_estimate, direction) %>%
  column_to_rownames("Pretty.names")

# Flip std_estimate sign for current amplitude rows to match flipped matrix
flipped_pretty_clean <- str_remove(to_flip, " †$")
ephys_row_annot <- ephys_row_annot %>%
  mutate(std_estimate = ifelse(
    rownames(.) %in% flipped_pretty_clean, -std_estimate, std_estimate
  ))
ephys_row_annot <- ephys_row_annot[
  str_remove(rownames(mat_plot), " †$"), , drop = FALSE
]

row_annot <- data.frame(
  `Effect size` = ephys_row_annot$std_estimate,
  Direction     = ephys_row_annot$direction,
  row.names     = rownames(mat_plot),
  check.names   = FALSE
)

# ── Shared diverging colour scale for effect sizes ────────────────────────────
# Row and column effect sizes are both standardised β (estimate / SD),
# comparable across modalities — use a single shared colour function
all_effect <- c(col_annot_full$`Effect size`, row_annot$`Effect size`)
max_abs    <- max(abs(all_effect), na.rm = TRUE)

col_fun_effect <- colorRamp2(
  c(-max_abs, 0, max_abs),
  c("#053059", "#F5F5F5", "#A6D278")
)

direction_cols <- c(
  "Higher in SCβ"  = "#A6D278",
  "Higher in beta" = "#053059",
  "Unchanged"      = "#C0CFDE"
)

col_fun_heatmap <- colorRamp2(
  c(0.3, 0.6, 0.9),
  scico::scico(3, palette = "devon", begin = 0.4)
)

# ── Build and save heatmap ────────────────────────────────────────────────────
col_ha <- HeatmapAnnotation(
  `Effect size` = col_annot_full$`Effect size`,
  Direction     = col_annot_full$Direction,
  col = list(`Effect size` = col_fun_effect, Direction = direction_cols),
  na_col = "grey90",
  annotation_name_side = "right",
  annotation_name_gp   = gpar(fontsize = 8),
  show_legend          = c(`Effect size` = TRUE, Direction = TRUE)
)

row_ha <- rowAnnotation(
  `Effect size` = row_annot$`Effect size`,
  Direction     = row_annot$Direction,
  col = list(`Effect size` = col_fun_effect, Direction = direction_cols),
  annotation_name_side = "top",
  annotation_name_gp   = gpar(fontsize = 8),
  # Suppress duplicate legend — row and col use same scale
  show_legend = FALSE
)

ht <- Heatmap(
  mat_plot,
  name                        = "Correlation\n(DIABLO)",
  col                         = col_fun_heatmap,
  top_annotation              = col_ha,
  left_annotation             = row_ha,
  cluster_rows                = TRUE,
  cluster_columns             = TRUE,
  clustering_distance_rows    = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows      = "ward.D2",
  clustering_method_columns   = "ward.D2",
  show_row_names              = TRUE,
  show_column_names           = TRUE,
  row_names_side              = "right",
  column_names_side           = "bottom",
  column_names_rot            = 45,
  row_names_gp                = gpar(fontsize = 10),
  column_names_gp             = gpar(fontsize = 8),
  column_title                = "Ephys \u00d7 Cross-block correlations (DIABLO)",
  column_title_gp             = gpar(fontsize = 12, fontface = "bold"),
  row_gap                     = unit(2, "mm"),
  border                      = FALSE,
  heatmap_legend_param        = list(title         = "Correlation\n(DIABLO)",
                                     legend_height = unit(4, "cm"))
)

svg("./for_manuscript/20260327_DIABLO_ephys_crossblock_heatmap_Unchanged.svg",
    width = 14, height = 4)
draw(ht,
     heatmap_legend_side    = "right",
     annotation_legend_side = "right",
     merge_legend           = TRUE)
dev.off()

cat("\nDone. Key outputs:\n")
cat("  diablo_50imps/feature_stability_comp1to4.csv\n")
cat("  SCENIC_AUC_lmer_beta_vs_SCbeta_JM.csv\n")
cat("  20260327_DIABLO_ephys_crossblock_heatmap_Unchanged.svg\n")
