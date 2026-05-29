# ==============================================================================
# 08_MOFA2.R
# ==============================================================================
# Unsupervised multi-omics factor analysis (MOFA2) across 50 imputation
# datasets. Extracts the primary separation factor (maximally separating
# SCβ/S6D1 from primary endocrine cells), aligns signs across runs, computes
# consensus factor scores, and generates UMAP embedding and lollipop plot
# of factor loadings.
#
# Corresponds to Methods: "Unsupervised multi-omics integration (MOFA2)"
#   "MOFA2 was run on single-cell data across 50 imputation datasets...
#    For each model, the factor with the maximum standardised effect size
#    separating SCβ/S6D1 cells from primary endocrine cells was selected
#    as the primary separation factor, after sign alignment. Consensus
#    factor scores were computed as the mean across 50 runs..."
#
# Input:
#   MOFA2_imputations/mofa_imp_*.hdf5
#     - One MOFA2 model per imputation dataset
#     - Generated externally using run_mofa() on each completed_list dataset
#   20260113_cells_plus_filtered_mi_SCENIC.rds
#     - cells.plus: for cell metadata
#
# Output:
#   sep_runs.rds
#     - Per-run separation factor info, loadings, and factor scores
#   MOFA2_consensus_factors_with_metadata.csv
#     - Consensus factor scores (mean across 50 runs) per cell with metadata
#     - Input to scripts/06_SCENIC_networks_control_energy.R (LASSO section)
#     - Input to python/04_control_energy_bootstrap.py
#   MOFA2_factor1_feature_weights.csv
#     - Feature weight stability metrics across 50 runs
#   for_manuscript/20260126_MOFA_factor1_violin_consensus.svg
#     - Violin plot of consensus separation factor by cell type/stage
#   for_manuscript/20260126_MOFA_umap_consensus.svg
#     - UMAP of consensus factor scores colored by cell type
#   for_manuscript/20260126_MOFA_factor2_loadings.svg
#     - Lollipop plot of factor loadings (ephys + SCENIC, top 25 per view)
#
# Dependencies:
#   MOFA2, tidyverse, matrixStats, data.table, uwot, ggpubr, ggplot2, scico
#
# Note on MOFA2 model training:
#   Model training (run_mofa()) is not included in this script — it is
#   computationally intensive and was run as a separate job. This script
#   assumes 50 trained .hdf5 model files exist in MOFA2_imputations/.
#   See MOFA2 documentation for training code.
# ==============================================================================

suppressPackageStartupMessages({
  library("MOFA2")
  library("tidyverse")
  library("matrixStats")
  library("data.table")
  library("uwot")
  library("ggpubr")
  library("scico")
})

# ═══════════════════════════════════════════════════════════════════════════════
# Section 0: Prepare input matrices and train MOFA2 models (50 imputations)
# ═══════════════════════════════════════════════════════════════════════════════
# This section trains one MOFA2 model per imputation dataset, saving each to
# MOFA2_imputations/mofa_imp_*.hdf5. This is computationally intensive
# (~hours total) and is run once; Section 1+ loads the saved models.
#
# Input views per model:
#   mRNA:   SCTransform scale.data (cells × top 4000 variable genes)
#           filtered by LOSO stability (unstable + weak genes removed)
#   ephys:  Imputed electrophysiology (cells × 13 parameters, sign-flipped
#           for inward currents so positive = more active)
#   scenic: AUC matrix (cells × regulons, NZV-filtered, scaled)
#
# S6D1 cells have no electrophysiology measurements — their ephys rows are
# padded with NA. MOFA2 handles missing values natively via expectation
# maximization.

imp_bundle     <- readRDS("20260120_imp_bundle.rds")
completed_list <- imp_bundle$completed_list

cells.plus <- readRDS("20260113_cells_plus_filtered_mi_SCENIC.rds")

# ── Inward current sign flip ──────────────────────────────────────────────────
# Flip before building MOFA input so positive = more electrophysiologically
# active (consistent with correlation and DIABLO analyses)
flip_vars_mofa <- c(
  "NormalizedPeakSodiumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedEarlyPeakCalciumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedLateCalciumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedHVA_.20mV_pA.pF",
  "NormalizedLVA_.60mV_pA.pF",
  "NormalizedHyperpolarizationactivatedcurrent.at.140mV_pA.pF",
  "CalciumIntegralNormalizedtoCellSize_pC.pF"
)

ephys_by_imp_mofa <- lapply(seq_along(completed_list), function(k) {
  d <- completed_list[[k]] %>%
    dplyr::filter(CellID %in% colnames(cells.plus))
  rownames(d) <- d$CellID
  d2 <- d[colnames(cells.plus), ephys_vars, drop = FALSE]
  rownames(d2) <- colnames(cells.plus)
  for (v in intersect(flip_vars_mofa, colnames(d2)))
    d2[, v] <- -d2[, v]
  d2$VoltageforSodiumPeakCurrent_mV <- as.numeric(d2$VoltageforSodiumPeakCurrent_mV)
  d2
})

# ── mRNA matrix preparation ───────────────────────────────────────────────────
# SCTransform scale.data; HGNC symbols; ribosomal genes removed;
# LOSO unstable + weak genes dropped; top 4000 by variance retained
mrna_mtx <- t(cells.plus[["SCT"]]$scale.data)

ensg  <- colnames(mrna_mtx)
symbs <- mapIds(EnsDb.Hsapiens.v86, keys = ensg, column = "SYMBOL",
                keytype = "GENEID", multiVals = "first")
symbs[is.na(symbs)]        <- names(symbs[is.na(symbs)])
symbs[duplicated(symbs)]   <- names(symbs[duplicated(symbs)])
colnames(mrna_mtx) <- as.character(symbs)
mrna_mtx <- mrna_mtx[, !str_detect(colnames(mrna_mtx), "^RP[LS]")]

# Remove LOSO-unstable, weak-effect genes
loso_summary <- readRDS("loso_summary.rds")
unstable_cut <- quantile(loso_summary$loso_rel_iqr, 0.90, na.rm = TRUE)
weak_cut     <- quantile(abs(loso_summary$signed_f_full), 0.50, na.rm = TRUE)
genes_drop   <- loso_summary %>%
  dplyr::filter(loso_rel_iqr >= unstable_cut,
                abs(signed_f_full) <= weak_cut) %>%
  pull(gene_name) %>% unique()
mrna_mtx <- mrna_mtx[, !colnames(mrna_mtx) %in% genes_drop, drop = FALSE]

# Keep top 4000 most variable genes
gene_var   <- apply(mrna_mtx, 2, var, na.rm = TRUE)
keep_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(4000, length(gene_var)))]
mrna_mtx   <- mrna_mtx[, keep_genes, drop = FALSE]
cat("mRNA features for MOFA:", ncol(mrna_mtx), "\n")

# ── SCENIC AUC matrix ─────────────────────────────────────────────────────────
auc_mtx_full   <- cells.plus[["AUC"]]$counts %>% t() %>% as.matrix()
nzv            <- nearZeroVar(auc_mtx_full, freqCut = 95/5, uniqueCut = 10)
auc_mtx_scaled <- scale(auc_mtx_full[, -nzv$Position])
cat("SCENIC features for MOFA:", ncol(auc_mtx_scaled), "\n")

# ── Row padding helper ────────────────────────────────────────────────────────
# Ensures all views have the same sample set.
# S6D1 cells are padded with NA in the ephys view.
.pad_rows_to <- function(mat, target_samples) {
  stopifnot(!is.null(rownames(mat)))
  out <- matrix(NA_real_, nrow = length(target_samples), ncol = ncol(mat),
                dimnames = list(target_samples, colnames(mat)))
  common     <- intersect(target_samples, rownames(mat))
  out[common, ] <- mat[common, , drop = FALSE]
  out
}

# ── MOFA2 training function ───────────────────────────────────────────────────
run_mofa_one_imp <- function(k, ephys_by_imp_mofa, mrna_mtx, auc_mtx_scaled,
                              out_dir = "./MOFA2_imputations",
                              num_factors = 10, seed = 42,
                              convergence_mode = "slow",
                              use_basilisk = TRUE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  mrna_samples   <- colnames(t(mrna_mtx))
  scenic_samples <- colnames(t(auc_mtx_scaled))
  target_samples <- intersect(mrna_samples, scenic_samples)
  if (length(target_samples) < 3)
    stop("Too few overlapping samples: ", length(target_samples))

  ephys_k        <- ephys_by_imp_mofa[[k]]
  ephys_scaled_k <- scale(as.matrix(ephys_k))
  ephys_scaled_k <- .pad_rows_to(ephys_scaled_k, target_samples)

  MOFA_data <- list(
    mRNA   = t(mrna_mtx)[, target_samples, drop = FALSE],
    ephys  = t(ephys_scaled_k),
    scenic = t(auc_mtx_scaled)[, target_samples, drop = FALSE]
  )

  stopifnot(
    identical(colnames(MOFA_data$mRNA),  colnames(MOFA_data$ephys)),
    identical(colnames(MOFA_data$mRNA),  colnames(MOFA_data$scenic))
  )

  MOFAobject <- create_mofa(MOFA_data)
  model_opts              <- get_default_model_options(MOFAobject)
  model_opts$num_factors  <- num_factors
  train_opts              <- get_default_training_options(MOFAobject)
  train_opts$convergence_mode <- convergence_mode
  train_opts$seed         <- seed + k   # varies per imputation; reproducible

  MOFAobject <- prepare_mofa(
    MOFAobject,
    data_options     = get_default_data_options(MOFAobject),
    model_options    = model_opts,
    training_options = train_opts
  )

  h5_file  <- file.path(out_dir, sprintf("mofa_imp_%02d.hdf5",  k))
  rds_file <- file.path(out_dir, sprintf("mofa_imp_%02d_diagnostics.rds", k))

  MOFAobject <- run_mofa(MOFAobject, outfile = h5_file,
                          use_basilisk = use_basilisk)

  diag <- list(
    k                  = k,
    h5_file            = h5_file,
    target_samples     = target_samples,
    variance_explained = MOFA2::calculate_variance_explained(MOFAobject),
    factors            = MOFA2::get_factors(MOFAobject, factors = "all")$group1,
    weights            = MOFA2::get_weights(MOFAobject, views = "all", factors = "all")
  )
  saveRDS(diag, rds_file)

  list(k = k, h5_file = h5_file, rds_file = rds_file,
       n_samples = length(target_samples), n_factors = num_factors)
}

# ── Train 50 MOFA2 models ─────────────────────────────────────────────────────
cat("Training MOFA2 models across 50 imputations...\n")
mofa_results <- lapply(seq_along(completed_list), function(k) {
  message("Running MOFA for imputation ", k, "/", length(completed_list))
  run_mofa_one_imp(
    k                  = k,
    ephys_by_imp_mofa  = ephys_by_imp_mofa,
    mrna_mtx           = mrna_mtx,
    auc_mtx_scaled     = auc_mtx_scaled,
    out_dir            = "./MOFA2_imputations",
    num_factors        = 10,
    seed               = 42,
    convergence_mode   = "slow",
    use_basilisk       = TRUE
  )
})

saveRDS(mofa_results,
        file = "./MOFA2_imputations/mofa_imputation_run_manifest.rds")
cat("All models trained. Manifest saved.\n")
cat("Model files: MOFA2_imputations/mofa_imp_01.hdf5 ... mofa_imp_50.hdf5\n\n")
# End of Section 0. Sections 1+ load the saved .hdf5 files.

# ── Colour palette ────────────────────────────────────────────────────────────
cluster_cols2 <- c(
  "#ff7db1","#e2ff50","#713be8","#669b00","#c734e9","#01a457","#e400b3",
  "#8cffb0","#010f92","#ffd26e","#015fde","#005a02","#ff6de1","#01deaf",
  "#9d0098","#4f6600","#a282ff","#867100","#170047","#b06100","#01b4e8",
  "#ff5662","#02cfcb","#ff3a85","#006c55","#ac004a","#baecff","#862200",
  "#006d95","#ff7d56","#003767","#ffcaa5","#33002d","#ffd6db","#101c00",
  "#f6bbff","#001e22","#ff9da8","#360c00","#5f3800","#001024","#d0f700"
)

# Cell type stage color scale (consistent across figures)
pop_colors_MOFA <-c(cluster_cols2[1:6])

# Block colors for lollipop plot
block_colours <- c(
  "ephys"  = scico(3, palette = "buda", end = 0.8)[1],
  "scenic" = scico(3, palette = "buda", end = 0.8)[3]
)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 1: Load data and metadata
# ═══════════════════════════════════════════════════════════════════════════════

# cells.plus already loaded in Section 0
combined_metadata <- cells.plus@meta.data[, c(2:3, 22, 53, 58, 60, 100, 169, 170)]
combined_metadata$sample <- colnames(cells.plus)
combined_metadata$CellID <- combined_metadata$sample

# Load electrophysiology pretty-name lookup (for lollipop plot labels)
ephys_names <- readxl::read_excel("250404_patch_key.xlsx",
                                   .name_repair = "universal") %>%
  dplyr::select(Pretty.names, name_in_xq_sheets)

# ── Load MOFA2 model files ────────────────────────────────────────────────────
models_dir  <- "MOFA2_imputations"
model_files <- list.files(models_dir, pattern = "^mofa_imp_\\d+\\.hdf5$",
                           full.names = TRUE)
stopifnot(length(model_files) == 50)
cat("MOFA2 models found:", length(model_files), "\n")

# ── Define SCβ/S6D1 vs endocrine group labels ─────────────────────────────────
# Separation factor = factor with maximum standardised effect size between
# SCβ/S6D1 cells and primary endocrine cells
grp <- combined_metadata %>%
  transmute(
    CellID,
    group_sep = ifelse(celltype %in% c("SCβ", "S6D1"), "SCB_S6D1", "Endocrine")
  )
cell_ids <- grp$CellID

# ═══════════════════════════════════════════════════════════════════════════════
# Section 2: Extract separation factor from each model
# ═══════════════════════════════════════════════════════════════════════════════
# For each model: find the factor that best separates SCβ/S6D1 from endocrine
# cells (by standardised Cohen's d-like score), flip sign if needed for
# consistency (SCB_S6D1 positive), and extract weights for mRNA/ephys/scenic.

score_factor_sep <- function(fvec, group_bin) {
  # group_bin: TRUE for SCB_S6D1, FALSE for Endocrine
  x1 <- fvec[group_bin]
  x0 <- fvec[!group_bin]
  if (length(x1) < 5 || length(x0) < 5) return(NA_real_)
  # Standardised effect size: equivalent to Cohen's d scaled by pooled SD
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / sd(fvec, na.rm = TRUE)
}

extract_sep_factor_and_weights <- function(model_file, grp_df,
                                           views = c("mRNA", "ephys", "scenic")) {
  m <- load_model(model_file)
  Z <- get_factors(m, factors = "all", groups = "all")[[1]]
  Z <- Z[grp_df$CellID, , drop = FALSE]

  group_bin <- grp_df$group_sep == "SCB_S6D1"
  scores    <- apply(Z, 2, score_factor_sep, group_bin = group_bin)
  j         <- which.max(abs(scores))
  best_factor <- colnames(Z)[j]
  best_score  <- scores[j]

  # Flip all factors if best score is negative (SCB_S6D1 should be positive)
  if (best_score < 0) {
    Z          <- -Z
    flip       <- TRUE
    best_score <- -best_score
  } else {
    flip <- FALSE
  }

  # Extract and optionally flip weights for the separation factor
  W_list <- lapply(views, function(v) {
    w <- get_weights(m, views = v, factors = best_factor, as.data.frame = TRUE)
    w$view <- v
    w
  }) %>% dplyr::bind_rows()

  if (flip) W_list$value <- -W_list$value

  list(
    file        = model_file,
    factor      = best_factor,
    sep_score   = best_score,
    flipped     = flip,
    Z_sep       = Z[, best_factor, drop = FALSE],
    Z_all       = Z,
    W_sep       = W_list
  )
}

cat("Extracting separation factors across 50 models...\n")
sep_runs <- lapply(model_files, extract_sep_factor_and_weights,
                   grp_df = grp, views = c("mRNA", "ephys", "scenic"))

saveRDS(sep_runs, file = "sep_runs.rds")

# ── Separation factor summary ─────────────────────────────────────────────────
sep_meta <- tibble(
  run       = seq_along(sep_runs),
  file      = vapply(sep_runs, `[[`, character(1), "file"),
  factor    = vapply(sep_runs, `[[`, character(1), "factor"),
  sep_score = vapply(sep_runs, `[[`, numeric(1),   "sep_score"),
  flipped   = vapply(sep_runs, `[[`, logical(1),   "flipped")
)

cat("Factor selection summary:\n")
print(sep_meta %>% count(factor) %>% arrange(desc(n)))
cat("Median separation score:", median(sep_meta$sep_score), "\n")
cat("Models flipped:", sum(sep_meta$flipped), "/ 50\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 3: Consensus factor scores
# ═══════════════════════════════════════════════════════════════════════════════

# ── Primary separation factor: mean across 50 runs ───────────────────────────
Z_mat <- do.call(cbind, lapply(sep_runs, function(x) x$Z_sep[, 1]))
rownames(Z_mat) <- cell_ids
colnames(Z_mat) <- paste0("run", seq_len(ncol(Z_mat)))

Z_consensus_sep   <- rowMeans(Z_mat, na.rm = TRUE)
Z_consensus_sep_z <- as.numeric(scale(Z_consensus_sep))

consensus_factors_df <- tibble(
  CellID  = cell_ids,
  F_sep   = Z_consensus_sep,
  F_sep_z = Z_consensus_sep_z
) %>%
  left_join(combined_metadata %>% mutate(CellID = sample), by = "CellID")

# ── All factors: mean across 50 runs (for UMAP) ──────────────────────────────
Z_list          <- lapply(sep_runs, `[[`, "Z_all")
common_factors  <- Reduce(intersect, lapply(Z_list, colnames))
Z_list          <- lapply(Z_list, function(z) z[, common_factors, drop = FALSE])
Z_arr           <- simplify2array(Z_list)

# Ensure cell order is consistent
if (!identical(dimnames(Z_arr)[[1]], grp$CellID))
  Z_arr <- aperm(Z_arr, c(2, 1, 3))

Z_consensus_all    <- apply(Z_arr, c(1, 2), mean, na.rm = TRUE)
Z_consensus_scaled <- scale(Z_consensus_all)

factors_df <- Z_consensus_all %>%
  as.data.frame() %>%
  mutate(CellID = rownames(Z_consensus_all))

factors_meta <- full_join(factors_df, combined_metadata, by = "CellID")
write.csv(factors_meta, file = "MOFA2_consensus_factors_with_metadata.csv",
          row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 4: Violin plot of consensus separation factor
# ═══════════════════════════════════════════════════════════════════════════════

ggplot(consensus_factors_df,
       aes(x = celltype_stage, y = F_sep, fill = celltype_stage)) +
  geom_violin(trim = FALSE) +
  geom_jitter(aes(fill = celltype_stage), colour = "grey10",
              alpha = 0.3, width = 0.05, shape = 21) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = 2) +
  ggpubr::theme_classic2() +
  theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = cluster_cols2) +
  scale_color_manual(values = cluster_cols2) +
  labs(y = "MOFA2 Latent Factor")

ggsave("./for_manuscript/20260126_MOFA_factor1_violin_consensus.svg",
       width = 9, height = 3.5)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 5: UMAP of consensus factor scores
# ═══════════════════════════════════════════════════════════════════════════════

set.seed(1)
um <- uwot::umap(
  X          = Z_consensus_scaled,
  n_neighbors = 30,
  min_dist    = 0.3,
  metric      = "euclidean",
  verbose     = TRUE
)

umap_df <- grp %>%
  mutate(UMAP1 = um[, 1], UMAP2 = um[, 2]) %>%
  left_join(combined_metadata, by = "CellID")

ggplot(umap_df, aes(UMAP1, UMAP2, fill = celltype)) +
  geom_point(shape = 21, size = 1.5, stroke = 0.05) +
  scale_fill_manual(values = pop_colors_MOFA) +
  theme_void()

ggsave("./for_manuscript/20260126_MOFA_umap_consensus.svg",
       width = 4.5, height = 3.5)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 6: Feature weight stability
# ═══════════════════════════════════════════════════════════════════════════════
# Summarize mean weight, sign stability, and mean absolute weight across runs.
# Only features with ≥90% sign consistency are shown in the lollipop plot.

all_W <- bind_rows(lapply(seq_along(sep_runs), function(i) {
  w      <- sep_runs[[i]]$W_sep
  w$run  <- i
  w
}))

W_stab <- all_W %>%
  group_by(view, feature) %>%
  summarise(
    n_runs              = n(),
    frac_pos            = mean(value > 0, na.rm = TRUE),
    frac_neg            = mean(value < 0, na.rm = TRUE),
    frac_consistent_sign = pmax(frac_pos, frac_neg),
    mean_value          = mean(value,     na.rm = TRUE),
    sd_value            = sd(value,       na.rm = TRUE),
    mean_abs            = mean(abs(value), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs))

write_csv(W_stab, "./for_manuscript/MOFA2_factor1_feature_weights.csv")

# Top 25 per view (ephys + SCENIC only; mRNA excluded from lollipop)
# with ≥90% sign consistency
W_stable_top <- W_stab %>%
  dplyr::filter(frac_consistent_sign >= 0.9) %>%
  group_by(view) %>%
  slice_max(order_by = mean_abs, n = 25) %>%
  ungroup()

top_features <- W_stable_top %>%
  transmute(feature, value = mean_value, view)

# Join pretty names for ephys features
top_features_plot <- top_features %>%
  left_join(W_stab %>% dplyr::select(feature, view, sd_value, frac_consistent_sign),
            by = c("feature", "view")) %>%
  left_join(ephys_names, by = c("feature" = "name_in_xq_sheets")) %>%
  dplyr::filter(view != "mRNA") %>%
  mutate(
    # Use Pretty.names for ephys; feature name for SCENIC
    Pretty.names = case_when(
      !is.na(Pretty.names) ~ Pretty.names,
      TRUE                 ~ feature
    ),
    Pretty.names = fct_reorder(Pretty.names, value),
    view         = factor(view, levels = c("ephys", "scenic"))
  )

# ═══════════════════════════════════════════════════════════════════════════════
# Section 7: Lollipop plot of MOFA2 factor loadings
# ═══════════════════════════════════════════════════════════════════════════════
# Positive loadings → SCβ/S6D1 state; negative loadings → endocrine state.
# Error bars = SD across 50 imputation runs.
# Features shown: top 25 per view (ephys + SCENIC), ≥90% sign consistency.

y_n <- nrow(top_features_plot)

p <- ggplot(top_features_plot,
            aes(x = value, y = Pretty.names, colour = view)) +
  # Zero line
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey50") +
  # Lollipop stems
  geom_segment(aes(x = 0, xend = value, yend = Pretty.names),
               linewidth = 0.6, alpha = 0.7) +
  # Points
  geom_point(aes(shape = view), size = 3.5) +
  # Population direction labels
  annotate("text",
           x = -0.35, y = y_n + 0.8,
           label = "\u2190 Endocrine cells",
           hjust = 0.5, size = 3.8,
           colour = "grey35", fontface = "italic") +
  annotate("text",
           x = 0.35, y = y_n + 0.8,
           label = "Beta-like cells \u2192",
           hjust = 0.5, size = 3.8,
           colour = "grey35", fontface = "italic") +
  scale_colour_manual(values = block_colours) +
  scale_shape_manual(values = c("ephys" = 16, "scenic" = 17)) +
  facet_grid(view ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "MOFA2 Factor loading",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(strip.text.y = element_text(angle = 0),
        panel.grid.minor = element_blank())

ggsave("./for_manuscript/20260126_MOFA_factor2_loadings.svg",
       p, width = 7, height = 9)

cat("\nDone. Key outputs:\n")
cat("  MOFA2_consensus_factors_with_metadata.csv\n")
cat("  MOFA2_factor1_feature_weights.csv\n")
cat("  for_manuscript/20260126_MOFA_factor1_violin_consensus.svg\n")
cat("  for_manuscript/20260126_MOFA_umap_consensus.svg\n")
cat("  for_manuscript/20260126_MOFA_factor2_loadings.svg\n")
