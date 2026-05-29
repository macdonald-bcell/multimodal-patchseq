# ==============================================================================
# 05_correlations.R
# ==============================================================================
# Bootstrapped Spearman correlations between gene expression and electrophysiology
# parameters across primary β-cells (meta-analytic, multi-study) and SCβ-cells
# (JM_patchSeq only). Computes Δz = z_SCβ − z_β as the primary metric of
# differential transcript–electrophysiology coupling.
#
# Corresponds to Methods: "Bootstrapped Spearman Correlations and Δz Analysis"
#   "...bootstrapped Spearman correlations were computed between log-normalized
#    gene expression and electrophysiology parameters. For primary β-cells,
#    within-study correlations were combined across five studies using a
#    DerSimonian-Laird random-effects meta-analytic framework on Fisher
#    z-transformed correlations. For SCβ-cells, correlations were computed
#    within the JM_patchSeq study only. Δz = z_SCβ − z_β was used as the
#    primary metric of differential coupling..."
#
# Input:
#   20260120_imp_bundle.rds
#     - imp_bundle$completed_list: 50 imputed electrophysiology datasets
#     - From scripts/03_imputation.R
#   20260113_cells_plus_filtered_mi_SCENIC.rds
#     - cells.plus: integrated Seurat object with AUC assay
#     - From scripts/02_integrate_patchseq_datasets.R + pySCENIC pipeline
#   250404_patch_key.xlsx
#     - Electrophysiology parameter pretty-name lookup table
#
# Output:
#   20260423_boot_store_REbeta_MI_sample_repaired.rds
#     - Raw bootstrap results (200 iterations × all ephys params × all genes)
#   for_manuscript/20260423_boot_meta_dz_all_params_gene_pairs.csv
#     - Summary table: dz, r_beta, r_scb, sign_stability per gene × param pair
#     - Primary analytical output; feeds scripts/07_DIABLO.R and figures
#   for_manuscript/fig4_*.pdf        - Main figure panels
#   for_manuscript/suppfig_*.pdf     - Supplementary figure panels
#
# Dependencies:
#   tidyverse, Seurat, Matrix, matrixStats, future, future.apply,
#   EnsDb.Hsapiens.v86, ensembldb, ComplexHeatmap, circlize, ggalt, ggtext,
#   patchwork, cowplot, viridis, scico, fgsea, msigdbr, tidytext, readxl
#
#
# Runtime note:
#   200 bootstrap iterations in parallel require ~20GB memory.
#   options(future.globals.maxSize = 20 * 1024^3) is set accordingly.
#   Adjust workers and memory limit to match your system.
# ==============================================================================

library("tidyverse")
library("Seurat")
library("Matrix")
library("matrixStats")
library("future")
library("future.apply")
library("EnsDb.Hsapiens.v86")
library("ensembldb")
library("ComplexHeatmap")
library("circlize")
library("ggalt")
library("ggtext")
library("patchwork")
library("cowplot")
library("viridis")
library("scico")
library("fgsea")
library("msigdbr")
library("tidytext")
library("readxl")

# ── Electrophysiology parameters ──────────────────────────────────────────────
# These are the 13 normalized parameters used throughout the correlation analysis
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

# ── Current amplitude sign convention ─────────────────────────────────────────
# These parameters are inward currents (negative by convention).
# Signs are flipped for plotting so more positive = more electrophysiologically
# active. Applied to z_scb, z_beta, dz, and r values before visualization.
flip_vars <- c(
  "NormalizedPeakSodiumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedEarlyPeakCalciumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedLateCalciumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedHVA_.20mV_pA.pF",
  "NormalizedLVA_.60mV_pA.pF",
  "NormalizedHyperpolarizationactivatedcurrent.at.140mV_pA.pF",
  "CalciumIntegralNormalizedtoCellSize_pC.pF"
)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 1: Load data
# ═══════════════════════════════════════════════════════════════════════════════

imp_bundle     <- readRDS("20260120_imp_bundle.rds")
completed_list <- imp_bundle$completed_list

cells.plus <- readRDS("20260113_cells_plus_filtered_mi_SCENIC.rds")

# Restrict to patched cells with valid cell size; exclude ERCC spike-ins and S6D1
cells.combined <- subset(cells.plus, celltype != "S6D1")
DefaultAssay(cells.combined) <- "RNA"

cells_to_keep <- rownames(cells.combined@meta.data)[
  !is.na(cells.combined$CellSize_pF) & cells.combined$CellSize_pF != 0
]
genes_to_keep <- rownames(cells.combined)[
  !str_detect(rownames(cells.combined), "^ERCC-\\d+")
]
cells.combined.patched <- subset(cells.combined,
                                  cells    = cells_to_keep,
                                  features = genes_to_keep)
cells.combined.patched.sub <- subset(cells.combined.patched,
                                      celltype %in% c("beta", "SCβ"))

# ═══════════════════════════════════════════════════════════════════════════════
# Section 2: Gene filtering
# ═══════════════════════════════════════════════════════════════════════════════

id_col       <- "CellID"
study_col    <- "Study"
celltype_col <- "celltype"

cells_beta_all <- cells.combined.patched.sub$CellID[
  cells.combined.patched.sub[[celltype_col]][, 1] == "beta"
]
cells_scb_jm <- cells.combined.patched.sub$CellID[
  cells.combined.patched.sub[[celltype_col]][, 1] == "SCβ"
]

meta_all <- cells.combined.patched.sub@meta.data
rownames(meta_all) <- meta_all$CellID

# Log-normalize RNA counts
DefaultAssay(cells.combined.patched.sub) <- "RNA"
cells.combined.patched.sub <- NormalizeData(cells.combined.patched.sub)
cells.combined.patched.sub <- JoinLayers(cells.combined.patched.sub)
rna <- cells.combined.patched.sub[["RNA"]]$data

cells_use        <- unique(c(cells_beta_all, cells_scb_jm))
colnames(rna)    <- cells.combined.patched.sub$CellID %>% unname()
rna              <- rna[, cells_use, drop = FALSE]

# Detection filter: gene must be detected in ≥10% of cells in BOTH populations
expr_min_pct <- 0.10
pct_beta     <- Matrix::rowMeans(rna[, cells_beta_all, drop = FALSE] > 0)
pct_scb      <- Matrix::rowMeans(rna[, cells_scb_jm,  drop = FALSE] > 0)
genes_det    <- c(
  names(pct_beta)[pct_beta >= expr_min_pct & pct_scb >= expr_min_pct],
  "eGFP"   # retain eGFP explicitly as a SCβ identity marker
)

# SD filter: remove bottom 10% lowest-variance genes in BOTH populations
# (gene must pass threshold in at least one population to be retained)
sd_beta   <- matrixStats::rowSds(as.matrix(rna[genes_det, cells_beta_all, drop = FALSE]))
sd_scb    <- matrixStats::rowSds(as.matrix(rna[genes_det, cells_scb_jm,   drop = FALSE]))
names(sd_beta) <- genes_det
names(sd_scb)  <- genes_det

cut_beta   <- quantile(sd_beta, probs = 0.10)
cut_scb    <- quantile(sd_scb,  probs = 0.10)
genes_keep <- names(sd_beta)[sd_beta > cut_beta | sd_scb > cut_scb]

message("Genes retained after detection + SD filters: ", length(genes_keep))

# Gene symbol lookup
tr2g <- tibble(
  gene      = row.names(rna),
  gene_name = mapIds(EnsDb.Hsapiens.v86,
                     keys      = row.names(rna),
                     column    = "SYMBOL",
                     keytype   = "GENEID",
                     multiVals = "first")
)

# Build dense cells × genes expression matrix
expr_mat <- t(as.matrix(rna[genes_keep, , drop = FALSE]))
storage.mode(expr_mat) <- "double"
stopifnot(identical(rownames(expr_mat), cells_use))
gene_names <- colnames(expr_mat)

# Align imputed ephys data to expression matrix row order
m            <- length(completed_list)
ephys_by_imp <- lapply(seq_len(m), function(k) {
  d              <- completed_list[[k]]
  d[[id_col]]    <- as.character(d[[id_col]])
  rownames(d)    <- d[[id_col]]
  stopifnot(all(ephys_vars %in% colnames(d)))
  d2             <- d[cells_use, c(id_col, ephys_vars), drop = FALSE]
  rownames(d2)   <- d2[[id_col]]
  d2             <- d2[, ephys_vars, drop = FALSE]
  d2             <- d2[rownames(expr_mat), , drop = FALSE]
  stopifnot(identical(rownames(d2), rownames(expr_mat)))
  d2
})

# ═══════════════════════════════════════════════════════════════════════════════
# Section 3: Core statistical functions
# ═══════════════════════════════════════════════════════════════════════════════

# ── DerSimonian-Laird random-effects meta-analysis on Fisher z ────────────────
# Used to combine within-study Spearman r values across primary β-cell studies.
# Returns pooled z (z_re), standard error, I² heterogeneity, and k (study count).
# Input: Rmat (studies × genes correlation matrix), Nvec (per-study sample sizes)
meta_RE_fisher <- function(Rmat, Nvec) {
  stopifnot(is.matrix(Rmat))
  G    <- ncol(Rmat)
  keep <- which(Nvec >= 4)
  if (length(keep) == 0) {
    return(list(z  = rep(NA_real_, G), se = rep(NA_real_, G),
                I2 = rep(NA_real_, G), k  = rep(0L, G)))
  }
  Ruse  <- Rmat[keep, , drop = FALSE]
  nuse  <- Nvec[keep]
  Rclip <- pmax(pmin(Ruse, 0.999999), -0.999999)
  Z     <- atanh(Rclip)
  w_fe  <- pmax(nuse - 3, 0)
  W     <- matrix(w_fe, nrow = length(w_fe), ncol = G)
  Z_na  <- is.na(Z)
  W[Z_na] <- 0
  den_fe <- colSums(W)
  den_fe[den_fe == 0] <- NA_real_
  z_fe     <- colSums(W * Z) / den_fe
  z_fe_mat <- matrix(z_fe, nrow = nrow(Z), ncol = G, byrow = TRUE)
  Q        <- colSums(W * (Z - z_fe_mat)^2)
  k        <- colSums(!Z_na)
  df       <- pmax(k - 1, 0)
  C        <- den_fe - (colSums(W^2) / den_fe)
  C[C <= 0] <- NA_real_
  tau2     <- pmax((Q - df) / C, 0)
  var_i    <- 1 / W
  var_i[W == 0] <- Inf
  tau2_mat <- matrix(tau2, nrow = nrow(Z), ncol = G, byrow = TRUE)
  W_re     <- 1 / (var_i + tau2_mat)
  W_re[Z_na] <- 0
  den_re   <- colSums(W_re)
  den_re[den_re == 0] <- NA_real_
  z_re     <- colSums(W_re * Z) / den_re
  se_re    <- sqrt(1 / den_re)
  I2       <- rep(0, G)
  okQ      <- is.finite(Q) & Q > 0 & (Q > df)
  I2[okQ]  <- pmax((Q[okQ] - df[okQ]) / Q[okQ], 0)
  list(z = z_re, se = se_re, I2 = I2, k = k)
}

# ── Spearman correlation of expression matrix against a single vector ─────────
# Returns a vector of correlations (one per gene). NA if n < min_n or SD = 0.
spearman_vec <- function(expr_mat, y, min_n = 4) {
  y  <- as.numeric(y)
  ok <- is.finite(y)
  if (sum(ok) < min_n) return(rep(NA_real_, ncol(expr_mat)))
  y2 <- y[ok]
  X2 <- expr_mat[ok, , drop = FALSE]
  if (stats::sd(y2) == 0) return(rep(NA_real_, ncol(expr_mat)))
  r  <- suppressWarnings(
    stats::cor(y2, X2, method = "spearman", use = "pairwise.complete.obs")
  )
  r  <- as.numeric(r)
  if (length(r) != ncol(expr_mat)) {
    r2 <- rep(NA_real_, ncol(expr_mat))
    r2[seq_len(min(length(r), length(r2)))] <- r[seq_len(min(length(r), length(r2)))]
    r <- r2
  }
  r
}

# ═══════════════════════════════════════════════════════════════════════════════
# Section 4: Bootstrap function (one iteration)
# ═══════════════════════════════════════════════════════════════════════════════
# Each iteration:
#   1. Randomly samples one imputation dataset from completed_list
#   2. Resamples β-cell pseudobulk samples with replacement WITHIN each study
#      (preserving study-level structure for the meta-analysis)
#   3. Resamples SCβ cells with replacement
#   4. Computes within-study Spearman r for β-cells, combines via DL meta-analysis
#   5. Computes Spearman r for SCβ cells
#   6. Computes Δz = z_SCβ − z_β (difference in Fisher z-transformed correlations)
#      and SE(Δz) = sqrt(SE_SCβ² + SE_β²)

beta_by_study <- split(cells_beta_all, meta_all[cells_beta_all, study_col])

run_one_bootstrap <- function(b) {
  tryCatch({
    # Sample imputation dataset randomly at each iteration
    k_imp     <- sample.int(length(ephys_by_imp), 1)
    ephys_mat <- ephys_by_imp[[k_imp]]

    # Resample within each β-cell study (preserves study composition)
    boot_beta <- unlist(lapply(beta_by_study, function(cells) {
      sample(cells, length(cells), replace = TRUE)
    }), use.names = FALSE)
    boot_scb  <- sample(cells_scb_jm, length(cells_scb_jm), replace = TRUE)

    expr_beta  <- expr_mat[boot_beta, , drop = FALSE]
    expr_scb   <- expr_mat[boot_scb,  , drop = FALSE]
    study_beta <- as.character(meta_all[boot_beta, study_col])
    gene_names <- colnames(expr_mat)

    out        <- vector("list", length(ephys_vars))
    names(out) <- ephys_vars

    for (v in ephys_vars) {
      studies     <- names(beta_by_study)
      Rmat        <- matrix(NA_real_, nrow = length(studies),
                            ncol  = length(gene_names),
                            dimnames = list(studies, gene_names))
      Nvec        <- integer(length(studies))
      names(Nvec) <- studies
      y_beta_all  <- ephys_mat[boot_beta, v]

      for (s in seq_along(studies)) {
        st      <- studies[s]
        ix      <- which(study_beta == st & is.finite(y_beta_all))
        Nvec[s] <- length(ix)
        if (Nvec[s] >= 4)
          Rmat[s, ] <- spearman_vec(expr_beta[ix, , drop = FALSE],
                                    y_beta_all[ix])
      }

      # DerSimonian-Laird meta-analytic estimate across β-cell studies
      beta_RE <- meta_RE_fisher(Rmat, Nvec)

      # SCβ correlations (JM_patchSeq only — single study, no meta-analysis needed)
      y_scb  <- ephys_mat[boot_scb, v]
      r_scb  <- spearman_vec(expr_scb, y_scb)
      r_scb  <- pmax(pmin(r_scb, 0.999999), -0.999999)
      z_scb  <- atanh(r_scb)
      n_scb  <- sum(is.finite(y_scb))
      se_scb <- if (n_scb >= 4) rep(1 / sqrt(n_scb - 3), length(z_scb)) else
                  rep(NA_real_, length(z_scb))

      # Δz: key output — difference in coupling strength between SCβ and β
      dz    <- z_scb - beta_RE$z
      se_dz <- sqrt(se_scb^2 + beta_RE$se^2)

      out[[v]] <- list(
        z_scb   = setNames(as.numeric(z_scb),      gene_names),
        z_beta  = setNames(as.numeric(beta_RE$z),  gene_names),
        dz      = setNames(as.numeric(dz),         gene_names),
        se_dz   = setNames(as.numeric(se_dz),      gene_names),
        I2_beta = setNames(as.numeric(beta_RE$I2), gene_names),
        k_beta  = setNames(as.numeric(beta_RE$k),  gene_names),
        imp     = k_imp
      )
    }
    out
  }, error = function(e) {
    structure(list(error = TRUE, message = conditionMessage(e)),
              class = "boot_error")
  })
}

# ═══════════════════════════════════════════════════════════════════════════════
# Section 5: Run 200 bootstrap iterations in parallel
# ═══════════════════════════════════════════════════════════════════════════════

plan(multicore, workers = max(1, parallel::detectCores() - 1))
options(future.globals.maxSize = 20 * 1024^3)
set.seed(42)

cat("Running 200 bootstrap iterations...\n")
boot_store <- future_lapply(seq_len(200), run_one_bootstrap,
                             future.seed = TRUE)

is_err <- vapply(boot_store, inherits, logical(1), what = "boot_error")
cat("Failed iterations:", sum(is_err), "/ 200\n")
if (any(is_err)) print(boot_store[is_err][1:min(5, sum(is_err))])

saveRDS(boot_store,
        "20260423_boot_store_REbeta_MI_sample_repaired.rds")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 6: Summarize bootstrap results
# ═══════════════════════════════════════════════════════════════════════════════

is_boot_ok <- function(x) {
  is.list(x) && !inherits(x, "boot_error") && !inherits(x, "error") &&
    (is.null(x$error) || !isTRUE(x$error))
}

# Extract a genes × bootstraps matrix for a given parameter and field
get_template_genes <- function(boot_store, param, field = "dz") {
  for (b in seq_along(boot_store)) {
    bs <- boot_store[[b]]
    if (!is_boot_ok(bs)) next
    obj <- bs[[param]]
    if (!is.list(obj)) next
    v <- obj[[field]]
    if (is.null(v)) next
    nm <- names(v)
    if (!is.null(nm) && length(nm) > 0) return(nm)
  }
  character(0)
}

safe_extract_vec <- function(bs, param, field, genes) {
  if (!is_boot_ok(bs)) return(rep(NA_real_, length(genes)))
  obj <- bs[[param]]
  if (!is.list(obj)) return(rep(NA_real_, length(genes)))
  v <- obj[[field]]
  if (is.null(v)) return(rep(NA_real_, length(genes)))
  if (!is.null(names(v))) {
    out <- v[genes]; out <- as.numeric(out); names(out) <- genes; return(out)
  }
  if (length(v) == length(genes)) {
    out <- as.numeric(v); names(out) <- genes; return(out)
  }
  rep(NA_real_, length(genes))
}

extract_gene_by_boot <- function(boot_store, param, field, genes = NULL) {
  if (is.null(genes))
    genes <- get_template_genes(boot_store, param, field)
  if (length(genes) == 0)
    stop("No template genes for param=", param, " field=", field)
  B <- length(boot_store)
  M <- matrix(NA_real_, nrow = length(genes), ncol = B,
               dimnames = list(genes, paste0("b", seq_len(B))))
  for (b in seq_len(B))
    M[, b] <- safe_extract_vec(boot_store[[b]], param, field, genes)
  as.matrix(M)
}

# Summarize Δz across bootstraps per gene × parameter pair
summarise_param <- function(boot_store, param, genes = NULL) {
  dz        <- extract_gene_by_boot(boot_store, param, "dz", genes)
  genes_use <- rownames(dz)
  se <- extract_gene_by_boot(boot_store, param, "se_dz",    genes_use)
  I2 <- extract_gene_by_boot(boot_store, param, "I2_beta",  genes_use)
  kb <- extract_gene_by_boot(boot_store, param, "k_beta",   genes_use)
  tibble(
    param = param, gene = genes_use,
    n_boot_finite_dz = rowSums(is.finite(dz)),
    dz               = rowMeans(dz, na.rm = TRUE),
    dz_lo            = apply(dz, 1, stats::quantile, probs = 0.025,
                              na.rm = TRUE, names = FALSE),
    dz_hi            = apply(dz, 1, stats::quantile, probs = 0.975,
                              na.rm = TRUE, names = FALSE),
    se_dz_mean       = rowMeans(se, na.rm = TRUE),
    I2_beta_median   = apply(I2, 1, median, na.rm = TRUE),
    k_beta_median    = apply(kb, 1, median, na.rm = TRUE),
    pct_kbeta_le1    = apply(kb, 1, function(x) mean(x <= 1, na.rm = TRUE)) * 100
  ) %>%
    mutate(abs_dz = abs(dz)) %>%
    arrange(desc(abs_dz))
}

all_param_summaries <- purrr::map(ephys_vars, ~ summarise_param(boot_store, .x))
names(all_param_summaries) <- ephys_vars
dz_all_tbl <- dplyr::bind_rows(all_param_summaries)

# Sign stability: fraction of bootstraps where Δz has the same sign as the median
get_dz_matrix <- function(boot_store, param) {
  dz_list <- lapply(boot_store, function(x) x[[param]]$dz)
  dz_mat  <- do.call(cbind, dz_list)
  rownames(dz_mat) <- names(dz_list[[1]])
  dz_mat
}

compute_sign_stability <- function(dz_mat) {
  tibble(
    gene = rownames(dz_mat),
    sign_stability = apply(dz_mat, 1, function(x) {
      mean(sign(x) == sign(median(x, na.rm = TRUE)), na.rm = TRUE)
    })
  )
}

sign_stability_tbl <- purrr::map(
  unique(dz_all_tbl$param),
  function(p) {
    get_dz_matrix(boot_store, p) %>%
      compute_sign_stability() %>%
      mutate(param = p)
  }
) %>% bind_rows()

# ═══════════════════════════════════════════════════════════════════════════════
# Section 7: Build and save summary table
# ═══════════════════════════════════════════════════════════════════════════════

# Mean z_scb and z_beta across bootstraps (back-transformed to r for reporting)
summarise_z_means <- function(boot_store, param) {
  z_scb  <- extract_gene_by_boot(boot_store, param, "z_scb")
  z_beta <- extract_gene_by_boot(boot_store, param, "z_beta")
  tibble(
    gene   = rownames(z_scb),
    param  = param,
    z_scb  = rowMeans(z_scb,  na.rm = TRUE),
    z_beta = rowMeans(z_beta, na.rm = TRUE),
    dz     = z_scb - z_beta
  )
}

z_means <- purrr::map(ephys_vars, ~ summarise_z_means(boot_store, .x)) %>%
  bind_rows() %>%
  left_join(tr2g, by = "gene") %>%
  mutate(SYMBOL = ifelse(is.na(gene_name), gene, gene_name))

# Primary output table: all gene × parameter pairs with Δz and correlation estimates
out_for_jasmine <- dz_all_tbl %>%
  dplyr::select(param, gene, n_boot_finite_dz, dz, dz_lo, dz_hi,
                se_dz_mean, I2_beta_median, k_beta_median,
                pct_kbeta_le1, abs_dz) %>%
  mutate(ci_width = dz_hi - dz_lo) %>%
  left_join(
    z_means %>% dplyr::select(param, gene, z_beta, z_scb, gene_name, SYMBOL),
    by = c("param", "gene")
  ) %>%
  mutate(r_beta = tanh(z_beta), r_scb = tanh(z_scb)) %>%
  left_join(sign_stability_tbl, by = c("param", "gene"))

write.csv(out_for_jasmine,
          "./for_manuscript/20260423_boot_meta_dz_all_params_gene_pairs.csv",
          row.names = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 8: Sign-flip current amplitude variables
# ═══════════════════════════════════════════════════════════════════════════════
# Inward current amplitudes are negative by electrophysiology convention.
# After flipping, more positive = more electrophysiologically active,
# consistent with the direction used for exocytosis and other parameters.
# Note: dz_lo and dz_hi are also swapped on flip to maintain correct CI order.

out_for_jasmine_flipped <- out_for_jasmine %>%
  mutate(
    flip   = param %in% flip_vars,
    r_scb  = ifelse(flip, -r_scb,  r_scb),
    r_beta = ifelse(flip, -r_beta, r_beta),
    z_scb  = ifelse(flip, -z_scb,  z_scb),
    z_beta = ifelse(flip, -z_beta, z_beta),
    dz     = ifelse(flip, -dz,     dz),
    dz_lo  = ifelse(flip, -dz_hi,  dz_lo),
    dz_hi  = ifelse(flip, -dz_lo,  dz_hi)
  ) %>%
  dplyr::select(-flip)

# ── Electrophysiology parameter pretty-name lookup ────────────────────────────
ephys_names <- read_excel("250404_patch_key.xlsx", .name_repair = "universal") %>%
  dplyr::select(Pretty.names, name_in_xq_sheets)

# ── CI width helper ───────────────────────────────────────────────────────────
summarise_diff_like <- function(boot_store, param) {
  dz <- extract_gene_by_boot(boot_store, param, "dz")
  tibble(
    gene = rownames(dz), param = param,
    mean_diff     = rowMeans(dz, na.rm = TRUE),
    lower         = apply(dz, 1, quantile, probs = 0.025, na.rm = TRUE),
    upper         = apply(dz, 1, quantile, probs = 0.975, na.rm = TRUE),
    ci_width      = upper - lower,
    n_boot_finite = rowSums(is.finite(dz))
  )
}

diff_like <- purrr::map(ephys_vars, ~ summarise_diff_like(boot_store, .x)) %>%
  bind_rows() %>%
  left_join(tr2g, by = "gene") %>%
  mutate(SYMBOL = ifelse(is.na(gene_name), gene, gene_name)) %>%
  left_join(sign_stability_tbl, by = c("gene", "param"))

# ═══════════════════════════════════════════════════════════════════════════════
# Section 9: Dumbbell/arrow plots of divergent and consistent gene–ephys pairs
# ═══════════════════════════════════════════════════════════════════════════════
# Arrow segments show the shift in Spearman r from β-cells (tail) to SCβ
# (arrowhead) for the most divergent (Fig 4) and most consistent gene–ephys pairs.
# Genes are filtered for: CI width ≤ 0.35, ≥180/200 finite bootstraps,
# sign stability ≥ 0.9, and exclusion of unannotated transcripts.

# Top divergent: ranked by |Δz|, filtered for stability
top_divergent <- out_for_jasmine_flipped %>%
  mutate(abs_dz = abs(dz)) %>%
  arrange(desc(abs_dz)) %>%
  slice_head(n = 25) %>%
  bind_rows(out_for_jasmine_flipped %>% dplyr::filter(gene == "eGFP")) %>%
  mutate(r_beta = tanh(z_beta), r_scb = tanh(z_scb)) %>%
  left_join(ephys_names, by = c("param" = "name_in_xq_sheets")) %>%
  mutate(Pretty.names = str_remove(Pretty.names, "Normalized ")) %>%
  inner_join(diff_like, by = c("gene", "SYMBOL", "gene_name", "param")) %>%
  dplyr::filter(!is.na(ci_width), ci_width <= 0.35,
                n_boot_finite >= 180, sign_stability >= 0.9) %>%
  mutate(Gene_label = case_when(
    SYMBOL != gene ~ paste0("<i>", SYMBOL, "</i>"),
    TRUE           ~ SYMBOL
  ))

# Top consistent: small |Δz|, both populations have meaningful correlation
top_consistent <- out_for_jasmine_flipped %>%
  mutate(abs_dz = abs(dz), r_beta = tanh(z_beta), r_scb = tanh(z_scb),
         abs_r_beta = abs(r_beta), abs_r_scb = abs(r_scb)) %>%
  dplyr::filter(!is.na(dz), abs_dz <= 0.05,
                abs_r_beta >= 0.15, abs_r_scb >= 0.15) %>%
  arrange(abs_dz, desc(pmax(abs_r_beta, abs_r_scb))) %>%
  slice_head(n = 25) %>%
  left_join(ephys_names, by = c("param" = "name_in_xq_sheets")) %>%
  mutate(Pretty.names = str_remove(Pretty.names, "Normalized ")) %>%
  inner_join(diff_like, by = c("gene", "SYMBOL", "gene_name", "param")) %>%
  dplyr::filter(!is.na(ci_width), ci_width <= 0.35, n_boot_finite >= 180) %>%
  mutate(Gene_label = case_when(
    SYMBOL != gene ~ paste0("<i>", SYMBOL, "</i>"),
    TRUE           ~ SYMBOL
  ))

# Shared color palette — one color per electrophysiology parameter
roma_colors <- scico(8, palette = "roma")
feature_colors <- c(
  "Ca²⁺ Influxed  (pC/pF)"     = roma_colors[1],
  "Cell Size (pF)"              = roma_colors[2],
  "Early  Ca²⁺ Current  (pA/pF)"  = roma_colors[3],
  "Early Exocytosis (fF/pF)"   = roma_colors[4],
  "Late Ca²⁺ Current  (pA/pF)" = roma_colors[5],
  "Na⁺ Current (pA/pF)"        = roma_colors[6],
  "Reversal Potential (mV)"    = roma_colors[7],
  "Total Exocytosis (fF/pF)"   = roma_colors[8]
)

# Divergent panel: arrow points from r_beta (tail) to r_scb (arrowhead)
pdiv <- ggplot(top_divergent,
               aes(y = Gene_label, x = r_beta, xend = r_scb,
                   group = interaction(SYMBOL, param),
                   color = Pretty.names)) +
  geom_segment(aes(yend = Gene_label), linewidth = 1,
               arrow = arrow(length = unit(0.15, "cm"))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_colour_manual(values = feature_colors, name = "Feature") +
  labs(x = "Spearman Correlation", y = "Transcript") +
  theme_classic() +
  theme(axis.text.y = element_markdown())

# Consistent panel
pcons <- ggplot(top_consistent,
                aes(y = Gene_label, x = r_beta, xend = r_scb,
                    color = Pretty.names)) +
  geom_segment(aes(yend = Gene_label), linewidth = 1,
               arrow = arrow(length = unit(0.15, "cm"))) +
  scale_colour_manual(values = feature_colors, name = "Feature") +
  scale_y_discrete(labels = top_consistent$Gene_label) +
  labs(x = "Spearman Correlation", y = "Transcript") +
  theme_classic() +
  theme(axis.text.y = element_markdown())

ggsave("./for_manuscript/fig4_A_consistent.pdf",
       pcons + theme(legend.position = "none"),
       width = 6, height = 7, units = "cm", device = cairo_pdf)
ggsave("./for_manuscript/fig4_B_divergent.pdf",
       pdiv + theme(legend.position = "none"),
       width = 6, height = 7, units = "cm", device = cairo_pdf)
legend_ab <- cowplot::get_legend(pdiv + theme(legend.position = "right"))
ggsave("./for_manuscript/fig4_AB_legend.pdf",
       cowplot::plot_grid(legend_ab),
       width = 4, height = 4, units = "cm", device = cairo_pdf)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 10: CHGB × Cell size scatter (example correlation panel)
# ═══════════════════════════════════════════════════════════════════════════════
# Illustrates the sign reversal for CHGB–cell size correlation between
# SCβ and primary β-cells, as an example of a strongly divergent gene–param pair.

cells_chgb <- rownames(cells.combined.patched.sub@meta.data)[
  !is.na(cells.combined.patched.sub@meta.data$CellSize_pF) &
    cells.combined.patched.sub@meta.data$celltype %in% c("beta", "SCβ")
]
r_vals_chgb <- out_for_jasmine %>%
  dplyr::filter(param == "CellSize_pF", SYMBOL == "CHGB") %>%
  dplyr::select(r_beta, r_scb)
annot_chgb <- data.frame(
  celltype_label = c("Primary β-cell", "SCβ"),
  label = c(paste0("r = ", round(r_vals_chgb$r_beta, 2)),
            paste0("r = ", round(r_vals_chgb$r_scb,  2))),
  x = c(Inf, Inf), y = c(Inf, Inf)
)
plot_chgb <- cells.combined.patched.sub@meta.data[cells_chgb, ] %>%
  mutate(
    CHGB = as.numeric(cells.combined.patched.sub[["RNA"]]$data[
      "ENSG00000089199", cells_chgb]),
    celltype_label = ifelse(celltype == "beta", "Primary β-cell", "SCβ")
  )
p_chgb <- ggplot(plot_chgb, aes(x = CellSize_pF, y = CHGB)) +
  geom_point(aes(color = celltype_label), alpha = 0.5, size = 1.2) +
  geom_smooth(aes(color = celltype_label), method = "lm", se = TRUE,
              linewidth = 0.8) +
  geom_text(data = annot_chgb, aes(x = x, y = y, label = label),
            hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic") +
  scale_color_manual(values = c("Primary β-cell" = "#053059",
                                "SCβ"            = "#A6D278"),
                     guide = "none") +
  facet_wrap(~ celltype_label, scales = "free_x") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9)) +
  labs(x = "Cell size (pF)",
       y = expression(italic(CHGB) ~ "expression (log norm)"))
ggsave("./for_manuscript/fig4_schematic_CHGB_cellsize.pdf", p_chgb,
       device = cairo_pdf, width = 8, height = 6.5, units = "cm")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 11: Exocytosis × gene expression heatmap
# ═══════════════════════════════════════════════════════════════════════════════
# Top 10 positive and negative Δz genes for total exocytosis in SCβ cells,
# ordered by exocytosis value. Expression values are z-scored within each gene
# across SCβ cells (not SCTransform Pearson residuals — this is a within-gene
# scaling across the SCβ subset for visualization purposes).

exo_var   <- "NormalizedTotalCapacitance_fF.pF"
stage_var <- "celltype_stage"
meta      <- cells.combined.patched.sub@meta.data

top_exo_genes <- out_for_jasmine %>%
  dplyr::filter(param == exo_var, !is.na(dz), !is.na(SYMBOL),
                !str_detect(SYMBOL, "^RP[SL]|^RP11-|^AC[0-9]|^AL[0-9]|^ENSG"),
                sign_stability >= 0.9) %>%
  arrange(desc(dz)) %>%
  mutate(r_beta = round(r_beta, 3), r_scb = round(r_scb, 3),
         dz     = round(dz, 3)) %>%
  { bind_rows(slice_head(., n = 10), slice_tail(., n = 10)) } %>%
  dplyr::select(SYMBOL, gene, dz, r_beta, r_scb, I2_beta_median)

gene_order  <- top_exo_genes$gene
gene_labels <- top_exo_genes$SYMBOL
anno_df     <- top_exo_genes %>% arrange(match(gene, gene_order))

cells_scb_exo  <- rownames(meta)[
  meta$celltype == "SCβ" & !is.na(meta[[exo_var]]) & meta[[exo_var]] > 0
]
exo_vals       <- meta[cells_scb_exo, exo_var]
cell_order_exo <- cells_scb_exo[order(exo_vals)]
exo_ordered    <- exo_vals[order(exo_vals)]
stage_ordered  <- meta[cell_order_exo, stage_var]

expr_exo        <- as.matrix(
  cells.combined.patched.sub[["RNA"]]$data[gene_order, cell_order_exo]
)
# z-score within each gene across SCβ cells (for visualization only)
expr_scaled_exo <- t(scale(t(expr_exo)))
rownames(expr_scaled_exo) <- gene_labels

col_exo   <- colorRamp2(quantile(exo_ordered, c(0, 0.5, 1), na.rm = TRUE),
                         magma(3))
col_stage <- c("SCB D23–35" = "#A6D278", "SCB D45–54" = "#327D89")
r_lim     <- ceiling(max(abs(c(anno_df$r_scb, anno_df$r_beta)),
                          na.rm = TRUE) * 10) / 10
col_r     <- colorRamp2(c(-r_lim, 0, r_lim), scico(3, palette = "cork"))
col_expr  <- colorRamp2(c(-2, 0, 2),
                         c(viridis(3)[1], viridis(3)[2], viridis(3)[3]))
row_split <- factor(c(rep("Higher in SCβ", 10), rep("Higher in β", 10)),
                    levels = c("Higher in SCβ", "Higher in β"))

ht <- Heatmap(
  expr_scaled_exo,
  name = "Expression\n(z-score)", col = col_expr,
  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
  row_labels = gene_labels, row_names_side = "left",
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  row_split = row_split, row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_gap = unit(2, "mm"),
  top_annotation = HeatmapAnnotation(
    Exocytosis = exo_ordered, Stage = stage_ordered,
    col = list(Exocytosis = col_exo, Stage = col_stage),
    annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8)
  ),
  right_annotation = rowAnnotation(
    `r SCβ` = anno_df$r_scb, `r β` = anno_df$r_beta,
    col = list(`r SCβ` = col_r, `r β` = col_r),
    annotation_name_side = "top", annotation_name_gp = gpar(fontsize = 8),
    show_legend = c(`r SCβ` = TRUE, `r β` = FALSE),
    annotation_legend_param = list(`r SCβ` = list(title = "Spearman r"))
  ),
  heatmap_legend_param = list(direction = "vertical")
)

cairo_pdf("./for_manuscript/fig4_C_exo_heatmap.pdf", width = 5, height = 4)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 20), "mm"))
dev.off()

# ═══════════════════════════════════════════════════════════════════════════════
# Section 12: Na⁺ current × Pancreas Beta Cells heatmap
# ═══════════════════════════════════════════════════════════════════════════════

hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

pancreas_beta_genes <- hallmark_sets[["HALLMARK_PANCREAS_BETA_CELLS"]]
na_var <- "NormalizedPeakSodiumCurrentAmplitudeat.10mV_pA.pF"

pancreas_beta_filtered <- out_for_jasmine_flipped %>%
  dplyr::filter(param == na_var, SYMBOL %in% pancreas_beta_genes,
                n_boot_finite_dz > 0, sign_stability >= 0.85) %>%
  dplyr::select(SYMBOL, gene, dz, r_beta, r_scb,
                sign_stability, I2_beta_median) %>%
  arrange(dz)

gene_order_na  <- pancreas_beta_filtered$gene
gene_labels_na <- pancreas_beta_filtered$SYMBOL
anno_df_na     <- pancreas_beta_filtered

cells_na         <- rownames(meta)[
  meta$celltype == "SCβ" & !is.na(meta[[na_var]])
]
na_vals          <- -meta[cells_na, na_var]  # flip: more positive = more active
cell_order_na    <- cells_na[order(na_vals)]
na_ordered       <- na_vals[order(na_vals)]
stage_ordered_na <- meta[cell_order_na, stage_var]

expr_na        <- as.matrix(
  cells.combined.patched.sub[["RNA"]]$data[gene_order_na, cell_order_na]
)
expr_scaled_na <- t(scale(t(expr_na)))
rownames(expr_scaled_na) <- gene_labels_na

col_na     <- colorRamp2(quantile(na_ordered, c(0, 0.5, 1), na.rm = TRUE),
                          magma(3))
r_lim_na   <- ceiling(max(abs(c(anno_df_na$r_scb, anno_df_na$r_beta)),
                            na.rm = TRUE) * 10) / 10
col_r_na   <- colorRamp2(c(-r_lim_na, 0, r_lim_na), scico(3, palette = "cork"))
col_expr_na <- colorRamp2(c(-2, 0, 2),
                           c(viridis(3)[1], viridis(3)[2], viridis(3)[3]))
row_split_na <- factor(
  ifelse(anno_df_na$dz < 0, "Lower Na⁺ activity", "Higher Na⁺ activity"),
  levels = c("Lower Na⁺ activity", "Higher Na⁺ activity")
)

ht_na <- Heatmap(
  expr_scaled_na,
  name = "Expression\n(z-score)", col = col_expr_na,
  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
  row_labels = gene_labels_na, row_names_side = "left",
  row_names_gp = gpar(fontsize = 9, fontface = "italic"),
  row_split = row_split_na, row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_gap = unit(2, "mm"),
  top_annotation = HeatmapAnnotation(
    `Na⁺ current` = na_ordered, Stage = stage_ordered_na,
    col = list(`Na⁺ current` = col_na, Stage = col_stage),
    annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8)
  ),
  right_annotation = rowAnnotation(
    `r SCβ` = anno_df_na$r_scb, `r β` = anno_df_na$r_beta,
    col = list(`r SCβ` = col_r_na, `r β` = col_r_na),
    annotation_name_side = "top", annotation_name_gp = gpar(fontsize = 8),
    show_legend = c(`r SCβ` = TRUE, `r β` = FALSE),
    annotation_legend_param = list(`r SCβ` = list(title = "Spearman r"))
  ),
  heatmap_legend_param = list(direction = "vertical")
)

cairo_pdf("./for_manuscript/fig4_D_na_heatmap.pdf", width = 5, height = 5)
draw(ht_na, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 20), "mm"))
dev.off()

# ═══════════════════════════════════════════════════════════════════════════════
# Section 13: GSEA on Δz-ranked genes
# ═══════════════════════════════════════════════════════════════════════════════
# Hallmark GSEA using Δz as the ranking metric, run separately for each
# electrophysiology parameter. Genes ranked by mean Δz (flipped for current
# amplitudes so positive = more active in SCβ).

pathways_H <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$ensembl_gene, f = .$gs_name)

run_gsea_for_param <- function(param_name, min_boots = 150) {
  ranks <- out_for_jasmine_flipped %>%
    dplyr::filter(param == param_name,
                  n_boot_finite_dz >= min_boots,
                  !is.na(dz), !is.na(gene)) %>%
    arrange(desc(dz)) %>%
    dplyr::select(gene, dz) %>%
    deframe()

  fgsea(
    pathways    = pathways_H,
    stats       = ranks,
    nPermSimple = 10000
  ) %>%
    as_tibble() %>%
    mutate(param = param_name) %>%
    arrange(padj)
}

gsea_results <- purrr::map(ephys_vars, run_gsea_for_param) %>%
  bind_rows()

saveRDS(gsea_results,
        "./for_manuscript/gsea_hallmark_dz_all_params.rds")

cat("Significant GSEA results (padj < 0.05):\n")
print(gsea_results %>%
        dplyr::filter(padj < 0.05) %>%
        dplyr::select(param, pathway, NES, padj) %>%
        arrange(padj))

# ═══════════════════════════════════════════════════════════════════════════════
# Section 14: MYC target and OXPHOS heatmap
# ═══════════════════════════════════════════════════════════════════════════════

myc_genes   <- hallmark_sets[["HALLMARK_MYC_TARGETS_V1"]]

myc_filtered <- out_for_jasmine_flipped %>%
  dplyr::filter(param == exo_var, SYMBOL %in% myc_genes,
                n_boot_finite_dz > 0, sign_stability >= 0.9) %>%
  dplyr::select(SYMBOL, gene, dz, r_beta, r_scb,
                sign_stability, I2_beta_median) %>%
  arrange(dz)

expr_myc        <- as.matrix(
  cells.combined.patched.sub[["RNA"]]$data[myc_filtered$gene, cell_order_exo]
)
expr_scaled_myc <- t(scale(t(expr_myc)))
rownames(expr_scaled_myc) <- myc_filtered$SYMBOL

ht_myc <- Heatmap(
  expr_scaled_myc,
  name = "Expression\n(z-score)",
  col  = colorRamp2(c(-2, 0, 2),
                     c(viridis(3)[1], viridis(3)[2], viridis(3)[3])),
  cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE,
  row_labels = myc_filtered$SYMBOL, row_names_side = "left",
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  row_split = factor(
    ifelse(myc_filtered$dz < 0, "Higher in β", "Higher in SCβ"),
    levels = c("Higher in β", "Higher in SCβ")
  ),
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_gap = unit(2, "mm"),
  top_annotation = HeatmapAnnotation(
    Exocytosis = exo_ordered, Stage = stage_ordered,
    col = list(Exocytosis = col_exo, Stage = col_stage),
    annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8)
  ),
  right_annotation = rowAnnotation(
    `r SCβ` = myc_filtered$r_scb, `r β` = myc_filtered$r_beta,
    col = list(`r SCβ` = col_r, `r β` = col_r),
    annotation_name_side = "top", annotation_name_gp = gpar(fontsize = 8),
    show_legend = c(`r SCβ` = TRUE, `r β` = FALSE),
    annotation_legend_param = list(`r SCβ` = list(title = "Spearman r"))
  ),
  heatmap_legend_param = list(direction = "vertical")
)

cairo_pdf("./for_manuscript/suppfig_heatmap_myc_targets.pdf",
          width = 7, height = 6)
draw(ht_myc, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 20), "mm"))
dev.off()

# ═══════════════════════════════════════════════════════════════════════════════
# Section 15: MYC + OXPHOS density plot
# ═══════════════════════════════════════════════════════════════════════════════
# Visualizes the left-shift of MYC target and OXPHOS gene correlations
# with exocytosis in SCβ vs. primary β-cells — the key finding motivating
# the "negative enrichment" interpretation of the Hallmark GSEA results.

oxphos_genes <- hallmark_sets[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]
density_df   <- out_for_jasmine_flipped %>%
  dplyr::filter(param == exo_var, n_boot_finite_dz > 0, !is.na(SYMBOL)) %>%
  mutate(gene_set = case_when(
    SYMBOL %in% myc_genes    ~ "MYC targets",
    SYMBOL %in% oxphos_genes ~ "Oxidative phosphorylation",
    TRUE                     ~ "Other genes"
  )) %>%
  pivot_longer(cols = c(r_scb, r_beta),
               names_to = "population", values_to = "r") %>%
  mutate(
    population = ifelse(population == "r_scb", "SCβ", "Primary β-cell"),
    gene_set   = factor(gene_set, levels = c("MYC targets",
                                              "Oxidative phosphorylation",
                                              "Other genes"))
  )

p_density <- ggplot(density_df, aes(x = r, color = gene_set, fill = gene_set)) +
  geom_density(alpha = 0.2, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  facet_wrap(~ population) +
  scale_fill_manual(values = c("MYC targets"              = "#FDE725",
                                "Oxidative phosphorylation" = "#35B779",
                                "Other genes"               = "#440154"),
                    name = NULL) +
  scale_color_manual(values = c("MYC targets"              = "#FDE725",
                                 "Oxidative phosphorylation" = "#35B779",
                                 "Other genes"               = "#440154"),
                     name = NULL) +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  labs(x = "Spearman r", y = "Density")

ggsave("./for_manuscript/suppfig_density_myc_oxphos.pdf", p_density,
       device = cairo_pdf, width = 9, height = 5, units = "cm")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 16: Supplementary panel — top correlations per parameter
# ═══════════════════════════════════════════════════════════════════════════════
# Error bar plot showing top 5 positive and negative correlates of four key
# electrophysiology parameters, separately for SCβ and primary β-cells.
# Bold labels indicate high-confidence associations (sign stability ≥ 0.9,
# CI width ≤ 0.3).

params_supp <- c(
  "CellSize_pF",
  "NormalizedTotalCapacitance_fF.pF",
  "NormalizedEarlyPeakCalciumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedPeakSodiumCurrentAmplitudeat.10mV_pA.pF"
)
params_supp_labels <- c(
  "CellSize_pF"                                              = "Cell Size (pF)",
  "NormalizedTotalCapacitance_fF.pF"                         = "Total Exocytosis (fF/pF)",
  "NormalizedEarlyPeakCalciumCurrentAmplitudeat.10mV_pA.pF"  = "Early Ca²⁺ Current (pA/pF)",
  "NormalizedPeakSodiumCurrentAmplitudeat.10mV_pA.pF"        = "Na⁺ Current (pA/pF)"
)

scb_summary <- purrr::map(params_supp, function(p) {
  mat <- extract_gene_by_boot(boot_store, p, "z_scb")
  if (p %in% flip_vars) mat <- -mat
  tibble(param = p, gene = rownames(mat),
         z_mean   = rowMeans(mat, na.rm = TRUE),
         z_lo     = apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE),
         z_hi     = apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE),
         n_finite = rowSums(is.finite(mat))) %>%
    mutate(r_mean = tanh(z_mean), r_lo = tanh(z_lo),
           r_hi = tanh(z_hi), population = "SCβ")
}) %>% bind_rows()

beta_summary <- purrr::map(params_supp, function(p) {
  mat <- extract_gene_by_boot(boot_store, p, "z_beta")
  if (p %in% flip_vars) mat <- -mat
  tibble(param = p, gene = rownames(mat),
         z_mean   = rowMeans(mat, na.rm = TRUE),
         z_lo     = apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE),
         z_hi     = apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE),
         n_finite = rowSums(is.finite(mat))) %>%
    mutate(r_mean = tanh(z_mean), r_lo = tanh(z_lo),
           r_hi = tanh(z_hi), population = "Primary β-cell")
}) %>% bind_rows()

combined_summary <- bind_rows(scb_summary, beta_summary) %>%
  left_join(tr2g, by = "gene") %>%
  mutate(SYMBOL = ifelse(is.na(gene_name), gene, gene_name)) %>%
  dplyr::filter(
    # Exclude unannotated transcripts and low-confidence estimates
    !str_detect(SYMBOL, "^RP[SL]|^RP11-|^AC[0-9]|^AL[0-9]|^ENSG"),
    n_finite >= 150
  ) %>%
  left_join(sign_stability_tbl,
            by = c("param", "gene")) %>%
  left_join(diff_like %>%
              dplyr::select(param, gene, ci_width) %>% distinct(),
            by = c("param", "gene")) %>%
  mutate(high_confidence = sign_stability >= 0.9 &
           ci_width <= 0.3 & !is.na(ci_width))

top_genes_supp <- combined_summary %>%
  mutate(param_label = factor(params_supp_labels[param],
                               levels = params_supp_labels),
         population  = factor(population,
                               levels = c("Primary β-cell", "SCβ"))) %>%
  group_by(param, population) %>%
  arrange(desc(r_mean)) %>% mutate(rank_pos = row_number()) %>%
  arrange(r_mean)        %>% mutate(rank_neg = row_number()) %>%
  dplyr::filter(rank_pos <= 5 | rank_neg <= 5) %>%
  ungroup()

p_supp_c <- ggplot(top_genes_supp,
                   aes(x = r_mean,
                       y = reorder_within(SYMBOL, r_mean,
                                           interaction(param_label, population)),
                       color = population)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = r_lo, xmax = r_hi), height = 0.2, linewidth = 0.4) +
  geom_point(size = 2) +
  geom_text(data = dplyr::filter(top_genes_supp,  high_confidence),
            aes(label = SYMBOL, x = r_lo - 0.015),
            hjust = 1, size = 2.2, fontface = "bold.italic",
            show.legend = FALSE) +
  geom_text(data = dplyr::filter(top_genes_supp, !high_confidence),
            aes(label = SYMBOL, x = r_lo - 0.015),
            hjust = 1, size = 2.2, fontface = "italic",
            show.legend = FALSE) +
  scale_y_reordered(labels = NULL) +
  scale_color_manual(values = c("Primary β-cell" = "#053059",
                                "SCβ"            = "#A6D278"),
                     name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.35, 0.1))) +
  facet_wrap(population ~ param_label, scales = "free", nrow = 2) +
  theme_classic(base_size = 9) +
  theme(strip.background  = element_blank(),
        strip.text        = element_text(face = "bold", size = 8),
        axis.text.y       = element_blank(),
        axis.ticks.y      = element_blank(),
        axis.title.y      = element_blank(),
        legend.position   = "bottom",
        panel.spacing.x   = unit(1.2, "lines"),
        panel.spacing.y   = unit(0.8, "lines")) +
  labs(x = "Spearman r (bootstrapped mean ± 95% CI)")

ggsave("./for_manuscript/suppfig_panel_c_divergent.pdf", p_supp_c,
       width = 16, height = 10, units = "cm", device = cairo_pdf)

message("All figures saved successfully.")
