# ==============================================================================
# 04_differential_expression.R
# ==============================================================================
# Bootstrapped pseudobulk differential expression between SCβ-cells and primary
# β-cells across all patch-seq studies, with leave-one-study-out (LOSO)
# stability analysis and triage tier assignment.
#
# Corresponds to Methods: "Differential Expression Analysis"
#   "Differential expression was assessed using a bootstrapped pseudobulk
#    approach. Donor-level pseudobulk counts were aggregated per cell type
#    and study using AggregateExpression(), and glmGamPoi was used to model
#    count data with deconvolution-based size factors. 100 bootstrap iterations
#    were run by resampling pseudobulk samples with replacement within each
#    cell type, and the mean signed F-statistic across iterations was used as
#    the ranking metric for GSEA. Stability was assessed by leave-one-study-out
#    (LOSO) resampling across β-cell studies..."
#
# Input:
#   20260113_cells_plus_filtered_mi_SCENIC.rds
#     - cells.plus: full integrated Seurat object with AUC assay
#     - From scripts/02_integrate_patchseq_datasets.R + pySCENIC pipeline
#
# Output:
#   bootstrap_de_list.rds        - Raw DE results per bootstrap iteration
#   bootstrap_gene_summary.rds   - Summarized bootstrap DE (mean logFC, sign stability)
#   bootstrap_fit_diagnostics.rds - QC diagnostics per iteration
#   fgsea_bootstrap_H.rds        - Hallmark GSEA on full bootstrap summary
#   fgsea_bootstrap_C2_REACTOME.rds - REACTOME GSEA
#   fgsea_bootstrap_C5_GOBP.rds  - GO:BP GSEA
#   loso_beta_study_results.rds  - Per-study LOSO results
#   loso_summary.rds             - LOSO stability metrics per gene
#   JM_only_beta_vs_scb_results.rds - JM_patchSeq-only sensitivity analysis
#   bootstrapped_DE_triage_table.csv - Final triage tier assignments per gene
#     (primary input for all downstream feature selection and annotation)
#
# Dependencies:
#   tidyverse, Seurat, BiocParallel, glmGamPoi, fgsea, msigdbr,
#   EnsDb.Hsapiens.v86, ensembldb
#
# Runtime note:
#   100 bootstrap iterations × glmGamPoi fits is compute-intensive.
#   MulticoreParam(workers = 7) is used for parallelization.
#   Adjust workers to match your system's available cores.
# ==============================================================================

suppressPackageStartupMessages({
  library("tidyverse")
  library("Seurat")
  library("BiocParallel")
  library("glmGamPoi")
  library("fgsea")
  library("msigdbr")
  library("EnsDb.Hsapiens.v86")
  library("ensembldb")
})

# ── Gene ID → symbol lookup table ────────────────────────────────────────────
# Used throughout to annotate DE results with HGNC symbols
tr2g <- tibble(
  gene      = keys(EnsDb.Hsapiens.v86, keytype = "GENEID"),
  gene_name = mapIds(EnsDb.Hsapiens.v86,
                     keys     = keys(EnsDb.Hsapiens.v86, keytype = "GENEID"),
                     column   = "SYMBOL",
                     keytype  = "GENEID",
                     multiVals = "first")
)

# ── Gene set databases (msigdbr, Ensembl IDs) ────────────────────────────────
# Using Ensembl IDs throughout since our count matrix uses Ensembl row names
pathways_H  <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  split(x = .$ensembl_gene, f = .$gs_name)

pathways_C2 <- msigdbr(species = "Homo sapiens", collection = "C2",
                        subcollection = "CP:REACTOME") %>%
  split(x = .$ensembl_gene, f = .$gs_name)

pathways_C5 <- msigdbr(species = "Homo sapiens", collection = "C5",
                        subcollection = "GO:BP") %>%
  split(x = .$ensembl_gene, f = .$gs_name)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 1: Build pseudobulk count matrix
# ═══════════════════════════════════════════════════════════════════════════════

cells.plus <- readRDS("20260113_cells_plus_filtered_mi_SCENIC.rds")

# Restrict to SCβ and primary β-cells; exclude S6D1 (not a primary comparison)
cells.combined     <- subset(cells.plus, celltype != "S6D1")
cells.combined.sub <- subset(cells.combined, celltype %in% c("beta", "SCβ"))

# Aggregate counts to pseudobulk level: one sample per Differentiation × Study × celltype
# Differentiation = donor batch for SCβ; donor ID for primary β-cells
bulk.expression <- AggregateExpression(
  cells.combined.sub,
  return.seurat = TRUE,
  slot          = "counts",
  assays        = "RNA",
  group.by      = c("Differentiation", "Study", "celltype")
)

bulk.counts <- LayerData(bulk.expression, layer = "counts", assay = "RNA") %>%
  as.matrix()
coldata     <- bulk.expression@meta.data
coldata     <- coldata[colnames(bulk.counts), , drop = FALSE]

cat("Pseudobulk samples:\n")
print(table(coldata$celltype, coldata$Study))

# ═══════════════════════════════════════════════════════════════════════════════
# Section 2: Bootstrap helper functions
# ═══════════════════════════════════════════════════════════════════════════════

# ── Compact per-iteration diagnostics ────────────────────────────────────────
.compute_diag <- function(counts_i, meta_i, sampled_beta, sampled_scb, elapsed_sec) {
  lib <- colSums(counts_i)
  tibble(
    n_samples_total  = ncol(counts_i),
    n_beta           = length(sampled_beta),
    n_scb            = length(sampled_scb),
    n_unique_beta    = length(unique(sampled_beta)),
    n_unique_scb     = length(unique(sampled_scb)),
    max_mult_beta    = max(table(sampled_beta)),
    max_mult_scb     = max(table(sampled_scb)),
    lib_min          = min(lib),
    lib_median       = median(lib),
    lib_max          = max(lib),
    lib_median_beta  = median(lib[sampled_beta], na.rm = TRUE),
    lib_median_scb   = median(lib[sampled_scb],  na.rm = TRUE),
    elapsed_sec      = elapsed_sec
  )
}

# ── One bootstrap iteration ───────────────────────────────────────────────────
# Resamples pseudobulk samples within each cell type with replacement,
# fits glmGamPoi, and extracts DE results for SCβ vs β.
# Returns NULL for DE on failure (logged in diagnostics).
.run_one_boot <- function(iter, bulk.counts, coldata,
                          samples_beta, samples_scb, seed = 123) {
  set.seed(seed + iter)
  t0 <- proc.time()[["elapsed"]]

  sampled_beta <- sample(samples_beta, length(samples_beta), replace = TRUE)
  sampled_scb  <- sample(samples_scb,  length(samples_scb),  replace = TRUE)
  sampled_all  <- c(sampled_beta, sampled_scb)

  counts_i <- bulk.counts[, sampled_all, drop = FALSE]
  meta_i   <- coldata[sampled_all, , drop = FALSE]
  colnames(counts_i) <- rownames(meta_i)

  out <- tryCatch({
    fit_i <- glm_gp(
      data           = counts_i,
      col_data       = meta_i,
      design         = ~ celltype,
      reference_level = "beta",
      size_factors   = "deconvolution",
      on_disk        = FALSE,
      verbose        = FALSE
    )
    de_i <- test_de(fit_i,
                    contrast = cond(celltype = "SCβ") - cond(celltype = "beta")) %>%
      as_tibble() %>%
      dplyr::rename(gene = name) %>%
      mutate(iter = iter)

    t1     <- proc.time()[["elapsed"]]
    diag_i <- .compute_diag(counts_i, meta_i, sampled_beta, sampled_scb,
                             elapsed_sec = t1 - t0) %>%
      mutate(iter = iter, ok = TRUE, err = NA_character_)

    list(de = de_i, diag = diag_i)
  }, error = function(e) {
    t1     <- proc.time()[["elapsed"]]
    diag_i <- .compute_diag(counts_i, meta_i, sampled_beta, sampled_scb,
                             elapsed_sec = t1 - t0) %>%
      mutate(iter = iter, ok = FALSE, err = conditionMessage(e))
    list(de = NULL, diag = diag_i)
  })
  out
}

# ── Summarize bootstrap DE results ───────────────────────────────────────────
# Computes per-gene mean logFC, sign stability (fraction of bootstraps where
# logFC > 0), and robust score (effect size penalized by instability).
# signed_f = mean_stat * sign(frac_pos - 0.5): effect size with consistent direction.
.summarize_boot_de <- function(de_list) {
  de_all <- bind_rows(de_list)
  if (nrow(de_all) == 0) stop("No successful DE results to summarize.")
  if (!"lfc" %in% colnames(de_all))
    stop("Expected 'lfc' column from glmGamPoi::test_de().")

  de_all %>%
    group_by(gene) %>%
    summarize(
      n_iter       = n(),
      mean_logFC   = mean(lfc,         na.rm = TRUE),
      sd_logFC     = sd(lfc,           na.rm = TRUE),
      frac_pos     = mean(lfc > 0,     na.rm = TRUE),
      mean_stat    = if ("f_statistic" %in% colnames(de_all))
                       mean(f_statistic, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      robust_score  = mean_logFC / (1 + sd_logFC),
      sign_stability = pmax(frac_pos, 1 - frac_pos),
      robust_score2  = robust_score * sign_stability
    )
}

# ── Run fgsea from a summary table ───────────────────────────────────────────
.run_fgsea <- function(summary_df, pathways, score_col = "signed_f",
                       nPermSimple = 10000) {
  stopifnot(score_col %in% colnames(summary_df))
  ranks <- summary_df %>%
    dplyr::select(gene, !!sym(score_col)) %>%
    drop_na() %>%
    arrange(desc(!!sym(score_col))) %>%
    deframe()
  fgsea(pathways = pathways, stats = ranks, nPermSimple = nPermSimple)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Section 3: Main bootstrap (100 iterations, parallelized)
# ═══════════════════════════════════════════════════════════════════════════════

n_boot <- 100
seed   <- 123

samples_beta <- rownames(coldata)[coldata$celltype == "beta"]
samples_scb  <- rownames(coldata)[coldata$celltype == "SCβ"]
stopifnot(length(samples_beta) > 0, length(samples_scb) > 0)

# MulticoreParam uses fork-based parallelism (Linux/Mac)
# Switch to SnowParam on Windows: SnowParam(workers = 7, type = "SOCK")
param <- MulticoreParam(workers = 7, RNGseed = seed)

cat("Running", n_boot, "bootstrap iterations...\n")
boot_out <- bplapply(seq_len(n_boot), function(i) {
  .run_one_boot(iter = i, bulk.counts = bulk.counts, coldata = coldata,
                samples_beta = samples_beta, samples_scb = samples_scb,
                seed = seed)
}, BPPARAM = param)

diag_df <- bind_rows(map(boot_out, "diag"))
de_list  <- map(boot_out, "de") %>% discard(is.null)

cat("Successful iterations:", length(de_list), "/", n_boot, "\n")
saveRDS(diag_df, file = "bootstrap_fit_diagnostics.rds")
saveRDS(de_list, file = "bootstrap_de_list.rds")

# ── Summarize and run GSEA ────────────────────────────────────────────────────
summary_boot <- .summarize_boot_de(de_list) %>%
  mutate(signed_f = mean_stat * sign(frac_pos - 0.5)) %>%
  left_join(tr2g, by = "gene")

saveRDS(summary_boot, file = "bootstrap_gene_summary.rds")

fgsea_H  <- .run_fgsea(summary_boot, pathways_H,  score_col = "signed_f")
fgsea_C2 <- .run_fgsea(summary_boot, pathways_C2, score_col = "signed_f")
fgsea_C5 <- .run_fgsea(summary_boot, pathways_C5, score_col = "signed_f")

saveRDS(fgsea_H,  file = "fgsea_bootstrap_H.rds")
saveRDS(fgsea_C2, file = "fgsea_bootstrap_C2_REACTOME.rds")
saveRDS(fgsea_C5, file = "fgsea_bootstrap_C5_GOBP.rds")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 4: LOSO stability analysis
# ═══════════════════════════════════════════════════════════════════════════════
# For each β-cell study, rerun bootstrapped DE with that study's β-cells
# excluded. This tests whether the DE signal is driven by any single study.
# loso_rel_iqr = IQR of signed_f across LOSO runs / |signed_f_full|:
# high values indicate genes whose effect estimates are sensitive to which
# β-cell studies are included.

beta_studies <- sort(unique(coldata$Study[coldata$celltype == "beta"]))
cat("Running LOSO across", length(beta_studies), "beta studies:",
    paste(beta_studies, collapse = ", "), "\n")

n_boot_loso <- 50

loso_results <- vector("list", length(beta_studies))
names(loso_results) <- beta_studies

for (st in beta_studies) {
  message("LOSO: excluding beta study = ", st)

  samples_beta_loso <- rownames(coldata)[
    coldata$celltype == "beta" & coldata$Study != st
  ]

  if (length(samples_beta_loso) < 2) {
    warning("Skipping LOSO for ", st, ": too few beta samples remain.")
    next
  }

  out_loso <- bplapply(seq_len(n_boot_loso), function(i) {
    .run_one_boot(
      iter         = i,
      bulk.counts  = bulk.counts,
      coldata      = coldata,
      samples_beta = samples_beta_loso,
      samples_scb  = samples_scb,
      seed         = seed + 10000 + match(st, beta_studies) * 1000
    )
  }, BPPARAM = param)

  de_loso  <- map(out_loso, "de") %>% discard(is.null)
  sum_loso <- .summarize_boot_de(de_loso) %>%
    mutate(signed_f = mean_stat * sign(frac_pos - 0.5)) %>%
    left_join(tr2g, by = "gene")

  loso_results[[st]] <- list(
    diagnostics  = bind_rows(map(out_loso, "diag")) %>%
                     mutate(loso_excluded_beta_study = st),
    gene_summary = sum_loso,
    fgsea_H      = .run_fgsea(sum_loso, pathways_H,  score_col = "signed_f"),
    fgsea_C2     = .run_fgsea(sum_loso, pathways_C2, score_col = "signed_f"),
    fgsea_C5     = .run_fgsea(sum_loso, pathways_C5, score_col = "signed_f")
  )
}

saveRDS(loso_results, file = "loso_beta_study_results.rds")

# ── Compute LOSO stability metrics ───────────────────────────────────────────
# For each gene: how much does signed_f vary across LOSO runs?
# loso_rel_iqr is used downstream to flag unstable genes for triage.
full_sf <- summary_boot %>%
  mutate(signed_f = mean_stat * sign(frac_pos - 0.5)) %>%
  dplyr::select(gene, gene_name, signed_f_full = signed_f)

loso_long <- imap(loso_results, function(x, excluded_study) {
  gs <- x$gene_summary
  if (!"signed_f" %in% names(gs))
    gs <- gs %>% mutate(signed_f = mean_stat * sign(frac_pos - 0.5))
  gs %>%
    dplyr::select(gene, signed_f) %>%
    mutate(excluded_beta_study = excluded_study)
}) %>% bind_rows()

loso_summary <- loso_long %>%
  group_by(gene) %>%
  summarize(
    loso_n             = n(),
    loso_mean_signed_f = mean(signed_f, na.rm = TRUE),
    loso_sd_signed_f   = sd(signed_f,   na.rm = TRUE),
    loso_iqr_signed_f  = IQR(signed_f,  na.rm = TRUE),
    loso_min_signed_f  = min(signed_f,  na.rm = TRUE),
    loso_max_signed_f  = max(signed_f,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    loso_long %>%
      left_join(full_sf, by = "gene") %>%
      group_by(gene) %>%
      summarize(
        loso_sign_flip_rate = mean(sign(signed_f) != sign(signed_f_full),
                                   na.rm = TRUE),
        .groups = "drop"
      ),
    by = "gene"
  ) %>%
  left_join(full_sf, by = "gene") %>%
  mutate(
    loso_rel_sd  = loso_sd_signed_f  / (abs(signed_f_full) + 1e-6),
    loso_rel_iqr = loso_iqr_signed_f / (abs(signed_f_full) + 1e-6)
  )

saveRDS(loso_summary, file = "loso_summary.rds")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 5: JM_patchSeq-only sensitivity analysis
# ═══════════════════════════════════════════════════════════════════════════════
# Reruns bootstrapped DE using only JM_patchSeq β-cells as the comparator.
# This avoids cross-study confounding entirely and is used to assess whether
# the full cross-study result is consistent with the within-study comparison.
# Results are used as a "JM support" flag in the triage table.

samples_beta_JM <- rownames(coldata)[
  coldata$celltype == "beta" & coldata$Study == "JM_patchSeq"
]
cat("JM_patchSeq beta samples:", length(samples_beta_JM), "\n")

if (length(samples_beta_JM) >= 2) {
  out_JM <- bplapply(seq_len(100), function(i) {
    .run_one_boot(
      iter         = i,
      bulk.counts  = bulk.counts,
      coldata      = coldata,
      samples_beta = samples_beta_JM,
      samples_scb  = samples_scb,
      seed         = seed + 20000
    )
  }, BPPARAM = param)

  de_JM  <- map(out_JM, "de") %>% discard(is.null)
  sum_JM <- .summarize_boot_de(de_JM) %>%
    mutate(signed_f = mean_stat * sign(frac_pos - 0.5)) %>%
    left_join(tr2g, by = "gene")

  JM_results <- list(
    diagnostics  = bind_rows(map(out_JM, "diag")),
    gene_summary = sum_JM,
    fgsea_H      = .run_fgsea(sum_JM, pathways_H,  score_col = "signed_f"),
    fgsea_C2     = .run_fgsea(sum_JM, pathways_C2, score_col = "signed_f"),
    fgsea_C5     = .run_fgsea(sum_JM, pathways_C5, score_col = "signed_f")
  )
  saveRDS(JM_results, file = "JM_only_beta_vs_scb_results.rds")
} else {
  warning("Not enough JM_patchSeq beta samples to run JM-only sensitivity.")
  JM_results <- NULL
}

# ═══════════════════════════════════════════════════════════════════════════════
# Section 6: Triage tier assignment
# ═══════════════════════════════════════════════════════════════════════════════
# Assigns each gene to a triage tier based on:
# - Effect size (strong = top 5% of |signed_f_full|)
# - LOSO stability (unstable = top 10% of loso_rel_iqr OR sign flip rate > 0)
# - JM support (strong effect in JM-only analysis, consistent direction)
#
# Tier definitions:
#   Tier1A: strong, stable, JM-supported      → primary candidates
#   Tier1B: strong, stable, JM data missing   → include with caveat
#   Tier2:  strong but unstable or JM-sensitive → interpret cautiously
#   Tier3:  weak but stable                   → low-confidence
#   Tier4:  weak and unstable                 → drop

summary_boot2 <- summary_boot %>%
  mutate(
    signed_f_full       = mean_stat * sign(frac_pos - 0.5),
    sign_stability_full = pmax(frac_pos, 1 - frac_pos)
  )

# Threshold definitions (quantile-based for robustness)
strong_cut   <- quantile(abs(summary_boot2$signed_f_full), 0.95, na.rm = TRUE)
weak_cut     <- quantile(abs(summary_boot2$signed_f_full), 0.50, na.rm = TRUE)
unstable_cut <- quantile(loso_summary$loso_rel_iqr,        0.90, na.rm = TRUE)

triage <- summary_boot2 %>%
  dplyr::select(gene, gene_name, mean_logFC, sd_logFC, frac_pos, mean_stat,
                signed_f_full, sign_stability_full) %>%
  left_join(
    loso_summary %>%
      dplyr::select(gene, loso_n, loso_mean_signed_f, loso_sd_signed_f,
                    loso_iqr_signed_f, loso_min_signed_f, loso_max_signed_f,
                    loso_sign_flip_rate, loso_rel_sd, loso_rel_iqr),
    by = "gene"
  ) %>%
  mutate(
    abs_signed_f   = abs(signed_f_full),
    abs_logFC      = abs(mean_logFC),
    logFC_clipped  = abs(mean_logFC) >= 9.999,
    # Flag unannotated or non-coding genes for awareness (not hard filtering)
    name_flag = case_when(
      is.na(gene_name) ~ "no_symbol",
      str_detect(gene_name,
                 "^RP11\\b|^CTD\\-|^AC0\\d|^AL\\d|^LINC\\b|^MIR\\b|^SNOR\\b|^RNU\\b|^RN7S\\b") ~
        "lncRNA_or_locus",
      str_detect(gene_name, "^ERCC\\b|^eGFP\\b") ~ "spikein_or_reporter",
      TRUE ~ "ok"
    ),
    is_strong   = abs_signed_f >= strong_cut,
    is_weak     = abs_signed_f <= weak_cut,
    is_unstable = (loso_rel_iqr >= unstable_cut) | (loso_sign_flip_rate > 0),
    triage_tier = case_when(
      is_strong & !is_unstable ~ "Tier1_strong_stable",
      is_strong &  is_unstable ~ "Tier2_strong_unstable",
      !is_strong & !is_unstable ~ "Tier3_weak_stable",
      TRUE ~ "Tier4_weak_unstable"
    )
  )

# ── Add JM-only sensitivity columns ──────────────────────────────────────────
if (!is.null(JM_results)) {
  jm_weak_cut      <- quantile(abs(JM_results$gene_summary$mean_logFC),
                                0.50, na.rm = TRUE)
  full_strong_cut2 <- quantile(abs(triage$signed_f_full), 0.95, na.rm = TRUE)

  sum_JM_triage <- JM_results$gene_summary %>%
    mutate(
      signed_f_JM       = mean_stat * sign(frac_pos - 0.5),
      sign_stability_JM = pmax(frac_pos, 1 - frac_pos)
    ) %>%
    dplyr::select(gene, signed_f_JM, mean_logFC_JM = mean_logFC,
                  frac_pos_JM = frac_pos, sign_stability_JM)

  triage2 <- triage %>%
    left_join(sum_JM_triage, by = "gene") %>%
    mutate(
      dir_full     = sign(signed_f_full),
      dir_JM       = sign(signed_f_JM),
      jm_dir_match = case_when(
        is.na(signed_f_JM)            ~ NA,
        dir_full == 0 | dir_JM == 0  ~ NA,
        TRUE                          ~ (dir_full == dir_JM)
      ),
      jm_supported = !is.na(signed_f_JM) &
        (abs(signed_f_JM) >= quantile(abs(signed_f_JM), 0.90, na.rm = TRUE)) &
        (sign_stability_JM >= 0.9) &
        (jm_dir_match %in% TRUE),
      strong_full        = abs(signed_f_full) >= full_strong_cut2,
      weak_in_JM         = !is.na(signed_f_JM) & abs(signed_f_JM) <= jm_weak_cut,
      jm_sensitivity_flag = strong_full & weak_in_JM,
      triage_tier2 = case_when(
        strong_full & !is_unstable & jm_supported  ~ "Tier1A_strong_stable_JM_supported",
        strong_full & !is_unstable & is.na(signed_f_JM) ~ "Tier1B_strong_stable_JM_missing",
        strong_full & (is_unstable | jm_sensitivity_flag |
                         (jm_dir_match %in% FALSE))  ~ "Tier2_strong_sensitive",
        !strong_full & !is_unstable                 ~ "Tier3_weak_stable",
        TRUE                                        ~ "Tier4_weak_unstable"
      )
    )
} else {
  # If JM-only analysis couldn't run, all strong stable genes get Tier1B
  triage2 <- triage %>%
    mutate(triage_tier2 = case_when(
      is_strong & !is_unstable  ~ "Tier1B_strong_stable_JM_missing",
      is_strong &  is_unstable  ~ "Tier2_strong_sensitive",
      !is_strong & !is_unstable ~ "Tier3_weak_stable",
      TRUE                      ~ "Tier4_weak_unstable"
    ))
}

write_csv(triage2, file = "bootstrapped_DE_triage_table.csv")

cat("\nTriage tier distribution:\n")
print(table(triage2$triage_tier2))

cat("\nDone. Key outputs:\n")
cat("  bootstrapped_DE_triage_table.csv  - primary downstream input\n")
cat("  loso_summary.rds                  - LOSO stability metrics\n")
cat("  fgsea_bootstrap_H.rds             - Hallmark GSEA results\n")
