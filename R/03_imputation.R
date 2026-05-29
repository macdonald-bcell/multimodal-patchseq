# ==============================================================================
# 03_imputation.R
# ==============================================================================
# Multiple imputation of missing electrophysiology values using MICE, with
# custom imputation methods tailored to the measurement characteristics of each
# parameter type. Produces 50 completed datasets propagated through all
# downstream analyses.
#
# Corresponds to Methods: "Multiple Imputations for Missing Electrophysiology
# Data"
#   "Missing electrophysiology values were handled using multiple imputation by
#    chained equations (MICE)... m = 50 imputed datasets. Imputation models were
#    tailored to the measurement characteristics of each parameter type.
#    Exocytosis parameters were imputed using a two-part model: logistic
#    regression for presence or absence... Inward current parameters were
#    imputed as magnitudes using a censored Gaussian model via Tobit
#    regression..."
#
# Input:
#   20241017_JM_islets_seurat_object_ensg_decon.RData
#     - JM_seurobj_ensg: JM study cells with metadata
#   ../human_islets_scRNAseq/20250708_patchseq_metadata.xlsx
#     - Patch-clamp metadata for all external studies (XQ, AG, HPAP, etc.)
#   ../human_islets_scRNAseq/20250708_pclamp_all.RData
#     - pclamp_all: Seurat object for all patched cells (external studies)
#
# Output:
#   20260120_imp.rds
#     - MICE mids object (50 imputed datasets, internal format)
#   20260120_imp_bundle.rds
#     - List: list(imp = imp, completed_list = completed_list)
#     - completed_list: 50 data frames with reconstructed signed values and
#       recomputed normalized parameters — primary input for all downstream
#       analyses
#   20260120_imputation_settings_bundle.rds
#     - All imputation settings, method assignments, predictor matrix, etc.
#     - For reproducibility: recreating the exact imputation from scratch
#   20251210_patchSeq_integrated_patched-only_decon_mi.RData
#     - pclamp_patched_all: Seurat object with Rubin-pooled imputed means and
#       SEs added to metadata as mi_mean_* and mi_se_* columns
#
# Dependencies:
#   tidyverse, openxlsx, Seurat, mice, AER, truncnorm, survival, matrixStats
#
# Runtime note:
#   MICE with m=50, maxit=20 is computationally intensive (~hours).
#   Save the imp object immediately after running and reload for all
#   downstream work.
# ==============================================================================

library("tidyverse")
library("openxlsx")
library("Seurat")
library("mice")
library("AER")
library("truncnorm")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 1: Load and combine metadata from all patch-seq studies
# ═══════════════════════════════════════════════════════════════════════════════
# The imputation is run on all studies combined so that cross-study information
# (e.g. shared cell types, glucose concentrations) can be used as predictors.
# Signed current conventions and column name harmonization are applied here.

metadata_patched <- openxlsx::read.xlsx(
  "../human_islets_scRNAseq/20250708_patchseq_metadata.xlsx"
)

load("20241017_JM_islets_seurat_object_ensg_decon.RData")
# JM_seurobj_ensg loaded

# ── Prepare JM metadata ───────────────────────────────────────────────────────

colnames(JM_seurobj_ensg) -> JM_seurobj_ensg$CellID

# Extract Differentiation batch from free-text ID field
diff <- str_extract(
  JM_seurobj_ensg$ID.From.Fclynne,
  "(?i)(\\d*\\s*diff(\\s*\\d+(/\\d+)?|\\s*\\d{6}|\\s*start\\s*date:\\d{8})|differentiation\\s*\\d+|\\d+\\s*differentiation)"
)
JM_seurobj_ensg$Differentiation <- diff
JM_seurobj_ensg <- subset(x = JM_seurobj_ensg, Group != "N/A" & !is.na(Group))

# Fill missing Differentiation for primary islet cells using donor/FL_ID
donors <- JM_seurobj_ensg$ID.From.Fclynne[
  is.na(JM_seurobj_ensg$Donor) & JM_seurobj_ensg$Group == "Human Primary"
]
names(donors) <- colnames(JM_seurobj_ensg)[
  is.na(JM_seurobj_ensg$Donor) & JM_seurobj_ensg$Group == "Human Primary"
]
JM_seurobj_ensg$Differentiation[names(donors)] <- donors
JM_seurobj_ensg$Donor[names(donors)]           <- donors

donors2 <- JM_seurobj_ensg$Donor[
  is.na(JM_seurobj_ensg$Differentiation) & JM_seurobj_ensg$Group == "Human Primary"
]
names(donors2) <- colnames(JM_seurobj_ensg)[
  is.na(JM_seurobj_ensg$Differentiation) & JM_seurobj_ensg$Group == "Human Primary"
]
JM_seurobj_ensg$Differentiation[names(donors2)] <- donors2

# Clean metadata
JM_meta_clean <- JM_seurobj_ensg@meta.data %>%
  dplyr::select(where(~ any(!is.na(.))))
JM_seurobj_ensg@meta.data <- JM_meta_clean

# Fix colon in label (causes issues in some file formats)
JM_seurobj_ensg$Differentiation[
  JM_seurobj_ensg$Differentiation == "35 Diff start date:20210905"
] <- "35 Diff start date 20210905"

JM_seurobj_ensg$Study <- "JM_patchSeq"
JM_seurobj_ensg       <- subset(JM_seurobj_ensg, GroupBin != "IPSC")
JM_metadata           <- JM_seurobj_ensg@meta.data

# ── Harmonize column names across studies ─────────────────────────────────────
# Some columns were renamed between recording sessions; standardize here

harmonize_col <- function(meta, target, source) {
  meta[[target]] <- NA
  meta[[target]][is.na(meta[[target]])] <- meta[[source]][is.na(meta[[target]])]
  meta
}

JM_metadata <- harmonize_col(JM_metadata,
  "NormalizedLVA_.60mV_pA.pF", "NormalizedLVA__.60mV_pA.pF")
JM_metadata <- harmonize_col(JM_metadata,
  "EarlyPeakCalciumCurrentAmplitudeat.10mV_pA", "EarlyPeakCaCurrentAmplitudeat.10mV_pA")
JM_metadata <- harmonize_col(JM_metadata,
  "NormalizedEarlyPeakCalciumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedEarlyPeakCaCurrentAmplitudeat.10mV_pA.pF")
JM_metadata <- harmonize_col(JM_metadata,
  "LateCalciumCurrentAmplitudeat.10mV_pA", "LateCaCurrentAmplitudeat.10mV_pA")
JM_metadata <- harmonize_col(JM_metadata,
  "NormalizedLateCalciumCurrentAmplitudeat.10mV_pA.pF",
  "NormalizedLateCaCurrentAmplitudeat.10mV_pA.pF")

# Fill missing glucose with 5 mM (standard recording solution)
JM_metadata$Glucose_mM[is.na(JM_metadata$Glucose_mM)] <- 5
JM_metadata$celltype[JM_metadata$celltype == "Unknown"] <- "SCβ"
JM_metadata$Patcher         <- "JM"
JM_metadata$Donor           <- JM_metadata$Differentiation
JM_metadata$DiabetesStatus  <- "Non-Diabetic"

# ── Merge metadata ────────────────────────────────────────────────────────────

common_cols        <- intersect(colnames(metadata_patched), colnames(JM_metadata))
metadata_patched   <- metadata_patched[, common_cols]
JM_metadata        <- JM_metadata[, common_cols]

df <- bind_rows(metadata_patched, JM_metadata)
df <- df %>% dplyr::filter(!is.na(CellSize_pF))

cat("Total cells for imputation:", nrow(df), "\n")
cat("Studies:", paste(unique(df$Study), collapse = ", "), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Section 2: Sign correction and data cleaning
# ═══════════════════════════════════════════════════════════════════════════════
# Inward current amplitudes (calcium, sodium) should be negative by convention.
# Positive values are confirmed typos and are flipped. Extreme outlier values
# and biologically impossible zeros are set to NA.

# Columns that must be ≤ 0 by electrophysiology convention
force_negative_cols <- c(
  "LateCalciumCurrentAmplitudeat.10mV_pA",
  "LVA_.60mV_pA",
  "HVA_.20mV_pA",
  "PeakSodiumCurrentAmplitudeat.10mV_pA",
  "EarlyPeakCalciumCurrentAmplitudeat.10mV_pA",
  "CalciumIntegralatFirstDepolarization_pC"
) %>% intersect(names(df))

df <- df %>%
  mutate(
    # Flip confirmed typos (positive values for inward currents)
    across(
      all_of(force_negative_cols),
      ~ if_else(is.na(.x) | .x <= 0, .x, -abs(.x))
    )
  ) %>%
  # Remove cell with CellSize_pF == 0 (should have been NA)
  dplyr::filter(is.na(CellSize_pF) | CellSize_pF != 0) %>%
  mutate(
    # Extreme LVA/HVA values are measurement artifacts
    LVA_.60mV_pA = if_else(LVA_.60mV_pA < -10000, NA_real_, LVA_.60mV_pA),
    HVA_.20mV_pA = if_else(HVA_.20mV_pA < -10000, NA_real_, HVA_.20mV_pA),
    # Zero half-inactivation voltage is not biologically meaningful
    HalfInactivationofSodiumCurrent_mV =
      na_if(HalfInactivationofSodiumCurrent_mV, 0)
  )

# ── Define parameter groups ───────────────────────────────────────────────────
# These groups determine imputation strategy (see Section 4)

exocytosis_cols <- c(
  "TotalCapacitance_fF",
  "FirstDepolarizationCapacitance_fF",
  "LateDepolarizationCapacitance"
) %>% intersect(names(df))

calcium_cols <- c(
  "EarlyPeakCalciumCurrentAmplitudeat.10mV_pA",
  "LateCalciumCurrentAmplitudeat.10mV_pA",
  "LVA_.60mV_pA",
  "HVA_.20mV_pA",
  "CalciumIntegralatFirstDepolarization_pC"
) %>% intersect(names(df))

sodium_cols <- c(
  "PeakSodiumCurrentAmplitudeat.10mV_pA"
) %>% intersect(names(df))

# Special parameters with non-standard imputation
v_sodium_peak_ord <- intersect("VoltageforSodiumPeakCurrent_mV", names(df))   # ordered factor
v_reversal_signed <- intersect("ReversalPotentialbyramp_mV",     names(df))   # continuous ±
v_half_inact_mv   <- intersect("HalfInactivationofSodiumCurrent_mV", names(df)) # continuous ±

# ═══════════════════════════════════════════════════════════════════════════════
# Section 3: Custom MICE imputation methods
# ═══════════════════════════════════════════════════════════════════════════════

# ── Helper: tiny epsilon from data scale ──────────────────────────────────────
tiny_eps <- function(x, fallback = 1e-12) {
  z <- abs(x[is.finite(x) & x != 0])
  if (length(z) >= 20) {
    out <- stats::quantile(z, 0.01, na.rm = TRUE) / 1000
    if (is.finite(out) && out > 0) return(out)
  }
  fallback
}

# ── Custom method: truncated positive imputation (Tobit/censored Gaussian) ───
# Used for inward current magnitudes. Fits a left-censored Gaussian (Tobit)
# model on observed positive values, falls back to truncated normal if fit fails.
mice.impute.trunc_pos <- function(y, ry, x, left = 1e-12, ...) {
  dots <- list(...)
  wy   <- if (!is.null(dots$wy)) dots$wy else !ry
  idx  <- (!ry) & isTRUE(wy)
  if (is.logical(wy) && length(wy) == length(y)) idx <- (!ry) & wy
  need <- sum(idx)
  if (need == 0L) return(numeric(0))

  xdf <- if (is.null(dim(x))) data.frame(x = x) else as.data.frame(x)
  dat <- data.frame(y = y, xdf, check.names = FALSE)
  obs <- dat[ry,  , drop = FALSE]
  mis <- dat[idx, , drop = FALSE]
  pos_y <- obs$y[is.finite(obs$y) & obs$y > left]

  # Standardize numeric predictors for numerical stability
  num_cols    <- vapply(obs, is.numeric, TRUE); num_cols[1] <- FALSE
  if (any(num_cols)) {
    cn     <- names(obs)[num_cols]
    center <- vapply(obs[cn], function(z) mean(z, na.rm = TRUE), numeric(1))
    scalev <- vapply(obs[cn], function(z) sd(z,   na.rm = TRUE), numeric(1))
    scalev[!is.finite(scalev) | scalev == 0] <- 1
    obs[cn] <- Map(function(z, m, s) (z - m) / s, obs[cn], center, scalev)
    mis[cn] <- Map(function(z, m, s) (z - m) / s, mis[cn], center, scalev)
  }

  mu_fallback <- if (length(pos_y) >= 3) mean(pos_y) else if (length(pos_y) >= 1) median(pos_y) else tiny_eps(y)
  if (!is.finite(mu_fallback) || mu_fallback <= 0) mu_fallback <- tiny_eps(y)
  sd_fallback <- if (length(pos_y) >= 2) stats::sd(pos_y) else max(1, 10 * tiny_eps(y))
  if (!is.finite(sd_fallback) || sd_fallback <= 0) sd_fallback <- max(1, 10 * tiny_eps(y))

  mu  <- rep(mu_fallback, need)
  sd0 <- sd_fallback

  # Attempt Tobit (censored Gaussian) fit
  fit <- tryCatch({
    y_cens <- survival::Surv(pmax(obs$y, left), obs$y > left, type = "left")
    suppressWarnings(survival::survreg(
      y_cens ~ ., data = obs, dist = "gaussian",
      control = survival::survreg.control(maxiter = 200, rel.tolerance = 1e-09)
    ))
  }, warning = function(w) NULL, error = function(e) NULL)

  if (!is.null(fit)) {
    pred <- tryCatch(stats::predict(fit, newdata = mis, type = "response"),
                     error = function(e) rep(NA_real_, need))
    sc   <- tryCatch(as.numeric(fit$scale), error = function(e) NA_real_)
    if (all(is.finite(pred)) && is.finite(sc) && sc > 0) {
      mu  <- pred
      sd0 <- sc
    }
  }

  mu[!is.finite(mu) | mu <= 0] <- tiny_eps(y)
  if (!is.finite(sd0) || sd0 <= 0) sd0 <- max(1, 10 * tiny_eps(y))
  a   <- max(left, tiny_eps(y))
  out <- truncnorm::rtruncnorm(n = need, a = a, b = Inf, mean = mu, sd = sd0)
  if (length(out) != need) stop("trunc_pos: length(out)=", length(out), " != need=", need)
  out
}

# ── Custom method: log1p PMM ──────────────────────────────────────────────────
# Used for exocytosis magnitudes and inward current magnitudes.
# Log-transforms positive values before PMM, back-transforms after.
# Handles NA draws with a parametric fallback.
mice.impute.log1p_pmm <- function(y, ry, x, wy = NULL, ...) {
  if (is.null(wy)) wy <- !ry
  need <- which(!ry & wy)
  if (!length(need)) return(numeric(0))

  z        <- y
  pos      <- is.finite(y) & (y > 0)
  z[pos]   <- log1p(y[pos])
  z[!pos]  <- NA_real_
  ry_z     <- ry & is.finite(z)

  imp_z <- mice:::mice.impute.pmm(z, ry = ry_z, x = x, wy = wy, ...)
  out   <- expm1(imp_z)

  if (anyNA(out)) {
    pos_y <- y[pos]
    if (length(pos_y) >= 5 && sd(log1p(pos_y)) > 0) {
      mu       <- mean(log1p(pos_y)); sg <- sd(log1p(pos_y))
      fallback <- expm1(rnorm(sum(is.na(out)), mean = mu, sd = sg))
    } else {
      eps      <- max(1e-12, quantile(abs(y[is.finite(y) & y != 0]), 0.01, na.rm = TRUE) / 1000)
      fallback <- rep(eps, sum(is.na(out)))
    }
    out[is.na(out)] <- fallback
  }

  pmax(out, .Machine$double.eps)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Section 4: Prepare imputation data frame
# ═══════════════════════════════════════════════════════════════════════════════
# Exocytosis parameters are decomposed into presence (binary) and magnitude
# (positive continuous) components for the two-part imputation model.
# Inward current parameters are decomposed into magnitude (positive, to impute)
# and reconstructed as signed negatives after imputation.

df_work <- df

# Exocytosis: decompose into presence (0/1) and magnitude (positive)
make_exo_parts <- function(x) {
  pres <- ifelse(is.na(x), NA_integer_,
                 ifelse(x > 0 | x == 0, 1L, 0L))
  mag  <- ifelse(is.na(x), NA_real_,
                 ifelse(x > 0, x, ifelse(x == 0, NA_real_, 0)))
  list(pres = pres, mag = mag)
}

for (v in exocytosis_cols) {
  pr <- make_exo_parts(df_work[[v]])
  df_work[[paste0(v, "__pres")]] <- pr$pres
  df_work[[paste0(v, "__mag")]]  <- pr$mag
}

exo_pres <- paste0(exocytosis_cols, "__pres")
exo_mag  <- paste0(exocytosis_cols, "__mag")

for (p in exo_pres[exo_pres %in% names(df_work)])
  df_work[[p]] <- factor(df_work[[p]], levels = c(0, 1))

# Inward currents: magnitude = |x| if x < 0; zeros/positives/NA → NA (to impute)
make_inward_mag <- function(x) ifelse(is.finite(x) & x < 0, -x, NA_real_)
for (v in c(calcium_cols, sodium_cols))
  df_work[[paste0(v, "__mag")]] <- make_inward_mag(df_work[[v]])

# Pre-clean: HalfInactivationofSodiumCurrent_mV exact 0 → NA
if (length(v_half_inact_mv))
  df_work[[v_half_inact_mv]] <- dplyr::na_if(df_work[[v_half_inact_mv]], 0)

# Build imputation frame: drop original signed columns (will be reconstructed)
drop_signed <- c(exocytosis_cols, calcium_cols, sodium_cols)
df_imp <- df_work[, setdiff(names(df_work), drop_signed)]
df_imp <- df_imp %>%
  mutate(cell_type = case_when(
    celltype == "alpha" ~ "alpha",
    celltype == "beta"  ~ "beta",
    celltype == "SCβ"   ~ "SCβ",
    celltype == "PP"    ~ "PP",
    TRUE                ~ "other"
  ))

inward_mag     <- paste0(c(calcium_cols, sodium_cols), "__mag") %>% intersect(names(df_imp))
specials       <- c(v_sodium_peak_ord, v_reversal_signed, v_half_inact_mv) %>% intersect(names(df_imp))
vars_to_impute <- unique(c(exo_pres, exo_mag, inward_mag, specials))

# Ensure presence is factor(0,1)
for (p in exo_pres) if (p %in% names(df_imp))
  df_imp[[p]] <- factor(df_imp[[p]], levels = c(0, 1))

# VoltageforSodiumPeakCurrent_mV must be ordered factor for polr
if (length(v_sodium_peak_ord) && !is.ordered(df_imp[[v_sodium_peak_ord]]))
  df_imp[[v_sodium_peak_ord]] <- factor(df_imp[[v_sodium_peak_ord]], ordered = TRUE)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 5: MICE configuration
# ═══════════════════════════════════════════════════════════════════════════════

# ── Method assignments ────────────────────────────────────────────────────────
ini     <- mice(df_imp, maxit = 0, print = FALSE)
meth    <- ini$method
meth[]  <- ""

meth[v_sodium_peak_ord] <- "polr"       # ordered logistic for voltage factor levels
meth[exo_pres]          <- "logreg"     # logistic regression for exocytosis presence
meth[exo_mag]           <- "log1p_pmm"  # log1p PMM for exocytosis magnitude
meth[inward_mag]        <- "log1p_pmm"  # log1p PMM for inward current magnitudes

# Hard-to-impute inward current magnitudes use standard PMM as fallback
hard_mag <- c(
  "LVA_.60mV_pA__mag", "HVA_.20mV_pA__mag",
  "EarlyPeakCalciumCurrentAmplitudeat.10mV_pA__mag",
  "LateCalciumCurrentAmplitudeat.10mV_pA__mag",
  "CalciumIntegralatFirstDepolarization_pC__mag"
)
meth[hard_mag] <- "pmm"

# Voltage and half-inactivation use midastouch (robust non-parametric matching)
meth["ReversalPotentialbyramp_mV"]          <- "midastouch"
meth["HalfInactivationofSodiumCurrent_mV"]  <- "midastouch"
meth["PeakSodiumCurrentAmplitudeat.10mV_pA__mag"] <- "midastouch"

# ── Predictor matrix ──────────────────────────────────────────────────────────
# Use a fixed, robust predictor set rather than automatic correlation-based
# selection, to ensure stability across imputation iterations.
# Predictors: Study, cell type, cell size, glucose, diabetes status, patcher

robust <- c("Study","cell_type","CellSize_pF","DiabetesStatus",
            "Glucose_mM","Patcher") %>% intersect(colnames(df_imp))

pred    <- quickpred(df_imp, mincor = 0.05, minpuc = 0.25,
                     include = intersect(c("Donor","Plate","Study","Patcher",
                                           "Glucose_mM","TimefromDispersion_days",
                                           "Fresh.or.Cryo","Age","Sex","BMI","HbA1c",
                                           "DiabetesStatus","cell_type","CellSize_pF"),
                                         names(df_imp)),
                     exclude = "CellID")
diag(pred) <- 0

# Override with fixed predictor set for all targets
pred[,] <- 0L
for (p in intersect(exo_pres, rownames(pred))) pred[p, robust] <- 1L
for (i in seq_along(exocytosis_cols)) {
  mag  <- paste0(exocytosis_cols[i], "__mag")
  pres <- paste0(exocytosis_cols[i], "__pres")
  if (mag %in% rownames(pred)) {
    pred[mag, robust] <- 1L
    if (pres %in% colnames(pred)) pred[mag, pres] <- 1L  # each magnitude uses only its own presence
  }
}

all_hard <- intersect(c(inward_mag, specials), rownames(pred))
for (v in all_hard) pred[v, robust] <- 1L

# ── Where matrix: which cells to impute ───────────────────────────────────────
where <- matrix(FALSE, nrow = nrow(df_imp), ncol = ncol(df_imp),
                dimnames = dimnames(ini$where))

for (p in exo_pres)  where[, p] <- is.na(df_imp[[p]])
for (i in seq_along(exocytosis_cols)) {
  pres <- exo_pres[i]; mag <- exo_mag[i]
  if (!all(c(pres, mag) %in% colnames(df_imp))) next
  pval <- as.character(df_imp[[pres]])
  where[, mag] <- is.na(df_imp[[mag]]) & (is.na(pval) | pval == "1")
}
for (m in inward_mag) where[, m] <- is.na(df_imp[[m]])
for (s in specials)   where[, s] <- is.na(df_imp[[s]])

stopifnot(is.logical(where), !anyNA(where))
stopifnot(all(colnames(df_imp)[colSums(where) > 0] %in% vars_to_impute))

# ── Post-processing constraints ───────────────────────────────────────────────
# Applied after each MICE draw to enforce biological sign constraints
post <- ini$post
tiny <- ".Machine$double.eps"
for (v in c(exo_mag, inward_mag)) {
  if (v %in% names(df_imp))
    post[v] <- paste0("imp[[j]][, i] <- pmax(imp[[j]][, i], ", tiny, ")")
}
post["ReversalPotentialbyramp_mV"]         <- "imp[[j]][,i] <- pmin(pmax(imp[[j]][,i], -120), 60)"
post["HalfInactivationofSodiumCurrent_mV"] <- "imp[[j]][,i] <- pmin(pmax(imp[[j]][,i], -120), 60)"
for (b in exocytosis_cols) {
  pres <- paste0(b, "__pres"); mag <- paste0(b, "__mag")
  post[mag] <- sprintf(
    'imp[[j]][,i] <- ifelse(as.character(data[where[, "%s"], "%s"])=="1",
                            pmax(imp[[j]][,i], .Machine$double.eps), 0)',
    mag, pres
  )
}

# ── Visit sequence ────────────────────────────────────────────────────────────
# Impute core inward currents first (most missing), then exocytosis, then specials
core_mag <- intersect(c("EarlyPeakCalciumCurrentAmplitudeat.10mV_pA__mag",
                        "LateCalciumCurrentAmplitudeat.10mV_pA__mag",
                        "LVA_.60mV_pA__mag","HVA_.20mV_pA__mag",
                        "PeakSodiumCurrentAmplitudeat.10mV_pA__mag",
                        "CalciumIntegralatFirstDepolarization_pC__mag",
                        "Hyperpolarizationactivatedcurrent.at.140mV_pA__mag"), names(df_imp))
visit <- c(core_mag,
           intersect(exo_pres, names(df_imp)),
           intersect(exo_mag,  names(df_imp)),
           intersect(c("ReversalPotentialbyramp_mV",
                       "HalfInactivationofSodiumCurrent_mV",
                       "VoltageforSodiumPeakCurrent_mV"), names(df_imp)))

ridge_val <- 1e-6  # ridge regularization for numerical stability

# ── QC: check for rows with no valid predictors ───────────────────────────────
check_empty_rows <- function(df, pred, var) {
  cols <- names(which(pred[var, ] == 1))
  miss <- which(is.na(df[[var]]))
  if (!length(miss)) return(0L)
  sum(rowSums(!is.na(df[miss, cols, drop = FALSE])) == 0)
}
empty_row_counts <- sapply(vars_to_impute, \(v) check_empty_rows(df_imp, pred, v))
if (any(empty_row_counts > 0)) {
  cat("⚠ Variables with rows having no non-missing predictors:\n")
  print(empty_row_counts[empty_row_counts > 0])
}

# ═══════════════════════════════════════════════════════════════════════════════
# Section 6: Run MICE
# ═══════════════════════════════════════════════════════════════════════════════

set.seed(42)
imp <- mice(
  df_imp,
  m               = 50,
  maxit           = 20,
  method          = meth,
  visitSequence   = visit,
  predictorMatrix = pred,
  where           = where,
  post            = post,
  ridge           = ridge_val,
  remove.collinear = TRUE,
  print           = TRUE
)

saveRDS(imp, file = "20260120_imp.rds")
cat("MICE complete. m =", imp$m, "imputations.\n")

imp <- readRDS("20260120_imp.rds")

# ── Quick convergence check ───────────────────────────────────────────────────
# Traces should show good mixing without trends
# plot(imp)  # uncomment to view interactively

# ── Verify coherence: no pres==1 with mag==NA ─────────────────────────────────
comp_test <- mice::complete(imp, 1)
pres_mag_checks <- sapply(exocytosis_cols, function(b) {
  pres <- paste0(b, "__pres"); mag <- paste0(b, "__mag")
  sum(as.character(comp_test[[pres]]) == "1" & is.na(comp_test[[mag]]), na.rm = TRUE)
})
cat("Pres==1 with mag==NA (should all be 0):\n"); print(pres_mag_checks)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 7: Reconstruct signed values and normalize
# ═══════════════════════════════════════════════════════════════════════════════
# After MICE imputation of magnitudes, reconstruct signed values:
# - Exocytosis: positive if pres==1, small negative if pres==0 ("no-event")
# - Inward currents: negate magnitude to restore signed negative values

# ── Helper: safe cell-size normalization ──────────────────────────────────────
safe_div <- function(num, den, use_abs_den = TRUE, absolute_floor = 0.5,
                     set_na_if_tiny = TRUE, clamp_if_tiny = FALSE) {
  out  <- rep(NA_real_, length(num))
  good <- is.finite(num) & is.finite(den)
  if (!any(good)) return(out)
  d <- den[good]
  if (use_abs_den) {
    d_abs <- abs(d)
    tiny  <- d_abs < absolute_floor
    d_eff <- sign(d) * ifelse(tiny & clamp_if_tiny, absolute_floor, d_abs)
  } else {
    tiny  <- (d >= 0 & d < absolute_floor) | (d < 0 & d > -absolute_floor)
    d_eff <- ifelse(d >= 0, pmax(d, absolute_floor), pmin(d, -absolute_floor))
    if (!clamp_if_tiny) d_eff[tiny] <- NA_real_
  }
  if (set_na_if_tiny && any(tiny)) {
    keep          <- !tiny
    out[good][keep] <- num[good][keep] / d_eff[keep]
  } else {
    out[good] <- num[good] / d_eff
  }
  out
}

# ── Reconstruct signed values from imputed magnitudes ─────────────────────────
reconstruct_signed <- function(comp, max_abs_no_event = 0.5,
                               nudge_min = 1e-6, nudge_max = 1e-3) {
  if (is.list(comp) && !is.null(comp$data)) comp <- comp$data
  out <- comp
  .eps <- function(x) {
    if (exists("tiny_eps", mode = "function")) tiny_eps(x) else {
      z <- abs(x[is.finite(x) & x != 0])
      if (length(z) >= 20) stats::quantile(z, 0.01, na.rm = TRUE) / 1000 else 1e-12
    }
  }

  # Exocytosis
  for (b in exocytosis_cols) {
    pres <- paste0(b, "__pres"); mag <- paste0(b, "__mag")
    if (!all(c(pres, mag) %in% names(out))) next
    val  <- rep(NA_real_, nrow(out))
    pchr <- as.character(out[[pres]])
    idx1 <- which(pchr == "1"); idx0 <- which(pchr == "0")
    idxNA <- which(is.na(pchr))

    if (length(idxNA)) {
      m <- out[[mag]][idxNA]; e <- .eps(m)
      use_mag <- which(is.finite(m) & m > 0)
      if (length(use_mag)) { m[use_mag] <- pmax(m[use_mag], e); val[idxNA[use_mag]] <- m[use_mag] }
      need_neg <- setdiff(seq_along(idxNA), use_mag)
      if (length(need_neg)) {
        lo <- max(max_abs_no_event * 0.1, nudge_min); hi <- max(max_abs_no_event, lo + 1e-6)
        val[idxNA[need_neg]] <- -runif(length(need_neg), lo, hi)
      }
    }
    if (length(idx1)) {
      m <- out[[mag]][idx1]; m[!is.finite(m)] <- NA_real_
      if (any(is.finite(m))) {
        e    <- .eps(m)
        nonpos <- which(is.finite(m) & m <= 0)
        if (length(nonpos)) m[nonpos] <- runif(length(nonpos), min = e, max = 2 * e)
        m <- pmax(m, e, na.rm = TRUE)
      }
      val[idx1] <- m
    }
    if (length(idx0)) {
      lo <- max(max_abs_no_event * 0.1, nudge_min); hi <- max(max_abs_no_event, lo + 1e-6)
      val[idx0] <- -runif(length(idx0), min = lo, max = hi)
    }
    out[[b]] <- val
  }

  # Inward currents: signed = −magnitude
  for (b in c(calcium_cols, sodium_cols)) {
    mcol <- paste0(b, "__mag")
    if (!mcol %in% names(out)) next
    s        <- -out[[mcol]]
    s[!is.finite(s)] <- NA_real_
    bad <- which(is.finite(s) & s >= 0)
    if (length(bad))
      s[bad] <- -pmax(runif(length(bad), min = nudge_min, max = nudge_max), nudge_min)
    out[[b]] <- s
  }
  out
}

# ── Recompute normalized parameters ───────────────────────────────────────────
# Normalization is recomputed from reconstructed signed values after each
# imputation rather than imputing normalized values directly, to ensure
# consistency between raw and normalized columns.
recompute_normalized_safe <- function(df, tiny_pf = 0.5, mode = c("na", "clamp")) {
  mode      <- match.arg(mode)
  use_na    <- (mode == "na")
  use_clamp <- (mode == "clamp")
  out       <- df

  sdiv <- function(n, d) safe_div(n, d, use_abs_den = TRUE, absolute_floor = tiny_pf,
                                  set_na_if_tiny = use_na, clamp_if_tiny = use_clamp)

  if (all(c("CalciumIntegralatFirstDepolarization_pC","CellSize_pF") %in% names(out)))
    out$CalciumIntegralNormalizedtoCellSize_pC.pF <-
      sdiv(out$CalciumIntegralatFirstDepolarization_pC, out$CellSize_pF)

  if (all(c("TotalCapacitance_fF","CellSize_pF") %in% names(out)))
    out$NormalizedTotalCapacitance_fF.pF <-
      sdiv(out$TotalCapacitance_fF, out$CellSize_pF)

  if (all(c("FirstDepolarizationCapacitance_fF","CellSize_pF") %in% names(out)))
    out$NormalizedFirstDepolarizationCapacitance_fF.pF <-
      sdiv(out$FirstDepolarizationCapacitance_fF, out$CellSize_pF)

  if (all(c("LateDepolarizationCapacitance","CellSize_pF") %in% names(out)))
    out$NormalizedLateDepolarizationCapacitance <-
      sdiv(out$LateDepolarizationCapacitance, out$CellSize_pF)

  inward <- c("EarlyPeakCalciumCurrentAmplitudeat.10mV_pA",
              "LateCalciumCurrentAmplitudeat.10mV_pA",
              "LVA_.60mV_pA","HVA_.20mV_pA",
              "PeakSodiumCurrentAmplitudeat.10mV_pA")
  nmaps  <- c("NormalizedEarlyPeakCalciumCurrentAmplitudeat.10mV_pA.pF",
              "NormalizedLateCalciumCurrentAmplitudeat.10mV_pA.pF",
              "NormalizedLVA_.60mV_pA.pF","NormalizedHVA_.20mV_pA.pF",
              "NormalizedPeakSodiumCurrentAmplitudeat.10mV_pA.pF",
              "NormalizedHyperpolarizationactivatedcurrent.at.140mV_pA.pF")
  for (i in seq_along(inward))
    if (all(c(inward[i], "CellSize_pF") %in% names(out)))
      out[[nmaps[i]]] <- sdiv(out[[inward[i]]], out$CellSize_pF)

  out
}

# ── Build 50 completed datasets ───────────────────────────────────────────────
completed_list <- lapply(seq_len(imp$m), function(k) {
  mice::complete(imp, k) %>%
    reconstruct_signed(max_abs_no_event = 0.5) %>%
    recompute_normalized_safe(mode = "na")
})

# ── QC: verify sign constraints ───────────────────────────────────────────────
check_signs <- function(comp) {
  bad_exo <- sum(unlist(lapply(exocytosis_cols, \(v) sum(comp[[v]] < 0, na.rm = TRUE))))
  bad_ca  <- sum(unlist(lapply(calcium_cols,    \(v) sum(comp[[v]] >= 0, na.rm = TRUE))))
  bad_na  <- sum(unlist(lapply(sodium_cols,     \(v) sum(comp[[v]] >= 0, na.rm = TRUE))))
  cat("Exocytosis negatives:", bad_exo,
      "| Ca non-negatives:", bad_ca,
      "| Na non-negatives:", bad_na, "\n")
}
lapply(completed_list, check_signs)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 8: Rubin's rules pooling for Seurat metadata
# ═══════════════════════════════════════════════════════════════════════════════
# Pool imputed values using Rubin's rules to get a single estimate per cell
# for attaching to Seurat metadata. This is used for QC and visualization;
# the full completed_list is used for all statistical analyses.

vars_to_pool <- c(
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

vars_to_pool_raw <- unique(
  c(exocytosis_cols, calcium_cols, sodium_cols,
    "ReversalPotentialbyramp_mV", "HalfInactivationofSodiumCurrent_mV",
    "VoltageforSodiumPeakCurrent_mV")[
    c(exocytosis_cols, calcium_cols, sodium_cols,
      "ReversalPotentialbyramp_mV", "HalfInactivationofSodiumCurrent_mV",
      "VoltageforSodiumPeakCurrent_mV") %in% names(df)]
)

id_col <- "CellID"
stopifnot(id_col %in% names(completed_list[[1]]))

# Stack all imputations
long <- dplyr::bind_rows(lapply(seq_along(completed_list), function(k) {
  d      <- completed_list[[k]][, c(id_col, vars_to_pool), drop = FALSE]
  d$.imp <- k; d
}))

long <- long %>%
  mutate(VoltageforSodiumPeakCurrent_mV =
           as.numeric(as.character(VoltageforSodiumPeakCurrent_mV)))

# Track which values were originally missing or zero
was_na   <- df[, c(id_col, vars_to_pool_raw), drop = FALSE] %>%
  mutate(across(all_of(vars_to_pool_raw), ~is.na(.x))) %>%
  rename_with(~paste0(.x, "_was_imputed"), all_of(vars_to_pool_raw))

was_zero <- df[, c(id_col, vars_to_pool_raw), drop = FALSE] %>%
  mutate(across(all_of(vars_to_pool_raw), ~!is.na(.x) & .x == 0)) %>%
  rename_with(~paste0(.x, "_was_zero"), all_of(vars_to_pool_raw))

# Pool: mean across imputations + Rubin SE
pooled_num <- long %>%
  tidyr::pivot_longer(cols = all_of(vars_to_pool), names_to = "var", values_to = "val") %>%
  dplyr::select(.imp, !!sym(id_col), var, val) %>%
  group_by(var, !!sym(id_col)) %>%
  summarise(
    m         = dplyr::n(),
    q_bar     = mean(val),           # pooled mean
    b         = stats::var(val),     # between-imputation variance
    T_var     = (1 + 1/m) * b,       # Rubin total variance
    se_pooled = sqrt(T_var),
    .groups   = "drop"
  )

# ── QC: SE distributions by imputation status ─────────────────────────────────
# Imputed values should have larger SE than observed; constructed zeros should
# have SE = 0. Flag any unexpected patterns.
pn <- pooled_num %>%
  dplyr::filter(var != "CapacitanceNormalizedtoCalcium_fF.pC")

se_summary <- pn %>%
  group_by(var) %>%
  summarise(
    n_cells     = n(),
    med_se      = median(se_pooled, na.rm = TRUE),
    p95_se      = quantile(se_pooled, 0.95, na.rm = TRUE),
    max_se      = max(se_pooled, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  arrange(desc(med_se))
cat("Pooled SE summary by variable:\n")
print(se_summary, n = Inf)

# ═══════════════════════════════════════════════════════════════════════════════
# Section 9: Attach pooled estimates to Seurat object and save
# ═══════════════════════════════════════════════════════════════════════════════

load("../human_islets_scRNAseq/20250708_pclamp_all.RData")
pclamp_patched_all <- subset(
  pclamp_all,
  cells = pclamp_all$CellID[!is.na(pclamp_all$CellSize_pF) &
                               pclamp_all$CellSize_pF != 0]
)
ids_seurat <- as.character(Cells(pclamp_patched_all))

# Build wide tables: one column per variable, prefixed mi_mean_ or mi_se_
pn_filt <- pn %>% dplyr::filter(CellID %in% ids_seurat)

pn_mean_wide <- pn_filt %>%
  dplyr::select(CellID, var, q_bar) %>%
  tidyr::pivot_wider(names_from = var, values_from = q_bar, names_prefix = "mi_mean_")

pn_se_wide <- pn_filt %>%
  dplyr::select(CellID, var, se_pooled) %>%
  tidyr::pivot_wider(names_from = var, values_from = se_pooled, names_prefix = "mi_se_")

# Attach to Seurat metadata
md <- pclamp_patched_all@meta.data %>%
  left_join(pn_mean_wide, by = "CellID") %>%
  left_join(pn_se_wide,   by = "CellID") %>%
  as.data.frame()
rownames(md) <- md$CellID
pclamp_patched_all@meta.data <- md

# Save Seurat object with imputed metadata
save(pclamp_patched_all,
     file = "20251210_patchSeq_integrated_patched-only_decon_mi.RData")

# Save imp bundle (primary output for all downstream analyses)
imp_bundle <- list(imp = imp, completed_list = completed_list)
saveRDS(imp_bundle, file = "20260120_imp_bundle.rds")

# Save imputation settings for reproducibility
imputation_settings <- list(
  timestamp        = Sys.time(),
  seed             = 42,
  m                = imp$m,
  maxit            = 20,
  ridge            = ridge_val,
  remove.collinear = TRUE,
  methods          = meth,
  visitSequence    = visit,
  predictorMatrix  = pred,
  where            = where,
  post             = post,
  specials = list(
    exocytosis_cols   = exocytosis_cols,
    calcium_cols      = calcium_cols,
    sodium_cols       = sodium_cols,
    special_ord       = v_sodium_peak_ord,
    special_signed0   = v_reversal_signed
  ),
  notes = list(
    ReversalPotentialbyramp_mV         = "midastouch, bounded [-120, 60]",
    HalfInactivationofSodiumCurrent_mV = "midastouch, bounded [-120, 60]",
    exocytosis_pipeline                = "two-part: logreg(__pres) + log1p_pmm(__mag) → reconstruct_signed()",
    inward_currents                    = "log1p_pmm on __mag > 0, then store as −__mag"
  ),
  sessionInfo = utils::sessionInfo()
)
saveRDS(imputation_settings, file = "20260120_imputation_settings_bundle.rds")

cat("\nDone. Outputs:\n")
cat("  20260120_imp.rds                              - MICE mids object\n")
cat("  20260120_imp_bundle.rds                       - imp + completed_list (50 datasets)\n")
cat("  20260120_imputation_settings_bundle.rds       - reproducibility record\n")
cat("  20251210_patchSeq_integrated_patched-only_decon_mi.RData - Seurat with mi_ metadata\n")