# ---------------------------------------------------------------------------------------------------------------
# Repository copy. Run from the repository root. Inputs: data/ (covariate-annotated Seurat objects from
# Full_merge_RDS_script_with_covariates.R) and results/<analysis output folder> (rename each timestamped output
# folder to the name used below, e.g. CX_TE_Control_<date> -> results/CX_TE_Control). R 4.4.3.
# ---------------------------------------------------------------------------------------------------------------
# ============================================================
# SCRIPT 5: PD vs CONTROL — INNATE IMMUNE SUB-PANEL COMPARISON
# Astrocytes vs Neurons, Covariate-Adjusted
# ============================================================
#
# PURPOSE:
#   Before examining dsRNA_Sensing × Immune coupling, establish
#   whether innate immune sub-panel scores differ between PD and
#   Control at the donor level, separately for astrocytes and
#   neurons.  This "baseline landscape" analysis answers:
#     Q1: Which immune pathways are up/down in PD astrocytes?
#     Q2: Which are up/down in PD neurons?
#     Q3: Are PD-vs-Control differences cell-type-specific?
#
# STATISTICAL APPROACH:
#   PRIMARY:   Covariate-adjusted linear model
#              lm(panel_score ~ age + pmi + sex + condition)
#              Reports β_PD, SE, t, p, BH-adjusted p
#   SECONDARY: Wilcoxon rank-sum (unadjusted, for robustness)
#   EFFECT SIZE: Hedges' g (bias-corrected Cohen's d)
#   INTERACTION: Cell-type × condition interaction test
#              lm(score ~ age + pmi + sex + condition * cell_type)
#
# IMMUNE SUB-PANELS (9 — same as Script 4):
#   ISG_Core, IFN_STAT, NF_kB, Complement, OAS_RNaseL,
#   Astro_React, Inflammasome, TLR_dsRNA, TLR_DAMP
#
# INPUTS:
#   - Astrocyte Seurat RDS (with covariates)
#   - Neuron Seurat RDS (with covariates)
#   - (Optional) 7-Module AllResults RDS for module scores
#
# OUTPUTS:
#   Tables:
#     - PD_vs_Control_Immune_Panels_Summary.xlsx
#       Sheets: Astrocytes, Neurons, Interaction, GeneLevelDE,
#               PanelSizeQC, Methods
#   Figures:
#     - Fig1_Violin_PanelScores_ByCondition.pdf
#     - Fig2_Heatmap_EffectSize_Astro_vs_Neuron.pdf
#     - Fig3_ForestPlot_CovAdj_Betas.pdf
#     - Fig4_Interaction_DotPlot.pdf
#     - Fig5_GeneLevel_Volcano.pdf
#
# ============================================================

rm(list = ls()); gc()

# =============================================================================
# SECTION 0: LIBRARIES
# =============================================================================

cat("========== LOADING LIBRARIES ==========\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(ggrepel)
  library(RColorBrewer)
  library(scales)
})

set.seed(2026)

# ---- Strip S4 Rle / list columns to base R types (Seurat v5 fix) ----
strip_rle <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (cn in colnames(df)) {
    x <- df[[cn]]
    if (!is.atomic(x)) {
      df[[cn]] <- tryCatch(as.vector(x),
                           error = function(e) as.character(x))
    }
  }
  df
}

# =============================================================================
# SECTION 1: CONFIGURATION
# =============================================================================

cat("\n========== CONFIGURATION ==========\n")

# --- Region: change to "SN" or "CX" per run ---
REGION <- "CX"

# --- Covariate control ---
COVARIATES_ENABLED <- TRUE
COV_AGE_COL <- "age"
COV_PMI_COL <- "pmi"
COV_SEX_COL <- "sex"
covariate_cols_to_use <- c()

# --- Interactive = FALSE → use hardcoded paths below ---
INTERACTIVE <- FALSE
if (!INTERACTIVE) {
  if (REGION == "SN") {
    astro_file  <- "data/SN_Astrocytes.rds"
    neuron_file <- "data/SN_Neurons_Total.rds"
    output_base <- "results"
  } else {
    astro_file  <- "data/CX_Astrocytes_withCov.rds"
    neuron_file <- "data/CX_Neurons_Total_withCov.rds"
    output_base <- "results"
  }
} else {
  library(tcltk)
  cat("Please select your ASTROCYTE RDS file...\n")
  astro_file <- tk_choose.files(
    caption = paste0("Select ", REGION, " ASTROCYTE RDS file"),
    filters = matrix(c("RDS Files", ".rds"), 1, 2)
  )
  if (length(astro_file) == 0 || is.na(astro_file)) stop("No astrocyte file selected.")

  cat("Please select your NEURON RDS file...\n")
  neuron_file <- tk_choose.files(
    caption = paste0("Select ", REGION, " NEURON RDS file"),
    filters = matrix(c("RDS Files", ".rds"), 1, 2)
  )
  if (length(neuron_file) == 0 || is.na(neuron_file)) stop("No neuron file selected.")

  cat("Please select OUTPUT folder...\n")
  output_base <- tk_choose.dir(caption = "Select OUTPUT folder")
  if (is.na(output_base)) stop("No folder selected.")
}

timestamp   <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- file.path(output_base,
                         paste0(REGION, "_PD_vs_Ctrl_ImmunePanels_", timestamp))
tbl_dir <- file.path(results_dir, "tables")
fig_dir <- file.path(results_dir, "figures")
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("Region:", REGION, "\n")
cat("Output directory:", results_dir, "\n\n")

# =============================================================================
# SECTION 2: IMMUNE SUB-PANEL DEFINITIONS (9 panels — matches Script 4)
# =============================================================================

cat("========== DEFINING IMMUNE SUB-PANELS ==========\n")

immune_subpanels <- list(
  ISG_Core      = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
                    "IFITM1", "IFITM2", "IFITM3", "RSAD2", "BST2", "USP18"),
  IFN_STAT      = c("STAT1", "STAT2", "IRF1", "IRF3", "IRF7", "IRF9",
                    "JAK1", "TYK2"),
  NF_kB         = c("NFKB1", "NFKB2", "RELA", "RELB", "NFKBIA",
                    "IL6", "TNF", "CCL2", "CCL5", "CXCL10"),
  Complement    = c("C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B",
                    "SERPING1", "CFB", "CFH"),
  OAS_RNaseL    = c("OAS1", "OAS2", "OAS3", "OASL", "RNASEL"),
  Astro_React   = c("GFAP", "VIM", "SERPINA3", "LCN2", "CHI3L1", "OSMR"),
  Inflammasome  = c("NLRP3", "NLRP1", "AIM2", "PYCARD", "CASP1", "CASP4",
                    "IL18", "IL1B", "GSDMD"),
  TLR_dsRNA     = c("TLR3", "TICAM1", "TICAM2"),
  TLR_DAMP      = c("TLR2", "TLR4", "MYD88", "TIRAP", "IRAK1", "IRAK4", "TRAF6")
)

panel_names <- names(immune_subpanels)
cat("Panels:", paste(panel_names, collapse = ", "), "\n")
cat("Total unique genes:", length(unique(unlist(immune_subpanels))), "\n\n")

# =============================================================================
# SECTION 3: LOAD SEURAT OBJECTS
# =============================================================================

cat("========== LOADING SEURAT OBJECTS ==========\n")

astro  <- readRDS(astro_file)
neuron <- readRDS(neuron_file)
DefaultAssay(astro)  <- "RNA"
DefaultAssay(neuron) <- "RNA"
if (length(Layers(astro))  > 1) astro  <- JoinLayers(astro)
if (length(Layers(neuron)) > 1) neuron <- JoinLayers(neuron)

cat("Astrocytes:", ncol(astro),  "cells\n")
cat("Neurons:   ", ncol(neuron), "cells\n\n")

# --- Metadata column detection ---
detect_col <- function(meta, candidates, label) {
  for (nm in candidates) if (nm %in% colnames(meta)) return(nm)
  warning("Missing column: ", label); NULL
}
donor_possible <- c("donor_id", "Donor_ID", "donor", "Donor", "sample_id")
cond_possible  <- c("condition", "Condition", "disease", "Disease",
                    "Status", "Group", "Diagnosis")

donor_col_astro  <- detect_col(astro@meta.data,  donor_possible, "donor (astro)")
cond_col_astro   <- detect_col(astro@meta.data,  cond_possible,  "condition (astro)")
donor_col_neuron <- detect_col(neuron@meta.data, donor_possible, "donor (neuron)")
cond_col_neuron  <- detect_col(neuron@meta.data, cond_possible,  "condition (neuron)")

cat("Donor column (astro):",  donor_col_astro,  "\n")
cat("Cond  column (astro):",  cond_col_astro,   "\n")
cat("Donor column (neuron):", donor_col_neuron, "\n")
cat("Cond  column (neuron):", cond_col_neuron,  "\n\n")

# --- Covariate detection ---
if (COVARIATES_ENABLED) {
  find_col <- function(meta, candidates) {
    for (nm in candidates) if (nm %in% colnames(meta) &&
                               sum(!is.na(meta[[nm]])) > 0) return(nm)
    NULL
  }
  age_col <- find_col(astro@meta.data, c(COV_AGE_COL, "age", "Age", "AGE", "age_death"))
  pmi_col <- find_col(astro@meta.data, c(COV_PMI_COL, "pmi", "PMI", "post_mortem_interval"))
  sex_col <- find_col(astro@meta.data, c(COV_SEX_COL, "sex", "Sex", "gender"))
  for (col in c(age_col, pmi_col, sex_col)) {
    if (!is.null(col)) covariate_cols_to_use <- c(covariate_cols_to_use, col)
  }
  cat("Active covariates:",
      if (length(covariate_cols_to_use) > 0)
        paste(covariate_cols_to_use, collapse = ", ") else "NONE",
      "\n\n")
}

# =============================================================================
# SECTION 4: HELPER FUNCTIONS
# =============================================================================

cat("========== DEFINING HELPER FUNCTIONS ==========\n")

# ---- Background-subtracted donor-level score (identical to Script 4) ----
donor_bg_subtracted_score <- function(obj, genes, donor_col, cond_col,
                                      n_ctrl = 100, n_bins = 25) {
  avail <- genes[genes %in% rownames(obj)]
  if (length(avail) < 2) return(NULL)
  expr_all <- GetAssayData(obj, layer = "data")
  cell_means <- Matrix::rowMeans(expr_all)
  bins <- cut(cell_means, breaks = n_bins, labels = FALSE)
  names(bins) <- rownames(expr_all)
  ctrl_pool <- c()
  for (g in avail) {
    b <- bins[g]; if (is.na(b)) next
    pool <- setdiff(names(bins)[which(bins == b)], avail)
    if (length(pool) >= n_ctrl) {
      ctrl_pool <- c(ctrl_pool, sample(pool, n_ctrl))
    } else if (length(pool) > 0) {
      ctrl_pool <- c(ctrl_pool, pool)
    }
  }
  ctrl_pool <- unique(ctrl_pool)
  if (length(ctrl_pool) == 0) return(NULL)
  panel_mean <- Matrix::colMeans(expr_all[avail,     , drop = FALSE])
  ctrl_mean  <- Matrix::colMeans(expr_all[ctrl_pool, , drop = FALSE])
  cell_score <- as.numeric(panel_mean - ctrl_mean)
  meta <- strip_rle(as.data.frame(obj@meta.data))
  agg_df <- data.frame(
    donor_id  = as.character(meta[[donor_col]]),
    condition = as.character(meta[[cond_col]]),
    score_raw = cell_score,
    stringsAsFactors = FALSE
  )
  donor_df <- aggregate(score_raw ~ donor_id + condition, data = agg_df,
                        FUN = mean, na.rm = TRUE)
  colnames(donor_df)[colnames(donor_df) == "score_raw"] <- "score"
  donor_df <- as.data.frame(donor_df, stringsAsFactors = FALSE)
  list(donor_df = donor_df, genes_used = avail)
}

# ---- Pull covariates for donors from a Seurat object ----
get_donor_cov_df <- function(obj, donor_col, donor_ids, cov_cols) {
  if (length(cov_cols) == 0) return(NULL)
  meta <- strip_rle(as.data.frame(obj@meta.data))
  keep <- cov_cols[cov_cols %in% colnames(meta)]
  if (length(keep) == 0) return(NULL)
  meta$donor_key <- as.character(meta[[donor_col]])
  unique_donors <- unique(meta$donor_key)
  donor_cov_list <- lapply(unique_donors, function(d) {
    sub <- meta[meta$donor_key == d, keep, drop = FALSE]
    vals <- sapply(keep, function(cn) {
      x <- sub[[cn]]
      if (is.numeric(x)) mean(x, na.rm = TRUE)
      else as.character(x[1])
    })
    vals
  })
  donor_meta <- as.data.frame(do.call(rbind, donor_cov_list),
                              stringsAsFactors = FALSE)
  colnames(donor_meta) <- keep
  donor_meta$donor_key <- unique_donors
  for (cn in keep) {
    num_try <- suppressWarnings(as.numeric(donor_meta[[cn]]))
    if (!all(is.na(num_try))) donor_meta[[cn]] <- num_try
  }
  donor_meta <- donor_meta[match(donor_ids, donor_meta$donor_key), , drop = FALSE]
  cov_df <- donor_meta[, keep, drop = FALSE]
  # Encode categorical → 0/1
  for (cn in colnames(cov_df)) {
    if (!is.numeric(cov_df[[cn]])) {
      cov_df[[cn]] <- as.numeric(factor(cov_df[[cn]])) - 1
    }
    if (any(is.na(cov_df[[cn]])))
      cov_df[[cn]][is.na(cov_df[[cn]])] <- mean(cov_df[[cn]], na.rm = TRUE)
  }
  keep2 <- sapply(cov_df, function(x) sd(x, na.rm = TRUE) > 0)
  cov_df[, keep2, drop = FALSE]
}

# ---- Hedges' g (bias-corrected Cohen's d) ----
hedges_g <- function(x_pd, x_ctrl) {
  n1 <- length(x_pd); n2 <- length(x_ctrl)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  m1 <- mean(x_pd, na.rm = TRUE); m2 <- mean(x_ctrl, na.rm = TRUE)
  s1 <- sd(x_pd, na.rm = TRUE);   s2 <- sd(x_ctrl, na.rm = TRUE)
  s_pooled <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (s_pooled == 0) return(NA_real_)
  d <- (m1 - m2) / s_pooled
  # Hedges correction factor
  J <- 1 - 3 / (4 * (n1 + n2 - 2) - 1)
  d * J
}

# ---- Significance stars ----
sig_star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*", ""))))
}

# =============================================================================
# SECTION 5: COMPUTE IMMUNE PANEL DONOR SCORES (bg-subtracted)
# =============================================================================

cat("\n========== COMPUTING IMMUNE PANEL DONOR SCORES ==========\n")

panel_scores <- list()
panel_qc     <- list()

for (ct in list(
  list(name = "Astrocytes", obj = astro,  d = donor_col_astro,  c = cond_col_astro),
  list(name = "Neurons",    obj = neuron, d = donor_col_neuron, c = cond_col_neuron)
)) {
  cat("--- ", ct$name, " ---\n")
  panel_scores[[ct$name]] <- list()
  for (pn in panel_names) {
    res <- donor_bg_subtracted_score(ct$obj, immune_subpanels[[pn]], ct$d, ct$c)
    if (is.null(res)) {
      cat(sprintf("  %-14s : SKIPPED (insufficient genes)\n", pn)); next
    }
    panel_scores[[ct$name]][[pn]] <- res$donor_df
    panel_qc[[paste0(ct$name, "_", pn)]] <- data.frame(
      cell_type   = ct$name,
      panel       = pn,
      total_genes = length(immune_subpanels[[pn]]),
      genes_found = length(res$genes_used),
      genes_used  = paste(res$genes_used, collapse = ", "),
      genes_miss  = paste(setdiff(immune_subpanels[[pn]], res$genes_used), collapse = ", "),
      n_donors    = nrow(res$donor_df),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %-14s : %d/%d genes, %d donors\n",
                pn, length(res$genes_used), length(immune_subpanels[[pn]]),
                nrow(res$donor_df)))
  }
}

panel_qc_df <- do.call(rbind, panel_qc)
rownames(panel_qc_df) <- NULL

# =============================================================================
# SECTION 6: PD vs CONTROL COMPARISON — COVARIATE-ADJUSTED LM + WILCOXON
# =============================================================================

cat("\n========== PD vs CONTROL: PANEL-LEVEL COMPARISONS ==========\n")

comparison_long <- list()

for (ct_name in c("Astrocytes", "Neurons")) {
  obj       <- if (ct_name == "Astrocytes") astro else neuron
  donor_col <- if (ct_name == "Astrocytes") donor_col_astro else donor_col_neuron
  ct_panels <- panel_scores[[ct_name]]
  if (length(ct_panels) == 0) next

  cat(sprintf("--- %s (%d panels) ---\n", ct_name, length(ct_panels)))

  for (pn in names(ct_panels)) {
    panel_df <- ct_panels[[pn]]
    if (is.null(panel_df) || nrow(panel_df) < 6) next

    panel_df$donor_id  <- trimws(as.character(panel_df$donor_id))
    panel_df$condition <- trimws(as.character(panel_df$condition))

    # Standardise condition labels to PD / Control
    panel_df$cond_std <- NA_character_
    panel_df$cond_std[grepl("PD|Parkinson|Disease", panel_df$condition,
                            ignore.case = TRUE)] <- "PD"
    panel_df$cond_std[grepl("Control|Ctrl|Healthy|HC",
                            panel_df$condition, ignore.case = TRUE)] <- "Control"
    if (any(is.na(panel_df$cond_std))) {
      # Fallback: keep original labels if standard mapping fails
      unmapped <- unique(panel_df$condition[is.na(panel_df$cond_std)])
      cat(sprintf("  [%s] WARNING: unmapped conditions: %s\n",
                  pn, paste(unmapped, collapse = ", ")))
      panel_df$cond_std[is.na(panel_df$cond_std)] <-
        panel_df$condition[is.na(panel_df$cond_std)]
    }

    # Split
    score_pd   <- panel_df$score[panel_df$cond_std == "PD"]
    score_ctrl <- panel_df$score[panel_df$cond_std == "Control"]
    n_pd   <- length(score_pd)
    n_ctrl <- length(score_ctrl)

    if (n_pd < 3 || n_ctrl < 3) {
      cat(sprintf("  %-14s : PD=%d, Ctrl=%d — too few, skipping\n",
                  pn, n_pd, n_ctrl)); next
    }

    # ---- Descriptive stats ----
    mean_pd   <- mean(score_pd,   na.rm = TRUE)
    mean_ctrl <- mean(score_ctrl, na.rm = TRUE)
    fc        <- mean_pd - mean_ctrl
    g         <- hedges_g(score_pd, score_ctrl)

    # ---- Wilcoxon (unadjusted) ----
    wt <- tryCatch(
      wilcox.test(score_pd, score_ctrl, exact = FALSE),
      error = function(e) list(statistic = NA, p.value = NA)
    )

    # ---- Covariate-adjusted linear model ----
    beta_cov <- NA_real_; se_cov <- NA_real_
    t_cov    <- NA_real_; p_cov <- NA_real_

    if (COVARIATES_ENABLED && length(covariate_cols_to_use) > 0) {
      cov_df <- get_donor_cov_df(obj, donor_col, panel_df$donor_id,
                                 covariate_cols_to_use)
      if (!is.null(cov_df) && nrow(cov_df) == nrow(panel_df)) {
        lm_df <- cbind(
          data.frame(score     = panel_df$score,
                     condition = factor(panel_df$cond_std, levels = c("Control", "PD"))),
          cov_df
        )
        cov_nms <- colnames(cov_df)
        fmla <- as.formula(paste0("score ~ ",
                                  paste(cov_nms, collapse = " + "),
                                  " + condition"))
        fit <- tryCatch(lm(fmla, data = lm_df), error = function(e) NULL)
        if (!is.null(fit)) {
          s <- summary(fit)$coefficients
          if ("conditionPD" %in% rownames(s)) {
            beta_cov <- s["conditionPD", "Estimate"]
            se_cov   <- s["conditionPD", "Std. Error"]
            t_cov    <- s["conditionPD", "t value"]
            p_cov    <- s["conditionPD", "Pr(>|t|)"]
          }
        }
      }
    }

    # If covariates failed or disabled, fall back to simple lm
    if (is.na(p_cov)) {
      lm_simple <- lm(score ~ condition,
                       data = data.frame(
                         score     = panel_df$score,
                         condition = factor(panel_df$cond_std,
                                            levels = c("Control", "PD"))
                       ))
      s2 <- summary(lm_simple)$coefficients
      if ("conditionPD" %in% rownames(s2)) {
        beta_cov <- s2["conditionPD", "Estimate"]
        se_cov   <- s2["conditionPD", "Std. Error"]
        t_cov    <- s2["conditionPD", "t value"]
        p_cov    <- s2["conditionPD", "Pr(>|t|)"]
      }
    }

    comparison_long[[length(comparison_long) + 1]] <- data.frame(
      cell_type  = ct_name,
      panel      = pn,
      n_PD       = n_pd,
      n_Control  = n_ctrl,
      mean_PD    = round(mean_pd, 6),
      mean_Ctrl  = round(mean_ctrl, 6),
      diff       = round(fc, 6),
      hedges_g   = round(g, 4),
      direction  = ifelse(fc > 0, "UP_in_PD", "DOWN_in_PD"),
      beta_lm    = round(beta_cov, 6),
      se_lm      = round(se_cov, 6),
      t_lm       = round(t_cov, 3),
      p_lm       = p_cov,
      p_wilcox   = wt$p.value,
      stringsAsFactors = FALSE
    )

    cat(sprintf("  %-14s : PD=%d Ctrl=%d  g=%.3f  β=%.4f  p_lm=%.4g  p_wilcox=%.4g  %s\n",
                pn, n_pd, n_ctrl, g, beta_cov, p_cov, wt$p.value,
                ifelse(fc > 0, "↑PD", "↓PD")))
  }
}

# --- Combine & BH-correct ---
comp_df <- do.call(rbind, comparison_long)
comp_df <- strip_rle(comp_df)

comp_df$p_BH_lm     <- NA_real_
comp_df$p_BH_wilcox <- NA_real_
comp_df$star_lm     <- ""
comp_df$star_wilcox <- ""

for (ct in unique(comp_df$cell_type)) {
  idx <- comp_df$cell_type == ct
  comp_df$p_BH_lm[idx]     <- p.adjust(comp_df$p_lm[idx],     method = "BH")
  comp_df$p_BH_wilcox[idx] <- p.adjust(comp_df$p_wilcox[idx], method = "BH")
  comp_df$star_lm[idx]     <- sig_star(comp_df$p_BH_lm[idx])
  comp_df$star_wilcox[idx] <- sig_star(comp_df$p_BH_wilcox[idx])
}

cat("\n=== Summary of FDR-significant panels (lm, BH<0.05) ===\n")
sig_rows <- comp_df[!is.na(comp_df$p_BH_lm) & comp_df$p_BH_lm < 0.05, ]
if (nrow(sig_rows) > 0) {
  print(sig_rows[, c("cell_type", "panel", "hedges_g", "direction",
                      "beta_lm", "p_BH_lm", "star_lm")])
} else {
  cat("  No panels reach FDR<0.05 — consistent with coupling-not-level narrative.\n")
}

# =============================================================================
# SECTION 7: INTERACTION TEST — CONDITION × CELL TYPE
# =============================================================================

cat("\n========== INTERACTION: CONDITION × CELL TYPE ==========\n")

# For each panel, fit combined model:
#   lm(score ~ covariates + condition * cell_type)
# The condition:cell_type interaction term tests whether PD-Control
# differs significantly between astrocytes and neurons.

interaction_long <- list()

for (pn in panel_names) {
  astro_df  <- panel_scores[["Astrocytes"]][[pn]]
  neuron_df <- panel_scores[["Neurons"]][[pn]]
  if (is.null(astro_df) || is.null(neuron_df)) next

  # Standardise and label
  astro_df$cond_std <- NA_character_
  astro_df$cond_std[grepl("PD|Parkinson|Disease", astro_df$condition, ignore.case = TRUE)] <- "PD"
  astro_df$cond_std[grepl("Control|Ctrl|Healthy|HC", astro_df$condition, ignore.case = TRUE)] <- "Control"
  astro_df$cell_type <- "Astrocytes"

  neuron_df$cond_std <- NA_character_
  neuron_df$cond_std[grepl("PD|Parkinson|Disease", neuron_df$condition, ignore.case = TRUE)] <- "PD"
  neuron_df$cond_std[grepl("Control|Ctrl|Healthy|HC", neuron_df$condition, ignore.case = TRUE)] <- "Control"
  neuron_df$cell_type <- "Neurons"

  combined <- rbind(
    astro_df[, c("donor_id", "score", "cond_std", "cell_type")],
    neuron_df[, c("donor_id", "score", "cond_std", "cell_type")]
  )
  combined <- combined[!is.na(combined$cond_std), ]
  combined$condition <- factor(combined$cond_std, levels = c("Control", "PD"))
  combined$cell_type <- factor(combined$cell_type, levels = c("Neurons", "Astrocytes"))

  if (nrow(combined) < 10) next

  # Get covariates (use astro object — same donors)
  cov_df_a <- get_donor_cov_df(astro, donor_col_astro,
                                astro_df$donor_id, covariate_cols_to_use)
  cov_df_n <- get_donor_cov_df(neuron, donor_col_neuron,
                                neuron_df$donor_id, covariate_cols_to_use)

  # Attach covariates by donor_id
  if (!is.null(cov_df_a) && nrow(cov_df_a) == nrow(astro_df)) {
    cov_a <- cbind(data.frame(donor_id = astro_df$donor_id), cov_df_a)
  } else { cov_a <- NULL }
  if (!is.null(cov_df_n) && nrow(cov_df_n) == nrow(neuron_df)) {
    cov_n <- cbind(data.frame(donor_id = neuron_df$donor_id), cov_df_n)
  } else { cov_n <- NULL }

  if (!is.null(cov_a) && !is.null(cov_n)) {
    cov_both <- rbind(cov_a, cov_n)
    lm_df <- merge(combined, cov_both, by = "donor_id", all.x = TRUE)
    cov_nms <- colnames(cov_df_a)
    fmla <- as.formula(paste0("score ~ ",
                              paste(cov_nms, collapse = " + "),
                              " + condition * cell_type"))
  } else {
    lm_df <- combined
    fmla <- as.formula("score ~ condition * cell_type")
  }

  fit <- tryCatch(lm(fmla, data = lm_df), error = function(e) NULL)
  if (is.null(fit)) next

  s <- summary(fit)$coefficients
  interaction_term <- "conditionPD:cell_typeAstrocytes"

  if (interaction_term %in% rownames(s)) {
    interaction_long[[length(interaction_long) + 1]] <- data.frame(
      panel          = pn,
      beta_interact  = round(s[interaction_term, "Estimate"], 6),
      se_interact    = round(s[interaction_term, "Std. Error"], 6),
      t_interact     = round(s[interaction_term, "t value"], 3),
      p_interact     = s[interaction_term, "Pr(>|t|)"],
      # Also pull main effects for reference
      beta_PD        = round(s["conditionPD", "Estimate"], 6),
      p_PD           = s["conditionPD", "Pr(>|t|)"],
      beta_astro     = round(s["cell_typeAstrocytes", "Estimate"], 6),
      p_astro        = s["cell_typeAstrocytes", "Pr(>|t|)"],
      n_total        = nrow(lm_df),
      stringsAsFactors = FALSE
    )
  }
}

if (length(interaction_long) > 0) {
  interact_df <- do.call(rbind, interaction_long)
  interact_df <- strip_rle(interact_df)
  interact_df$p_BH_interact <- p.adjust(interact_df$p_interact, method = "BH")
  interact_df$star           <- sig_star(interact_df$p_BH_interact)
  interact_df <- interact_df[order(interact_df$p_BH_interact), ]

  cat("Interaction results (top 5):\n")
  print(head(interact_df[, c("panel", "beta_interact", "t_interact",
                              "p_interact", "p_BH_interact", "star")], 5))
} else {
  interact_df <- data.frame()
  cat("WARNING: no interaction tests computed.\n")
}

# =============================================================================
# SECTION 8: GENE-LEVEL DESeq2-STYLE PD vs CONTROL (per cell type)
# =============================================================================

cat("\n========== GENE-LEVEL PD vs CONTROL (all immune panel genes) ==========\n")

# Simple donor-level mean expression → Wilcoxon + effect size per gene
# (Full DESeq2 is in Script 1; this is a focused supplementary view)

all_immune_genes <- unique(unlist(immune_subpanels))

gene_de_long <- list()

for (ct in list(
  list(name = "Astrocytes", obj = astro,  d = donor_col_astro,  c = cond_col_astro),
  list(name = "Neurons",    obj = neuron, d = donor_col_neuron, c = cond_col_neuron)
)) {
  genes_avail <- intersect(all_immune_genes, rownames(ct$obj))
  if (length(genes_avail) == 0) next

  cat(sprintf("--- %s: %d/%d immune genes found ---\n",
              ct$name, length(genes_avail), length(all_immune_genes)))

  expr <- GetAssayData(ct$obj, layer = "data")[genes_avail, , drop = FALSE]
  meta <- strip_rle(as.data.frame(ct$obj@meta.data))
  meta$donor_key <- as.character(meta[[ct$d]])
  meta$cond_key  <- as.character(meta[[ct$c]])

  # Donor-level mean per gene
  donors <- unique(meta$donor_key)
  gene_donor_mat <- matrix(NA, nrow = length(genes_avail), ncol = length(donors),
                           dimnames = list(genes_avail, donors))
  for (d in donors) {
    cells <- which(meta$donor_key == d)
    if (length(cells) > 0) {
      gene_donor_mat[, d] <- Matrix::rowMeans(expr[, cells, drop = FALSE])
    }
  }

  cond_lookup <- unique(meta[, c("donor_key", "cond_key")])
  donor_cond  <- cond_lookup$cond_key[match(donors, cond_lookup$donor_key)]

  # Standardise
  cond_std <- rep(NA_character_, length(donor_cond))
  cond_std[grepl("PD|Parkinson|Disease", donor_cond, ignore.case = TRUE)] <- "PD"
  cond_std[grepl("Control|Ctrl|Healthy|HC", donor_cond, ignore.case = TRUE)] <- "Control"

  idx_pd   <- which(cond_std == "PD")
  idx_ctrl <- which(cond_std == "Control")

  for (g in genes_avail) {
    x_pd   <- gene_donor_mat[g, idx_pd]
    x_ctrl <- gene_donor_mat[g, idx_ctrl]
    x_pd   <- x_pd[!is.na(x_pd)]
    x_ctrl <- x_ctrl[!is.na(x_ctrl)]
    if (length(x_pd) < 3 || length(x_ctrl) < 3) next

    wt <- tryCatch(wilcox.test(x_pd, x_ctrl, exact = FALSE),
                   error = function(e) list(p.value = NA))
    g_es <- hedges_g(x_pd, x_ctrl)

    # Identify which panel(s) this gene belongs to
    panels_for_gene <- names(immune_subpanels)[sapply(immune_subpanels,
                                                       function(pg) g %in% pg)]

    gene_de_long[[length(gene_de_long) + 1]] <- data.frame(
      cell_type  = ct$name,
      gene       = g,
      panel      = paste(panels_for_gene, collapse = "; "),
      mean_PD    = round(mean(x_pd), 5),
      mean_Ctrl  = round(mean(x_ctrl), 5),
      log2FC     = round(log2((mean(x_pd) + 1e-6) / (mean(x_ctrl) + 1e-6)), 4),
      hedges_g   = round(g_es, 4),
      p_wilcox   = wt$p.value,
      n_PD       = length(x_pd),
      n_Ctrl     = length(x_ctrl),
      stringsAsFactors = FALSE
    )
  }
}

if (length(gene_de_long) > 0) {
  gene_de_df <- do.call(rbind, gene_de_long)
  gene_de_df <- strip_rle(gene_de_df)
  gene_de_df$p_BH <- NA_real_
  for (ct in unique(gene_de_df$cell_type)) {
    idx <- gene_de_df$cell_type == ct
    gene_de_df$p_BH[idx] <- p.adjust(gene_de_df$p_wilcox[idx], method = "BH")
  }
  gene_de_df$star <- sig_star(gene_de_df$p_BH)
  gene_de_df <- gene_de_df[order(gene_de_df$cell_type, gene_de_df$p_BH), ]
  cat("Gene-level DE rows:", nrow(gene_de_df), "\n")
  cat("Sig at FDR<0.05:", sum(gene_de_df$p_BH < 0.05, na.rm = TRUE), "\n")
} else {
  gene_de_df <- data.frame()
  cat("WARNING: no gene-level DE computed.\n")
}

# =============================================================================
# SECTION 9: FIGURES
# =============================================================================

cat("\n========== GENERATING FIGURES ==========\n")

# ---- Collect all panel scores into a long data.frame for plotting ----
plot_long <- list()
for (ct_name in c("Astrocytes", "Neurons")) {
  for (pn in names(panel_scores[[ct_name]])) {
    df <- panel_scores[[ct_name]][[pn]]
    df$cond_std <- NA_character_
    df$cond_std[grepl("PD|Parkinson|Disease", df$condition, ignore.case = TRUE)] <- "PD"
    df$cond_std[grepl("Control|Ctrl|Healthy|HC", df$condition, ignore.case = TRUE)] <- "Control"
    plot_long[[length(plot_long) + 1]] <- data.frame(
      cell_type = ct_name,
      panel     = pn,
      donor_id  = df$donor_id,
      condition = df$cond_std,
      score     = df$score,
      stringsAsFactors = FALSE
    )
  }
}
plot_df <- do.call(rbind, plot_long)
plot_df <- plot_df[!is.na(plot_df$condition), ]
plot_df$panel     <- factor(plot_df$panel, levels = panel_names)
plot_df$condition <- factor(plot_df$condition, levels = c("Control", "PD"))

# ---- FIGURE 1: Violin + boxplot of panel scores by condition ----
cat("--- Figure 1: Violin plots ---\n")

# Merge significance annotations
annot_df <- comp_df[, c("cell_type", "panel", "p_BH_lm", "star_lm")]
annot_df$panel <- factor(annot_df$panel, levels = panel_names)
# Y-position for annotation
y_maxes <- aggregate(score ~ cell_type + panel, data = plot_df, FUN = max, na.rm = TRUE)
colnames(y_maxes)[3] <- "ymax"
annot_df <- merge(annot_df, y_maxes, by = c("cell_type", "panel"))
annot_df$y_pos <- annot_df$ymax * 1.1

fig1 <- ggplot(plot_df, aes(x = condition, y = score, fill = condition)) +
  geom_violin(alpha = 0.45, width = 0.8, color = NA) +
  geom_boxplot(width = 0.2, outlier.size = 1, alpha = 0.8) +
  geom_jitter(width = 0.08, size = 1.2, alpha = 0.5, color = "grey30") +
  geom_text(data = annot_df[annot_df$star_lm != "", ],
            aes(x = 1.5, y = y_pos, label = star_lm),
            size = 5, inherit.aes = FALSE, color = "black") +
  facet_grid(cell_type ~ panel, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#4DAF4A", PD = "#E41A1C")) +
  labs(title = paste0("Innate Immune Sub-Panel Scores: PD vs Control (", REGION, ")"),
       subtitle = "Background-subtracted donor-level scores | * FDR<0.05 ** FDR<0.01 *** FDR<0.001",
       x = "", y = "Panel score (bg-subtracted)") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "Fig1_Violin_PanelScores_ByCondition.pdf"),
       fig1, width = 16, height = 7)
ggsave(file.path(fig_dir, "Fig1_Violin_PanelScores_ByCondition.png"),
       fig1, width = 16, height = 7, dpi = 300)

# ---- FIGURE 2: Heatmap of Hedges' g effect sizes ----
cat("--- Figure 2: Effect size heatmap ---\n")

hm_df <- comp_df[, c("cell_type", "panel", "hedges_g", "star_lm")]
hm_df$panel <- factor(hm_df$panel, levels = rev(panel_names))

g_limit <- max(abs(hm_df$hedges_g), na.rm = TRUE)
g_limit <- min(g_limit, 2)  # cap for colour scale

fig2 <- ggplot(hm_df, aes(x = cell_type, y = panel, fill = hedges_g)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f%s", hedges_g, star_lm)),
            size = 3.8, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-g_limit, g_limit),
                       name = "Hedges' g\n(+PD > Ctrl)") +
  labs(title = paste0("Effect Size of PD vs Control (", REGION, ")"),
       subtitle = "Hedges' g | * FDR<0.05 (cov-adjusted lm)",
       x = "", y = "Immune sub-panel") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 10))

ggsave(file.path(fig_dir, "Fig2_Heatmap_EffectSize_Astro_vs_Neuron.pdf"),
       fig2, width = 6, height = 6)
ggsave(file.path(fig_dir, "Fig2_Heatmap_EffectSize_Astro_vs_Neuron.png"),
       fig2, width = 6, height = 6, dpi = 300)

# ---- FIGURE 3: Forest plot of covariate-adjusted betas ----
cat("--- Figure 3: Forest plot ---\n")

forest_df <- comp_df[, c("cell_type", "panel", "beta_lm", "se_lm",
                          "p_BH_lm", "star_lm")]
forest_df$ci_lo <- forest_df$beta_lm - 1.96 * forest_df$se_lm
forest_df$ci_hi <- forest_df$beta_lm + 1.96 * forest_df$se_lm
forest_df$sig   <- ifelse(!is.na(forest_df$p_BH_lm) & forest_df$p_BH_lm < 0.05,
                          "FDR<0.05", "n.s.")
forest_df$panel <- factor(forest_df$panel, levels = rev(panel_names))

fig3 <- ggplot(forest_df, aes(x = beta_lm, y = panel, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3,
                 color = "grey50") +
  geom_point(size = 3) +
  geom_text(aes(label = star_lm), hjust = -0.5, size = 4, show.legend = FALSE) +
  facet_wrap(~ cell_type, ncol = 2) +
  scale_color_manual(values = c("FDR<0.05" = "#E41A1C", "n.s." = "grey55")) +
  labs(title = paste0("PD vs Control: Covariate-Adjusted Beta (", REGION, ")"),
       subtitle = "lm(score ~ age + pmi + sex + condition) | 95% CI",
       x = "β (PD − Control)", y = "Immune sub-panel", color = "") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "Fig3_ForestPlot_CovAdj_Betas.pdf"),
       fig3, width = 11, height = 6)
ggsave(file.path(fig_dir, "Fig3_ForestPlot_CovAdj_Betas.png"),
       fig3, width = 11, height = 6, dpi = 300)

# ---- FIGURE 4: Interaction dot plot ----
cat("--- Figure 4: Interaction plot ---\n")

if (nrow(interact_df) > 0) {
  int_plot <- interact_df
  int_plot$panel <- factor(int_plot$panel, levels = rev(panel_names))
  int_plot$sig   <- ifelse(!is.na(int_plot$p_BH_interact) &
                             int_plot$p_BH_interact < 0.05,
                           "FDR<0.05", "n.s.")

  fig4 <- ggplot(int_plot, aes(x = beta_interact, y = panel, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(size = 4) +
    geom_text(aes(label = sprintf("p=%.3g%s",
                                  p_BH_interact, star)),
              hjust = -0.15, size = 3.2, show.legend = FALSE) +
    scale_color_manual(values = c("FDR<0.05" = "#E41A1C", "n.s." = "grey55")) +
    labs(title = paste0("Condition × Cell-Type Interaction (", REGION, ")"),
         subtitle = "β_interact > 0 ⇒ PD effect larger in astrocytes vs neurons",
         x = "Interaction β (conditionPD : cell_typeAstrocytes)",
         y = "Immune sub-panel", color = "") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  ggsave(file.path(fig_dir, "Fig4_Interaction_DotPlot.pdf"),
         fig4, width = 8, height = 5.5)
  ggsave(file.path(fig_dir, "Fig4_Interaction_DotPlot.png"),
         fig4, width = 8, height = 5.5, dpi = 300)
}

# ---- FIGURE 5: Gene-level volcano ----
cat("--- Figure 5: Gene-level volcano ---\n")

if (nrow(gene_de_df) > 0) {
  vol_df <- gene_de_df
  vol_df$neg_log10_p <- -log10(vol_df$p_BH + 1e-300)
  vol_df$sig <- ifelse(!is.na(vol_df$p_BH) & vol_df$p_BH < 0.05 &
                         abs(vol_df$hedges_g) > 0.3,
                       "Significant", "n.s.")
  vol_df$label <- ifelse(vol_df$sig == "Significant", vol_df$gene, "")

  fig5 <- ggplot(vol_df, aes(x = hedges_g, y = neg_log10_p,
                              color = sig, label = label)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_text_repel(size = 2.8, max.overlaps = 20, segment.size = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-0.3, 0.3), linetype = "dotted", color = "grey70") +
    facet_wrap(~ cell_type, scales = "free") +
    scale_color_manual(values = c("Significant" = "#E41A1C", "n.s." = "grey60")) +
    labs(title = paste0("Gene-Level PD vs Control (Immune Panel Genes, ", REGION, ")"),
         subtitle = "x = Hedges' g | y = −log10(FDR) | labelled: FDR<0.05 & |g|>0.3",
         x = "Hedges' g (PD − Control)", y = "-log10(FDR p)", color = "") +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(fig_dir, "Fig5_GeneLevel_Volcano.pdf"),
         fig5, width = 12, height = 6)
  ggsave(file.path(fig_dir, "Fig5_GeneLevel_Volcano.png"),
         fig5, width = 12, height = 6, dpi = 300)
}

# =============================================================================
# SECTION 10: SAVE EXCEL WORKBOOK
# =============================================================================

cat("\n========== SAVING WORKBOOK ==========\n")

wb <- createWorkbook()

# Sheet 1: Astrocytes
astro_comp <- comp_df[comp_df$cell_type == "Astrocytes", ]
astro_comp <- astro_comp[order(astro_comp$p_BH_lm), ]
addWorksheet(wb, "Astrocytes")
writeData(wb, "Astrocytes", astro_comp)

# Sheet 2: Neurons
neuron_comp <- comp_df[comp_df$cell_type == "Neurons", ]
neuron_comp <- neuron_comp[order(neuron_comp$p_BH_lm), ]
addWorksheet(wb, "Neurons")
writeData(wb, "Neurons", neuron_comp)

# Sheet 3: Interaction
if (nrow(interact_df) > 0) {
  addWorksheet(wb, "Interaction")
  writeData(wb, "Interaction", interact_df)
}

# Sheet 4: Gene-level DE
if (nrow(gene_de_df) > 0) {
  addWorksheet(wb, "GeneLevelDE")
  writeData(wb, "GeneLevelDE", gene_de_df)
}

# Sheet 5: Panel QC
addWorksheet(wb, "PanelSizeQC")
writeData(wb, "PanelSizeQC", panel_qc_df)

# Sheet 6: Methods
methods_note <- data.frame(
  Section = c("Purpose", "Panel_Scoring", "Primary_Test", "Secondary_Test",
              "Effect_Size", "Interaction", "GeneLevelDE",
              "Multiple_Testing", "Interpretation_Note"),
  Description = c(
    paste0("PD vs Control comparison of 9 innate immune sub-panel scores in ",
           "astrocytes and neurons. Region: ", REGION, "."),
    paste0("Background-subtracted scoring (AddModuleScore-style): panel mean minus ",
           "mean of expression-matched control genes (100 controls/gene, 25 bins). ",
           "Aggregated to donor-level mean."),
    paste0("Covariate-adjusted linear model: lm(score ~ ",
           if (COVARIATES_ENABLED && length(covariate_cols_to_use) > 0)
             paste(covariate_cols_to_use, collapse = " + ")
           else "no covariates",
           " + condition). Reports beta, SE, t, p for conditionPD."),
    "Wilcoxon rank-sum (unadjusted), provided as robustness check.",
    paste0("Hedges' g (bias-corrected Cohen's d). Positive = higher in PD. ",
           "Thresholds: |g|<0.2 negligible, 0.2-0.5 small, 0.5-0.8 medium, >0.8 large."),
    paste0("Condition x cell_type interaction: lm(score ~ covariates + condition * cell_type). ",
           "Tests whether PD-Control difference is cell-type-specific."),
    paste0("Per-gene Wilcoxon on donor-level mean expression for all immune panel genes. ",
           "Supplementary to DESeq2 (Script 1)."),
    "BH correction applied across 9 panels within each cell type.",
    paste0("NULL results are expected and informative: if panel-level scores ",
           "do not differ by condition, this supports the manuscript's coupling ",
           "narrative — the signal is in co-regulation structure, not marginal levels.")
  )
)
addWorksheet(wb, "Methods")
writeData(wb, "Methods", methods_note)

xlsx_filename <- paste0(REGION, "_PD_vs_Control_Immune_Panels.xlsx")
old_wd <- getwd()
setwd(results_dir)
saveWorkbook(wb, xlsx_filename, overwrite = TRUE)
setwd(old_wd)
cat("Workbook saved.\n")

# =============================================================================
# SECTION 11: SAVE RDS
# =============================================================================

saveRDS(
  list(
    region             = REGION,
    comparison_df      = comp_df,
    interaction_df     = interact_df,
    gene_de_df         = gene_de_df,
    panel_scores       = panel_scores,
    panel_qc           = panel_qc_df,
    immune_subpanels   = immune_subpanels,
    covariate_cols_used = covariate_cols_to_use
  ),
  file.path(results_dir, paste0(REGION, "_PD_vs_Ctrl_ImmunePanels_Results.rds"))
)
cat("RDS saved.\n")

# =============================================================================
# SECTION 12: SESSION SUMMARY
# =============================================================================

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Region:", REGION, "\n")
cat("Output directory:\n  ", results_dir, "\n\n")

cat("=== Tables (xlsx) ===\n")
cat("  Sheet 'Astrocytes'   — 9 panels, PD vs Ctrl, cov-adjusted\n")
cat("  Sheet 'Neurons'      — same for neurons\n")
cat("  Sheet 'Interaction'  — condition × cell_type interaction β\n")
cat("  Sheet 'GeneLevelDE'  — per-gene Wilcoxon + Hedges' g\n")
cat("  Sheet 'PanelSizeQC'  — genes found/missing per panel\n")
cat("  Sheet 'Methods'      — documentation\n\n")

cat("=== Figures ===\n")
cat("  Fig1 — Violin: panel scores by condition (facet: panel × cell_type)\n")
cat("  Fig2 — Heatmap: Hedges' g effect sizes (astro vs neuron)\n")
cat("  Fig3 — Forest: cov-adjusted β ± 95% CI\n")
cat("  Fig4 — Dot: condition × cell_type interaction\n")
cat("  Fig5 — Volcano: gene-level FDR vs effect size\n\n")

cat("=== Quick summary ===\n")
for (ct in c("Astrocytes", "Neurons")) {
  sub <- comp_df[comp_df$cell_type == ct, ]
  n_sig <- sum(!is.na(sub$p_BH_lm) & sub$p_BH_lm < 0.05)
  cat(sprintf("  %s: %d/%d panels FDR<0.05\n", ct, n_sig, nrow(sub)))
}
if (nrow(interact_df) > 0) {
  n_int_sig <- sum(!is.na(interact_df$p_BH_interact) &
                     interact_df$p_BH_interact < 0.05)
  cat(sprintf("  Interaction: %d/%d panels with cell-type-specific PD effect (FDR<0.05)\n",
              n_int_sig, nrow(interact_df)))
}

cat("\n=== Narrative guidance ===\n")
cat("If most panels are NULL: 'Innate immune pathway activity does not differ\n")
cat("  significantly between PD and control at the donor level, consistent with\n")
cat("  the coupling-not-level hypothesis — the disease signal resides in the\n")
cat("  co-regulation structure between dsRNA sensing and immune execution.'\n\n")
cat("If some panels are significant: these become the motivating context for\n")
cat("  the coupling analysis (Script 4), e.g., 'Complement is elevated in PD\n")
cat("  astrocytes; we next ask whether this elevation is coupled to dsRNA sensing.'\n")

cat("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")
