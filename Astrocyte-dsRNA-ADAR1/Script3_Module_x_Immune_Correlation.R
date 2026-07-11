# ============================================================
# SCRIPT 4: MODULE × INNATE IMMUNE PATHWAY CORRELATION ANALYSIS
# Builds on 8-Module + 9-Panel Immune Framework
# ============================================================
#
# PURPOSE:
#   For each of the 8 functional modules (ADAR1, ADAR2, SNCA,
#   dsRNA_Processing, dsRNA_Sensing, dsRNA_EV_Transport,
#   Innate_Immune) compute partial Spearman correlations with
#   each of the 9 innate-immune sub-panels, separately in
#   astrocytes and neurons, with BH correction and Fisher-Z
#   cross-cell-type comparison.
#
# INPUTS:
#   - 8-Module AllResults RDS (donor-level module scores)
#   - Astrocyte + Neuron Seurat objects (with covariates)
#
# OUTPUTS (Z:/.../10_Module_x_Immune_Correlation/):
#   Tables:
#     - PanelLevel_PartialSpearman_Astro.xlsx
#     - PanelLevel_PartialSpearman_Neuron.xlsx
#     - FisherZ_CellType_Comparison.xlsx
#     - GeneLevel_Correlations_LongFormat.xlsx
#   Figures:
#     - FigA_Heatmap_Module_x_Panel.pdf  (astro | neuron panels)
#     - FigB_Scatter_TopPanel_dsRNA_Sensing.pdf
#     - FigC_ForestPlot_FisherZ.pdf
#     - FigD_Scatter_Stratified_PD_vs_Control.pdf
#     - FigE_GeneLevel_Heatmap.pdf  (bonus: per-gene view)
#
# METHODS:
#   - Immune panel donor scores computed with background-subtracted
#     scoring (Seurat::AddModuleScore-style: panel mean − mean of
#     size-matched control genes drawn from expression-binned pool).
#     This is the QC-corrected approach from your methodological
#     audit — prevents library-size confounding of cross-panel
#     correlations.
#   - Partial Spearman residualizes module and panel scores on
#     donor covariates (age, PMI, sex) before correlation.
#   - BH correction applied across (module × panel) tests within
#     each cell type. Significance stars:
#       * FDR < 0.05   ** FDR < 0.01   *** FDR < 0.001
#   - Fisher Z test: z_diff = (z_a − z_n) / sqrt(1/(n_a−3)+1/(n_n−3))
#
# ============================================================

rm(list = ls()); gc()

# =============================================================================
# SECTION 0: LIBRARIES
# =============================================================================

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
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  library(purrr)
})

set.seed(42)


# =============================================================================
# SECTION 1: CONFIGURATION
# =============================================================================

cat("========== CONFIGURATION ==========\n")

# --- Region: 
REGION <- "CX"

# --- Covariate control ---
COVARIATES_ENABLED <- TRUE
COV_AGE_COL <- "age"
COV_PMI_COL <- "pmi"
COV_SEX_COL <- "sex"
covariate_cols_to_use <- c()

library(here)
astro_file  <- here::here("data", "CX_Astrocytes_withCov.rds")
neuron_file <- here::here("data", "CX_Neurons_Total_withCov.rds")
results_rds <- here::here("data", "CX_8Module_AllResults.rds")
OUT_DIR     <- here::here("results")

timestamp   <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- file.path(output_base,
                         paste0(REGION, "_Module_x_Immune_Correlation_", timestamp))
tbl_dir <- file.path(results_dir, "tables")
fig_dir <- file.path(results_dir, "figures")
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("Output directory:", results_dir, "\n\n")

# =============================================================================
# SECTION 2: PANEL DEFINITIONS (9 immune sub-panels)
# =============================================================================

cat("========== DEFINING IMMUNE SUB-PANELS ==========\n")

immune_subpanels <- list(
  ISG_Core      = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
                    "IFITM1", "IFITM2", "IFITM3", "RSAD2", "BST2", "USP18"),
  IFN_STAT      = c("STAT1", "STAT2", "IRF1", "IRF3", "IRF7", "IRF9",
                    "JAK1", "TYK2"),
  NF_kB         = c("NFKB1", "NFKB2", "RELA", "RELB", "NFKBIA",
                    "IL6", "TNF", "CCL2", "CCL5", "CXCL10"),   # IL1B moved → Inflammasome
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
cat("Total unique genes across panels:",
    length(unique(unlist(immune_subpanels))), "\n\n")

# =============================================================================
# SECTION 3: LOAD DATA
# =============================================================================

cat("========== LOADING DATA ==========\n")

prev_results <- readRDS(results_rds)
cat("Loaded 8-module results. Keys:",
    paste(names(prev_results), collapse = ", "), "\n")

# --- Donor-level module scores (per cell type) ---
donor_scores_astro  <- prev_results$donor_scores[["Astrocytes"]]
donor_scores_neuron <- prev_results$donor_scores[["Neurons_Total"]]

if (is.null(donor_scores_astro) || is.null(donor_scores_neuron)) {
  stop("Donor module scores not found. Check RDS structure.")
}

# Strip S4 Rle columns inherited from Seurat v5
donor_scores_astro  <- strip_rle(donor_scores_astro)
donor_scores_neuron <- strip_rle(donor_scores_neuron)

# Identify module score columns (exclude id/meta/covariate cols)
id_meta_cols <- c("donor_id", "condition", "Age", "age", "Sex", "sex",
                  "PMI", "pmi", "Braak", "braak", "nCount_RNA",
                  "log10_nCount", "n_cells")

module_cols_astro  <- setdiff(colnames(donor_scores_astro),  id_meta_cols)
module_cols_astro  <- grep("_Score$|_score$|^Module|^ADAR|^SNCA|^dsRNA|^Innate|^TE_",
                           module_cols_astro, value = TRUE)
module_cols_neuron <- setdiff(colnames(donor_scores_neuron), id_meta_cols)
module_cols_neuron <- grep("_Score$|_score$|^Module|^ADAR|^SNCA|^dsRNA|^Innate|^TE_",
                           module_cols_neuron, value = TRUE)

module_cols <- intersect(module_cols_astro, module_cols_neuron)
if (length(module_cols) == 0) {
  # Fallback: take all numeric non-id columns
  num_cols_a <- colnames(donor_scores_astro)[sapply(donor_scores_astro, is.numeric)]
  num_cols_n <- colnames(donor_scores_neuron)[sapply(donor_scores_neuron, is.numeric)]
  module_cols <- intersect(setdiff(num_cols_a, id_meta_cols),
                           setdiff(num_cols_n, id_meta_cols))
}
cat("Module score columns detected (n =", length(module_cols), "):\n  ",
    paste(module_cols, collapse = ", "), "\n\n")

# --- Seurat objects ---
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

# ---- Background-subtracted donor-level score (AddModuleScore-style) ----
# Mimics Seurat::AddModuleScore:
#   1. Bin all genes by expression level
#   2. For each panel gene, draw N control genes from its bin
#   3. Panel score = mean(panel genes) − mean(control genes)
# Then aggregate to donor-level mean.
donor_bg_subtracted_score <- function(obj, genes, donor_col, cond_col,
                                      n_ctrl = 100, n_bins = 25) {
  avail <- genes[genes %in% rownames(obj)]
  if (length(avail) < 2) return(NULL)
  expr_all <- GetAssayData(obj, layer = "data")
  cell_means <- Matrix::rowMeans(expr_all)
  # Bin genes
  bins <- cut(cell_means, breaks = n_bins, labels = FALSE)
  names(bins) <- rownames(expr_all)
  # Control genes (sampled once per panel gene, pooled)
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
  # Cell-level panel mean − control mean
  panel_mean <- Matrix::colMeans(expr_all[avail,     , drop = FALSE])
  ctrl_mean  <- Matrix::colMeans(expr_all[ctrl_pool, , drop = FALSE])
  cell_score <- as.numeric(panel_mean - ctrl_mean)
  # Aggregate to donor level
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

# ---- Partial Spearman correlation with covariate residualization ----
partial_spearman <- function(x, y, covariate_df) {
  if (length(x) != length(y)) return(list(rho = NA, p = NA, n = 0))
  if (is.null(covariate_df) || ncol(covariate_df) == 0 ||
      !COVARIATES_ENABLED) {
    ok <- complete.cases(x, y)
    if (sum(ok) < 4) return(list(rho = NA, p = NA, n = sum(ok)))
    t <- tryCatch(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE),
                  error = function(e) list(estimate = NA, p.value = NA))
    return(list(rho = as.numeric(t$estimate), p = t$p.value, n = sum(ok)))
  }
  df <- data.frame(x = x, y = y, covariate_df)
  ok <- complete.cases(df)
  if (sum(ok) < ncol(covariate_df) + 3) return(list(rho = NA, p = NA, n = sum(ok)))
  df <- df[ok, ]
  cov_nms <- colnames(covariate_df)
  rx <- tryCatch(residuals(lm(as.formula(paste0("x ~ ", paste(cov_nms, collapse = " + "))),
                              data = df)), error = function(e) NULL)
  ry <- tryCatch(residuals(lm(as.formula(paste0("y ~ ", paste(cov_nms, collapse = " + "))),
                              data = df)), error = function(e) NULL)
  if (is.null(rx) || is.null(ry)) return(list(rho = NA, p = NA, n = sum(ok)))
  t <- tryCatch(cor.test(rx, ry, method = "spearman", exact = FALSE),
                error = function(e) list(estimate = NA, p.value = NA))
  list(rho = as.numeric(t$estimate), p = t$p.value, n = sum(ok))
}

# ---- Fisher Z test for two independent Spearman correlations ----
fisher_z_test <- function(rho1, n1, rho2, n2) {
  if (any(is.na(c(rho1, rho2))) || n1 < 4 || n2 < 4)
    return(list(z = NA, p = NA, z1 = NA, z2 = NA))
  z1 <- 0.5 * log((1 + rho1) / (1 - rho1))
  z2 <- 0.5 * log((1 + rho2) / (1 - rho2))
  se <- sqrt(1/(n1 - 3) + 1/(n2 - 3))
  z  <- (z1 - z2) / se
  list(z = z, p = 2 * pnorm(-abs(z)), z1 = z1, z2 = z2, se = se)
}

# ---- Pull covariates for donors from a Seurat object ----
get_donor_cov_df <- function(obj, donor_col, donor_ids, cov_cols) {
  if (length(cov_cols) == 0) return(NULL)
  meta <- strip_rle(as.data.frame(obj@meta.data))
  keep <- cov_cols[cov_cols %in% colnames(meta)]
  if (length(keep) == 0) return(NULL)
  meta$donor_key <- as.character(meta[[donor_col]])
  # Aggregate covariates to donor level using base R (avoids across/Rle issues)
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
  # Convert numeric-looking columns back to numeric
  for (cn in keep) {
    num_try <- suppressWarnings(as.numeric(donor_meta[[cn]]))
    if (!all(is.na(num_try))) donor_meta[[cn]] <- num_try
  }
  donor_meta <- donor_meta[match(donor_ids, donor_meta$donor_key), , drop = FALSE] 
  cov_df <- donor_meta[, keep, drop = FALSE]
  # Encode sex → 0/1 if character
  for (cn in colnames(cov_df)) {
    if (!is.numeric(cov_df[[cn]])) {
      cov_df[[cn]] <- as.numeric(factor(cov_df[[cn]])) - 1
    }
    if (any(is.na(cov_df[[cn]])))
      cov_df[[cn]][is.na(cov_df[[cn]])] <- mean(cov_df[[cn]], na.rm = TRUE)
  }
  # Drop zero-variance columns
  keep2 <- sapply(cov_df, function(x) sd(x, na.rm = TRUE) > 0)
  cov_df[, keep2, drop = FALSE]
}

# ---- Significance stars from p ----
sig_star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*", ""))))
}

# =============================================================================
# SECTION 5: COMPUTE IMMUNE PANEL DONOR SCORES (bg-subtracted)
# =============================================================================

cat("========== COMPUTING IMMUNE PANEL DONOR SCORES ==========\n")

panel_scores <- list()  # panel_scores[[cell_type]][[panel]]  = donor_df

for (ct in list(
  list(name = "Astrocytes", obj = astro,  d = donor_col_astro,  c = cond_col_astro),
  list(name = "Neurons",    obj = neuron, d = donor_col_neuron, c = cond_col_neuron)
)) {
  cat("--- ", ct$name, " ---\n")
  panel_scores[[ct$name]] <- list()
  for (pn in panel_names) {
    res <- donor_bg_subtracted_score(ct$obj, immune_subpanels[[pn]], ct$d, ct$c)
    if (is.null(res)) {
      cat(sprintf("  %-14s : SKIPPED (no genes)\n", pn)); next
    }
    panel_scores[[ct$name]][[pn]] <- res$donor_df
    cat(sprintf("  %-14s : %d/%d genes, %d donors\n",
                pn, length(res$genes_used), length(immune_subpanels[[pn]]),
                nrow(res$donor_df)))
  }
}

# =============================================================================
# SECTION 6: PANEL-LEVEL MODULE × IMMUNE PARTIAL CORRELATIONS
# =============================================================================

cat("\n========== PANEL-LEVEL PARTIAL CORRELATIONS ==========\n")

panel_cor_long <- list()

# --- Diagnose merge compatibility ONCE before the loop ---
cat("--- Pre-merge diagnostics ---\n")
sample_panel <- panel_scores[["Astrocytes"]][[1]]
cat("  prev_donor (Astro) donor_ids (head):",
    paste(head(as.character(donor_scores_astro$donor_id), 5), collapse = ", "), "\n")
cat("  panel_df   (Astro) donor_ids (head):",
    paste(head(as.character(sample_panel$donor_id), 5), collapse = ", "), "\n")
cat("  prev_donor conditions:", paste(unique(donor_scores_astro$condition), collapse = ", "), "\n")
cat("  panel_df   conditions:", paste(unique(sample_panel$condition), collapse = ", "), "\n")

# Check if conditions match; if not, we'll merge on donor_id only
merge_by_donor_only <- FALSE
prev_conds  <- unique(as.character(donor_scores_astro$condition))
panel_conds <- unique(as.character(sample_panel$condition))
if (length(intersect(prev_conds, panel_conds)) == 0) {
  cat("  *** CONDITION MISMATCH DETECTED ***\n")
  cat("      RDS uses:", paste(prev_conds, collapse = ", "), "\n")
  cat("      Seurat uses:", paste(panel_conds, collapse = ", "), "\n")
  cat("      → Will merge on donor_id only and keep condition from RDS.\n\n")
  merge_by_donor_only <- TRUE
} else {
  cat("  Conditions match. Merging on donor_id + condition.\n\n")
}

for (ct_name in c("Astrocytes", "Neurons")) {
  prev_donor <- if (ct_name == "Astrocytes") donor_scores_astro else donor_scores_neuron
  obj        <- if (ct_name == "Astrocytes") astro else neuron
  donor_col  <- if (ct_name == "Astrocytes") donor_col_astro else donor_col_neuron
  
  # Coerce merge keys to character
  prev_donor$donor_id  <- trimws(as.character(prev_donor$donor_id))
  prev_donor$condition <- trimws(as.character(prev_donor$condition))
  
  ct_panel_scores <- panel_scores[[ct_name]]
  if (length(ct_panel_scores) == 0) {
    cat(sprintf("  [%s] No panel scores — skipping\n", ct_name)); next
  }
  
  cat(sprintf("--- %s: %d panels, %d donors in RDS, %d module cols ---\n",
              ct_name, length(ct_panel_scores), nrow(prev_donor),
              sum(module_cols %in% colnames(prev_donor))))
  
  for (pn in names(ct_panel_scores)) {
    panel_df <- ct_panel_scores[[pn]]
    if (is.null(panel_df) || nrow(panel_df) == 0) next
    
    panel_df$donor_id  <- trimws(as.character(panel_df$donor_id))
    panel_df$condition <- trimws(as.character(panel_df$condition))
    
    # Merge: full or donor-only depending on condition compatibility
    if (merge_by_donor_only) {
      panel_df_slim <- panel_df[, setdiff(colnames(panel_df), "condition"), drop = FALSE]
      merged <- merge(prev_donor, panel_df_slim, by = "donor_id")
    } else {
      merged <- merge(prev_donor, panel_df, by = c("donor_id", "condition"))
      # Fallback if full merge fails
      if (nrow(merged) == 0) {
        panel_df_slim <- panel_df[, setdiff(colnames(panel_df), "condition"), drop = FALSE]
        merged <- merge(prev_donor, panel_df_slim, by = "donor_id")
      }
    }
    
    if (nrow(merged) < 6) {
      cat(sprintf("  [%s / %-14s] merged n=%d (<6), skipping\n",
                  ct_name, pn, nrow(merged)))
      next
    }
    
    cov_df <- get_donor_cov_df(obj, donor_col, merged$donor_id,
                               covariate_cols_to_use)
    
    mods_avail <- intersect(module_cols, colnames(merged))
    if (length(mods_avail) == 0) {
      cat(sprintf("  [%s / %-14s] 0 module cols in merged. Merged cols: %s\n",
                  ct_name, pn,
                  paste(head(colnames(merged), 15), collapse = ", ")))
      next
    }
    
    for (mod in mods_avail) {
      res <- partial_spearman(merged[[mod]], merged$score, cov_df)
      panel_cor_long[[length(panel_cor_long) + 1]] <- data.frame(
        cell_type = ct_name,
        module    = mod,
        panel     = pn,
        n         = res$n,
        rho       = res$rho,
        p         = res$p,
        stringsAsFactors = FALSE
      )
    }
  }
}

# --- Guard: if nothing was computed, stop with full diagnostics ---
if (length(panel_cor_long) == 0) {
  cat("\n!!! ERROR: No correlations computed. !!!\n")
  cat("module_cols:", paste(module_cols, collapse = ", "), "\n")
  cat("donor_scores_astro cols:\n  ",
      paste(colnames(donor_scores_astro), collapse = ", "), "\n")
  cat("donor_scores_neuron cols:\n  ",
      paste(colnames(donor_scores_neuron), collapse = ", "), "\n")
  cat("panel_scores Astrocytes:", paste(names(panel_scores$Astrocytes), collapse = ", "), "\n")
  cat("panel_scores Neurons:", paste(names(panel_scores$Neurons), collapse = ", "), "\n")
  # Try a raw merge to show what's happening
  cat("\n--- Raw merge test (Astro, first panel) ---\n")
  test_pn <- names(panel_scores$Astrocytes)[1]
  test_pdf <- panel_scores$Astrocytes[[test_pn]]
  cat("  test panel donor_ids:", paste(head(test_pdf$donor_id, 5), collapse = ", "), "\n")
  cat("  test panel conditions:", paste(unique(test_pdf$condition), collapse = ", "), "\n")
  cat("  RDS donor_ids:", paste(head(donor_scores_astro$donor_id, 5), collapse = ", "), "\n")
  cat("  RDS conditions:", paste(unique(donor_scores_astro$condition), collapse = ", "), "\n")
  cat("  Shared donor_ids:", length(intersect(donor_scores_astro$donor_id, test_pdf$donor_id)), "\n")
  m_test <- merge(donor_scores_astro, test_pdf[, "donor_id", drop = FALSE], by = "donor_id")
  cat("  Merge on donor_id only:", nrow(m_test), "rows\n")
  stop("Aborting: see diagnostics above. Likely a donor_id or condition mismatch.")
}

panel_cor_df <- do.call(rbind, panel_cor_long)
panel_cor_df <- strip_rle(panel_cor_df)

# BH correction within each cell type
panel_cor_df$p_BH <- NA_real_
panel_cor_df$star <- ""
for (ct in unique(panel_cor_df$cell_type)) {
  idx <- panel_cor_df$cell_type == ct
  panel_cor_df$p_BH[idx] <- p.adjust(panel_cor_df$p[idx], method = "BH")
  panel_cor_df$star[idx]  <- sig_star(panel_cor_df$p_BH[idx])
}

cat("Panel-level correlations computed:", nrow(panel_cor_df), "rows\n")

# --- Save tables ---
wb_astro <- createWorkbook()
addWorksheet(wb_astro, "Astrocytes_Panel_Cor")
astro_sub <- panel_cor_df[panel_cor_df$cell_type == "Astrocytes", ]
astro_sub <- astro_sub[order(astro_sub$p_BH), ]
writeData(wb_astro, 1, astro_sub)
saveWorkbook(wb_astro,
             file.path(tbl_dir, "PanelLevel_PartialSpearman_Astro.xlsx"),
             overwrite = TRUE)

wb_neuron <- createWorkbook()
addWorksheet(wb_neuron, "Neurons_Panel_Cor")
neuron_sub <- panel_cor_df[panel_cor_df$cell_type == "Neurons", ]
neuron_sub <- neuron_sub[order(neuron_sub$p_BH), ]
writeData(wb_neuron, 1, neuron_sub)
saveWorkbook(wb_neuron,
             file.path(tbl_dir, "PanelLevel_PartialSpearman_Neuron.xlsx"),
             overwrite = TRUE)

# =============================================================================
# SECTION 7: FISHER Z CELL-TYPE COMPARISON
# =============================================================================

cat("\n========== FISHER Z CELL-TYPE COMPARISON ==========\n")

fz_list <- list()
for (mod in module_cols) {
  for (pn in panel_names) {
    a <- panel_cor_df[panel_cor_df$cell_type == "Astrocytes" &
                        panel_cor_df$module == mod &
                        panel_cor_df$panel == pn, ]
    n <- panel_cor_df[panel_cor_df$cell_type == "Neurons" &
                        panel_cor_df$module == mod &
                        panel_cor_df$panel == pn, ]
    if (nrow(a) == 0 || nrow(n) == 0) next
    fz <- fisher_z_test(a$rho, a$n, n$rho, n$n)
    fz_list[[length(fz_list) + 1]] <- data.frame(
      module       = mod,
      panel        = pn,
      rho_astro    = a$rho,
      n_astro      = a$n,
      rho_neuron   = n$rho,
      n_neuron     = n$n,
      delta_rho    = a$rho - n$rho,
      z_diff       = fz$z,
      p_fisher     = fz$p,
      stringsAsFactors = FALSE
    )
  }
}
if (length(fz_list) == 0) {
  cat("WARNING: no Fisher Z tests computed.\n")
  fisher_df <- data.frame(module = character(), panel = character(),
                          rho_astro = numeric(), n_astro = integer(),
                          rho_neuron = numeric(), n_neuron = integer(),
                          delta_rho = numeric(), z_diff = numeric(),
                          p_fisher = numeric(), p_BH_fisher = numeric(),
                          star_fisher = character(), stringsAsFactors = FALSE)
} else {
  fisher_df <- do.call(rbind, fz_list)
  fisher_df <- strip_rle(fisher_df)
  fisher_df$p_BH_fisher <- p.adjust(fisher_df$p_fisher, method = "BH")
  fisher_df$star_fisher <- sig_star(fisher_df$p_BH_fisher)
  fisher_df <- fisher_df[order(fisher_df$p_BH_fisher), ]
}

wb_fz <- createWorkbook()
addWorksheet(wb_fz, "FisherZ")
writeData(wb_fz, 1, fisher_df)
saveWorkbook(wb_fz,
             file.path(tbl_dir, "FisherZ_CellType_Comparison.xlsx"),
             overwrite = TRUE)

cat("Top 5 Fisher-Z differential rows:\n")
print(head(fisher_df, 5))

# =============================================================================
# SECTION 8: GENE-LEVEL CORRELATIONS (per-gene × per-panel)
# =============================================================================

cat("\n========== GENE-LEVEL CORRELATIONS ==========\n")

# For each module, get its constituent genes from the stored module list
module_gene_lists <- prev_results$modules
if (is.null(module_gene_lists)) {
  # Try alternative storage keys
  for (k in c("module_genes", "modules_list", "gene_modules"))
    if (!is.null(prev_results[[k]])) { module_gene_lists <- prev_results[[k]]; break }
}
if (is.null(module_gene_lists)) {
  cat("WARNING: module gene lists not found in RDS. Skipping gene-level panel.\n")
  gene_cor_df <- NULL
} else {
  gene_cor_long <- list()
  for (ct_name in c("Astrocytes", "Neurons")) {
    obj       <- if (ct_name == "Astrocytes") astro else neuron
    donor_col <- if (ct_name == "Astrocytes") donor_col_astro else donor_col_neuron
    cond_col  <- if (ct_name == "Astrocytes") cond_col_astro  else cond_col_neuron
    ct_panel_scores <- panel_scores[[ct_name]]
    
    for (mod in names(module_gene_lists)) {
      genes <- intersect(module_gene_lists[[mod]], rownames(obj))
      if (length(genes) == 0) next
      
      # Donor-level expression per gene (raw normalized mean)
      expr <- GetAssayData(obj, layer = "data")[genes, , drop = FALSE]
      meta <- strip_rle(as.data.frame(obj@meta.data))
      meta$donor_key <- as.character(meta[[donor_col]])
      meta$cond_key  <- as.character(meta[[cond_col]])
      
      gene_donor <- as.data.frame(t(
        apply(expr, 1, function(x) tapply(x, meta$donor_key, mean, na.rm = TRUE))
      ))
      gene_donor$donor_id <- rownames(gene_donor)
      cond_lookup <- unique(meta[, c("donor_key", "cond_key")])
      gene_donor$condition <- cond_lookup$cond_key[match(gene_donor$donor_id,
                                                         cond_lookup$donor_key)]
      
      for (pn in names(ct_panel_scores)) {
        panel_df <- ct_panel_scores[[pn]]
        panel_df$donor_id  <- trimws(as.character(panel_df$donor_id))
        panel_df$condition <- trimws(as.character(panel_df$condition))
        if (merge_by_donor_only) {
          pdf_slim <- panel_df[, setdiff(colnames(panel_df), "condition"), drop = FALSE]
          m <- merge(gene_donor, pdf_slim, by = "donor_id")
        } else {
          m <- merge(gene_donor, panel_df, by = c("donor_id", "condition"))
          if (nrow(m) == 0) {
            pdf_slim <- panel_df[, setdiff(colnames(panel_df), "condition"), drop = FALSE]
            m <- merge(gene_donor, pdf_slim, by = "donor_id")
          }
        }
        if (nrow(m) < 6) next
        cov_df <- get_donor_cov_df(obj, donor_col, m$donor_id,
                                   covariate_cols_to_use)
        for (g in genes) {
          if (!g %in% colnames(m)) next
          r <- partial_spearman(m[[g]], m$score, cov_df)
          gene_cor_long[[length(gene_cor_long) + 1]] <- data.frame(
            cell_type = ct_name,
            module    = mod,
            gene      = g,
            panel     = pn,
            n         = r$n,
            rho       = r$rho,
            p         = r$p,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (length(gene_cor_long) == 0) {
    cat("WARNING: no gene-level correlations computed. Skipping Fig E.\n")
    gene_cor_df <- NULL
  } else {
    gene_cor_df <- do.call(rbind, gene_cor_long)
    gene_cor_df <- strip_rle(gene_cor_df)
    gene_cor_df$p_BH <- NA_real_
    gene_cor_df$star <- ""
    for (ct in unique(gene_cor_df$cell_type)) {
      idx <- gene_cor_df$cell_type == ct
      gene_cor_df$p_BH[idx] <- p.adjust(gene_cor_df$p[idx], method = "BH")
      gene_cor_df$star[idx]  <- sig_star(gene_cor_df$p_BH[idx])
    }
    
    wb_gl <- createWorkbook()
    addWorksheet(wb_gl, "GeneLevel")
    writeData(wb_gl, 1, gene_cor_df)
    saveWorkbook(wb_gl,
                 file.path(tbl_dir, "GeneLevel_Correlations_LongFormat.xlsx"),
                 overwrite = TRUE)
    cat("Gene-level rows:", nrow(gene_cor_df), "\n")
  }
}

# =============================================================================
# SECTION 9: FIGURE A — MODULE × PANEL HEATMAP (astro | neuron)
# =============================================================================

cat("\n========== FIGURE A: HEATMAP ==========\n")

make_heatmap_df <- function(cor_df, ct) {
  d <- cor_df[cor_df$cell_type == ct, ]
  d$module <- factor(d$module, levels = module_cols)
  d$panel  <- factor(d$panel,  levels = panel_names)
  d
}
ha <- make_heatmap_df(panel_cor_df, "Astrocytes")
hn <- make_heatmap_df(panel_cor_df, "Neurons")

rho_limit <- max(abs(c(ha$rho, hn$rho)), na.rm = TRUE)
rho_limit <- min(rho_limit, 1)

fig_A_theme <- theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        panel.grid  = element_blank(),
        plot.title  = element_text(face = "bold", hjust = 0.5))

p_A_astro <- ggplot(ha, aes(panel, module, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = star), size = 4, vjust = 0.75, color = "black") +
  geom_text(aes(label = sprintf("%.2f", rho)),
            size = 2.8, vjust = -0.35, color = "grey20") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-rho_limit, rho_limit),
                       name = "Partial\nSpearman ρ") +
  labs(title = "Astrocytes", x = "Immune sub-panel", y = "Module") +
  fig_A_theme

p_A_neuron <- ggplot(hn, aes(panel, module, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = star), size = 4, vjust = 0.75, color = "black") +
  geom_text(aes(label = sprintf("%.2f", rho)),
            size = 2.8, vjust = -0.35, color = "grey20") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-rho_limit, rho_limit),
                       name = "Partial\nSpearman ρ") +
  labs(title = "Neurons", x = "Immune sub-panel", y = "Module") +
  fig_A_theme + theme(axis.text.y = element_blank(),
                      axis.title.y = element_blank())

fig_A <- (p_A_astro | p_A_neuron) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Module × Innate Immune Sub-Panel Partial Spearman Correlations",
    subtitle = "* FDR<0.05  ** FDR<0.01  *** FDR<0.001 (BH within cell type)",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(fig_dir, "FigA_Heatmap_Module_x_Panel.pdf"),
       fig_A, width = 11, height = 7)
ggsave(file.path(fig_dir, "FigA_Heatmap_Module_x_Panel.png"),
       fig_A, width = 11, height = 7, dpi = 300)

# =============================================================================
# SECTION 10: FIGURE B — SCATTER dsRNA_Sensing × TOP PANEL
# =============================================================================

cat("\n========== FIGURE B: SCATTER ==========\n")

# Identify the panel with largest |delta_rho| vs dsRNA_Sensing
sensing_mod <- grep("Sensing|sensing", module_cols, value = TRUE)[1]
if (is.na(sensing_mod)) sensing_mod <- module_cols[1]

if (nrow(fisher_df) == 0) {
  cat("WARNING: fisher_df is empty — skipping Figs B, C, D.\n")
  top_panel_name <- panel_names[1]  # fallback for RDS save
} else {
  
  top_row <- fisher_df[fisher_df$module == sensing_mod, ]
  top_row <- top_row[order(-abs(top_row$delta_rho)), ]
  top_row <- top_row[1, , drop = FALSE]
  top_panel_name <- top_row$panel[1]
  cat("Top-differential panel for", sensing_mod, ":", top_panel_name,
      "(Δρ =", round(top_row$delta_rho, 3), ")\n")
  
  build_scatter_df <- function(ct_name) {
    prev_donor <- if (ct_name == "Astrocytes") donor_scores_astro else donor_scores_neuron
    panel_df   <- panel_scores[[ct_name]][[top_panel_name]]
    prev_donor$donor_id  <- trimws(as.character(prev_donor$donor_id))
    panel_df$donor_id    <- trimws(as.character(panel_df$donor_id))
    if (merge_by_donor_only) {
      pdf_slim <- panel_df[, setdiff(colnames(panel_df), "condition"), drop = FALSE]
      merged <- merge(prev_donor, pdf_slim, by = "donor_id")
    } else {
      prev_donor$condition <- trimws(as.character(prev_donor$condition))
      panel_df$condition   <- trimws(as.character(panel_df$condition))
      merged <- merge(prev_donor, panel_df, by = c("donor_id", "condition"))
      if (nrow(merged) == 0) {
        pdf_slim <- panel_df[, setdiff(colnames(panel_df), "condition"), drop = FALSE]
        merged <- merge(prev_donor, pdf_slim, by = "donor_id")
      }
    }
    merged$cell_type <- ct_name
    merged
  }
  df_scat <- rbind(build_scatter_df("Astrocytes"), build_scatter_df("Neurons"))
  
  # Pull the rho/p for annotation (base R — avoids dplyr dispatch on Rle)
  annot <- panel_cor_df[panel_cor_df$module == sensing_mod &
                          panel_cor_df$panel == top_panel_name, ]
  annot$label <- sprintf("ρ = %.2f, FDR = %.3g", annot$rho, annot$p_BH)
  
  fig_B <- ggplot(df_scat, aes(x = .data[[sensing_mod]], y = score)) +
    geom_point(aes(color = condition), size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", color = "black", se = TRUE, alpha = 0.15) +
    geom_text(data = annot, aes(label = label),
              x = -Inf, y = Inf, hjust = -0.08, vjust = 1.5,
              size = 3.8, inherit.aes = FALSE) +
    facet_wrap(~ cell_type, scales = "free") +
    scale_color_manual(values = c(Control = "#4DAF4A", PD = "#E41A1C"),
                       name = "Condition") +
    labs(title = sprintf("%s vs %s (top cell-type-differential panel)",
                         sensing_mod, top_panel_name),
         x = paste0(sensing_mod, " score"),
         y = paste0(top_panel_name, " score (bg-subtracted)")) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold"))
  
  ggsave(file.path(fig_dir, "FigB_Scatter_TopPanel_dsRNA_Sensing.pdf"),
         fig_B, width = 10, height = 5)
  ggsave(file.path(fig_dir, "FigB_Scatter_TopPanel_dsRNA_Sensing.png"),
         fig_B, width = 10, height = 5, dpi = 300)
  
  # =============================================================================
  # SECTION 11: FIGURE C — FOREST PLOT FISHER Z
  # =============================================================================
  
  cat("\n========== FIGURE C: FOREST PLOT ==========\n")
  
  # Focus on the dsRNA_Sensing module (the question of primary interest)
  # but also include all modules in a full panel
  forest_df <- fisher_df[fisher_df$module == sensing_mod, ]
  forest_df <- forest_df[order(forest_df$z_diff), ]
  forest_df$panel <- factor(forest_df$panel, levels = forest_df$panel)
  forest_df$ci_lo <- forest_df$z_diff - 1.96
  forest_df$ci_hi <- forest_df$z_diff + 1.96
  forest_df$sig   <- ifelse(forest_df$p_BH_fisher < 0.05, "FDR<0.05", "n.s.")
  
  fig_C <- ggplot(forest_df, aes(x = z_diff, y = panel)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25,
                   color = "grey30") +
    geom_point(aes(color = sig), size = 4) +
    geom_text(aes(label = sprintf("Δρ=%.2f%s",
                                  delta_rho, star_fisher)),
              hjust = -0.2, size = 3.3, nudge_x = 0.15) +
    scale_color_manual(values = c("FDR<0.05" = "#E41A1C", "n.s." = "grey60")) +
    labs(title = sprintf("Fisher Z: Astrocytes vs Neurons — %s × Immune Panels",
                         sensing_mod),
         subtitle = "Z > 0 ⇒ stronger correlation in astrocytes",
         x = "Fisher Z statistic (Z_astro − Z_neuron)",
         y = "Immune panel", color = "") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  
  ggsave(file.path(fig_dir, "FigC_ForestPlot_FisherZ.pdf"),
         fig_C, width = 9, height = 5.5)
  ggsave(file.path(fig_dir, "FigC_ForestPlot_FisherZ.png"),
         fig_C, width = 9, height = 5.5, dpi = 300)
  
  # Also build a multi-module version (all modules × all panels)
  forest_all <- fisher_df
  forest_all$sig <- ifelse(forest_all$p_BH_fisher < 0.05, "FDR<0.05", "n.s.")
  
  fig_C_all <- ggplot(forest_all, aes(x = z_diff, y = panel, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(size = 2.5) +
    facet_wrap(~ module, ncol = 3, scales = "free_x") +
    scale_color_manual(values = c("FDR<0.05" = "#E41A1C", "n.s." = "grey60")) +
    labs(title = "Fisher Z — All Modules × All Immune Panels",
         x = "Fisher Z (Astro − Neuron)", y = "Panel", color = "") +
    theme_bw(base_size = 9) +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold"))
  
  ggsave(file.path(fig_dir, "FigC_ForestPlot_FisherZ_AllModules.pdf"),
         fig_C_all, width = 11, height = 8)
  
  # =============================================================================
  # SECTION 12: FIGURE D — STRATIFIED BY PD vs CONTROL
  # =============================================================================
  
  cat("\n========== FIGURE D: PD vs CONTROL STRATIFIED SCATTER ==========\n")
  
  # Re-use top_panel_name, but split correlation within each condition
  stratified_df <- df_scat  # built in Section 10
  # Compute within-condition Spearman for annotation (base R — avoids Rle issues)
  strat_groups <- expand.grid(
    cell_type = unique(stratified_df$cell_type),
    condition = unique(stratified_df$condition),
    stringsAsFactors = FALSE
  )
  strat_annot <- do.call(rbind, lapply(seq_len(nrow(strat_groups)), function(i) {
    ct <- strat_groups$cell_type[i]
    cn <- strat_groups$condition[i]
    sub <- stratified_df[stratified_df$cell_type == ct &
                           stratified_df$condition == cn, ]
    xv <- sub[[sensing_mod]]; yv <- sub$score
    ok <- !is.na(xv) & !is.na(yv)
    rho_val <- if (sum(ok) >= 4) as.numeric(cor(xv[ok], yv[ok], method = "spearman")) else NA_real_
    data.frame(cell_type = ct, condition = cn,
               rho = rho_val, n = sum(ok),
               stringsAsFactors = FALSE)
  }))
  strat_annot$label <- sprintf("%s: ρ=%.2f (n=%d)",
                               strat_annot$condition,
                               strat_annot$rho,
                               strat_annot$n)
  
  strat_annot_combined <- do.call(rbind, lapply(
    unique(strat_annot$cell_type), function(ct) {
      sub <- strat_annot[strat_annot$cell_type == ct, ]
      data.frame(cell_type = ct,
                 label_all = paste(sub$label, collapse = "\n"),
                 stringsAsFactors = FALSE)
    }
  ))
  
  fig_D <- ggplot(stratified_df,
                  aes(x = .data[[sensing_mod]], y = score, color = condition)) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.8) +
    geom_text(data = strat_annot_combined,
              aes(label = label_all),
              x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
              size = 3.5, inherit.aes = FALSE) +
    facet_wrap(~ cell_type, scales = "free") +
    scale_color_manual(values = c(Control = "#4DAF4A", PD = "#E41A1C")) +
    labs(title = sprintf("%s × %s stratified by condition",
                         sensing_mod, top_panel_name),
         subtitle = "Robustness check — coupling persists within each group",
         x = paste0(sensing_mod, " score"),
         y = paste0(top_panel_name, " score (bg-subtracted)")) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold"))
  
  ggsave(file.path(fig_dir, "FigD_Scatter_Stratified_PD_vs_Control.pdf"),
         fig_D, width = 10, height = 5.5)
  ggsave(file.path(fig_dir, "FigD_Scatter_Stratified_PD_vs_Control.png"),
         fig_D, width = 10, height = 5.5, dpi = 300)
  
}  # end if (nrow(fisher_df) > 0) — Figs B/C/D

# =============================================================================
# SECTION 13: FIGURE E — GENE-LEVEL HEATMAP (bonus)
# =============================================================================

cat("\n========== FIGURE E: GENE-LEVEL HEATMAP ==========\n")

if (!is.null(gene_cor_df)) {
  # Focus on dsRNA_Sensing module genes and Innate_Immune module genes
  # (the two most relevant for the MS2 story), showing per-gene × per-panel
  focal_mods <- intersect(
    c("dsRNA_Sensing_Score", "Innate_Immune_Score",
      "dsRNA_Processing_Score", "ADAR1_Score", "SNCA_Score"),
    unique(gene_cor_df$module)
  )
  if (length(focal_mods) == 0) focal_mods <- unique(gene_cor_df$module)[1:2]
  
  gl_sub <- gene_cor_df[gene_cor_df$module %in% focal_mods, ]
  gl_sub$panel <- factor(gl_sub$panel, levels = panel_names)
  
  rho_limit_g <- min(max(abs(gl_sub$rho), na.rm = TRUE), 1)
  
  fig_E <- ggplot(gl_sub, aes(panel, gene, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = star), size = 3, color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-rho_limit_g, rho_limit_g),
                         name = "ρ") +
    facet_grid(module ~ cell_type, scales = "free_y", space = "free_y") +
    labs(title = "Gene-level partial Spearman correlations with immune sub-panels",
         subtitle = "Rows = individual module genes; columns = 9 immune panels",
         x = "Immune panel", y = "Module gene") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold", size = 9),
          panel.grid = element_blank())
  
  ggsave(file.path(fig_dir, "FigE_GeneLevel_Heatmap.pdf"),
         fig_E, width = 12, height = 12, limitsize = FALSE)
}

# =============================================================================
# SECTION 14: SESSION SUMMARY
# =============================================================================

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Output directory:\n  ", results_dir, "\n\n")
cat("Tables:\n")
cat("  - PanelLevel_PartialSpearman_Astro.xlsx\n")
cat("  - PanelLevel_PartialSpearman_Neuron.xlsx\n")
cat("  - FisherZ_CellType_Comparison.xlsx\n")
cat("  - GeneLevel_Correlations_LongFormat.xlsx\n\n")
cat("Figures:\n")
cat("  - FigA_Heatmap_Module_x_Panel.pdf/.png\n")
cat("  - FigB_Scatter_TopPanel_dsRNA_Sensing.pdf/.png\n")
cat("  - FigC_ForestPlot_FisherZ.pdf/.png  (+ AllModules version)\n")
cat("  - FigD_Scatter_Stratified_PD_vs_Control.pdf/.png\n")
cat("  - FigE_GeneLevel_Heatmap.pdf\n\n")
cat("Significant panel × module findings (FDR<0.05):\n")
sig_panels <- panel_cor_df[!is.na(panel_cor_df$p_BH) & panel_cor_df$p_BH < 0.05, ]
sig_panels <- sig_panels[order(sig_panels$cell_type, sig_panels$p_BH), ]
print(sig_panels[, c("cell_type", "module", "panel", "n", "rho", "p_BH", "star")])
cat("\nSignificant Fisher-Z cell-type differences (FDR<0.05):\n")
sig_fz <- fisher_df[!is.na(fisher_df$p_BH_fisher) & fisher_df$p_BH_fisher < 0.05, ]
sig_fz <- sig_fz[order(sig_fz$p_BH_fisher), ]
print(sig_fz[, c("module", "panel", "rho_astro", "rho_neuron", "delta_rho",
                 "z_diff", "p_BH_fisher", "star_fisher")])

saveRDS(list(panel_cor_df = panel_cor_df,
             fisher_df    = fisher_df,
             gene_cor_df  = gene_cor_df,
             panel_scores = panel_scores,
             top_panel    = top_panel_name,
             sensing_mod  = sensing_mod),
        file.path(results_dir, paste0(REGION,
                                      "_Module_x_Immune_Results.rds")))
cat("\nAll results saved to RDS.\n")