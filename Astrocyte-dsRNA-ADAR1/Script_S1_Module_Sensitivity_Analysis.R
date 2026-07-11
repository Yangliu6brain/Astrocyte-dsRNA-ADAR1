# ============================================================
# SCRIPT S1: MODULE SENSITIVITY ANALYSIS
# Supplementary Analysis
# ============================================================
#
# PURPOSE:
#   Address three questions for 8 framework:
#     Q1: How sensitive are key results to individual gene inclusion?
#     Q2: Do results replicate with canonical (KEGG/Reactome) gene sets?
#     Q3: Are observed couplings distinguishable from chance?
#
# ANALYSIS PARTS:
#   PART A: Leave-One-Out (LOO) — Complement panel
#           Remove each gene one at a time from the 9-gene Complement
#           panel; re-score, re-correlate, re-run Fisher Z.
#
#   PART B: Leave-One-Out (LOO) — dsRNA_Sensing module
#           Same procedure for the 8-gene dsRNA_Sensing module.
#
#   PART C: Leave-One-Out (LOO) — Innate_Immune module
#           Same procedure for the 21-gene Innate_Immune module.
#
#   PART D: Canonical Pathway Substitution
#           Replace custom panels with KEGG/Reactome/GO-derived gene sets.
#           Tests: (i) KEGG Complement & Coagulation (hsa04610 subset),
#                  (ii) Reactome Complement Cascade (R-HSA-166658),
#                  (iii) GO:0006958 (complement activation),
#                  (iv) Hallmark Complement (MSigDB).
#
#   PART E: Permutation Null Distribution
#           For each key module, draw N random gene sets of matching size
#           from the expressed-gene background, compute correlations,
#           and derive empirical p-values.
#
#   PART F: Module Overlap with Canonical Pathways
#           Jaccard index + hypergeometric enrichment of each custom
#           module against KEGG/Reactome/GO terms → shows systematic
#           derivation, not cherry-picking.
#
# INPUTS:
#   - 8-Module AllResults RDS (donor-level module scores from Script 1)
#   - Astrocyte + Neuron Seurat objects (with covariates)
#
# OUTPUTS:
#   Tables:
#     - Sensitivity_LOO_Complement.xlsx
#     - Sensitivity_LOO_dsRNA_Sensing.xlsx
#     - Sensitivity_LOO_Innate_Immune.xlsx
#     - Sensitivity_Canonical_Substitution.xlsx
#     - Sensitivity_Permutation_Results.xlsx
#     - Sensitivity_Pathway_Overlap.xlsx
#   Figures:
#     - FigS1_LOO_Complement_ForestPlot.pdf
#     - FigS2_LOO_dsRNA_Sensing_ForestPlot.pdf
#     - FigS3_Canonical_vs_Custom_BarPlot.pdf
#     - FigS4_Permutation_Null_Histogram.pdf
#     - FigS5_Jaccard_Heatmap.pdf
#
# METHODS:
#   All scoring, correlation, and testing procedures are identical
#   to Script 4 (background-subtracted AddModuleScore-style scoring,
#   partial Spearman with covariate residualization, Fisher Z for
#   cross-cell-type comparison).  Only the gene sets vary.
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
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

set.seed(2026)

# =============================================================================
# SECTION 1: CONFIGURATION — EDIT THIS SECTION
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("=== MODULE SENSITIVITY ANALYSIS (Supplementary)          ===\n")
cat("============================================================\n\n")

REGION <- "CX"

# --- Covariate settings (must match Script 1/4) ---
COVARIATES_ENABLED <- TRUE
COV_AGE_COL <- "age"
COV_PMI_COL <- "pmi"
COV_SEX_COL <- "sex"

# --- Permutation settings ---
N_PERMUTATIONS <- 1000   # Number of random gene sets for null distribution

library(here)
astro_file  <- here::here("data", "CX_Astrocytes_withCov.rds")
neuron_file <- here::here("data", "CX_Neurons_Total_withCov.rds")
results_rds <- here::here("data", "CX_8Module_AllResults.rds")
OUT_DIR     <- here::here("results")

timestamp   <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- file.path(output_base,
                         paste0(REGION, "_Sensitivity_", timestamp))
tbl_dir <- file.path(results_dir, "tables")
fig_dir <- file.path(results_dir, "figures")
dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("Region:", REGION, "\n")
cat("Output directory:", results_dir, "\n")
cat("Permutations:", N_PERMUTATIONS, "\n\n")

# =============================================================================
# SECTION 2: HELPER FUNCTIONS (identical to Script 4)
# =============================================================================

cat("========== DEFINING HELPER FUNCTIONS ==========\n")

# ---- Strip S4 Rle / list columns to base R types ----
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

# ---- Background-subtracted donor-level score (AddModuleScore-style) ----
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

# ---- Partial Spearman correlation with covariate residualization ----
partial_spearman <- function(x, y, covariate_df) {
  if (length(x) != length(y)) return(list(rho = NA, p = NA, n = 0))
  if (is.null(covariate_df) || ncol(covariate_df) == 0) {
    ok <- complete.cases(data.frame(x = x, y = y))
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
  unique_donors <- unique(meta$donor_key)
  donor_cov_list <- lapply(unique_donors, function(d) {
    sub <- meta[meta$donor_key == d, keep, drop = FALSE]
    vals <- sapply(keep, function(cn) {
      x <- sub[[cn]]
      if (is.numeric(x)) mean(x, na.rm = TRUE) else as.character(x[1])
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

# ---- Significance stars from p ----
sig_star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*", ""))))
}

# ---- Detect metadata column by trying several candidate names ----
detect_column <- function(meta, possible_names) {
  for (nm in possible_names) {
    if (nm %in% colnames(meta)) return(nm)
  }
  return(NULL)
}

# =============================================================================
# MASTER FUNCTION: Run one sensitivity iteration
# =============================================================================

run_one_sensitivity <- function(target_genes,         # genes for the panel being varied
                                ref_module_col,       # column name in RDS, e.g. "dsRNA_Sensing_Score"
                                donor_scores_a,       # RDS donor_scores for astrocytes (data.frame)
                                donor_scores_n,       # RDS donor_scores for neurons (data.frame)
                                astro_obj, neuron_obj,
                                donor_col_a, cond_col_a,
                                donor_col_n, cond_col_n,
                                cov_cols,
                                merge_donor_only = FALSE) {

  na_result <- data.frame(
    rho_astro = NA, p_astro = NA, n_astro = NA,
    rho_neuron = NA, p_neuron = NA, n_neuron = NA,
    delta_rho = NA, z_diff = NA, p_fisher = NA,
    genes_used_astro = NA, genes_used_neuron = NA,
    stringsAsFactors = FALSE
  )

  # --- Compute bg-subtracted panel scores (the thing being varied) ---
  panel_a <- donor_bg_subtracted_score(astro_obj,  target_genes, donor_col_a, cond_col_a)
  panel_n <- donor_bg_subtracted_score(neuron_obj, target_genes, donor_col_n, cond_col_n)

  if (is.null(panel_a) || is.null(panel_n)) return(na_result)

  # --- Merge panel scores with RDS module scores (exactly as Script 4) ---
  safe_merge_rds <- function(rds_df, panel_df) {
    rds_df$donor_id  <- trimws(as.character(rds_df$donor_id))
    rds_df$condition <- trimws(as.character(rds_df$condition))
    panel_df$donor_id  <- trimws(as.character(panel_df$donor_id))
    panel_df$condition <- trimws(as.character(panel_df$condition))

    if (merge_donor_only) {
      pdf_slim <- panel_df[, setdiff(colnames(panel_df), "condition"), drop = FALSE]
      m <- merge(rds_df, pdf_slim, by = "donor_id")
    } else {
      m <- merge(rds_df, panel_df, by = c("donor_id", "condition"))
      if (nrow(m) == 0) {
        # Fallback: condition labels differ between RDS and Seurat
        pdf_slim <- panel_df[, setdiff(colnames(panel_df), "condition"), drop = FALSE]
        m <- merge(rds_df, pdf_slim, by = "donor_id")
      }
    }
    m
  }

  merge_a <- safe_merge_rds(donor_scores_a, panel_a$donor_df)
  merge_n <- safe_merge_rds(donor_scores_n, panel_n$donor_df)

  if (nrow(merge_a) < 5 || nrow(merge_n) < 5) return(na_result)

  # Check that the module column exists in merged data
  if (!(ref_module_col %in% colnames(merge_a)) ||
      !(ref_module_col %in% colnames(merge_n))) {
    warning("Module column '", ref_module_col, "' not found in merged data.")
    return(na_result)
  }

  # --- Covariates ---
  cov_a <- get_donor_cov_df(astro_obj,  donor_col_a, merge_a$donor_id, cov_cols)
  cov_n <- get_donor_cov_df(neuron_obj, donor_col_n, merge_n$donor_id, cov_cols)

  # --- Partial Spearman: module_score ~ panel_score | covariates ---
  # "score" column = bg-subtracted panel score (from donor_bg_subtracted_score)
  # ref_module_col = pre-computed module score from RDS (Script 1 mean expression)
  ps_a <- partial_spearman(merge_a[[ref_module_col]], merge_a$score, cov_a)
  ps_n <- partial_spearman(merge_n[[ref_module_col]], merge_n$score, cov_n)

  # --- Fisher Z ---
  fz <- fisher_z_test(ps_a$rho, ps_a$n, ps_n$rho, ps_n$n)

  data.frame(
    rho_astro  = ps_a$rho,
    p_astro    = ps_a$p,
    n_astro    = ps_a$n,
    rho_neuron = ps_n$rho,
    p_neuron   = ps_n$p,
    n_neuron   = ps_n$n,
    delta_rho  = ifelse(is.na(ps_a$rho) | is.na(ps_n$rho), NA,
                        ps_a$rho - ps_n$rho),
    z_diff     = fz$z,
    p_fisher   = fz$p,
    genes_used_astro  = paste(panel_a$genes_used, collapse = ", "),
    genes_used_neuron = paste(panel_n$genes_used, collapse = ", "),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# SECTION 3: LOAD DATA
# =============================================================================

cat("========== LOADING DATA ==========\n")

prev_results <- readRDS(results_rds)
cat("Loaded 7-module results.\n")

# --- Seurat objects ---
astro  <- readRDS(astro_file)
neuron <- readRDS(neuron_file)
DefaultAssay(astro)  <- "RNA"
DefaultAssay(neuron) <- "RNA"
if (length(Layers(astro))  > 1) astro  <- JoinLayers(astro)
if (length(Layers(neuron)) > 1) neuron <- JoinLayers(neuron)

cat("Astrocytes:", ncol(astro),  "cells\n")
cat("Neurons:   ", ncol(neuron), "cells\n\n")

# --- Detect donor and condition columns ---
# MANUAL OVERRIDE: If auto-detection fails, set these directly.
# Set to NULL to use auto-detection; set to a string to force that column name.
MANUAL_DONOR_COL <- NULL    # e.g. "donor_id" or "subject_id"
MANUAL_COND_COL  <- NULL    # e.g. "condition" or "Status"

meta_a <- strip_rle(as.data.frame(astro@meta.data))
meta_n <- strip_rle(as.data.frame(neuron@meta.data))

# Print available columns for diagnostics
cat("--- Astrocyte metadata columns ---\n  ",
    paste(colnames(meta_a), collapse = ", "), "\n")
cat("--- Neuron metadata columns ---\n  ",
    paste(colnames(meta_n), collapse = ", "), "\n\n")

donor_possible <- c("donor_id", "Donor_ID", "donor", "Donor",
                     "sample_id", "Sample_ID", "subject_id", "individual",
                     "Sample", "sample", "orig.ident")
cond_possible  <- c("condition", "Condition", "disease", "Disease",
                     "Status", "Group", "group", "Diagnosis", "diagnosis",
                     "pathology", "Pathology")

if (!is.null(MANUAL_DONOR_COL)) {
  donor_col_astro  <- MANUAL_DONOR_COL
  donor_col_neuron <- MANUAL_DONOR_COL
  cat("Using manual donor column:", MANUAL_DONOR_COL, "\n")
} else {
  donor_col_astro  <- detect_column(meta_a, donor_possible)
  donor_col_neuron <- detect_column(meta_n, donor_possible)
}

if (!is.null(MANUAL_COND_COL)) {
  cond_col_astro  <- MANUAL_COND_COL
  cond_col_neuron <- MANUAL_COND_COL
  cat("Using manual condition column:", MANUAL_COND_COL, "\n")
} else {
  cond_col_astro  <- detect_column(meta_a, cond_possible)
  cond_col_neuron <- detect_column(meta_n, cond_possible)
}

# --- Validate: stop early with helpful message if columns not found ---
if (is.null(donor_col_astro) || is.null(donor_col_neuron)) {
  stop("Could not detect donor column.\n",
       "  Astrocyte columns: ", paste(colnames(meta_a), collapse = ", "), "\n",
       "  Neuron columns:    ", paste(colnames(meta_n), collapse = ", "), "\n",
       "  → Set MANUAL_DONOR_COL at top of Section 3 to the correct name.")
}
if (is.null(cond_col_astro) || is.null(cond_col_neuron)) {
  stop("Could not detect condition column.\n",
       "  Astrocyte columns: ", paste(colnames(meta_a), collapse = ", "), "\n",
       "  Neuron columns:    ", paste(colnames(meta_n), collapse = ", "), "\n",
       "  → Set MANUAL_COND_COL at top of Section 3 to the correct name.")
}

cat("Donor columns:   astro =", donor_col_astro, ", neuron =", donor_col_neuron, "\n")
cat("Condition columns: astro =", cond_col_astro, ", neuron =", cond_col_neuron, "\n")

# --- Covariates ---
covariate_cols <- c()
if (COVARIATES_ENABLED) {
  if (COV_AGE_COL %in% colnames(meta_a)) covariate_cols <- c(covariate_cols, COV_AGE_COL)
  if (COV_PMI_COL %in% colnames(meta_a)) covariate_cols <- c(covariate_cols, COV_PMI_COL)
  if (COV_SEX_COL %in% colnames(meta_a)) covariate_cols <- c(covariate_cols, COV_SEX_COL)
}
cat("Covariates:", paste(covariate_cols, collapse = ", "), "\n\n")

# --- Get expressed genes (background pool for permutations) ---
expressed_genes_astro  <- rownames(astro)[Matrix::rowMeans(GetAssayData(astro, layer = "data") > 0) > 0.01]
expressed_genes_neuron <- rownames(neuron)[Matrix::rowMeans(GetAssayData(neuron, layer = "data") > 0) > 0.01]
expressed_genes_shared <- intersect(expressed_genes_astro, expressed_genes_neuron)
cat("Expressed genes: astro =", length(expressed_genes_astro),
    ", neuron =", length(expressed_genes_neuron),
    ", shared =", length(expressed_genes_shared), "\n\n")


# =============================================================================
# SECTION 4: DEFINE GENE SETS
# =============================================================================

cat("========== DEFINING GENE SETS ==========\n")

# --- Custom modules (as in Script 1) ---
module_dsRNA_Sensing <- c("DDX58", "IFIH1", "DHX58", "TLR3",
                          "MAVS", "TICAM1", "TBK1", "IKBKE")

module_Innate_Immune <- c("JAK1", "TYK2", "STAT1", "STAT2", "IRF9", "RELA",
                          "ISG15", "USP18", "OAS1", "OAS2", "OAS3", "OASL",
                          "RNASEL", "EIF2AK2", "PRKRA",
                          "MX1", "IFIT1", "IFIT2", "IFIT3",
                          "C3", "C4A")

panel_Complement <- c("C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B",
                      "SERPING1", "CFB", "CFH")

# --- Canonical pathway gene sets for substitution (Part D) ---
# KEGG hsa04610: Complement and Coagulation Cascades (complement arm only)
kegg_complement <- c(
  "C1QA", "C1QB", "C1QC", "C1R", "C1S",
  "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
  "CFB", "CFD", "CFH", "CFI", "CFP",
  "MASP1", "MASP2", "MBL2",
  "SERPING1", "CD55", "CD59", "CR1", "CR2",
  "C3AR1", "C5AR1", "C5AR2",
  "ITGAM", "ITGAX", "ITGB2"
)

# Reactome R-HSA-166658: Complement Cascade
reactome_complement <- c(
  "C1QA", "C1QB", "C1QC", "C1R", "C1S",
  "C2", "C3", "C4A", "C4B", "C4BPA", "C4BPB",
  "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
  "CFB", "CFD", "CFH", "CFHR1", "CFHR2", "CFHR3", "CFHR4", "CFHR5",
  "CFI", "CFP",
  "MASP1", "MASP2", "MBL2", "FCN1", "FCN2", "FCN3",
  "SERPING1", "CD46", "CD55", "CD59",
  "CR1", "CR2", "C3AR1", "C5AR1"
)

# GO:0006958 complement activation
go_complement_activation <- c(
  "C1QA", "C1QB", "C1QC", "C1R", "C1S",
  "C2", "C3", "C4A", "C4B", "C5",
  "CFB", "CFD", "CFP",
  "MASP1", "MASP2", "MBL2", "FCN1", "FCN2",
  "CR1", "CD55", "SERPING1", "CFH", "CFI"
)

# MSigDB Hallmark Complement (representative subset)
hallmark_complement <- c(
  "C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2", "C3", "C4A", "C4B",
  "C5", "C6", "C7", "C8A", "C8B", "C9",
  "CFB", "CFD", "CFH", "CFI", "CFP",
  "SERPING1", "CD55", "CD59", "CR1",
  "MBL2", "MASP1", "MASP2",
  "C3AR1", "C5AR1"
)

# Canonical dsRNA sensing sets for comparison
kegg_rig_i <- c(
  "DDX58", "IFIH1", "DHX58",         # sensors
  "MAVS", "TRAF3", "TRAF6",          # adaptors
  "TBK1", "IKBKE", "IRF3", "IRF7",   # kinases/TFs
  "NFKB1", "NFKBIA",                 # NF-kB downstream
  "TRIM25", "ISG15", "RIPK1"         # regulators
)

# Reactome / GO cytosolic dsRNA sensing
reactome_dsRNA_sensing <- c(
  "DDX58", "IFIH1", "DHX58",
  "MAVS", "TRAF3", "TRAF6", "TANK",
  "TBK1", "IKBKE", "CHUK", "IKBKB",
  "IRF3", "IRF7", "NFKB1", "RELA",
  "TRIM25", "RNF135", "USP15"
)

canonical_sets <- list(
  Complement_Custom    = panel_Complement,
  Complement_KEGG      = kegg_complement,
  Complement_Reactome  = reactome_complement,
  Complement_GO        = go_complement_activation,
  Complement_Hallmark  = hallmark_complement,
  dsRNA_Sensing_Custom = module_dsRNA_Sensing,
  dsRNA_Sensing_KEGG_RIGI = kegg_rig_i,
  dsRNA_Sensing_Reactome  = reactome_dsRNA_sensing
)

# Print availability
cat("\n--- Gene Set Availability in Data ---\n")
for (nm in names(canonical_sets)) {
  genes <- canonical_sets[[nm]]
  avail <- sum(genes %in% expressed_genes_shared)
  cat(sprintf("  %-28s: %d / %d genes expressed\n", nm, avail, length(genes)))
}
cat("\n")



# =============================================================================
# SECTION 5: LOAD RDS MODULE SCORES
# =============================================================================

cat("========== LOADING RDS MODULE SCORES ==========\n")

donor_scores_astro  <- prev_results$donor_scores[["Astrocytes"]]
donor_scores_neuron <- prev_results$donor_scores[["Neurons_Total"]]

if (is.null(donor_scores_astro) || is.null(donor_scores_neuron)) {
  stop("Donor module scores not found in RDS. Check: names(prev_results$donor_scores)")
}

donor_scores_astro  <- strip_rle(donor_scores_astro)
donor_scores_neuron <- strip_rle(donor_scores_neuron)

# --- Detect module score columns ---
id_meta_cols <- c("donor_id", "condition", "Age", "age", "Sex", "sex",
                  "PMI", "pmi", "Braak", "braak", "nCount_RNA",
                  "log10_nCount", "n_cells")
module_cols <- setdiff(colnames(donor_scores_astro), id_meta_cols)
module_cols <- grep("_Score$|_score$|^Module|^ADAR|^SNCA|^dsRNA|^Innate|^TE_",
                    module_cols, value = TRUE)
module_cols <- intersect(module_cols, colnames(donor_scores_neuron))

cat("Module score columns from RDS (n =", length(module_cols), "):\n  ",
    paste(module_cols, collapse = ", "), "\n")

# --- Identify the key module columns we'll use as references ---
sensing_col  <- grep("Sensing", module_cols, value = TRUE)[1]
immune_col   <- grep("Innate_Immune|Innate", module_cols, value = TRUE)[1]

if (is.na(sensing_col)) stop("Cannot find dsRNA_Sensing_Score column in RDS.")
if (is.na(immune_col))  stop("Cannot find Innate_Immune_Score column in RDS.")

cat("Reference module for Complement LOO:    ", sensing_col, "\n")
cat("Reference module for dsRNA_Sensing LOO:  Complement panel (bg-subtracted)\n")
cat("Reference module for Innate_Immune LOO: ", sensing_col, "\n")

# --- Detect merge strategy (same logic as Script 4) ---
sample_panel <- donor_bg_subtracted_score(astro, panel_Complement,
                                          donor_col_astro, cond_col_astro)
rds_conds   <- unique(trimws(as.character(donor_scores_astro$condition)))
panel_conds <- unique(trimws(as.character(sample_panel$donor_df$condition)))
MERGE_DONOR_ONLY <- length(intersect(rds_conds, panel_conds)) == 0

if (MERGE_DONOR_ONLY) {
  cat("\n*** CONDITION MISMATCH: RDS uses [", paste(rds_conds, collapse=", "),
      "], Seurat uses [", paste(panel_conds, collapse=", "),
      "] → merging on donor_id only ***\n\n")
} else {
  cat("Conditions match between RDS and Seurat. Merging on donor_id + condition.\n\n")
}

# --- Baseline validation: reproduce Script 4's dsRNA_Sensing × Complement ---
cat("--- Baseline validation (should match Script 4) ---\n")
val_result <- run_one_sensitivity(
  target_genes    = panel_Complement,
  ref_module_col  = sensing_col,
  donor_scores_a  = donor_scores_astro,
  donor_scores_n  = donor_scores_neuron,
  astro_obj = astro, neuron_obj = neuron,
  donor_col_a = donor_col_astro, cond_col_a = cond_col_astro,
  donor_col_n = donor_col_neuron, cond_col_n = cond_col_neuron,
  cov_cols = covariate_cols,
  merge_donor_only = MERGE_DONOR_ONLY
)
cat(sprintf("  %s × Complement:\n", sensing_col))
cat(sprintf("    rho_astro  = %.6f  (Script 4 reference: 0.550136)\n", val_result$rho_astro))
cat(sprintf("    rho_neuron = %.6f  (Script 4 reference: -0.597101)\n", val_result$rho_neuron))
cat(sprintf("    Fisher Z p = %.6e  (Script 4 reference: 5.035e-09)\n", val_result$p_fisher))
if (!is.na(val_result$rho_astro) && abs(val_result$rho_astro - 0.550) > 0.05) {
  cat("  *** WARNING: rho_astro deviates from Script 4. Check scoring method. ***\n")
}
cat("\n")


# Function to compute Script-1-style module scores (simple colMeans, no bg subtraction)
script1_donor_score <- function(obj, genes, donor_col, cond_col) {
  avail <- genes[genes %in% rownames(obj)]
  if (length(avail) < 1) return(NULL)
  expr <- GetAssayData(obj, layer = "data")
  if (length(avail) == 1) {
    cell_score <- as.numeric(expr[avail, ])
  } else {
    cell_score <- Matrix::colMeans(expr[avail, , drop = FALSE])
  }
  meta <- strip_rle(as.data.frame(obj@meta.data))
  agg_df <- data.frame(
    donor_id  = as.character(meta[[donor_col]]),
    condition = as.character(meta[[cond_col]]),
    score     = as.numeric(cell_score),
    stringsAsFactors = FALSE
  )
  donor_df <- aggregate(score ~ donor_id + condition, data = agg_df,
                        FUN = mean, na.rm = TRUE)
  donor_df <- as.data.frame(donor_df, stringsAsFactors = FALSE)
  list(donor_df = donor_df, genes_used = avail)
}

cat("Reference scores ready.\n\n")


# =============================================================================
# PART A: LEAVE-ONE-OUT — COMPLEMENT PANEL
# =============================================================================
# Design: Fix RDS dsRNA_Sensing_Score, vary the Complement panel genes
# This directly mirrors Script 4's framework.

cat("========== PART A: LOO — COMPLEMENT ==========\n")

loo_complement_results <- list()

# --- Full model (baseline) ---
cat("  Running full Complement panel (baseline)...\n")
full_result <- run_one_sensitivity(
  target_genes    = panel_Complement,
  ref_module_col  = sensing_col,
  donor_scores_a  = donor_scores_astro,
  donor_scores_n  = donor_scores_neuron,
  astro_obj = astro, neuron_obj = neuron,
  donor_col_a = donor_col_astro, cond_col_a = cond_col_astro,
  donor_col_n = donor_col_neuron, cond_col_n = cond_col_neuron,
  cov_cols = covariate_cols,
  merge_donor_only = MERGE_DONOR_ONLY
)
full_result$removed_gene <- "None (full)"
full_result$n_genes      <- length(panel_Complement)
loo_complement_results[["Full"]] <- full_result

# --- LOO iterations ---
for (gene in panel_Complement) {
  reduced_set <- setdiff(panel_Complement, gene)
  cat(sprintf("  LOO: removing %-10s (%d genes remaining)...\n", gene, length(reduced_set)))
  res <- run_one_sensitivity(
    target_genes    = reduced_set,
    ref_module_col  = sensing_col,
    donor_scores_a  = donor_scores_astro,
    donor_scores_n  = donor_scores_neuron,
    astro_obj = astro, neuron_obj = neuron,
    donor_col_a = donor_col_astro, cond_col_a = cond_col_astro,
    donor_col_n = donor_col_neuron, cond_col_n = cond_col_neuron,
    cov_cols = covariate_cols,
    merge_donor_only = MERGE_DONOR_ONLY
  )
  res$removed_gene <- gene
  res$n_genes      <- length(reduced_set)
  loo_complement_results[[gene]] <- res
}

loo_compl_df <- do.call(rbind, loo_complement_results)
rownames(loo_compl_df) <- NULL

# --- Summary statistics ---
loo_compl_df$survives_005    <- !is.na(loo_compl_df$p_fisher) & loo_compl_df$p_fisher < 0.05
loo_compl_df$rho_change_pct  <- round(100 * (loo_compl_df$rho_astro -
                                              loo_compl_df$rho_astro[1]) /
                                       abs(loo_compl_df$rho_astro[1]), 1)

cat("\n--- LOO Complement Summary ---\n")
cat("Full model rho_astro:", round(loo_compl_df$rho_astro[1], 3),
    "  p_fisher:", signif(loo_compl_df$p_fisher[1], 3), "\n")
cat("Survived (p_fisher < 0.05):",
    sum(loo_compl_df$survives_005[-1], na.rm = TRUE), "/",
    nrow(loo_compl_df) - 1, "removals\n\n")

# --- Write table ---
wb_loo_c <- createWorkbook()
addWorksheet(wb_loo_c, "LOO_Complement")
writeData(wb_loo_c, 1, loo_compl_df)
addWorksheet(wb_loo_c, "Methods")
writeData(wb_loo_c, "Methods", data.frame(
  Item = c("Target_Panel", "Reference_Module", "Panel_Scoring", "Module_Scoring",
           "Correlation", "Comparison", "Covariates", "Interpretation"),
  Description = c(
    "Complement sub-panel (9 genes); each gene removed one at a time",
    paste0(sensing_col, " — held fixed from RDS (Script 1)"),
    "Background-subtracted AddModuleScore-style (100 ctrl genes/bin, 25 bins)",
    "Simple mean normalized expression (colMeans, Script 1) — from pre-computed RDS",
    "Partial Spearman (residualized on covariates)",
    "Fisher Z test: astrocyte rho vs neuron rho",
    paste(covariate_cols, collapse = ", "),
    "If p_fisher < 0.05 for all removals, no single gene drives the coupling"
  )
))
saveWorkbook(wb_loo_c,
             file.path(tbl_dir, "Sensitivity_LOO_Complement.xlsx"),
             overwrite = TRUE)


# =============================================================================
# PART B: LEAVE-ONE-OUT — dsRNA_Sensing MODULE
# =============================================================================

cat("========== PART B: LOO — dsRNA_Sensing ==========\n")

# Helper: Re-score module with reduced gene set, inject into donor_scores copy,
#         then run sensitivity against Complement panel.
run_module_loo <- function(module_genes, module_col_name,
                           panel_genes, astro_obj, neuron_obj,
                           donor_scores_a_orig, donor_scores_n_orig,
                           donor_col_a, cond_col_a,
                           donor_col_n, cond_col_n,
                           cov_cols, merge_donor_only) {

  # Re-score via Script 1 method (colMeans, no bg subtraction)
  new_a <- script1_donor_score(astro_obj,  module_genes, donor_col_a, cond_col_a)
  new_n <- script1_donor_score(neuron_obj, module_genes, donor_col_n, cond_col_n)

  if (is.null(new_a) || is.null(new_n)) {
    return(data.frame(rho_astro = NA, p_astro = NA, n_astro = NA,
                      rho_neuron = NA, p_neuron = NA, n_neuron = NA,
                      delta_rho = NA, z_diff = NA, p_fisher = NA,
                      genes_used_astro = NA, genes_used_neuron = NA,
                      stringsAsFactors = FALSE))
  }

  # Inject new module scores into copies of donor_scores
  ds_a <- donor_scores_a_orig
  ds_n <- donor_scores_n_orig
  ds_a$donor_id <- trimws(as.character(ds_a$donor_id))
  ds_n$donor_id <- trimws(as.character(ds_n$donor_id))
  new_a$donor_df$donor_id <- trimws(as.character(new_a$donor_df$donor_id))
  new_n$donor_df$donor_id <- trimws(as.character(new_n$donor_df$donor_id))

  # Replace module column
  ds_a[[module_col_name]] <- new_a$donor_df$score[match(ds_a$donor_id, new_a$donor_df$donor_id)]
  ds_n[[module_col_name]] <- new_n$donor_df$score[match(ds_n$donor_id, new_n$donor_df$donor_id)]

  # Run sensitivity with the modified donor_scores
  run_one_sensitivity(
    target_genes     = panel_genes,
    ref_module_col   = module_col_name,
    donor_scores_a   = ds_a,
    donor_scores_n   = ds_n,
    astro_obj = astro_obj, neuron_obj = neuron_obj,
    donor_col_a = donor_col_a, cond_col_a = cond_col_a,
    donor_col_n = donor_col_n, cond_col_n = cond_col_n,
    cov_cols = cov_cols,
    merge_donor_only = merge_donor_only
  )
}

loo_sensing_results <- list()

# --- Full model (baseline) ---
cat("  Running full dsRNA_Sensing module (baseline)...\n")
full_s <- run_module_loo(
  module_genes = module_dsRNA_Sensing,
  module_col_name = sensing_col,
  panel_genes = panel_Complement,
  astro_obj = astro, neuron_obj = neuron,
  donor_scores_a_orig = donor_scores_astro,
  donor_scores_n_orig = donor_scores_neuron,
  donor_col_a = donor_col_astro, cond_col_a = cond_col_astro,
  donor_col_n = donor_col_neuron, cond_col_n = cond_col_neuron,
  cov_cols = covariate_cols, merge_donor_only = MERGE_DONOR_ONLY
)
full_s$removed_gene <- "None (full)"
full_s$n_genes      <- length(module_dsRNA_Sensing)
loo_sensing_results[["Full"]] <- full_s

for (gene in module_dsRNA_Sensing) {
  reduced_set <- setdiff(module_dsRNA_Sensing, gene)
  cat(sprintf("  LOO: removing %-10s (%d genes remaining)...\n", gene, length(reduced_set)))
  res <- run_module_loo(
    module_genes = reduced_set,
    module_col_name = sensing_col,
    panel_genes = panel_Complement,
    astro_obj = astro, neuron_obj = neuron,
    donor_scores_a_orig = donor_scores_astro,
    donor_scores_n_orig = donor_scores_neuron,
    donor_col_a = donor_col_astro, cond_col_a = cond_col_astro,
    donor_col_n = donor_col_neuron, cond_col_n = cond_col_neuron,
    cov_cols = covariate_cols, merge_donor_only = MERGE_DONOR_ONLY
  )
  res$removed_gene <- gene
  res$n_genes      <- length(reduced_set)
  loo_sensing_results[[gene]] <- res
}

loo_sensing_df <- do.call(rbind, loo_sensing_results)
rownames(loo_sensing_df) <- NULL
loo_sensing_df$survives_005   <- !is.na(loo_sensing_df$p_fisher) & loo_sensing_df$p_fisher < 0.05
loo_sensing_df$rho_change_pct <- round(100 * (loo_sensing_df$rho_astro -
                                               loo_sensing_df$rho_astro[1]) /
                                        abs(loo_sensing_df$rho_astro[1]), 1)

cat("\n--- LOO dsRNA_Sensing Summary ---\n")
cat("Full model rho_astro:", round(loo_sensing_df$rho_astro[1], 3),
    "  p_fisher:", signif(loo_sensing_df$p_fisher[1], 3), "\n")
cat("Survived:", sum(loo_sensing_df$survives_005[-1], na.rm = TRUE),
    "/", nrow(loo_sensing_df) - 1, "\n\n")

wb_loo_s <- createWorkbook()
addWorksheet(wb_loo_s, "LOO_dsRNA_Sensing")
writeData(wb_loo_s, 1, loo_sensing_df)
saveWorkbook(wb_loo_s,
             file.path(tbl_dir, "Sensitivity_LOO_dsRNA_Sensing.xlsx"),
             overwrite = TRUE)


# =============================================================================
# PART C: LEAVE-ONE-OUT — INNATE_IMMUNE MODULE
# =============================================================================


cat("========== PART C: LOO — Innate_Immune ==========\n")

loo_immune_results <- list()

full_i <- run_module_loo(
  module_genes = module_Innate_Immune,
  module_col_name = immune_col,
  panel_genes = panel_Complement,
  astro_obj = astro, neuron_obj = neuron,
  donor_scores_a_orig = donor_scores_astro,
  donor_scores_n_orig = donor_scores_neuron,
  donor_col_a = donor_col_astro, cond_col_a = cond_col_astro,
  donor_col_n = donor_col_neuron, cond_col_n = cond_col_neuron,
  cov_cols = covariate_cols, merge_donor_only = MERGE_DONOR_ONLY
)
full_i$removed_gene <- "None (full)"
full_i$n_genes      <- length(module_Innate_Immune)
loo_immune_results[["Full"]] <- full_i

for (gene in module_Innate_Immune) {
  reduced_set <- setdiff(module_Innate_Immune, gene)
  cat(sprintf("  LOO: removing %-10s (%d genes remaining)...\n", gene, length(reduced_set)))
  res <- run_module_loo(
    module_genes = reduced_set,
    module_col_name = immune_col,
    panel_genes = panel_Complement,
    astro_obj = astro, neuron_obj = neuron,
    donor_scores_a_orig = donor_scores_astro,
    donor_scores_n_orig = donor_scores_neuron,
    donor_col_a = donor_col_astro, cond_col_a = cond_col_astro,
    donor_col_n = donor_col_neuron, cond_col_n = cond_col_neuron,
    cov_cols = covariate_cols, merge_donor_only = MERGE_DONOR_ONLY
  )
  res$removed_gene <- gene
  res$n_genes      <- length(reduced_set)
  loo_immune_results[[gene]] <- res
}

loo_immune_df <- do.call(rbind, loo_immune_results)
rownames(loo_immune_df) <- NULL
loo_immune_df$survives_005   <- !is.na(loo_immune_df$p_fisher) & loo_immune_df$p_fisher < 0.05
loo_immune_df$rho_change_pct <- round(100 * (loo_immune_df$rho_astro -
                                              loo_immune_df$rho_astro[1]) /
                                       abs(loo_immune_df$rho_astro[1]), 1)

cat("\n--- LOO Innate_Immune Summary ---\n")
cat("Full model rho_astro:", round(loo_immune_df$rho_astro[1], 3),
    "  p_fisher:", signif(loo_immune_df$p_fisher[1], 3), "\n")
cat("Survived:", sum(loo_immune_df$survives_005[-1], na.rm = TRUE),
    "/", nrow(loo_immune_df) - 1, "\n\n")

wb_loo_i <- createWorkbook()
addWorksheet(wb_loo_i, "LOO_Innate_Immune")
writeData(wb_loo_i, 1, loo_immune_df)
saveWorkbook(wb_loo_i,
             file.path(tbl_dir, "Sensitivity_LOO_Innate_Immune.xlsx"),
             overwrite = TRUE)


# =============================================================================
# PART D: CANONICAL PATHWAY SUBSTITUTION
# =============================================================================


cat("========== PART D: CANONICAL PATHWAY SUBSTITUTION ==========\n")
cat("Testing coupling with canonical pathway gene sets...\n\n")

canonical_results <- list()

for (set_name in names(canonical_sets)) {
  genes <- canonical_sets[[set_name]]
  avail <- genes[genes %in% expressed_genes_shared]
  cat(sprintf("  %-28s : %d / %d expressed ... ", set_name, length(avail), length(genes)))

  if (length(avail) < 3) {
    cat("SKIPPED (< 3 genes)\n")
    next
  }

  if (grepl("^Complement", set_name)) {
    # Complement variant: bg-subtracted panel, fixed RDS dsRNA_Sensing_Score
    res <- run_one_sensitivity(
      target_genes    = avail,
      ref_module_col  = sensing_col,
      donor_scores_a  = donor_scores_astro,
      donor_scores_n  = donor_scores_neuron,
      astro_obj = astro, neuron_obj = neuron,
      donor_col_a = donor_col_astro, cond_col_a = cond_col_astro,
      donor_col_n = donor_col_neuron, cond_col_n = cond_col_neuron,
      cov_cols = covariate_cols,
      merge_donor_only = MERGE_DONOR_ONLY
    )
  } else {
    # dsRNA_Sensing variant: re-score module with canonical genes,
    # correlate against Complement panel
    res <- run_module_loo(
      module_genes = avail,
      module_col_name = sensing_col,
      panel_genes = panel_Complement,
      astro_obj = astro, neuron_obj = neuron,
      donor_scores_a_orig = donor_scores_astro,
      donor_scores_n_orig = donor_scores_neuron,
      donor_col_a = donor_col_astro, cond_col_a = cond_col_astro,
      donor_col_n = donor_col_neuron, cond_col_n = cond_col_neuron,
      cov_cols = covariate_cols, merge_donor_only = MERGE_DONOR_ONLY
    )
  }

  res$gene_set     <- set_name
  res$n_genes_total <- length(genes)
  res$n_genes_avail <- length(avail)
  res$source       <- ifelse(grepl("Custom", set_name), "Custom",
                      ifelse(grepl("KEGG", set_name), "KEGG",
                      ifelse(grepl("Reactome", set_name), "Reactome",
                      ifelse(grepl("GO", set_name), "GO", "Hallmark"))))
  canonical_results[[set_name]] <- res
  cat(sprintf("rho_a = %.3f, p_fisher = %s\n",
              round(res$rho_astro, 3),
              ifelse(is.na(res$p_fisher), "NA", signif(res$p_fisher, 3))))
}

canonical_df <- do.call(rbind, canonical_results)
rownames(canonical_df) <- NULL
canonical_df$sig_star_fisher <- sig_star(canonical_df$p_fisher)

wb_canon <- createWorkbook()
addWorksheet(wb_canon, "Canonical_Substitution")
writeData(wb_canon, 1, canonical_df)
addWorksheet(wb_canon, "Methods")
writeData(wb_canon, "Methods", data.frame(
  Item = c("Design", "Panel_Scoring", "Module_Scoring",
           "KEGG_Source", "Reactome_Source", "GO_Source",
           "Hallmark_Source", "Interpretation"),
  Description = c(
    "Replace custom gene sets with canonical pathway annotations; run identical analysis pipeline",
    "Complement variants: bg-subtracted (matches Script 4 panel scoring)",
    "dsRNA_Sensing variants: simple colMeans (matches Script 1 module scoring)",
    "hsa04610 Complement and Coagulation Cascades (complement arm)",
    "R-HSA-166658 Complement Cascade",
    "GO:0006958 complement activation",
    "MSigDB Hallmark Complement",
    "If canonical sets reproduce the coupling, results are not an artifact of gene selection"
  )
))
saveWorkbook(wb_canon,
             file.path(tbl_dir, "Sensitivity_Canonical_Substitution.xlsx"),
             overwrite = TRUE)


# =============================================================================
# PART E: PERMUTATION NULL DISTRIBUTION
# =============================================================================


cat("========== PART E: PERMUTATION NULL ==========\n")
cat(sprintf("Running %d random gene sets (size-matched to Complement = %d genes)...\n",
            N_PERMUTATIONS, length(panel_Complement)))

# Exclude all custom module genes from the random pool to avoid overlap
all_module_genes <- unique(c(module_dsRNA_Sensing, module_Innate_Immune,
                             panel_Complement, "ADAR", "ADARB1", "SNCA"))
random_pool <- setdiff(expressed_genes_shared, all_module_genes)
target_size <- length(panel_Complement)

perm_rho_astro  <- numeric(N_PERMUTATIONS)
perm_rho_neuron <- numeric(N_PERMUTATIONS)
perm_z_diff     <- numeric(N_PERMUTATIONS)
perm_p_fisher   <- numeric(N_PERMUTATIONS)

cat("  Progress: ")
for (i in seq_len(N_PERMUTATIONS)) {
  if (i %% 100 == 0) cat(i, ".. ")

  rand_genes <- sample(random_pool, target_size, replace = FALSE)
  res <- run_one_sensitivity(
    target_genes    = rand_genes,
    ref_module_col  = sensing_col,
    donor_scores_a  = donor_scores_astro,
    donor_scores_n  = donor_scores_neuron,
    astro_obj = astro, neuron_obj = neuron,
    donor_col_a = donor_col_astro, cond_col_a = cond_col_astro,
    donor_col_n = donor_col_neuron, cond_col_n = cond_col_neuron,
    cov_cols = covariate_cols,
    merge_donor_only = MERGE_DONOR_ONLY
  )
  perm_rho_astro[i]  <- res$rho_astro
  perm_rho_neuron[i] <- res$rho_neuron
  perm_z_diff[i]     <- res$z_diff
  perm_p_fisher[i]   <- res$p_fisher
}
cat("Done.\n\n")

# --- Empirical p-values ---
observed_rho   <- loo_compl_df$rho_astro[1]    # Full model
observed_z     <- loo_compl_df$z_diff[1]

emp_p_rho  <- mean(abs(perm_rho_astro) >= abs(observed_rho), na.rm = TRUE)
emp_p_z    <- mean(abs(perm_z_diff) >= abs(observed_z), na.rm = TRUE)

cat("--- Permutation Results ---\n")
cat(sprintf("Observed rho_astro = %.3f  |  Empirical p = %.4f  (%.0f / %d permutations as extreme)\n",
            observed_rho, emp_p_rho, emp_p_rho * N_PERMUTATIONS, N_PERMUTATIONS))
cat(sprintf("Observed z_diff    = %.3f  |  Empirical p = %.4f\n",
            observed_z, emp_p_z))
cat(sprintf("Null distribution: rho_astro mean = %.3f, SD = %.3f\n",
            mean(perm_rho_astro, na.rm = TRUE), sd(perm_rho_astro, na.rm = TRUE)))
cat(sprintf("Null distribution: z_diff mean    = %.3f, SD = %.3f\n\n",
            mean(perm_z_diff, na.rm = TRUE), sd(perm_z_diff, na.rm = TRUE)))

# --- Write permutation table ---
perm_summary <- data.frame(
  Metric = c("Observed_rho_astro", "Null_mean_rho", "Null_SD_rho", "Empirical_p_rho",
             "Observed_z_diff", "Null_mean_z", "Null_SD_z", "Empirical_p_z",
             "N_permutations", "Random_pool_size", "Target_gene_set_size"),
  Value = c(observed_rho, mean(perm_rho_astro, na.rm = TRUE),
            sd(perm_rho_astro, na.rm = TRUE), emp_p_rho,
            observed_z, mean(perm_z_diff, na.rm = TRUE),
            sd(perm_z_diff, na.rm = TRUE), emp_p_z,
            N_PERMUTATIONS, length(random_pool), target_size)
)
wb_perm <- createWorkbook()
addWorksheet(wb_perm, "Summary")
writeData(wb_perm, 1, perm_summary)
addWorksheet(wb_perm, "NullDistribution")
writeData(wb_perm, "NullDistribution", data.frame(
  iteration = 1:N_PERMUTATIONS,
  rho_astro = perm_rho_astro,
  rho_neuron = perm_rho_neuron,
  z_diff = perm_z_diff,
  p_fisher = perm_p_fisher
))
saveWorkbook(wb_perm,
             file.path(tbl_dir, "Sensitivity_Permutation_Results.xlsx"),
             overwrite = TRUE)


# =============================================================================
# PART F: MODULE OVERLAP WITH CANONICAL PATHWAYS
# =============================================================================

cat("========== PART F: PATHWAY OVERLAP ANALYSIS ==========\n")

jaccard <- function(a, b) {
  int <- length(intersect(a, b))
  uni <- length(union(a, b))
  if (uni == 0) return(0)
  int / uni
}

hypergeom_test <- function(custom_genes, canonical_genes, background_n) {
  k <- length(intersect(custom_genes, canonical_genes))
  m <- length(canonical_genes)
  n <- background_n - m
  q <- length(custom_genes)
  p <- phyper(k - 1, m, n, q, lower.tail = FALSE)
  list(overlap = k, custom_n = q, canonical_n = m, bg_n = background_n, p = p)
}

custom_modules_for_overlap <- list(
  Complement_Panel  = panel_Complement,
  dsRNA_Sensing     = module_dsRNA_Sensing,
  Innate_Immune     = module_Innate_Immune
)

canonical_for_overlap <- list(
  KEGG_Complement      = kegg_complement,
  Reactome_Complement  = reactome_complement,
  GO_Complement        = go_complement_activation,
  KEGG_RIG_I           = kegg_rig_i,
  Reactome_dsRNA       = reactome_dsRNA_sensing
)

bg_n <- length(expressed_genes_shared)

overlap_results <- list()
idx <- 1
for (cust_name in names(custom_modules_for_overlap)) {
  cust_genes <- custom_modules_for_overlap[[cust_name]]
  for (canon_name in names(canonical_for_overlap)) {
    canon_genes <- canonical_for_overlap[[canon_name]]
    j <- jaccard(cust_genes, canon_genes)
    ht <- hypergeom_test(cust_genes, canon_genes, bg_n)
    overlap_results[[idx]] <- data.frame(
      Custom_Module  = cust_name,
      Canonical_Set  = canon_name,
      Jaccard        = round(j, 3),
      Overlap_Genes  = ht$overlap,
      Custom_Size    = ht$custom_n,
      Canonical_Size = ht$canonical_n,
      Hypergeom_p    = signif(ht$p, 4),
      Shared_Genes   = paste(intersect(cust_genes, canon_genes), collapse = ", "),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}

overlap_df <- do.call(rbind, overlap_results)
rownames(overlap_df) <- NULL

cat("--- Pathway Overlap Summary ---\n")
print(overlap_df[, c("Custom_Module", "Canonical_Set", "Jaccard",
                      "Overlap_Genes", "Custom_Size", "Canonical_Size", "Hypergeom_p")])
cat("\n")

wb_overlap <- createWorkbook()
addWorksheet(wb_overlap, "Overlap")
writeData(wb_overlap, 1, overlap_df)
saveWorkbook(wb_overlap,
             file.path(tbl_dir, "Sensitivity_Pathway_Overlap.xlsx"),
             overwrite = TRUE)


# =============================================================================
# SECTION 6: FIGURES
# =============================================================================

cat("========== GENERATING FIGURES ==========\n")

# ---- Fig S1: LOO Complement Forest Plot ----
cat("  FigS1: LOO Complement forest plot...\n")

plot_loo_forest <- function(loo_df, title_text, filename) {
  df <- loo_df
  df$label <- ifelse(df$removed_gene == "None (full)", "Full model",
                     paste0("−", df$removed_gene))
  df$label <- factor(df$label, levels = rev(df$label))
  df$is_full <- df$removed_gene == "None (full)"

  p <- ggplot(df, aes(x = rho_astro, y = label)) +
    geom_vline(xintercept = df$rho_astro[df$removed_gene == "None (full)"],
               linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_point(aes(color = survives_005, shape = is_full), size = 3) +
    geom_text(aes(label = sprintf("%.3f (p=%.3f)", rho_astro,
                                  ifelse(is.na(p_fisher), 1, p_fisher))),
              hjust = -0.1, size = 2.8) +
    scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D6604D"),
                       labels = c("TRUE" = "Survives (p<0.05)", "FALSE" = "Lost"),
                       name = "Fisher Z") +
    scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
    labs(title = title_text,
         subtitle = "Astrocyte rho (partial Spearman, covariate-adjusted)",
         x = expression(rho[astrocyte]),
         y = "Gene removed") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank()) +
    xlim(min(df$rho_astro, na.rm = TRUE) - 0.1,
         max(df$rho_astro, na.rm = TRUE) + 0.3)

  ggsave(file.path(fig_dir, filename), p,
         width = 7, height = max(4, nrow(df) * 0.35 + 1.5), dpi = 300)
}

plot_loo_forest(loo_compl_df,
                "LOO Sensitivity: Complement × dsRNA_Sensing Coupling",
                "FigS1_LOO_Complement_ForestPlot.pdf")

# ---- Fig S2: LOO dsRNA_Sensing Forest Plot ----
cat("  FigS2: LOO dsRNA_Sensing forest plot...\n")
plot_loo_forest(loo_sensing_df,
                "LOO Sensitivity: dsRNA_Sensing × Complement Coupling",
                "FigS2_LOO_dsRNA_Sensing_ForestPlot.pdf")

# ---- Fig S2b: LOO Innate_Immune Forest Plot ----
cat("  FigS2b: LOO Innate_Immune forest plot...\n")
plot_loo_forest(loo_immune_df,
                "LOO Sensitivity: Innate_Immune × dsRNA_Sensing Coupling",
                "FigS2b_LOO_Innate_Immune_ForestPlot.pdf")


# ---- Fig S3: Canonical vs Custom Bar Plot ----
cat("  FigS3: Canonical substitution bar plot...\n")

if (nrow(canonical_df) > 0) {
  canonical_df$gene_set <- factor(canonical_df$gene_set,
                                  levels = rev(canonical_df$gene_set))
  p3 <- ggplot(canonical_df, aes(x = rho_astro, y = gene_set, fill = source)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_vline(xintercept = 0, color = "grey50") +
    geom_text(aes(label = sprintf("p=%.3f", ifelse(is.na(p_fisher), 1, p_fisher))),
              hjust = ifelse(canonical_df$rho_astro >= 0, -0.1, 1.1), size = 2.8) +
    scale_fill_brewer(palette = "Set2", name = "Source") +
    labs(title = "Canonical Pathway Substitution",
         subtitle = "Astrocyte rho with dsRNA_Sensing (or Complement) module",
         x = expression(rho[astrocyte]), y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

  ggsave(file.path(fig_dir, "FigS3_Canonical_vs_Custom_BarPlot.pdf"), p3,
         width = 8, height = max(4, nrow(canonical_df) * 0.5 + 2), dpi = 300)
}


# ---- Fig S4: Permutation Null Histogram ----
cat("  FigS4: Permutation null histogram...\n")

perm_df <- data.frame(rho = perm_rho_astro)

p4a <- ggplot(perm_df, aes(x = rho)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey40", alpha = 0.8) +
  geom_vline(xintercept = observed_rho, color = "#D6604D", linewidth = 1.2,
             linetype = "solid") +
  annotate("text", x = observed_rho, y = Inf,
           label = sprintf("Observed\nrho = %.3f\nemp. p = %.4f",
                           observed_rho, emp_p_rho),
           vjust = 1.5, hjust = -0.1, color = "#D6604D", size = 3.5) +
  labs(title = "Permutation Null: Complement × dsRNA_Sensing",
       subtitle = sprintf("%d random gene sets (size = %d)", N_PERMUTATIONS, target_size),
       x = expression(rho[astrocyte] ~ "(random panel × dsRNA_Sensing)"),
       y = "Count") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p4b <- ggplot(data.frame(z = perm_z_diff), aes(x = z)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey40", alpha = 0.8) +
  geom_vline(xintercept = observed_z, color = "#D6604D", linewidth = 1.2) +
  annotate("text", x = observed_z, y = Inf,
           label = sprintf("Observed\nz = %.3f\nemp. p = %.4f",
                           observed_z, emp_p_z),
           vjust = 1.5, hjust = -0.1, color = "#D6604D", size = 3.5) +
  labs(title = "Permutation Null: Fisher Z (Astro vs Neuron)",
       x = expression(z[diff] ~ "(random panel)"),
       y = "Count") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "FigS4_Permutation_Null_Histogram.pdf"),
       p4a / p4b, width = 7, height = 9, dpi = 300)


# ---- Fig S5: Jaccard Heatmap ----
cat("  FigS5: Pathway overlap heatmap...\n")

if (nrow(overlap_df) > 0) {
  jaccard_wide <- overlap_df %>%
    select(Custom_Module, Canonical_Set, Jaccard) %>%
    pivot_wider(names_from = Canonical_Set, values_from = Jaccard,
                values_fill = 0) %>%
    column_to_rownames("Custom_Module") %>%
    as.matrix()

  pdf(file.path(fig_dir, "FigS5_Jaccard_Heatmap.pdf"), width = 8, height = 5)
  pheatmap(jaccard_wide,
           display_numbers = TRUE,
           number_format = "%.2f",
           color = colorRampPalette(c("white", "#2166AC"))(50),
           main = "Jaccard Overlap: Custom Modules vs Canonical Pathways",
           fontsize = 10,
           angle_col = 45,
           cluster_rows = FALSE,
           cluster_cols = FALSE)
  dev.off()
}


# =============================================================================
# SECTION 7: GRAND SUMMARY
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("=== SENSITIVITY ANALYSIS COMPLETE                        ===\n")
cat("============================================================\n\n")

cat("--- Results Summary ---\n\n")

cat("PART A (LOO Complement):\n")
cat(sprintf("  Full model:   rho = %.3f, Fisher Z p = %s\n",
            loo_compl_df$rho_astro[1],
            signif(loo_compl_df$p_fisher[1], 3)))
cat(sprintf("  Survived:     %d / %d gene removals (p < 0.05)\n",
            sum(loo_compl_df$survives_005[-1], na.rm = TRUE),
            nrow(loo_compl_df) - 1))
cat(sprintf("  Max rho change: %.1f%%\n\n",
            max(abs(loo_compl_df$rho_change_pct[-1]), na.rm = TRUE)))

cat("PART B (LOO dsRNA_Sensing):\n")
cat(sprintf("  Full model:   rho = %.3f, Fisher Z p = %s\n",
            loo_sensing_df$rho_astro[1],
            signif(loo_sensing_df$p_fisher[1], 3)))
cat(sprintf("  Survived:     %d / %d gene removals\n\n",
            sum(loo_sensing_df$survives_005[-1], na.rm = TRUE),
            nrow(loo_sensing_df) - 1))

cat("PART C (LOO Innate_Immune):\n")
cat(sprintf("  Full model:   rho = %.3f, Fisher Z p = %s\n",
            loo_immune_df$rho_astro[1],
            signif(loo_immune_df$p_fisher[1], 3)))
cat(sprintf("  Survived:     %d / %d gene removals\n\n",
            sum(loo_immune_df$survives_005[-1], na.rm = TRUE),
            nrow(loo_immune_df) - 1))

cat("PART D (Canonical Substitution):\n")
if (nrow(canonical_df) > 0) {
  for (i in seq_len(nrow(canonical_df))) {
    cat(sprintf("  %-28s: rho = %.3f, Fisher Z p = %s %s\n",
                canonical_df$gene_set[i],
                canonical_df$rho_astro[i],
                signif(canonical_df$p_fisher[i], 3),
                canonical_df$sig_star_fisher[i]))
  }
}

cat(sprintf("\nPART E (Permutation, n = %d):\n", N_PERMUTATIONS))
cat(sprintf("  Empirical p (rho):     %.4f\n", emp_p_rho))
cat(sprintf("  Empirical p (z_diff):  %.4f\n\n", emp_p_z))

cat("PART F (Pathway Overlap):\n")
cat("  See Jaccard heatmap and hypergeometric tests in tables.\n\n")

cat("Output directory:", results_dir, "\n")
cat("============================================================\n")

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "sessionInfo.txt"))

cat("\nDone.\n")
