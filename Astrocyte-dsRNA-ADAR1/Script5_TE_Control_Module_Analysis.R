# ---------------------------------------------------------------------------------------------------------------
# Repository copy. Run from the repository root. Inputs: data/ (covariate-annotated Seurat objects from
# Full_merge_RDS_script_with_covariates.R) and results/<analysis output folder> (rename each timestamped output
# folder to the name used below, e.g. CX_TE_Control_<date> -> results/CX_TE_Control). R 4.4.3.
# ---------------------------------------------------------------------------------------------------------------
# ============================================================
# SCRIPT 5: TE_CONTROL MODULE — CONSOLIDATED ANALYSIS
# Adds TE_Control (Module 8) to the full analytical pipeline
# ============================================================
#
# PURPOSE:
#   This script runs every analysis that was previously done for
#   the 7-module framework, but specifically for the 8th module
#   (TE_Control: MOV10, TRIM28, SETDB1, DNMT1, DICER1, AGO2,
#   L1TD1, PIWIL1). It consolidates the essential TE_Control
#   parts from Scripts 1–4 + SNCA Trajectory into one run.
#
# ANALYSES INCLUDED:
#   PART A: DESeq2 pseudobulk PD vs Control (astro + neuron subtypes)
#   PART B: Module score comparisons — PD vs Control (covariate-adjusted)
#   PART C: Cell-type comparisons — Astrocyte vs Neuron (paired Wilcoxon)
#   PART D: Cross-module correlations (TE_Control vs all 7 other modules)
#           — ALL, PD-only, Control-only; partial Spearman with covariates
#   PART E: Module × Immune sub-panel correlations (Script 4 logic)
#           — TE_Control vs 9 immune panels + Fisher Z astro vs neuron
#   PART F: SNCA trajectory coupling (TE_Control along SNCA axis)
#   PART G: Summary workbook
#
# INPUTS:
#   - Astrocyte Seurat RDS (with covariates)
#   - Neuron Seurat RDS (with covariates)
#
# OUTPUTS:
#   - Excel workbook: {REGION}_TE_Control_Comprehensive.xlsx
#   - CSV tables per analysis
#   - PDF/PNG figures (heatmaps, scatters, boxplots)
#   - RDS with all results
#
# NOTE: This is a STANDALONE script. It does NOT require
#       the 7-Module AllResults RDS. It recomputes all 8 module
#       scores from scratch so that TE_Control is fully integrated.
#
# ============================================================

rm(list = ls()); gc()

# =============================================================================
# SECTION 0: Libraries
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
  library(ggpubr)
  library(corrplot)
  library(scales)
})

suppressPackageStartupMessages(library(DESeq2))
cat("DESeq2 loaded.\n")

set.seed(42)

# =============================================================================
# SECTION 1: USER SETTINGS
# =============================================================================

cat("\n============================================================\n")
cat("=== SCRIPT 5: TE_CONTROL MODULE ANALYSIS ===\n")
cat("============================================================\n\n")

# ---- REGION ----
REGION <- "CX"   # Change to "SN" for substantia nigra

# ---- COVARIATE SETTINGS ----
COVARIATES_ENABLED <- TRUE
USE_COV_AGE <- TRUE
USE_COV_PMI <- TRUE
USE_COV_SEX <- TRUE
COV_AGE_COL <- "age"
COV_PMI_COL <- "pmi"
COV_SEX_COL <- "sex"

# ---- INPUT FILES ----
INPUT_METHOD <- "path"   # "dialog" or "path"

if (INPUT_METHOD == "path") {
  if (REGION == "CX") {
    astro_file  <- "data/CX_Astrocytes_withCov.rds"
    neuron_file <- "data/CX_Neurons_Total_withCov.rds"
    output_base <- "results"
  } else {
    astro_file  <- "data/SN_Astrocytes.rds"
    neuron_file <- "data/SN_Neurons_Total.rds"
    output_base <- "results"
  }
} else {
  library(tcltk)
  astro_file  <- tk_choose.files(caption = "Select ASTROCYTE RDS",
                                 filters = matrix(c("RDS", ".rds"), 1, 2))
  neuron_file <- tk_choose.files(caption = "Select NEURON RDS",
                                 filters = matrix(c("RDS", ".rds"), 1, 2))
  output_base <- tk_choose.dir(caption = "Select OUTPUT folder")
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
results_main <- file.path(output_base,
                          paste0(REGION, "_TE_Control_", timestamp))
dir.create(results_main, showWarnings = FALSE, recursive = TRUE)
cat("Output:", results_main, "\n\n")

# =============================================================================
# SECTION 2: MODULE DEFINITIONS — ALL 8 MODULES
# =============================================================================

cat("========== DEFINING 8 MODULES ==========\n")

single_gene_modules <- list(
  ADAR1 = "ADAR",
  ADAR2 = "ADARB1",
  SNCA  = "SNCA"
)

pathway_modules <- list(
  dsRNA_Processing = c(
    "ILF2", "ILF3", "STRBP", "ZFR", "DHX9",
    "DDX39B", "SFPQ", "NONO"
  ),
  dsRNA_EV_Transport = c(
    "SIDT1", "SIDT2", "LAMP2", "RAB27A", "RAB27B", "RAB35",
    "SMPD3", "CD63", "CD9", "CD81", "PDCD6IP", "TSG101",
    "SDCBP", "MVP"
  ),
  dsRNA_Sensing = c(
    "DDX58", "IFIH1", "DHX58",
    "TLR3",
    "MAVS", "TICAM1",
    "TBK1", "IKBKE"
  ),
  Innate_Immune = c(
    "JAK1", "TYK2", "STAT1", "STAT2", "IRF9", "RELA",
    "ISG15", "USP18", "OAS1", "OAS2", "OAS3", "OASL",
    "RNASEL", "EIF2AK2", "PRKRA",
    "MX1", "IFIT1", "IFIT2", "IFIT3",
    "C3", "C4A"
  ),
  # MODULE 8: TE CONTROL — THE STAR OF THIS SCRIPT
  TE_Control = c(
    "MOV10", "TRIM28", "SETDB1", "DNMT1",
    "DICER1", "AGO2", "L1TD1", "PIWIL1"
  )
)

module_gene_sets <- c(single_gene_modules, pathway_modules)
all_genes <- unique(unlist(module_gene_sets))

cat("\n=== 8-MODULE SUMMARY ===\n")
for (mod in names(module_gene_sets)) {
  cat(sprintf("  %s: %d gene(s)\n", mod, length(module_gene_sets[[mod]])))
}
cat(sprintf("Total unique genes: %d\n\n", length(all_genes)))

# =============================================================================
# SECTION 3: IMMUNE SUB-PANEL DEFINITIONS (for Part E)
# =============================================================================

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

# =============================================================================
# SECTION 4: HELPER FUNCTIONS
# =============================================================================

cat("========== DEFINING HELPER FUNCTIONS ==========\n")

# --- Strip S4 Rle columns (Seurat v5 safety) ---
strip_rle <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (cn in colnames(df)) {
    x <- df[[cn]]
    if (!is.atomic(x)) {
      df[[cn]] <- tryCatch(as.vector(x), error = function(e) as.character(x))
    }
  }
  df
}

# --- Metadata column auto-detection ---
detect_column <- function(meta, possible_names, description) {
  for (name in possible_names) {
    if (name %in% colnames(meta)) {
      cat("Found", description, "column:", name, "\n")
      return(name)
    }
  }
  cat("WARNING: Could not find", description, "column\n")
  return(NULL)
}

# --- Gene availability ---
check_gene_availability <- function(seurat_obj, genes) {
  present <- rownames(seurat_obj)
  available <- genes[genes %in% present]
  missing <- genes[!genes %in% present]
  list(available = available, missing = missing,
       pct_available = round(100 * length(available) / length(genes), 1))
}

# --- Module filtering ---
filter_modules <- function(seurat_obj, module_list, min_genes = 1) {
  present <- rownames(seurat_obj)
  filtered <- lapply(module_list, function(genes) unique(genes[genes %in% present]))
  qc <- tibble(
    Module = names(filtered),
    Total_Genes = sapply(module_list, length),
    Found_Genes = sapply(filtered, length),
    Percent = round(100 * Found_Genes / Total_Genes, 1),
    Genes_Found = sapply(filtered, paste, collapse = ", "),
    Genes_Missing = sapply(names(module_list), function(m) {
      miss <- setdiff(module_list[[m]], filtered[[m]])
      if (length(miss) > 0) paste(miss, collapse = ", ") else "None"
    })
  )
  valid <- names(filtered)[sapply(filtered, length) >= min_genes]
  cat("  Valid modules (>=", min_genes, "genes):", length(valid), "/", length(filtered), "\n")
  list(filtered = filtered[valid], qc = qc)
}

# --- Module scoring: raw colMeans ---
calculate_module_scores <- function(seurat_obj, module_list) {
  expr <- GetAssayData(seurat_obj, layer = "data")
  scores <- sapply(names(module_list), function(mod) {
    genes <- module_list[[mod]]
    if (length(genes) == 0) return(rep(NA, ncol(seurat_obj)))
    genes <- genes[genes %in% rownames(expr)]
    if (length(genes) == 0) return(rep(NA, ncol(seurat_obj)))
    if (length(genes) == 1) return(as.numeric(expr[genes, ]))
    colMeans(expr[genes, , drop = FALSE])
  })
  colnames(scores) <- paste0(names(module_list), "_Score")
  as.data.frame(scores)
}

# --- Background-subtracted module scoring (Seurat AddModuleScore-style) ---
calculate_module_scores_bgsub <- function(seurat_obj, module_list, n_background = 100) {
  expr <- GetAssayData(seurat_obj, layer = "data")
  gene_means <- Matrix::rowMeans(expr)
  n_bins <- 25
  bins <- cut(gene_means, breaks = n_bins, labels = FALSE)
  names(bins) <- rownames(expr)
  
  scores <- sapply(names(module_list), function(mod) {
    genes <- module_list[[mod]]
    genes <- genes[genes %in% rownames(expr)]
    if (length(genes) == 0) return(rep(NA, ncol(seurat_obj)))
    fg <- colMeans(expr[genes, , drop = FALSE])
    ctrl_genes <- c()
    for (g in genes) {
      b <- bins[g]; if (is.na(b)) next
      candidates <- names(bins)[bins == b]
      candidates <- setdiff(candidates, genes)
      n_sample <- min(ceiling(n_background / length(genes)), length(candidates))
      if (n_sample > 0) ctrl_genes <- c(ctrl_genes, sample(candidates, n_sample))
    }
    ctrl_genes <- unique(ctrl_genes)
    if (length(ctrl_genes) == 0) return(fg)
    bg <- colMeans(expr[ctrl_genes, , drop = FALSE])
    as.numeric(fg - bg)
  })
  colnames(scores) <- paste0(names(module_list), "_Score")
  as.data.frame(scores)
}

# --- Aggregate cell-level to donor-level ---
aggregate_scores_to_donor <- function(data_df, score_cols, donor_col,
                                       condition_col = NULL, covariate_cols = NULL) {
  data_df <- as.data.frame(data_df)
  data_df$donor_id <- as.character(data_df[[donor_col]])
  donors <- unique(data_df$donor_id)
  
  agg_list <- lapply(donors, function(d) {
    subset_df <- data_df[data_df$donor_id == d, , drop = FALSE]
    result <- data.frame(donor_id = d, n_cells = nrow(subset_df), stringsAsFactors = FALSE)
    for (col in score_cols) {
      result[[col]] <- if (col %in% colnames(subset_df)) mean(subset_df[[col]], na.rm = TRUE) else NA
    }
    if (!is.null(condition_col) && condition_col %in% colnames(subset_df)) {
      result[[condition_col]] <- as.character(subset_df[[condition_col]][1])
    }
    if (!is.null(covariate_cols)) {
      for (cov in covariate_cols) {
        if (cov %in% colnames(subset_df)) {
          vals <- subset_df[[cov]]
          result[[cov]] <- if (is.numeric(vals)) mean(vals, na.rm = TRUE) else as.character(vals[1])
        }
      }
    }
    result
  })
  do.call(rbind, agg_list)
}

# --- Pseudobulk creation ---
create_pseudobulk <- function(seurat_obj, donor_col, condition_col, genes = NULL,
                              covariate_cols = NULL) {
  counts <- GetAssayData(seurat_obj, layer = "counts")
  if (!is.null(genes)) {
    available_genes <- genes[genes %in% rownames(counts)]
    counts <- counts[available_genes, , drop = FALSE]
  }
  meta <- as.data.frame(seurat_obj@meta.data)
  meta[[donor_col]] <- as.character(meta[[donor_col]])
  meta[[condition_col]] <- as.character(meta[[condition_col]])
  donors <- unique(meta[[donor_col]])
  
  pb_counts <- matrix(0, nrow = nrow(counts), ncol = length(donors),
                      dimnames = list(rownames(counts), donors))
  for (donor in donors) {
    cells <- rownames(meta)[meta[[donor_col]] == donor]
    if (length(cells) > 0)
      pb_counts[, donor] <- Matrix::rowSums(counts[, cells, drop = FALSE])
  }
  
  sample_meta <- data.frame(donor_id = donors, stringsAsFactors = FALSE)
  names(sample_meta)[1] <- donor_col
  sample_meta[[condition_col]] <- sapply(donors, function(d)
    meta[[condition_col]][meta[[donor_col]] == d][1])
  rownames(sample_meta) <- sample_meta[[donor_col]]
  sample_meta <- sample_meta[colnames(pb_counts), , drop = FALSE]
  cell_counts <- table(meta[[donor_col]])
  sample_meta$n_cells <- as.numeric(cell_counts[rownames(sample_meta)])
  
  if (!is.null(covariate_cols)) {
    for (cov in covariate_cols) {
      if (cov %in% colnames(meta)) {
        sample_meta[[cov]] <- sapply(donors, function(d) {
          vals <- meta[[cov]][meta[[donor_col]] == d]
          vals <- vals[!is.na(vals)]
          if (length(vals) == 0) return(NA)
          if (is.numeric(vals)) mean(vals) else vals[1]
        })
      }
    }
  }
  
  list(counts = pb_counts, metadata = sample_meta,
       n_donors = length(donors), n_genes = nrow(pb_counts))
}

# --- DESeq2 pseudobulk ---
run_deseq2_pseudobulk <- function(pb_data, condition_col, reference = "Control",
                                   test_level = "PD", covariate_cols = NULL) {
  pb_data$metadata[[condition_col]] <- factor(pb_data$metadata[[condition_col]],
                                              levels = c(reference, test_level))
  n_ref  <- sum(pb_data$metadata[[condition_col]] == reference)
  n_test <- sum(pb_data$metadata[[condition_col]] == test_level)
  cat(sprintf("  DESeq2: %d %s vs %d %s\n", n_test, test_level, n_ref, reference))
  if (n_ref < 2 || n_test < 2) { warning("Insufficient samples"); return(NULL) }
  
  active_covariates <- c()
  if (!is.null(covariate_cols)) {
    for (cov in covariate_cols) {
      if (!cov %in% colnames(pb_data$metadata)) next
      vals <- pb_data$metadata[[cov]]
      if (is.numeric(vals) && any(is.na(vals))) {
        pb_data$metadata[[cov]][is.na(vals)] <- mean(vals, na.rm = TRUE)
        vals <- pb_data$metadata[[cov]]
      }
      if (all(is.na(vals))) next
      if (is.numeric(vals) && sd(vals, na.rm = TRUE) == 0) next
      if (!is.numeric(vals)) {
        pb_data$metadata[[cov]] <- as.factor(pb_data$metadata[[cov]])
        if (nlevels(pb_data$metadata[[cov]]) < 2) next
      }
      active_covariates <- c(active_covariates, cov)
    }
  }
  
  design_formula <- if (length(active_covariates) > 0) {
    as.formula(paste0("~ ", paste(active_covariates, collapse = " + "), " + ", condition_col))
  } else {
    as.formula(paste0("~", condition_col))
  }
  cat("  Design:", deparse(design_formula), "\n")
  
  dds <- tryCatch(
    DESeqDataSetFromMatrix(countData = round(pb_data$counts),
                           colData = pb_data$metadata, design = design_formula),
    error = function(e) { warning(e$message); NULL })
  if (is.null(dds)) return(NULL)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  if (nrow(dds) == 0) return(NULL)
  cat(sprintf("  Genes passing filter: %d\n", nrow(dds)))
  
  dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) { warning(e$message); NULL })
  if (is.null(dds)) return(NULL)
  
  res <- results(dds, contrast = c(condition_col, test_level, reference))
  res_df <- as.data.frame(res) %>%
    rownames_to_column("Gene") %>%
    arrange(padj) %>%
    mutate(
      significant = padj < 0.05 & !is.na(padj),
      direction = case_when(
        log2FoldChange > 0 & significant ~ "Up_in_PD",
        log2FoldChange < 0 & significant ~ "Down_in_PD",
        TRUE ~ "NS"),
      FoldChange = 2^log2FoldChange)
  
  list(results = res_df, dds = dds,
       normalized_counts = counts(dds, normalized = TRUE),
       n_significant = sum(res_df$significant, na.rm = TRUE))
}

# --- Covariate-adjusted group comparison (lm residuals) ---
compare_groups_covariate_adjusted <- function(data, cols, group_col, group1, group2,
                                               covariate_cols = NULL) {
  active_covs <- covariate_cols[covariate_cols %in% colnames(data)]
  active_covs <- active_covs[sapply(active_covs, function(cv) {
    vals <- data[[cv]]
    !all(is.na(vals)) && (is.numeric(vals) && sd(vals, na.rm = TRUE) > 0 ||
                          !is.numeric(vals) && length(unique(na.omit(vals))) > 1)
  })]
  
  results <- lapply(cols, function(col) {
    if (!(col %in% colnames(data))) return(NULL)
    formula_str <- paste0("`", col, "` ~ ",
                          if (length(active_covs) > 0) paste(active_covs, collapse = " + ") else "1",
                          " + ", group_col)
    lm_data <- data[data[[group_col]] %in% c(group1, group2), ]
    lm_data[[group_col]] <- factor(lm_data[[group_col]], levels = c(group2, group1))
    
    for (cv in active_covs) {
      if (is.numeric(lm_data[[cv]]) && any(is.na(lm_data[[cv]])))
        lm_data[[cv]][is.na(lm_data[[cv]])] <- mean(lm_data[[cv]], na.rm = TRUE)
    }
    
    complete_rows <- complete.cases(lm_data[, c(col, active_covs, group_col)])
    lm_data <- lm_data[complete_rows, ]
    if (nrow(lm_data) < (length(active_covs) + 3)) return(NULL)
    
    fit <- tryCatch(lm(as.formula(formula_str), data = lm_data), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    
    sm  <- summary(fit)
    cof <- as.data.frame(sm$coefficients)
    condition_row <- grep(group_col, rownames(cof), value = TRUE)[1]
    if (is.na(condition_row)) return(NULL)
    
    g1_raw <- data[[col]][data[[group_col]] == group1]
    g2_raw <- data[[col]][data[[group_col]] == group2]
    mean1 <- mean(g1_raw, na.rm = TRUE); mean2 <- mean(g2_raw, na.rm = TRUE)
    log2fc <- if (mean2 > 0 && mean1 > 0) log2(mean1/mean2) else NA
    
    data.frame(
      Variable = col, Comparison = paste(group1, "vs", group2),
      Covariates = paste(active_covs, collapse = "+"),
      N_G1 = sum(lm_data[[group_col]] == group1),
      N_G2 = sum(lm_data[[group_col]] == group2),
      Mean_G1_raw = mean1, Mean_G2_raw = mean2, Log2FC_raw = log2fc,
      lm_Beta_G1 = cof[condition_row, "Estimate"],
      lm_SE = cof[condition_row, "Std. Error"],
      lm_t = cof[condition_row, "t value"],
      lm_p = cof[condition_row, "Pr(>|t|)"],
      lm_R2 = sm$r.squared,
      stringsAsFactors = FALSE)
  })
  
  results_clean <- results[!sapply(results, is.null)]
  if (length(results_clean) == 0) return(NULL)
  result_df <- dplyr::bind_rows(results_clean)
  result_df$lm_p_adj <- p.adjust(result_df$lm_p, method = "BH")
  result_df$Direction <- case_when(
    is.na(result_df$lm_Beta_G1) | is.na(result_df$lm_p_adj) ~ "NA",
    result_df$lm_Beta_G1 > 0 & result_df$lm_p_adj < 0.05 ~ paste0("Up_in_", group1),
    result_df$lm_Beta_G1 < 0 & result_df$lm_p_adj < 0.05 ~ paste0("Down_in_", group1),
    TRUE ~ "NS")
  return(result_df)
}

# --- Partial Spearman with covariate residualization ---
partial_spearman <- function(x, y, covariate_df = NULL) {
  if (is.null(covariate_df) || ncol(covariate_df) == 0) {
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

# --- Fisher Z test ---
fisher_z_test <- function(rho1, n1, rho2, n2) {
  if (any(is.na(c(rho1, rho2))) || n1 < 4 || n2 < 4)
    return(list(z = NA, p = NA))
  z1 <- 0.5 * log((1 + rho1) / (1 - rho1))
  z2 <- 0.5 * log((1 + rho2) / (1 - rho2))
  se <- sqrt(1/(n1 - 3) + 1/(n2 - 3))
  z  <- (z1 - z2) / se
  list(z = z, p = 2 * pnorm(-abs(z)), z1 = z1, z2 = z2)
}

# --- Background-subtracted donor-level panel score (for immune panels) ---
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
    if (length(pool) >= n_ctrl) ctrl_pool <- c(ctrl_pool, sample(pool, n_ctrl))
    else if (length(pool) > 0) ctrl_pool <- c(ctrl_pool, pool)
  }
  ctrl_pool <- unique(ctrl_pool)
  if (length(ctrl_pool) == 0) return(NULL)
  panel_mean <- Matrix::colMeans(expr_all[avail, , drop = FALSE])
  ctrl_mean  <- Matrix::colMeans(expr_all[ctrl_pool, , drop = FALSE])
  cell_score <- as.numeric(panel_mean - ctrl_mean)
  meta <- strip_rle(as.data.frame(obj@meta.data))
  agg_df <- data.frame(donor_id = as.character(meta[[donor_col]]),
                       condition = as.character(meta[[cond_col]]),
                       score_raw = cell_score, stringsAsFactors = FALSE)
  donor_df <- aggregate(score_raw ~ donor_id + condition, data = agg_df,
                        FUN = mean, na.rm = TRUE)
  colnames(donor_df)[colnames(donor_df) == "score_raw"] <- "score"
  list(donor_df = as.data.frame(donor_df, stringsAsFactors = FALSE),
       genes_used = avail)
}

# --- Get donor-level covariates from Seurat ---
get_donor_cov_df <- function(obj, donor_col, donor_ids, cov_cols) {
  if (length(cov_cols) == 0) return(NULL)
  meta <- strip_rle(as.data.frame(obj@meta.data))
  keep <- cov_cols[cov_cols %in% colnames(meta)]
  if (length(keep) == 0) return(NULL)
  meta$donor_key <- as.character(meta[[donor_col]])
  unique_donors <- unique(meta$donor_key)
  donor_cov_list <- lapply(unique_donors, function(d) {
    sub <- meta[meta$donor_key == d, keep, drop = FALSE]
    sapply(keep, function(cn) {
      x <- sub[[cn]]
      if (is.numeric(x)) mean(x, na.rm = TRUE) else as.character(x[1])
    })
  })
  donor_meta <- as.data.frame(do.call(rbind, donor_cov_list), stringsAsFactors = FALSE)
  colnames(donor_meta) <- keep
  donor_meta$donor_key <- unique_donors
  for (cn in keep) {
    num_try <- suppressWarnings(as.numeric(donor_meta[[cn]]))
    if (!all(is.na(num_try))) donor_meta[[cn]] <- num_try
  }
  donor_meta <- donor_meta[match(donor_ids, donor_meta$donor_key), , drop = FALSE]
  cov_df <- donor_meta[, keep, drop = FALSE]
  for (cn in colnames(cov_df)) {
    if (!is.numeric(cov_df[[cn]])) cov_df[[cn]] <- as.numeric(factor(cov_df[[cn]])) - 1
    if (any(is.na(cov_df[[cn]])))
      cov_df[[cn]][is.na(cov_df[[cn]])] <- mean(cov_df[[cn]], na.rm = TRUE)
  }
  keep2 <- sapply(cov_df, function(x) sd(x, na.rm = TRUE) > 0)
  cov_df[, keep2, drop = FALSE]
}

sig_star <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

# =============================================================================
# SECTION 5: LOAD DATA & DETECT METADATA
# =============================================================================

cat("\n========== LOADING DATA ==========\n")

astro  <- readRDS(astro_file);  cat("  Astrocytes:", ncol(astro), "cells\n")
neuron <- readRDS(neuron_file); cat("  Neurons:", ncol(neuron), "cells\n")

DefaultAssay(astro)  <- "RNA"
DefaultAssay(neuron) <- "RNA"
if (length(Layers(astro))  > 1) astro  <- JoinLayers(astro)
if (length(Layers(neuron)) > 1) neuron <- JoinLayers(neuron)

astro_meta  <- astro@meta.data
neuron_meta <- neuron@meta.data

donor_possible <- c("donor_id", "Donor_ID", "donor", "Donor", "sample_id",
                    "Sample_ID", "subject_id", "Individual")
condition_possible <- c("condition", "Condition", "disease", "Disease",
                        "Status", "status", "Group", "group", "Diagnosis")
celltype_possible <- c("cell_type", "celltype", "CellType", "cell_type_broad",
                       "annotation", "subtype", "Subtype")

donor_col_astro  <- detect_column(astro_meta, donor_possible, "donor (astro)")
cond_col_astro   <- detect_column(astro_meta, condition_possible, "condition (astro)")
donor_col_neuron <- detect_column(neuron_meta, donor_possible, "donor (neuron)")
cond_col_neuron  <- detect_column(neuron_meta, condition_possible, "condition (neuron)")
ct_col_neuron    <- detect_column(neuron_meta, celltype_possible, "cell type (neuron)")

# --- Covariate detection ---
covariate_cols_to_use <- c()
if (COVARIATES_ENABLED) {
  for (cov_info in list(
    list(use = USE_COV_AGE, col = COV_AGE_COL,
         candidates = c("age", "Age", "AGE", "age_death", "AoD")),
    list(use = USE_COV_PMI, col = COV_PMI_COL,
         candidates = c("pmi", "PMI", "post_mortem_interval")),
    list(use = USE_COV_SEX, col = COV_SEX_COL,
         candidates = c("sex", "Sex", "gender", "Gender"))
  )) {
    if (!cov_info$use) next
    candidates <- unique(c(cov_info$col, cov_info$candidates))
    found <- detect_column(astro_meta, candidates, cov_info$col)
    if (!is.null(found)) covariate_cols_to_use <- c(covariate_cols_to_use, found)
  }
  cat(sprintf("Active covariates: %s\n\n",
              if (length(covariate_cols_to_use) > 0) paste(covariate_cols_to_use, collapse = ", ") else "none"))
}

# --- Neuron subtype splitting ---
neurons_exc <- NULL; neurons_inh <- NULL
if (!is.null(ct_col_neuron)) {
  cell_types_found <- unique(neuron_meta[[ct_col_neuron]])
  exc_pattern <- cell_types_found[grepl("excit|glut|pyramidal|projection", cell_types_found, ignore.case = TRUE)]
  if (length(exc_pattern) > 0) {
    neurons_exc <- tryCatch(
      subset(neuron, cells = rownames(neuron_meta)[neuron_meta[[ct_col_neuron]] %in% exc_pattern]),
      error = function(e) NULL)
    if (!is.null(neurons_exc) && ncol(neurons_exc) < 50) neurons_exc <- NULL
  }
  inh_pattern <- cell_types_found[grepl("inhib|gaba|inter", cell_types_found, ignore.case = TRUE)]
  if (length(inh_pattern) > 0) {
    neurons_inh <- tryCatch(
      subset(neuron, cells = rownames(neuron_meta)[neuron_meta[[ct_col_neuron]] %in% inh_pattern]),
      error = function(e) NULL)
    if (!is.null(neurons_inh) && ncol(neurons_inh) < 50) neurons_inh <- NULL
  }
}

cell_type_list <- list(
  Astrocytes    = list(obj = astro,  donor_col = donor_col_astro,  cond_col = cond_col_astro),
  Neurons_Total = list(obj = neuron, donor_col = donor_col_neuron, cond_col = cond_col_neuron)
)
if (!is.null(neurons_exc))
  cell_type_list$Neurons_Excitatory <- list(obj = neurons_exc, donor_col = donor_col_neuron, cond_col = cond_col_neuron)
if (!is.null(neurons_inh))
  cell_type_list$Neurons_Inhibitory <- list(obj = neurons_inh, donor_col = donor_col_neuron, cond_col = cond_col_neuron)

cat("Cell types:", paste(names(cell_type_list), collapse = ", "), "\n")

# --- Gene availability for TE_Control specifically ---
cat("\n=== TE_CONTROL GENE AVAILABILITY ===\n")
te_genes <- module_gene_sets[["TE_Control"]]
for (ct_name in names(cell_type_list)) {
  obj <- cell_type_list[[ct_name]]$obj
  avail <- te_genes[te_genes %in% rownames(obj)]
  miss  <- te_genes[!te_genes %in% rownames(obj)]
  cat(sprintf("  %s: %d/%d TE_Control genes available\n", ct_name, length(avail), length(te_genes)))
  if (length(miss) > 0) cat("    Missing:", paste(miss, collapse = ", "), "\n")
}

# Common genes across cell types
common_genes <- Reduce(intersect, lapply(cell_type_list, function(x) rownames(x$obj)))
common_genes <- intersect(common_genes, all_genes)

# --- Create output directories ---
for (subdir in c("A_DESeq2", "B_ModScores", "C_CellType", "D_CrossModule",
                 "E_Immune", "F_SNCA_Trajectory", "tables", "figures")) {
  dir.create(file.path(results_main, subdir), showWarnings = FALSE, recursive = TRUE)
}

# =============================================================================
# SECTION 6: COMPUTE ALL MODULE SCORES (8 modules, donor-level)
# =============================================================================

cat("\n========== COMPUTING 8-MODULE SCORES ==========\n")

donor_scores_all <- list()

for (ct_name in names(cell_type_list)) {
  cat("\n--- Module Scores:", ct_name, "---\n")
  ct_info <- cell_type_list[[ct_name]]
  ct_mod <- filter_modules(ct_info$obj, module_gene_sets, min_genes = 1)
  scores <- calculate_module_scores(ct_info$obj, ct_mod$filtered)
  cell_df <- cbind(as.data.frame(ct_info$obj@meta.data), scores)
  cell_df$condition <- as.character(cell_df[[ct_info$cond_col]])
  cell_df$donor_id  <- as.character(cell_df[[ct_info$donor_col]])
  
  score_cols <- colnames(scores)
  donor_df <- aggregate_scores_to_donor(cell_df, score_cols, "donor_id", "condition",
                                         covariate_cols_to_use)
  
  # Make numeric covariates numeric
  for (cv in covariate_cols_to_use) {
    if (cv %in% colnames(donor_df) && !is.numeric(donor_df[[cv]]))
      donor_df[[cv]] <- {  # v3 fix: this line turned the donor-level sex ("Male"/"Female") into NA, so sex silently dropped out
        num <- suppressWarnings(as.numeric(donor_df[[cv]]))   # of Parts B, D and F. Numbers stored as text are still converted;
        if (all(is.na(num) == is.na(donor_df[[cv]]))) num else {   # a two-level text covariate is coded 0/1, as get_donor_cov_df() does
          stopifnot(nlevels(factor(donor_df[[cv]])) <= 2)          # for Part E (sex: Female = 0, Male = 1)
          as.numeric(factor(donor_df[[cv]])) - 1 } }
  }
  
  donor_scores_all[[ct_name]] <- donor_df
  cat("  Donors:", nrow(donor_df), "| Score cols:", sum(grepl("_Score", names(donor_df))), "\n")
}

score_cols_all <- names(donor_scores_all[["Astrocytes"]])[grepl("_Score", names(donor_scores_all[["Astrocytes"]]))]
cat("\nModule score columns:", paste(score_cols_all, collapse = ", "), "\n")

# #####################################################################
# PART A: DESeq2 PSEUDOBULK — PD vs CONTROL (ALL 8 MODULES)
# #####################################################################

cat("\n============================================================\n")
cat("=== PART A: DESeq2 PSEUDOBULK (8 MODULES) ===\n")
cat("============================================================\n")

deseq_results <- list()

for (ct_name in names(cell_type_list)) {
  cat("\n--- DESeq2:", ct_name, "---\n")
  ct_info <- cell_type_list[[ct_name]]
  
  pb <- create_pseudobulk(ct_info$obj, ct_info$donor_col, ct_info$cond_col,
                          genes = common_genes,
                          covariate_cols = if (COVARIATES_ENABLED) covariate_cols_to_use else NULL)
  
  deseq_res <- run_deseq2_pseudobulk(pb, ct_info$cond_col,
                                     covariate_cols = if (COVARIATES_ENABLED) covariate_cols_to_use else NULL)
  deseq_results[[ct_name]] <- deseq_res
  
  if (!is.null(deseq_res)) {
    # Extract TE_Control gene results
    te_res <- deseq_res$results %>% filter(Gene %in% te_genes)
    cat("  TE_Control DESeq2 results:\n")
    if (nrow(te_res) > 0) {
      print(te_res[, c("Gene", "log2FoldChange", "padj", "direction")], row.names = FALSE)
    } else {
      cat("    No TE_Control genes passed DESeq2 filter\n")
    }
    write.csv(deseq_res$results,
              file.path(results_main, "A_DESeq2",
                        paste0(ct_name, "_DESeq2_8Module.csv")),
              row.names = FALSE)
  }
}

# Module-level summary
deseq_module_summary <- do.call(rbind, lapply(names(deseq_results), function(ct) {
  if (is.null(deseq_results[[ct]])) return(NULL)
  do.call(rbind, lapply(names(module_gene_sets), function(mod) {
    genes <- module_gene_sets[[mod]]
    mod_res <- deseq_results[[ct]]$results %>% filter(Gene %in% genes)
    if (nrow(mod_res) == 0) return(NULL)
    data.frame(Region = REGION, Cell_Type = ct, Module = mod,
               N_Tested = nrow(mod_res),
               N_Sig = sum(mod_res$significant, na.rm = TRUE),
               N_Up = sum(mod_res$direction == "Up_in_PD", na.rm = TRUE),
               N_Down = sum(mod_res$direction == "Down_in_PD", na.rm = TRUE),
               Mean_L2FC = mean(mod_res$log2FoldChange, na.rm = TRUE),
               Min_padj = min(mod_res$padj, na.rm = TRUE),
               Sig_Genes = paste(mod_res$Gene[mod_res$significant], collapse = ", "),
               stringsAsFactors = FALSE)
  }))
}))

write.csv(deseq_module_summary,
          file.path(results_main, "A_DESeq2", "DESeq2_8Module_Summary.csv"),
          row.names = FALSE)
cat("\nPart A complete.\n")

# #####################################################################
# PART B: MODULE SCORE COMPARISONS — PD vs CONTROL (covariate-adjusted)
# #####################################################################

cat("\n============================================================\n")
cat("=== PART B: PD vs CONTROL MODULE SCORES (8 MODULES) ===\n")
cat("============================================================\n")

all_comparisons <- list()

for (ct_name in names(cell_type_list)) {
  cat("\n--- PD vs Control:", ct_name, "---\n")
  donor_df <- donor_scores_all[[ct_name]]
  
  comp <- compare_groups_covariate_adjusted(
    donor_df, score_cols_all, "condition", "PD", "Control",
    covariate_cols = covariate_cols_to_use)
  
  if (!is.null(comp)) {
    comp$Cell_Type <- ct_name
    comp$Region <- REGION
    all_comparisons[[ct_name]] <- comp
    
    # Print TE_Control row
    te_row <- comp[comp$Variable == "TE_Control_Score", ]
    if (nrow(te_row) > 0) {
      cat(sprintf("  TE_Control: Beta=%.4f, p=%.4f, p_adj=%.4f, Direction=%s\n",
                  te_row$lm_Beta_G1, te_row$lm_p, te_row$lm_p_adj, te_row$Direction))
    }
  }
  
  # Boxplot for TE_Control
  if ("TE_Control_Score" %in% colnames(donor_df)) {
    p <- ggplot(donor_df, aes(x = condition, y = TE_Control_Score, fill = condition)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
      geom_jitter(width = 0.15, size = 3, alpha = 0.8) +
      scale_fill_manual(values = c("Control" = "#4DAF4A", "PD" = "#E41A1C")) +
      labs(title = paste0(REGION, " ", ct_name, ": TE_Control (Donor Level)"),
           x = "", y = "TE_Control Score") +
      theme_minimal() + theme(legend.position = "none")
    ggsave(file.path(results_main, "B_ModScores",
                     paste0(ct_name, "_TE_Control_Boxplot.png")),
           p, width = 5, height = 6, dpi = 300)
  }
}

master_comparisons <- dplyr::bind_rows(all_comparisons)
write.csv(master_comparisons,
          file.path(results_main, "B_ModScores",
                    paste0(REGION, "_8Module_PDvsControl_CovAdj.csv")),
          row.names = FALSE)
cat("\nPart B complete.\n")

# #####################################################################
# PART C: CELL-TYPE COMPARISONS — ASTRO vs NEURON (paired Wilcoxon)
# #####################################################################

cat("\n============================================================\n")
cat("=== PART C: ASTROCYTE vs NEURON (8 MODULES) ===\n")
cat("============================================================\n")

neuron_types <- names(cell_type_list)[grepl("Neuron", names(cell_type_list))]
celltype_results <- list()

for (n_type in neuron_types) {
  cat("\n--- Astrocytes vs", n_type, "---\n")
  astro_donor  <- donor_scores_all[["Astrocytes"]]
  neuron_donor <- donor_scores_all[[n_type]]
  score_cols <- intersect(names(astro_donor)[grepl("_Score", names(astro_donor))],
                          names(neuron_donor)[grepl("_Score", names(neuron_donor))])
  
  for (group_name in c("ALL", "PD", "Control")) {
    a_df <- if (group_name == "ALL") astro_donor else astro_donor[astro_donor$condition == group_name, ]
    n_df <- if (group_name == "ALL") neuron_donor else neuron_donor[neuron_donor$condition == group_name, ]
    common_donors <- intersect(a_df$donor_id, n_df$donor_id)
    if (length(common_donors) < 3) next
    
    a_matched <- a_df[match(sort(common_donors), a_df$donor_id), ]
    n_matched <- n_df[match(sort(common_donors), n_df$donor_id), ]
    
    comp_list <- lapply(score_cols, function(mod) {
      a_vals <- a_matched[[mod]]; n_vals <- n_matched[[mod]]
      test <- tryCatch(wilcox.test(a_vals, n_vals, paired = TRUE),
                       error = function(e) list(p.value = NA))
      mean_a <- mean(a_vals, na.rm = TRUE); mean_n <- mean(n_vals, na.rm = TRUE)
      log2fc <- if (mean_n > 0 && mean_a > 0) log2(mean_a / mean_n) else NA
      data.frame(Module = mod, Comparison = paste0("Astro_vs_", n_type),
                 Group = group_name, N_Pairs = length(common_donors),
                 Mean_Astro = mean_a, Mean_Neuron = mean_n,
                 Log2FC = log2fc, Paired_p = test$p.value,
                 stringsAsFactors = FALSE)
    })
    comp_df <- do.call(rbind, comp_list)
    comp_df$p_adj <- p.adjust(comp_df$Paired_p, method = "BH")
    comp_df$Region <- REGION
    celltype_results[[paste(n_type, group_name)]] <- comp_df
  }
}

master_celltype <- do.call(rbind, celltype_results)
write.csv(master_celltype,
          file.path(results_main, "C_CellType",
                    paste0(REGION, "_8Module_CellType_Comparisons.csv")),
          row.names = FALSE)

# Print TE_Control results
cat("\n=== TE_Control: Astrocyte vs Neuron ===\n")
te_ct <- master_celltype[grepl("TE_Control", master_celltype$Module), ]
print(te_ct[, c("Module", "Comparison", "Group", "Mean_Astro", "Mean_Neuron", "Log2FC", "p_adj")],
      row.names = FALSE)
cat("\nPart C complete.\n")

# #####################################################################
# PART D: CROSS-MODULE CORRELATIONS (TE_Control vs ALL other modules)
# #####################################################################

cat("\n============================================================\n")
cat("=== PART D: CROSS-MODULE CORRELATIONS (8 MODULES) ===\n")
cat("============================================================\n")

all_cor_results <- list()

for (ct_name in names(cell_type_list)) {
  cat("\n--- Correlations:", ct_name, "---\n")
  donor_df <- donor_scores_all[[ct_name]]
  score_cols <- names(donor_df)[grepl("_Score", names(donor_df))]
  
  for (group_name in c("ALL", "PD", "Control")) {
    sub_df <- if (group_name == "ALL") donor_df else donor_df[donor_df$condition == group_name, ]
    cat(sprintf("  %s (%s): n=%d\n", ct_name, group_name, nrow(sub_df)))
    
    # Residualize scores against covariates
    cov_df <- NULL
    if (COVARIATES_ENABLED && length(covariate_cols_to_use) > 0) {
      active_covs <- covariate_cols_to_use[covariate_cols_to_use %in% colnames(sub_df)]
      active_covs <- active_covs[sapply(active_covs, function(cv) {
        vals <- sub_df[[cv]]
        !all(is.na(vals)) && is.numeric(vals) && sd(vals, na.rm = TRUE) > 0
      })]
      if (length(active_covs) > 0) {
        # Impute
        for (cv in active_covs) {
          if (any(is.na(sub_df[[cv]])))
            sub_df[[cv]][is.na(sub_df[[cv]])] <- mean(sub_df[[cv]], na.rm = TRUE)
        }
        # Residualize each score
        resid_data <- sub_df
        for (col in score_cols) {
          if (!col %in% colnames(sub_df)) next
          form <- as.formula(paste0("`", col, "` ~ ", paste(active_covs, collapse = " + ")))
          complete <- complete.cases(sub_df[, c(col, active_covs)])
          if (sum(complete) < length(active_covs) + 3) next
          fit <- tryCatch(lm(form, data = sub_df[complete, ]), error = function(e) NULL)
          if (!is.null(fit)) {
            resid_data[[col]] <- NA
            resid_data[[col]][complete] <- residuals(fit)
          }
        }
        resid_clean <- resid_data[, score_cols, drop = FALSE]
        resid_clean <- resid_clean[complete.cases(resid_clean), ]
      } else {
        resid_clean <- sub_df[, score_cols, drop = FALSE]
        resid_clean <- resid_clean[complete.cases(resid_clean), ]
      }
    } else {
      resid_clean <- sub_df[, score_cols, drop = FALSE]
      resid_clean <- resid_clean[complete.cases(resid_clean), ]
    }
    
    if (nrow(resid_clean) < 5) next
    
    # Pairwise correlations
    pairs <- list()
    for (i in 1:(length(score_cols) - 1)) {
      for (j in (i + 1):length(score_cols)) {
        if (!score_cols[i] %in% colnames(resid_clean) || !score_cols[j] %in% colnames(resid_clean)) next
        test <- tryCatch(
          cor.test(resid_clean[[score_cols[i]]], resid_clean[[score_cols[j]]], method = "spearman"),
          error = function(e) list(estimate = NA, p.value = NA))
        pairs[[length(pairs) + 1]] <- data.frame(
          Var1 = score_cols[i], Var2 = score_cols[j],
          rho = as.numeric(test$estimate), p_value = test$p.value,
          n = nrow(resid_clean), Cell_Type = ct_name, Group = group_name,
          Region = REGION, stringsAsFactors = FALSE)
      }
    }
    pairs_df <- do.call(rbind, pairs)
    if (!is.null(pairs_df) && nrow(pairs_df) > 0) {
      pairs_df$p_adj <- p.adjust(pairs_df$p_value, method = "BH")
      all_cor_results[[paste(ct_name, group_name)]] <- pairs_df
    }
    
    # Heatmap
    cor_matrix <- cor(resid_clean, method = "spearman", use = "pairwise.complete.obs")
    png(file.path(results_main, "D_CrossModule",
                  paste0(ct_name, "_8Module_CorHeatmap_", group_name, ".png")),
        width = 10, height = 10, units = "in", res = 300)
    pheatmap(cor_matrix,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             breaks = seq(-1, 1, length.out = 101),
             main = paste0(REGION, " ", ct_name, " (", group_name, ") — 8 Modules"),
             fontsize = 8, display_numbers = TRUE,
             number_format = "%.2f", number_color = "black")
    dev.off()
  }
}

master_correlations <- dplyr::bind_rows(all_cor_results)
write.csv(master_correlations,
          file.path(results_main, "D_CrossModule",
                    paste0(REGION, "_8Module_AllCorrelations.csv")),
          row.names = FALSE)

# Print TE_Control correlation highlights
cat("\n=== TE_Control Cross-Module Correlations (FDR < 0.05) ===\n")
te_cors <- master_correlations[grepl("TE_Control", master_correlations$Var1) |
                               grepl("TE_Control", master_correlations$Var2), ]
te_cors_sig <- te_cors[!is.na(te_cors$p_adj) & te_cors$p_adj < 0.05, ]
if (nrow(te_cors_sig) > 0) {
  print(te_cors_sig[order(te_cors_sig$p_adj),
                    c("Var1", "Var2", "rho", "p_adj", "Cell_Type", "Group")],
        row.names = FALSE)
} else {
  cat("  None at FDR < 0.05\n")
}
cat("\n=== ALL TE_Control correlations ===\n")
print(te_cors[order(te_cors$Cell_Type, te_cors$Group, te_cors$p_adj),
              c("Var1", "Var2", "rho", "p_value", "p_adj", "Cell_Type", "Group")],
      row.names = FALSE)
cat("\nPart D complete.\n")

# #####################################################################
# PART E: TE_CONTROL × IMMUNE SUB-PANEL CORRELATIONS + FISHER Z
# #####################################################################

cat("\n============================================================\n")
cat("=== PART E: TE_CONTROL × IMMUNE PANELS (Script 4 logic) ===\n")
cat("============================================================\n")

# Compute immune panel donor scores (bg-subtracted)
panel_scores <- list()
for (ct in list(
  list(name = "Astrocytes", obj = astro,  d = donor_col_astro,  c = cond_col_astro),
  list(name = "Neurons",    obj = neuron, d = donor_col_neuron, c = cond_col_neuron)
)) {
  cat("--- ", ct$name, "---\n")
  panel_scores[[ct$name]] <- list()
  for (pn in panel_names) {
    res <- donor_bg_subtracted_score(ct$obj, immune_subpanels[[pn]], ct$d, ct$c)
    if (!is.null(res)) {
      panel_scores[[ct$name]][[pn]] <- res$donor_df
      cat(sprintf("  %-14s: %d/%d genes, %d donors\n",
                  pn, length(res$genes_used), length(immune_subpanels[[pn]]),
                  nrow(res$donor_df)))
    }
  }
}

# Correlate ALL 8 module scores × 9 immune panels
panel_cor_long <- list()

for (ct_name in c("Astrocytes", "Neurons")) {
  prev_donor <- if (ct_name == "Astrocytes") donor_scores_all[["Astrocytes"]] else donor_scores_all[["Neurons_Total"]]
  obj <- if (ct_name == "Astrocytes") astro else neuron
  donor_col <- if (ct_name == "Astrocytes") donor_col_astro else donor_col_neuron
  ct_panel_scores <- panel_scores[[ct_name]]
  if (length(ct_panel_scores) == 0) next
  
  prev_donor$donor_id  <- trimws(as.character(prev_donor$donor_id))
  prev_donor$condition <- trimws(as.character(prev_donor$condition))
  
  for (pn in names(ct_panel_scores)) {
    panel_df <- ct_panel_scores[[pn]]
    if (is.null(panel_df)) next
    panel_df$donor_id <- trimws(as.character(panel_df$donor_id))
    
    # Merge on donor_id only (safest approach)
    pdf_slim <- panel_df[, setdiff(colnames(panel_df), "condition"), drop = FALSE]
    merged <- merge(prev_donor, pdf_slim, by = "donor_id")
    if (nrow(merged) < 6) next
    
    cov_df <- get_donor_cov_df(obj, donor_col, merged$donor_id, covariate_cols_to_use)
    
    mods_avail <- intersect(score_cols_all, colnames(merged))
    for (mod in mods_avail) {
      res <- partial_spearman(merged[[mod]], merged$score, cov_df)
      panel_cor_long[[length(panel_cor_long) + 1]] <- data.frame(
        cell_type = ct_name, module = mod, panel = pn,
        n = res$n, rho = res$rho, p = res$p,
        stringsAsFactors = FALSE)
    }
  }
}

if (length(panel_cor_long) > 0) {
  panel_cor_df <- do.call(rbind, panel_cor_long)
  panel_cor_df$p_BH <- NA_real_
  for (ct in unique(panel_cor_df$cell_type)) {
    idx <- panel_cor_df$cell_type == ct
    panel_cor_df$p_BH[idx] <- p.adjust(panel_cor_df$p[idx], method = "BH")
  }
  panel_cor_df$star <- sig_star(panel_cor_df$p_BH)
  
  write.csv(panel_cor_df,
            file.path(results_main, "E_Immune",
                      paste0(REGION, "_8Module_x_ImmunePanel.csv")),
            row.names = FALSE)
  
  # --- Fisher Z: Astrocyte vs Neuron ---
  fz_list <- list()
  module_cols_fz <- unique(panel_cor_df$module)
  for (mod in module_cols_fz) {
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
        module = mod, panel = pn,
        rho_astro = a$rho, n_astro = a$n,
        rho_neuron = n$rho, n_neuron = n$n,
        delta_rho = a$rho - n$rho,
        z_diff = fz$z, p_fisher = fz$p,
        stringsAsFactors = FALSE)
    }
  }
  fisher_df <- if (length(fz_list) > 0) {
    fd <- do.call(rbind, fz_list)
    fd$p_BH_fisher <- p.adjust(fd$p_fisher, method = "BH")
    fd$star_fisher <- sig_star(fd$p_BH_fisher)
    fd[order(fd$p_BH_fisher), ]
  } else data.frame()
  
  write.csv(fisher_df,
            file.path(results_main, "E_Immune",
                      paste0(REGION, "_8Module_FisherZ_CellType.csv")),
            row.names = FALSE)
  
  # Print TE_Control immune results
  cat("\n=== TE_Control × Immune Panels ===\n")
  te_immune <- panel_cor_df[grepl("TE_Control", panel_cor_df$module), ]
  print(te_immune[order(te_immune$cell_type, te_immune$p_BH),
                  c("cell_type", "panel", "rho", "p_BH", "star")],
        row.names = FALSE)
  
  cat("\n=== TE_Control Fisher Z (Astro vs Neuron) ===\n")
  te_fz <- fisher_df[grepl("TE_Control", fisher_df$module), ]
  if (nrow(te_fz) > 0) {
    print(te_fz[, c("panel", "rho_astro", "rho_neuron", "delta_rho", "p_BH_fisher", "star_fisher")],
          row.names = FALSE)
  }
  
  # --- Heatmap: 8 modules × 9 panels (astro | neuron) ---
  rho_limit <- max(abs(panel_cor_df$rho), na.rm = TRUE)
  rho_limit <- min(rho_limit, 1)
  
  fig_theme <- theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = element_text(size = 9),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  for (ct_label in c("Astrocytes", "Neurons")) {
    d <- panel_cor_df[panel_cor_df$cell_type == ct_label, ]
    d$module <- factor(d$module, levels = score_cols_all)
    d$panel  <- factor(d$panel, levels = panel_names)
    
    p <- ggplot(d, aes(panel, module, fill = rho)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = star), size = 4, vjust = 0.75) +
      geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5, vjust = -0.35, color = "grey20") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, limits = c(-rho_limit, rho_limit)) +
      labs(title = paste0(REGION, " ", ct_label, " — 8 Modules × Immune Panels"),
           x = "Immune sub-panel", y = "Module") + fig_theme
    
    ggsave(file.path(results_main, "E_Immune",
                     paste0(REGION, "_", ct_label, "_8Module_x_Immune_Heatmap.png")),
           p, width = 10, height = 8, dpi = 300)
  }
} else {
  cat("WARNING: No immune panel correlations computed.\n")
  panel_cor_df <- data.frame()
  fisher_df <- data.frame()
}

cat("\nPart E complete.\n")

# #####################################################################
# PART F: SNCA TRAJECTORY COUPLING — TE_Control along SNCA axis
# #####################################################################

cat("\n============================================================\n")
cat("=== PART F: SNCA TRAJECTORY — TE_CONTROL COUPLING ===\n")
cat("============================================================\n")

snca_cor_results <- list()
snca_col <- "SNCA_Score"
target_mods <- setdiff(score_cols_all, snca_col)

for (ct_name in c("Astrocytes", "Neurons_Total")) {
  ct_label <- gsub("_Total", "", ct_name)
  donor_df <- donor_scores_all[[ct_name]]
  if (!snca_col %in% colnames(donor_df)) next
  
  for (group_name in c("ALL", "PD", "Control")) {
    sub_df <- if (group_name == "ALL") donor_df else donor_df[donor_df$condition == group_name, ]
    if (nrow(sub_df) < 5) next
    
    # Get covariates for residualization
    cov_df <- NULL
    if (COVARIATES_ENABLED && length(covariate_cols_to_use) > 0) {
      active_covs <- covariate_cols_to_use[covariate_cols_to_use %in% colnames(sub_df)]
      active_covs <- active_covs[sapply(active_covs, function(cv) {
        vals <- sub_df[[cv]]
        !all(is.na(vals)) && is.numeric(vals) && sd(vals, na.rm = TRUE) > 0
      })]
      if (length(active_covs) > 0) {
        for (cv in active_covs) {
          if (any(is.na(sub_df[[cv]])))
            sub_df[[cv]][is.na(sub_df[[cv]])] <- mean(sub_df[[cv]], na.rm = TRUE)
        }
        cov_df <- sub_df[, active_covs, drop = FALSE]
      }
    }
    
    for (mod in target_mods) {
      if (!mod %in% colnames(sub_df)) next
      res <- partial_spearman(sub_df[[snca_col]], sub_df[[mod]], cov_df)
      snca_cor_results[[length(snca_cor_results) + 1]] <- data.frame(
        Cell_Type = ct_label, Group = group_name,
        Module = mod, n = res$n, rho = res$rho, p = res$p,
        stringsAsFactors = FALSE)
    }
  }
}

if (length(snca_cor_results) > 0) {
  snca_cor_df <- do.call(rbind, snca_cor_results)
  snca_cor_df$p_adj <- p.adjust(snca_cor_df$p, method = "BH")
  snca_cor_df$star <- sig_star(snca_cor_df$p_adj)
  snca_cor_df$Region <- REGION
  
  write.csv(snca_cor_df,
            file.path(results_main, "F_SNCA_Trajectory",
                      paste0(REGION, "_SNCA_vs_8Modules.csv")),
            row.names = FALSE)
  
  cat("\n=== SNCA × TE_Control Trajectory Correlations ===\n")
  te_snca <- snca_cor_df[grepl("TE_Control", snca_cor_df$Module), ]
  print(te_snca[, c("Cell_Type", "Group", "rho", "p", "p_adj", "star")],
        row.names = FALSE)
  
  # Scatter: SNCA vs TE_Control
  for (ct_name in c("Astrocytes", "Neurons_Total")) {
    ct_label <- gsub("_Total", "", ct_name)
    donor_df <- donor_scores_all[[ct_name]]
    if (!all(c(snca_col, "TE_Control_Score") %in% colnames(donor_df))) next
    
    test <- tryCatch(cor.test(donor_df[[snca_col]], donor_df[["TE_Control_Score"]],
                              method = "spearman"),
                     error = function(e) list(estimate = NA, p.value = NA))
    
    p <- ggplot(donor_df, aes(x = .data[[snca_col]], y = TE_Control_Score, color = condition)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
      scale_color_manual(values = c("Control" = "#4DAF4A", "PD" = "#E41A1C")) +
      labs(title = paste0(REGION, " ", ct_label, ": SNCA vs TE_Control"),
           subtitle = sprintf("rho=%.3f, p=%.4f", test$estimate, test$p.value),
           x = "SNCA Score", y = "TE_Control Score") +
      theme_minimal() + theme(legend.position = "bottom")
    
    ggsave(file.path(results_main, "F_SNCA_Trajectory",
                     paste0(REGION, "_", ct_label, "_SNCA_vs_TE_Control.png")),
           p, width = 7, height = 6, dpi = 300)
  }
} else {
  snca_cor_df <- data.frame()
}
cat("\nPart F complete.\n")

# #####################################################################
# PART G: SUMMARY WORKBOOK + RDS
# #####################################################################

cat("\n============================================================\n")
cat("=== PART G: SUMMARY WORKBOOK ===\n")
cat("============================================================\n")

wb <- createWorkbook()

# Module definitions
addWorksheet(wb, "Module_Definitions")
writeData(wb, 1, data.frame(
  Module = names(module_gene_sets),
  N_Genes = sapply(module_gene_sets, length),
  Genes = sapply(module_gene_sets, paste, collapse = ", ")))

# DESeq2 summary
if (!is.null(deseq_module_summary) && nrow(deseq_module_summary) > 0) {
  addWorksheet(wb, "DESeq2_Summary")
  writeData(wb, "DESeq2_Summary", deseq_module_summary)
}

# PD vs Control
if (!is.null(master_comparisons) && nrow(master_comparisons) > 0) {
  addWorksheet(wb, "PDvsControl_CovAdj")
  writeData(wb, "PDvsControl_CovAdj", master_comparisons)
}

# Cell-type comparisons
if (!is.null(master_celltype) && nrow(master_celltype) > 0) {
  addWorksheet(wb, "Astro_vs_Neuron")
  writeData(wb, "Astro_vs_Neuron", master_celltype)
}

# Cross-module correlations
if (!is.null(master_correlations) && nrow(master_correlations) > 0) {
  addWorksheet(wb, "CrossModule_Correlations")
  writeData(wb, "CrossModule_Correlations", master_correlations)
}

# Immune panel correlations
if (nrow(panel_cor_df) > 0) {
  addWorksheet(wb, "Module_x_ImmunePanel")
  writeData(wb, "Module_x_ImmunePanel", panel_cor_df)
}

# Fisher Z
if (nrow(fisher_df) > 0) {
  addWorksheet(wb, "FisherZ_CellType")
  writeData(wb, "FisherZ_CellType", fisher_df)
}

# SNCA trajectory
if (nrow(snca_cor_df) > 0) {
  addWorksheet(wb, "SNCA_Trajectory")
  writeData(wb, "SNCA_Trajectory", snca_cor_df)
}

# Methods
addWorksheet(wb, "Methods")
writeData(wb, "Methods", data.frame(
  Section = c("Script", "Region", "Modules", "Covariates",
              "DESeq2", "Comparisons", "Correlations",
              "Immune", "SNCA_Trajectory"),
  Description = c(
    "Script 5: TE_Control consolidated analysis (8-module framework)",
    REGION,
    "8 modules: ADAR1, ADAR2, SNCA, dsRNA_Processing, dsRNA_EV_Transport, dsRNA_Sensing, Innate_Immune, TE_Control",
    paste("Covariates:", if (COVARIATES_ENABLED) paste(covariate_cols_to_use, collapse = ", ") else "none"),
    "DESeq2 pseudobulk with covariate-adjusted design",
    "Donor-level covariate-adjusted lm (PD vs Control) + paired Wilcoxon (Astro vs Neuron)",
    "Partial Spearman correlations on covariate-residualized scores, BH-corrected",
    "8 module × 9 immune sub-panel partial Spearman + Fisher Z cell-type comparison",
    "SNCA pseudobulk as continuous axis, TE_Control coupling with covariates")))

xlsx_file <- paste0(REGION, "_TE_Control_Comprehensive.xlsx")
old_wd <- getwd()
setwd(results_main)
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
setwd(old_wd)

# Save RDS
saveRDS(list(
  region = REGION,
  modules = module_gene_sets,
  cell_types = names(cell_type_list),
  deseq_results = deseq_results,
  deseq_module_summary = deseq_module_summary,
  donor_scores = donor_scores_all,
  master_comparisons = master_comparisons,
  master_celltype = master_celltype,
  master_correlations = master_correlations,
  panel_cor_df = panel_cor_df,
  fisher_df = fisher_df,
  snca_trajectory = snca_cor_df
), file.path(results_main, paste0(REGION, "_8Module_AllResults.rds")))

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n============================================================\n")
cat("=== SCRIPT 5 COMPLETE ===\n")
cat("============================================================\n")
cat("Region:", REGION, "\n")
cat("Output:", results_main, "\n")
cat("\n=== FOLDER STRUCTURE ===\n")
cat("  A_DESeq2/           — DESeq2 PD vs Control (8 modules)\n")
cat("  B_ModScores/        — Module score comparisons + boxplots\n")
cat("  C_CellType/         — Astrocyte vs Neuron paired comparisons\n")
cat("  D_CrossModule/      — Cross-module correlation heatmaps\n")
cat("  E_Immune/           — Module × immune panel heatmaps + Fisher Z\n")
cat("  F_SNCA_Trajectory/  — TE_Control coupling along SNCA axis\n")
cat("\n=== KEY TE_CONTROL RESULTS ===\n")

# Quick TE_Control dashboard
cat("\n  [DESeq2] Per-gene results:\n")
for (ct in names(deseq_results)) {
  if (!is.null(deseq_results[[ct]])) {
    te_res <- deseq_results[[ct]]$results %>% filter(Gene %in% te_genes)
    if (nrow(te_res) > 0) {
      cat(sprintf("    %s: ", ct))
      sig_g <- te_res$Gene[te_res$significant]
      cat(if (length(sig_g) > 0) paste(sig_g, collapse = ", ") else "none significant", "\n")
    }
  }
}

cat("\n  [PD vs Control] Module score:\n")
if (!is.null(master_comparisons)) {
  te_comp <- master_comparisons[master_comparisons$Variable == "TE_Control_Score", ]
  for (i in seq_len(nrow(te_comp))) {
    cat(sprintf("    %s: Beta=%.4f, p_adj=%.4f, Direction=%s\n",
                te_comp$Cell_Type[i], te_comp$lm_Beta_G1[i],
                te_comp$lm_p_adj[i], te_comp$Direction[i]))
  }
}

cat("\n  [Cross-module] TE_Control × dsRNA_Sensing (strongest known pair):\n")
te_sens <- master_correlations[
  (grepl("TE_Control", master_correlations$Var1) & grepl("Sensing", master_correlations$Var2)) |
  (grepl("TE_Control", master_correlations$Var2) & grepl("Sensing", master_correlations$Var1)), ]
if (nrow(te_sens) > 0) {
  for (i in seq_len(nrow(te_sens))) {
    cat(sprintf("    %s/%s: rho=%.3f, p_adj=%.4f\n",
                te_sens$Cell_Type[i], te_sens$Group[i],
                te_sens$rho[i], te_sens$p_adj[i]))
  }
}

cat("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")
