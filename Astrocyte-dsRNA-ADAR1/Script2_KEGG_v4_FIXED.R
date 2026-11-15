# ---------------------------------------------------------------------------------------------------------------
# Repository copy. Run from the repository root. Inputs: data/ (covariate-annotated Seurat objects from
# Full_merge_RDS_script_with_covariates.R) and results/<analysis output folder> (rename each timestamped output
# folder to the name used below, e.g. CX_TE_Control_<date> -> results/CX_TE_Control). R 4.4.3.
# ---------------------------------------------------------------------------------------------------------------
# ============================================================
# UNIFIED MODULE-DRIVEN GENE RANKING + fGSEA KEGG ENRICHMENT
# All 8 modules x Astrocytes vs Neurons (ALL donors only)
# ============================================================
#
# PURPOSE:
#   1. For each of 8 module scores, rank ALL expressed genes by
#      Spearman correlation (cell-level, covariate-residualised).
#   2. Summary table: key genes' ranks across all 8 modules.
#   3. fGSEA on KEGG for each module x cell-type ranking,
#      run SEPARATELY for KEGG_LEGACY and KEGG_MEDICUS.
#   4. Compare astrocyte vs neuron pathway enrichment.
#
# REVISION NOTES (v4):
#   - ADDED TE_Control module (8th module).
#   - FIXED KEGG pooling. Previously CP:KEGG_LEGACY and
#     CP:KEGG_MEDICUS were concatenated into one gene set object and
#     BH-corrected together. MEDICUS network sets are heavily nested
#     within canonical pathways and near-duplicates of one another
#     (pairwise Jaccard 0.85-0.98 for electron transport sets;
#     0.83-0.98 for 26S proteasome sets), which inflates the number
#     of tests and distorts BH correction. LEGACY and MEDICUS are now
#     run as two independent fGSEA analyses with independent BH
#     correction.
#   - ADDED expression filter to module scoring, matching the
#     min_frac filter already used for gene ranking, so that genes
#     too sparse to enter the ranking also do not contribute to the
#     module score. Dropped genes are logged.
#
# DESIGN:
#   - One unified loop (no Arms, no PD/Control split).
#   - ALL donors only.
#   - Covariate adjustment via residualisation.
#
# METHODS (for manuscript):
#   In astrocytes and neurons, all expressed genes (detected in
#   >=1% of cells) were ranked by Spearman correlation with each
#   of eight covariate-residualised module scores (ADAR1, ADAR2,
#   SNCA, dsRNA_Processing, dsRNA_EV_Transport, dsRNA_Sensing,
#   Innate_Immune, TE_Control). fGSEA (Korotkevich et al., 2021)
#   was run on each ranked list against KEGG pathway gene sets
#   from MSigDB, with C2:CP:KEGG_LEGACY (n = 186) and
#   C2:CP:KEGG_MEDICUS (n = 658) analysed and BH-adjusted
#   separately.
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
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(fgsea)
  library(msigdbr)
})

set.seed(42)

# =============================================================================
# SECTION 1: USER SETTINGS
# =============================================================================

cat("\n============================================================\n")
cat("=== UNIFIED MODULE RANKING + fGSEA KEGG ===\n")
cat("============================================================\n\n")

REGION <- "CX"   # "CX" or "SN"

# --- Covariate settings ---
COVARIATES_ENABLED <- TRUE
COV_AGE_COL <- "age"
COV_PMI_COL <- "pmi"
COV_SEX_COL <- "sex"

# --- File paths ---
INPUT_METHOD <- "path"   # "dialog" or "path"

if (INPUT_METHOD == "path") {
  results_rds <- "results/CX_TE_Control/CX_8Module_AllResults.rds"

  
  if (REGION == "CX") {
    astro_file  <- "data/CX_Astrocytes_withCov.rds"
    neuron_file <- "data/CX_Neurons_Total_withCov.rds"
    output_base <- "results"
  } else {
    astro_file  <- "data/03_astrocytes_SN.rds"
    neuron_file <- "data/03_neurons_SN.rds"
    output_base <- "results"
  }
} else {
  library(tcltk)
  cat("Select the 8-Module AllResults RDS file...\n")
  results_rds <- tk_choose.files(caption = "Select 8Module_AllResults.rds",
                                 filters = matrix(c("RDS Files", ".rds"), 1, 2))
  if (length(results_rds) == 0 || is.na(results_rds)) stop("No RDS file selected.")
  cat(paste0("Select ", REGION, " ASTROCYTE RDS...\n"))
  astro_file <- tk_choose.files(caption = paste0("Select ", REGION, " Astrocyte RDS"),
                                filters = matrix(c("RDS Files", ".rds"), 1, 2))
  if (length(astro_file) == 0 || is.na(astro_file)) stop("No file selected.")
  cat(paste0("Select ", REGION, " NEURON RDS...\n"))
  neuron_file <- tk_choose.files(caption = paste0("Select ", REGION, " Neuron RDS"),
                                 filters = matrix(c("RDS Files", ".rds"), 1, 2))
  if (length(neuron_file) == 0 || is.na(neuron_file)) stop("No file selected.")
  cat("Select OUTPUT folder...\n")
  output_base <- tk_choose.dir(caption = "Select OUTPUT folder")
  if (is.na(output_base)) stop("No folder selected.")
}

# fGSEA parameters
MIN_GENE_SET_SIZE <- 15
MAX_GENE_SET_SIZE <- 500
FGSEA_NPERM       <- 10000
PADJ_CUTOFF       <- 0.05
MIN_EXPR_FRAC     <- 0.01

# Key genes to track across all module rankings
KEY_GENES <- c(
  # Sensors
  "IFIH1", "DDX58", "DHX58", "TLR3",
  # Adaptors / kinases
  "MAVS", "TICAM1", "TBK1", "IKBKE",
  # ADAR / editing
  "ADAR", "ADARB1",
  # SNCA
  "SNCA",
  # Processing
  "ILF2", "ILF3", "DHX9", "STRBP", "ZFR", "DDX39B", "SFPQ", "NONO",
  # Transport / EV
  "SIDT1", "SIDT2", "LAMP2", "RAB27A", "RAB27B", "CD63", "CD9", "CD81",
  "SMPD3", "RAB35", "PDCD6IP", "TSG101", "SDCBP", "MVP",
  # Innate immune
  "JAK1", "TYK2", "STAT1", "STAT2", "IRF9", "RELA",
  "ISG15", "USP18", "OAS1", "OAS2", "OAS3", "OASL",
  "RNASEL", "EIF2AK2", "PRKRA",
  "MX1", "IFIT1", "IFIT2", "IFIT3",
  "C3", "C4A",
  # Complement / reactivity extras
  "C1QA", "C1QB", "C1QC",
  "GFAP", "VIM"
)

# Create output directory
timestamp   <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- file.path(output_base, paste0(REGION, "_Unified_fGSEA_", timestamp))
dir.create(file.path(results_dir, "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", results_dir, "\n\n")

# =============================================================================
# SECTION 2: MODULE DEFINITIONS (V4 -- 8 modules)
# =============================================================================

cat("========== DEFINING 8 MODULES (V4) ==========\n")

module_gene_sets <- list(
  ADAR1 = "ADAR",
  ADAR2 = "ADARB1",
  SNCA  = "SNCA",
  
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
  # MODULE 8: TE CONTROL
  # NOTE: L1TD1 and PIWIL1 are defined here for completeness but are
  # expressed in only ~0.2% of cortical nuclei and are therefore
  # dropped by the expression filter in calculate_module_scores()
  # (see SECTION 5). The effective module is the remaining 6 genes.
  # This must be stated in the manuscript Methods.
  TE_Control = c(
    "MOV10", "TRIM28", "SETDB1", "DNMT1",
    "DICER1", "AGO2", "L1TD1", "PIWIL1"
  )
)

module_names <- names(module_gene_sets)
cat("Modules (", length(module_names), "):",
    paste(module_names, collapse = ", "), "\n\n")

# =============================================================================
# SECTION 3: LOAD KEGG PATHWAYS
# =============================================================================

cat("========== LOADING KEGG PATHWAYS ==========\n")

# --- Load LEGACY and MEDICUS SEPARATELY (do not concatenate) ---
# The previous version broke out of the first loop then appended the remaining
# subcollections, silently pooling them into one object that was BH-corrected
# as a single family of tests. See revision notes at the top of this file.

load_kegg <- function(subcat, tag) {
  df <- tryCatch(
    msigdbr(species = "Homo sapiens", category = "C2", subcategory = subcat),
    error = function(e) NULL
  )
  if (!is.null(df) && nrow(df) > 0) {
    cat("  Loaded", subcat, ":", length(unique(df$gs_name)), "gene sets\n")
    df$kegg_source <- tag
    return(df)
  }
  cat("  Not available:", subcat, "\n")
  NULL
}

kegg_legacy_df  <- load_kegg("CP:KEGG_LEGACY",  "LEGACY")
if (is.null(kegg_legacy_df)) kegg_legacy_df <- load_kegg("CP:KEGG", "LEGACY")
kegg_medicus_df <- load_kegg("CP:KEGG_MEDICUS", "MEDICUS")

if (is.null(kegg_legacy_df) && is.null(kegg_medicus_df)) {
  all_cp <- tryCatch(
    msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP"),
    error = function(e) NULL
  )
  if (!is.null(all_cp)) {
    kegg_all <- all_cp %>% filter(grepl("^KEGG_", gs_name))
    kegg_legacy_df  <- kegg_all %>% filter(!grepl("^KEGG_MEDICUS_", gs_name)) %>%
      mutate(kegg_source = "LEGACY")
    kegg_medicus_df <- kegg_all %>% filter(grepl("^KEGG_MEDICUS_", gs_name)) %>%
      mutate(kegg_source = "MEDICUS")
  }
}
if (is.null(kegg_legacy_df) || nrow(kegg_legacy_df) == 0) {
  stop("Could not load KEGG LEGACY pathways. Check msigdbr_collections().")
}

# Two SEPARATE pathway lists -- each gets its own fGSEA run and BH correction
kegg_legacy  <- split(kegg_legacy_df$gene_symbol,  kegg_legacy_df$gs_name)
kegg_medicus <- if (!is.null(kegg_medicus_df) && nrow(kegg_medicus_df) > 0) {
  split(kegg_medicus_df$gene_symbol, kegg_medicus_df$gs_name)
} else NULL

# Combined lookup, retained ONLY so that plotEnrichment() calls further down
# can resolve any pathway name. It is never passed to fgsea().
kegg_pathways <- c(kegg_legacy, if (!is.null(kegg_medicus)) kegg_medicus)

kegg_source_map <- c(
  setNames(rep("LEGACY",  length(kegg_legacy)),  names(kegg_legacy)),
  if (!is.null(kegg_medicus))
    setNames(rep("MEDICUS", length(kegg_medicus)), names(kegg_medicus))
)

cat("KEGG gene sets:\n")
cat("  LEGACY :", length(kegg_legacy), "\n")
cat("  MEDICUS:", ifelse(is.null(kegg_medicus), 0, length(kegg_medicus)), "\n\n")

# --- Redundancy diagnostics (for manuscript Methods / reviewer response) ---
if (!is.null(kegg_medicus) && "KEGG_OXIDATIVE_PHOSPHORYLATION" %in% names(kegg_legacy)) {
  etc_sets <- grep("MEDICUS.*ELECTRON_TRANSFER", names(kegg_medicus), value = TRUE)
  if (length(etc_sets) > 1) {
    oxphos <- kegg_legacy[["KEGG_OXIDATIVE_PHOSPHORYLATION"]]
    pct <- sapply(kegg_medicus[etc_sets],
                  function(g) 100 * length(intersect(g, oxphos)) / length(g))
    cat("MEDICUS electron-transport sets: % genes also in KEGG_OXIDATIVE_PHOSPHORYLATION\n")
    cat("  n =", length(etc_sets), "| range",
        round(min(pct), 1), "-", round(max(pct), 1),
        "| median", round(median(pct), 1), "\n\n")
  }
}


# =============================================================================
# SECTION 4: LOAD DATA
# =============================================================================

cat("========== LOADING DATA ==========\n")

cat("Loading 8-Module results RDS...\n")
prev_results <- readRDS(results_rds)
cat("  Objects:", paste(names(prev_results), collapse = ", "), "\n\n")

cat("Loading Astrocyte Seurat...\n")
astro <- readRDS(astro_file)
cat("  Astrocytes:", ncol(astro), "cells\n")

cat("Loading Neuron Seurat...\n")
neuron <- readRDS(neuron_file)
cat("  Neurons:", ncol(neuron), "cells\n\n")

# --- Detect donor & condition columns ---
find_col <- function(meta, candidates) {
  for (nm in candidates) if (nm %in% colnames(meta)) return(nm)
  return(NULL)
}

donor_col_astro  <- find_col(astro@meta.data,  c("donor_id","Donor","donor","subject_id","sample_id","orig.ident"))
donor_col_neuron <- find_col(neuron@meta.data, c("donor_id","Donor","donor","subject_id","sample_id","orig.ident"))
cond_col_astro   <- find_col(astro@meta.data,  c("condition","Condition","disease","Disease","group","Group","diagnosis"))
cond_col_neuron  <- find_col(neuron@meta.data, c("condition","Condition","disease","Disease","group","Group","diagnosis"))

cat("Donor columns: astro =", donor_col_astro, "| neuron =", donor_col_neuron, "\n")
cat("Condition columns: astro =", cond_col_astro, "| neuron =", cond_col_neuron, "\n\n")

# --- Detect covariates ---
covariate_cols_to_use <- c()
if (COVARIATES_ENABLED) {
  for (cv in c(COV_AGE_COL, COV_PMI_COL, COV_SEX_COL)) {
    if (!is.null(find_col(astro@meta.data, cv))) {
      covariate_cols_to_use <- c(covariate_cols_to_use, cv)
      cat("  Covariate found:", cv, "\n")
    }
  }
  if (length(covariate_cols_to_use) == 0) {
    cat("  WARNING: No covariates found. Rankings will be unadjusted.\n")
    COVARIATES_ENABLED <- FALSE
  } else {
    cat("  Active covariates:", paste(covariate_cols_to_use, collapse = ", "), "\n")
  }
}

# =============================================================================
# SECTION 5: HELPER FUNCTIONS
# =============================================================================

cat("\n========== DEFINING HELPER FUNCTIONS ==========\n")

calculate_module_scores <- function(seurat_obj, module_list, min_frac = 0.01) {
  # min_frac matches the filter used in compute_gene_rankings(), so that a
  # gene too sparse to enter the ranked list also does not contribute to
  # the module score it is being correlated against.
  expr <- GetAssayData(seurat_obj, layer = "data")
  gene_frac <- rowMeans(expr > 0)
  
  scores <- sapply(names(module_list), function(mod) {
    genes <- module_list[[mod]]
    genes <- genes[genes %in% rownames(expr)]
    
    dropped <- genes[gene_frac[genes] <  min_frac]
    genes   <- genes[gene_frac[genes] >= min_frac]
    
    if (length(dropped) > 0) {
      cat("    [", mod, "] DROPPED (<", min_frac * 100, "% of cells): ",
          paste0(dropped, " (", round(gene_frac[dropped] * 100, 2), "%)",
                 collapse = ", "), "\n", sep = "")
    }
    cat("    [", mod, "] scored on ", length(genes), " gene(s): ",
        paste(genes, collapse = ", "), "\n", sep = "")
    
    if (length(genes) == 0) return(rep(NA_real_, ncol(seurat_obj)))
    if (length(genes) == 1) return(as.numeric(expr[genes, ]))
    colMeans(expr[genes, , drop = FALSE])
  })
  colnames(scores) <- paste0(names(module_list), "_Score")
  as.data.frame(scores)
}

residualise_score <- function(seurat_obj, score_vec, donor_col, covariate_cols) {
  if (!COVARIATES_ENABLED || length(covariate_cols) == 0) return(score_vec)
  meta <- as.data.frame(seurat_obj@meta.data)
  donor_ids <- as.character(meta[[donor_col]])
  cov_df <- as.data.frame(lapply(covariate_cols, function(cv) {
    if (!cv %in% colnames(meta)) return(rep(NA, nrow(meta)))
    vals <- meta[[cv]]
    if (is.numeric(vals)) {
      donor_means <- tapply(vals, donor_ids, mean, na.rm = TRUE)
      ifelse(is.na(vals), donor_means[donor_ids], vals)
    } else {
      as.character(vals)
    }
  }))
  colnames(cov_df) <- covariate_cols
  valid_covs <- covariate_cols[sapply(covariate_cols, function(cv) {
    v <- cov_df[[cv]]
    if (all(is.na(v))) return(FALSE)
    if (is.numeric(v)) sd(v, na.rm = TRUE) > 0 else length(unique(na.omit(v))) > 1
  })]
  if (length(valid_covs) == 0) return(score_vec)
  cov_df <- cov_df[, valid_covs, drop = FALSE]
  for (cv in valid_covs) {
    if (!is.numeric(cov_df[[cv]])) cov_df[[cv]] <- as.factor(cov_df[[cv]])
    else if (any(is.na(cov_df[[cv]]))) cov_df[[cv]][is.na(cov_df[[cv]])] <- mean(cov_df[[cv]], na.rm = TRUE)
  }
  lm_df <- data.frame(score = score_vec, cov_df)
  complete_rows <- complete.cases(lm_df)
  form <- as.formula(paste0("score ~ ", paste(valid_covs, collapse = " + ")))
  fit <- tryCatch(lm(form, data = lm_df[complete_rows, ]), error = function(e) NULL)
  if (is.null(fit)) return(score_vec)
  resid_vec <- rep(NA, length(score_vec))
  resid_vec[complete_rows] <- residuals(fit)
  grand_mean <- mean(score_vec[complete_rows], na.rm = TRUE)
  resid_vec[!complete_rows] <- score_vec[!complete_rows] - grand_mean
  resid_vec
}

compute_gene_rankings <- function(seurat_obj, score_vec, min_frac = 0.01) {
  expr <- GetAssayData(seurat_obj, layer = "data")
  gene_frac <- rowMeans(expr > 0)
  keep_genes <- names(gene_frac[gene_frac >= min_frac])
  cat("    Genes passing filter (>=", min_frac * 100, "% cells):",
      length(keep_genes), "/", nrow(expr), "\n")
  expr_sub <- expr[keep_genes, , drop = FALSE]
  
  # 保持稀疏,分块计算,避免 2.6 GB 稠密矩阵
  stopifnot(length(score_vec) == ncol(expr))
  ok <- !is.na(score_vec)
  sv <- score_vec[ok]
  sv_rank <- rank(sv)
  
  chunk <- 2000
  starts <- seq(1, nrow(expr_sub), by = chunk)
  res_list <- vector("list", length(starts))
  
  for (k in seq_along(starts)) {
    idx <- starts[k]:min(starts[k] + chunk - 1, nrow(expr_sub))
    block <- as.matrix(expr_sub[idx, ok, drop = FALSE])
    res_list[[k]] <- apply(block, 1, function(g) {
      suppressWarnings(cor(rank(g), sv_rank, method = "pearson"))
    })
    rm(block); gc(verbose = FALSE)
  }
  rankings <- unlist(res_list)
  rankings <- rankings[!is.na(rankings)]
  sort(rankings, decreasing = TRUE)
}

run_fgsea_analysis <- function(gene_rankings, pathways,
                               min_size = 15, max_size = 500, nperm = 10000) {
  fgsea_res <- fgsea(pathways    = pathways,
                     stats       = gene_rankings,
                     minSize     = min_size,
                     maxSize     = max_size,
                     nPermSimple = nperm)
  fgsea_res <- fgsea_res %>%
    arrange(padj) %>%
    mutate(
      Direction = ifelse(NES > 0, "Positive", "Negative"),
      Sig_Level = case_when(
        padj < 0.001 ~ "***",
        padj < 0.01  ~ "**",
        padj < 0.05  ~ "*",
        TRUE ~ "NS"
      ),
      Pathway_Short = gsub("^KEGG_", "", pathway) %>%
        gsub("_", " ", .) %>%
        tools::toTitleCase()
    )
  as.data.frame(fgsea_res)
}

plot_fgsea_barplot <- function(fgsea_df, n_top = 20, title = "fGSEA Results") {
  sig_df <- fgsea_df %>%
    filter(padj < 0.1) %>%
    arrange(desc(abs(NES))) %>%
    head(n_top) %>%
    arrange(NES)
  if (nrow(sig_df) == 0) { cat("    No significant pathways for plot.\n"); return(NULL) }
  sig_df$Pathway_Short <- factor(sig_df$Pathway_Short, levels = sig_df$Pathway_Short)
  ggplot(sig_df, aes(x = NES, y = Pathway_Short, fill = Direction)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("Positive" = "#D73027", "Negative" = "#4575B4")) +
    geom_text(aes(label = paste0("q=", formatC(padj, format = "e", digits = 1))),
              hjust = ifelse(sig_df$NES > 0, -0.1, 1.1), size = 2.8) +
    labs(title = title,
         subtitle = paste0("Top ", nrow(sig_df), " pathways by |NES| (padj < 0.1)"),
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 13),
          axis.text.y = element_text(size = 8))
}

plot_nes_heatmap <- function(all_results_list, pathway_filter = NULL,
                             padj_cutoff = 0.05, title = "NES Heatmap") {
  sig_pathways <- unique(unlist(lapply(all_results_list, function(df) {
    df$pathway[df$padj < padj_cutoff]
  })))
  if (!is.null(pathway_filter)) sig_pathways <- intersect(sig_pathways, pathway_filter)
  if (length(sig_pathways) < 2) { cat("    <2 sig pathways for heatmap.\n"); return(NULL) }
  
  nes_mat  <- matrix(0, nrow = length(sig_pathways), ncol = length(all_results_list),
                     dimnames = list(gsub("^KEGG_", "", sig_pathways),
                                     names(all_results_list)))
  padj_mat <- matrix(1, nrow = length(sig_pathways), ncol = length(all_results_list),
                     dimnames = dimnames(nes_mat))
  for (nm in names(all_results_list)) {
    df <- all_results_list[[nm]]
    for (pw in sig_pathways) {
      row <- df[df$pathway == pw, ]
      if (nrow(row) == 1) {
        nes_mat[gsub("^KEGG_", "", pw), nm]  <- row$NES
        padj_mat[gsub("^KEGG_", "", pw), nm] <- row$padj
      }
    }
  }
  sig_labels <- matrix("", nrow = nrow(nes_mat), ncol = ncol(nes_mat))
  sig_labels[padj_mat < 0.001] <- "***"
  sig_labels[padj_mat >= 0.001 & padj_mat < 0.01] <- "**"
  sig_labels[padj_mat >= 0.01 & padj_mat < 0.05] <- "*"
  rownames(nes_mat) <- gsub("_", " ", rownames(nes_mat))
  rownames(sig_labels) <- rownames(nes_mat)
  pheatmap(nes_mat,
           color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
           breaks = seq(-max(abs(nes_mat), na.rm = TRUE),
                        max(abs(nes_mat), na.rm = TRUE), length.out = 101),
           display_numbers = sig_labels,
           number_color = "black",
           fontsize_number = 12, fontsize_row = 8, fontsize_col = 10,
           cluster_cols = FALSE, cluster_rows = TRUE, angle_col = 45,
           main = title)
}

cat("Helper functions defined.\n\n")

# =============================================================================
# SECTION 6: COMPUTE CELL-LEVEL MODULE SCORES
# =============================================================================

cat("========== COMPUTING CELL-LEVEL MODULE SCORES ==========\n")

filter_modules <- function(obj, module_list) {
  present <- rownames(obj)
  filtered <- lapply(module_list, function(genes) unique(genes[genes %in% present]))
  valid <- names(filtered)[sapply(filtered, length) >= 1]
  cat("  Valid modules:", length(valid), "/", length(filtered), "\n")
  filtered[valid]
}

astro_mods  <- filter_modules(astro, module_gene_sets)
neuron_mods <- filter_modules(neuron, module_gene_sets)

astro_scores  <- calculate_module_scores(astro, astro_mods)
neuron_scores <- calculate_module_scores(neuron, neuron_mods)

astro  <- AddMetaData(astro, astro_scores)
neuron <- AddMetaData(neuron, neuron_scores)

cat("Module scores added.\n")
cat("  Astro:", paste(colnames(astro_scores), collapse = ", "), "\n")
cat("  Neuron:", paste(colnames(neuron_scores), collapse = ", "), "\n\n")

# =============================================================================
# SECTION 7: MAIN LOOP -- GENE RANKING + fGSEA (ALL donors, 8 modules x 2 CT)
# =============================================================================

cat("\n============================================================\n")
cat("=== MAIN ANALYSIS: 8 MODULES x 2 CELL TYPES ===\n")
cat("============================================================\n\n")

all_rankings      <- list()
all_fgsea         <- list()   # LEGACY (primary)
all_fgsea_medicus <- list()   # MEDICUS (corroboration)
key_gene_rows     <- list()

cell_type_list <- list(
  Astrocytes = list(obj = astro,  donor_col = donor_col_astro),
  Neurons    = list(obj = neuron, donor_col = donor_col_neuron)
)

for (ct_name in names(cell_type_list)) {
  ct <- cell_type_list[[ct_name]]
  obj <- ct$obj
  dc  <- ct$donor_col
  
  available_scores <- grep("_Score$", colnames(obj@meta.data), value = TRUE)
  cat("=== Cell type:", ct_name, "| Available scores:", paste(available_scores, collapse = ", "), "===\n\n")
  
  for (mod_name in module_names) {
    score_col <- paste0(mod_name, "_Score")
    if (!score_col %in% colnames(obj@meta.data)) {
      cat("  Skipping", mod_name, "- score not found.\n\n")
      next
    }
    
    label <- paste0(ct_name, "_", mod_name)
    cat("--- ", label, " ---\n")
    
    score_vec <- obj@meta.data[[score_col]]
    if (all(is.na(score_vec)) || sd(score_vec, na.rm = TRUE) == 0) {
      cat("  No variance in score, skipping.\n\n")
      next
    }
    
    # Residualise
    if (COVARIATES_ENABLED && length(covariate_cols_to_use) > 0) {
      cat("  Residualising", score_col, "...\n")
      score_vec <- residualise_score(obj, score_vec, dc, covariate_cols_to_use)
    }
    
    # Rank all genes
    cat("  Ranking genes...\n")
    rankings <- compute_gene_rankings(obj, score_vec, min_frac = MIN_EXPR_FRAC)
    all_rankings[[label]] <- rankings
    
    # --- Key gene tracking ---
    total_genes <- length(rankings)
    gene_names_ranked <- names(rankings)
    for (kg in KEY_GENES) {
      if (kg %in% gene_names_ranked) {
        rank_pos <- which(gene_names_ranked == kg)
        rho_val  <- rankings[kg]
        key_gene_rows[[length(key_gene_rows) + 1]] <- data.frame(
          Gene      = kg,
          Cell_Type = ct_name,
          Module    = mod_name,
          Rank      = rank_pos,
          Total_Genes = total_genes,
          Percentile = round(100 * (1 - rank_pos / total_genes), 1),
          Spearman_r = round(rho_val, 4),
          stringsAsFactors = FALSE
        )
      }
    }
    
    # --- fGSEA: LEGACY and MEDICUS run and BH-corrected SEPARATELY ---
    annotate_res <- function(df, src) {
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$Region       <- REGION
      df$Cell_Type    <- ct_name
      df$Module       <- mod_name
      df$kegg_source  <- src
      df$Cov_Adjusted <- COVARIATES_ENABLED && length(covariate_cols_to_use) > 0
      df
    }
    
    cat("  Running fGSEA [LEGACY]...\n")
    res_legacy <- annotate_res(
      run_fgsea_analysis(rankings, kegg_legacy,
                         min_size = MIN_GENE_SET_SIZE,
                         max_size = MAX_GENE_SET_SIZE,
                         nperm    = FGSEA_NPERM),
      "LEGACY"
    )
    
    res_medicus <- NULL
    if (!is.null(kegg_medicus)) {
      cat("  Running fGSEA [MEDICUS]...\n")
      res_medicus <- annotate_res(
        run_fgsea_analysis(rankings, kegg_medicus,
                           min_size = MIN_GENE_SET_SIZE,
                           max_size = MAX_GENE_SET_SIZE,
                           nperm    = FGSEA_NPERM),
        "MEDICUS"
      )
    }
    
    # Primary result object = LEGACY (used for main figures)
    fgsea_res <- res_legacy
    
    all_fgsea[[label]]         <- res_legacy
    all_fgsea_medicus[[label]] <- res_medicus
    
    cat("  Significant [LEGACY] (padj <", PADJ_CUTOFF, "):",
        sum(res_legacy$padj < PADJ_CUTOFF, na.rm = TRUE), "\n")
    if (!is.null(res_medicus)) {
      cat("  Significant [MEDICUS] (padj <", PADJ_CUTOFF, "):",
          sum(res_medicus$padj < PADJ_CUTOFF, na.rm = TRUE), "\n")
    }
    
    
    # --- Save per-run tables ---
    write.csv(data.frame(Gene = names(rankings), Spearman_r = rankings),
              file.path(results_dir, "tables",
                        paste0(REGION, "_GeneRankings_", label, ".csv")),
              row.names = FALSE)
    
    write.csv(fgsea_res %>% select(-leadingEdge),
              file.path(results_dir, "tables",
                        paste0(REGION, "_fGSEA_FullKEGG_", label, ".csv")),
              row.names = FALSE)
    
    # --- Bar plot ---
    p <- plot_fgsea_barplot(fgsea_res, n_top = 25,
                            title = paste0(REGION, " ", ct_name,
                                           "\nAll KEGG ranked by cor(", mod_name, ")"))
    if (!is.null(p)) {
      ggsave(file.path(results_dir, "figures",
                       paste0(REGION, "_BarPlot_", label, ".png")),
             p, width = 10, height = 8, dpi = 300)
    }
    
    # --- Enrichment curves for top 4 ---
    top4 <- fgsea_res %>% filter(padj < PADJ_CUTOFF) %>% arrange(padj) %>% head(4)
    if (nrow(top4) > 0) {
      e_plots <- list()
      for (i in seq_len(nrow(top4))) {
        pw <- top4$pathway[i]
        if (pw %in% names(kegg_pathways)) {
          e_plots[[pw]] <- plotEnrichment(kegg_pathways[[pw]], rankings) +
            labs(title = paste0(gsub("^KEGG_", "", pw)),
                 subtitle = paste0(ct_name, " | cor(", mod_name, ")")) +
            theme_minimal(base_size = 12)
        }
      }
      if (length(e_plots) > 0) {
        ggsave(file.path(results_dir, "figures",
                         paste0(REGION, "_Enrichment_", label, ".png")),
               wrap_plots(e_plots, ncol = 2),
               width = 12, height = 3 * ceiling(length(e_plots) / 2), dpi = 300)
      }
    }
    
    cat("\n")
  }
}

# =============================================================================
# SECTION 8: KEY GENE SUMMARY TABLE
# =============================================================================

cat("========== KEY GENE SUMMARY TABLE ==========\n")

key_gene_df <- do.call(rbind, key_gene_rows)

if (!is.null(key_gene_df) && nrow(key_gene_df) > 0) {
  for (ct in unique(key_gene_df$Cell_Type)) {
    ct_sub <- key_gene_df %>% filter(Cell_Type == ct)
    
    rank_wide <- ct_sub %>%
      select(Gene, Module, Rank) %>%
      pivot_wider(names_from = Module, values_from = Rank, names_prefix = "Rank_")
    
    pctl_wide <- ct_sub %>%
      select(Gene, Module, Percentile) %>%
      pivot_wider(names_from = Module, values_from = Percentile, names_prefix = "Pctl_")
    
    rho_wide <- ct_sub %>%
      select(Gene, Module, Spearman_r) %>%
      pivot_wider(names_from = Module, values_from = Spearman_r, names_prefix = "Rho_")
    
    summary_wide <- rank_wide %>%
      left_join(pctl_wide, by = "Gene") %>%
      left_join(rho_wide, by = "Gene")
    
    write.csv(summary_wide,
              file.path(results_dir, "tables",
                        paste0(REGION, "_KeyGene_Summary_", ct, ".csv")),
              row.names = FALSE)
    cat("  Saved key gene summary for", ct, ":", nrow(summary_wide), "genes\n")
  }
  
  write.csv(key_gene_df,
            file.path(results_dir, "tables",
                      paste0(REGION, "_KeyGene_AllModules_Long.csv")),
            row.names = FALSE)
  cat("  Saved long-format key gene table.\n")
} else {
  cat("  WARNING: No key gene data collected.\n")
}

# =============================================================================
# SECTION 9: COMBINED fGSEA TABLE + NES HEATMAPS
# =============================================================================

cat("\n========== COMBINED RESULTS & HEATMAPS ==========\n")

fgsea_combined <- do.call(rbind, all_fgsea)
if (!is.null(fgsea_combined) && nrow(fgsea_combined) > 0) {
  fgsea_export <- fgsea_combined %>% select(-leadingEdge)
  write.csv(fgsea_export,
            file.path(results_dir, "tables",
                      paste0(REGION, "_fGSEA_AllModules_Combined.csv")),
            row.names = FALSE)
  cat("Combined fGSEA table:", nrow(fgsea_export), "rows saved.\n")
}

# --- NES heatmap: Astrocytes ---
astro_fgsea <- all_fgsea[grep("^Astrocytes_", names(all_fgsea))]
if (length(astro_fgsea) >= 2) {
  png(file.path(results_dir, "figures",
                paste0(REGION, "_NES_Heatmap_Astrocytes.png")),
      width = 14, height = 14, units = "in", res = 300)
  tryCatch(
    plot_nes_heatmap(astro_fgsea, padj_cutoff = PADJ_CUTOFF,
                     title = paste0(REGION, " Astrocytes: KEGG NES across 8 Module Rankings")),
    error = function(e) {
      plot.new(); text(0.5, 0.5, paste("Error:", e$message), cex = 0.8)
    }
  )
  dev.off()
  cat("  Astrocyte NES heatmap saved.\n")
}

# --- NES heatmap: Neurons ---
neuron_fgsea <- all_fgsea[grep("^Neurons_", names(all_fgsea))]
if (length(neuron_fgsea) >= 2) {
  png(file.path(results_dir, "figures",
                paste0(REGION, "_NES_Heatmap_Neurons.png")),
      width = 14, height = 14, units = "in", res = 300)
  tryCatch(
    plot_nes_heatmap(neuron_fgsea, padj_cutoff = PADJ_CUTOFF,
                     title = paste0(REGION, " Neurons: KEGG NES across 8 Module Rankings")),
    error = function(e) {
      plot.new(); text(0.5, 0.5, paste("Error:", e$message), cex = 0.8)
    }
  )
  dev.off()
  cat("  Neuron NES heatmap saved.\n")
}

# --- NES heatmap: Astro vs Neuron combined ---
if (length(all_fgsea) >= 4) {
  png(file.path(results_dir, "figures",
                paste0(REGION, "_NES_Heatmap_AstroVsNeuron_AllModules.png")),
      width = 18, height = 16, units = "in", res = 300)
  tryCatch(
    plot_nes_heatmap(all_fgsea, padj_cutoff = PADJ_CUTOFF,
                     title = paste0(REGION, ": Astrocyte vs Neuron KEGG NES\nacross all 8 Module Rankings")),
    error = function(e) {
      plot.new(); text(0.5, 0.5, paste("Error:", e$message), cex = 0.8)
    }
  )
  dev.off()
  cat("  Combined Astro vs Neuron NES heatmap saved.\n")
}

# =============================================================================
# SECTION 10: EXCEL WORKBOOK
# =============================================================================

cat("\n========== SAVING EXCEL WORKBOOK ==========\n")

wb <- createWorkbook()

if (!is.null(key_gene_df) && nrow(key_gene_df) > 0) {
  addWorksheet(wb, "KeyGene_Summary")
  writeData(wb, "KeyGene_Summary", key_gene_df)
}

if (!is.null(fgsea_combined) && nrow(fgsea_combined) > 0) {
  sig_only <- fgsea_combined %>%
    filter(padj < PADJ_CUTOFF) %>%
    select(-leadingEdge) %>%
    arrange(Cell_Type, Module, padj)
  addWorksheet(wb, "fGSEA_Significant")
  writeData(wb, "fGSEA_Significant", sig_only)
  
  addWorksheet(wb, "fGSEA_All")
  writeData(wb, "fGSEA_All", fgsea_combined %>% select(-leadingEdge))
}

methods_df <- data.frame(
  Item = c("Analysis", "Modules", "Cell types", "Group", "Gene filter",
           "Score adjustment", "Ranking", "fGSEA", "KEGG source"),
  Description = c(
    "Unified module-driven gene ranking + fGSEA KEGG enrichment",
    paste(module_names, collapse = ", "),
    "Astrocytes, Neurons (Total)",
    "ALL donors (PD + Control combined)",
    paste0("Genes expressed in >= ", MIN_EXPR_FRAC * 100, "% of cells"),
    paste0("Module scores residualised against: ",
           if (COVARIATES_ENABLED) paste(covariate_cols_to_use, collapse = ", ") else "none"),
    "Per-gene Spearman correlation with residualised module score (cell-level; computed as Pearson correlation of ranks)",
    paste0("fgsea(nPermSimple=", FGSEA_NPERM, ", minSize=", MIN_GENE_SET_SIZE,
           ", maxSize=", MAX_GENE_SET_SIZE, "); BH p-adjustment"),
    "MSigDB C2:CP:KEGG_LEGACY (n=186) and C2:CP:KEGG_MEDICUS (n=658), analysed and BH-adjusted SEPARATELY"
  )
)
addWorksheet(wb, "Methods")
writeData(wb, "Methods", methods_df)

xlsx_file <- paste0(REGION, "_Unified_fGSEA_Summary.xlsx")
old_wd <- getwd()
setwd(results_dir)
saveWorkbook(wb, xlsx_file, overwrite = TRUE)
setwd(old_wd)
cat("Excel workbook saved:", xlsx_file, "\n")

# =============================================================================
# SECTION 11: SAVE RDS
# =============================================================================

cat("\n========== SAVING RDS ==========\n")

saveRDS(
  list(
    region         = REGION,
    modules        = module_gene_sets,
    key_gene_table = key_gene_df,
    gene_rankings  = all_rankings,
    fgsea_results  = all_fgsea,
    fgsea_results_medicus = all_fgsea_medicus,
    kegg_source_map       = kegg_source_map,
    fgsea_combined = fgsea_combined,
    parameters     = list(
      min_expr_frac  = MIN_EXPR_FRAC,
      fgsea_nperm    = FGSEA_NPERM,
      min_set_size   = MIN_GENE_SET_SIZE,
      max_set_size   = MAX_GENE_SET_SIZE,
      padj_cutoff    = PADJ_CUTOFF,
      covariates     = covariate_cols_to_use,
      cov_adjusted   = COVARIATES_ENABLED
    )
  ),
  file.path(results_dir, paste0(REGION, "_Unified_fGSEA_AllResults.rds"))
)
cat("RDS saved.\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n============================================================\n")
cat("=== ANALYSIS COMPLETE ===\n")
cat("============================================================\n")
cat("Region:", REGION, "\n")
cat("Output:", results_dir, "\n\n")

cat("=== OUTPUT FILES ===\n")
cat("  tables/\n")
cat("    *_GeneRankings_*.csv          Full gene rankings per module x cell type\n")
cat("    *_KeyGene_Summary_*.csv       Key genes rank/percentile/rho across modules\n")
cat("    *_KeyGene_AllModules_Long.csv  Long format for all key genes\n")
cat("    *_fGSEA_FullKEGG_*.csv        fGSEA results per module x cell type\n")
cat("    *_fGSEA_AllModules_Combined   All fGSEA results in one table\n")
cat("  figures/\n")
cat("    *_BarPlot_*.png               Top KEGG pathways per module x cell type\n")
cat("    *_Enrichment_*.png            Enrichment curves for top hits\n")
cat("    *_NES_Heatmap_*.png           NES heatmaps (astro, neuron, combined)\n")
cat("  *_Unified_fGSEA_Summary.xlsx    Excel workbook\n")
cat("  *_Unified_fGSEA_AllResults.rds  Full results object\n")

cat("\n=== INTERPRETATION GUIDE ===\n")
cat("1. Key Gene Summary: Where do IFIH1, DDX58, SIDT1/2 etc. rank\n")
cat("   across all 8 modules in astrocytes vs neurons?\n")
cat("   High percentile (>90) = strong co-variation with that module.\n")
cat("2. NES Heatmaps: Which KEGG pathways are enriched in astrocyte\n")
cat("   vs neuron rankings for each module? Cell-type-specific coupling\n")
cat("   visible when a pathway is significant in one CT but not the other.\n")
cat("3. Bar Plots: Dominant biological programmes co-regulated with\n")
cat("   each module (endosomal, IFN, complement, synaptic, etc.)\n")

cat("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")
