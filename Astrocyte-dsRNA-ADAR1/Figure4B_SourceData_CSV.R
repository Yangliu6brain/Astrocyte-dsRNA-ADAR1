# ---------------------------------------------------------------------------------------------------------------
# Repository copy. Run from the repository root. Inputs: data/ (covariate-annotated Seurat objects from
# Full_merge_RDS_script_with_covariates.R) and results/<analysis output folder> (rename each timestamped output
# folder to the name used below, e.g. CX_TE_Control_<date> -> results/CX_TE_Control). R 4.4.3.
# ---------------------------------------------------------------------------------------------------------------
# ============================================================
# MS2 FIGURE 4B SOURCE DATA — PD vs CONTROL MODULE SCORES
# Exports donor-level scores + covariate-adjusted lm() results
# Companion to Figure4_SourceData_CSV.R (A-panel exporter)
# ============================================================
#
# Model:
#   lm(score ~ diagnosis + age + PMI + sex)
#   — performed per module, per cell type
#
# Outputs (long format, one row per donor × module × cell type):
#   Fig4B_PDvsControl_ModuleScores_CX.csv
#     donor_id, cell_type, module, diagnosis,
#     raw_score, adjusted_score (residual + grand mean),
#     age, PMI, sex
#
#   Fig4B_LM_Results_CX.csv
#     module, cell_type,
#     beta_diagnosis, SE, t, p_raw, p_adj_BH,
#     n_PD, n_Ctrl, sig
#
# ============================================================

rm(list = ls()); gc()

# =============================================================================
# SECTION 0: LIBRARIES & SETTINGS
# =============================================================================

cat("========== MS2 FIG 4B SOURCE DATA EXPORT ==========\n\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

set.seed(42)

# -------------------------------------------------------
# >>> USER SETTINGS — EDIT THESE PATHS <<<
# -------------------------------------------------------

REGION <- "CX"

# Path to Script 1 output RDS (7/8-Module AllResults)
RESULTS_RDS <- "results/CX_TE_Control/CX_8Module_AllResults.rds"

# Paths to original Seurat objects (for donor-level covariates)
ASTRO_FILE  <- "data/CX_Astrocytes_withCov.rds"
NEURON_FILE <- "data/CX_Neurons_Total_withCov.rds"

# Output directory (reuse or create new timestamped dir)
OUTPUT_BASE <- "results"

timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
csv_dir   <- file.path(OUTPUT_BASE, paste0(REGION, "_Fig4B_SourceData_", timestamp))
dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
cat("CSV output directory:", csv_dir, "\n\n")


# =============================================================================
# SECTION 1: LOAD DATA
# =============================================================================

cat("========== LOADING DATA ==========\n")

prev <- readRDS(RESULTS_RDS)
donor_scores_astro  <- prev$donor_scores[["Astrocytes"]]
donor_scores_neuron <- prev$donor_scores[["Neurons_Total"]]

cat("  Astrocyte donors:", nrow(donor_scores_astro), "\n")
cat("  Neuron donors:",    nrow(donor_scores_neuron), "\n")

astro  <- readRDS(ASTRO_FILE)
neuron <- readRDS(NEURON_FILE)
if (length(Layers(astro))  > 1) astro  <- JoinLayers(astro)
if (length(Layers(neuron)) > 1) neuron <- JoinLayers(neuron)


# =============================================================================
# SECTION 2: EXTRACT DONOR-LEVEL COVARIATES
# =============================================================================

cat("\n========== EXTRACTING COVARIATES ==========\n")

# Flexible column detection
detect_col <- function(meta, candidates, label) {
  for (nm in candidates) if (nm %in% colnames(meta)) return(nm)
  stop(paste("Could not find", label))
}
donor_possible <- c("donor_id","Donor_ID","donor","Donor","sample_id","Individual")
cond_possible  <- c("condition","Condition","disease","Disease","Status","Group","Diagnosis","diagnosis")
age_possible   <- c("age","Age","AGE","age_years")
pmi_possible   <- c("PMI","pmi","post_mortem_interval","PMI_hours")
sex_possible   <- c("sex","Sex","SEX","gender","Gender")

get_donor_covs <- function(seurat_obj, label) {
  meta <- seurat_obj@meta.data
  donor_c <- detect_col(meta, donor_possible, paste0("donor (", label, ")"))
  cond_c  <- detect_col(meta, cond_possible,  paste0("condition (", label, ")"))
  age_c   <- detect_col(meta, age_possible,   paste0("age (", label, ")"))
  pmi_c   <- detect_col(meta, pmi_possible,   paste0("PMI (", label, ")"))
  sex_c   <- detect_col(meta, sex_possible,   paste0("sex (", label, ")"))

  covs <- meta %>%
    dplyr::select(all_of(c(donor_c, cond_c, age_c, pmi_c, sex_c))) %>%
    distinct() %>%
    setNames(c("donor_id", "diagnosis", "age", "PMI", "sex"))
  covs
}

covs_astro  <- get_donor_covs(astro,  "astro")
covs_neuron <- get_donor_covs(neuron, "neuron")

cat("  Astrocyte donor covariates:", nrow(covs_astro), "\n")
cat("  Neuron    donor covariates:", nrow(covs_neuron), "\n")

# Normalise diagnosis labels → {PD, Control}
normalise_dx <- function(x) {
  x <- as.character(x)
  ifelse(grepl("PD|Parkinson", x, ignore.case = TRUE), "PD",
  ifelse(grepl("Ctrl|Control|HC|Healthy", x, ignore.case = TRUE), "Control", x))
}
covs_astro$diagnosis  <- normalise_dx(covs_astro$diagnosis)
covs_neuron$diagnosis <- normalise_dx(covs_neuron$diagnosis)


# =============================================================================
# SECTION 3: BUILD LONG-FORMAT SCORE TABLE
# =============================================================================

cat("\n========== BUILDING SCORE TABLE ==========\n")

# Module score columns (drop TE_Control duplicates if present)
score_cols <- names(donor_scores_astro)[grepl("_Score$", names(donor_scores_astro))]
cat("  Modules found:", length(score_cols), "\n")
cat("  ", paste(score_cols, collapse = ", "), "\n")

build_long <- function(scores_df, covs_df, cell_label) {
  df <- scores_df %>%
    dplyr::select(donor_id, all_of(score_cols)) %>%
    pivot_longer(cols = all_of(score_cols),
                 names_to = "module", values_to = "raw_score") %>%
    mutate(module = gsub("_Score$", "", module),
           cell_type = cell_label) %>%
    left_join(covs_df, by = "donor_id")
  df
}

long_astro  <- build_long(donor_scores_astro,  covs_astro,  "Astrocyte")
long_neuron <- build_long(donor_scores_neuron, covs_neuron, "Neuron")
long_all    <- bind_rows(long_astro, long_neuron)

# Drop rows with missing covariates
before <- nrow(long_all)
# v3 fix: PMI is missing for every GSE253462 donor; instead of dropping those donors (n = 26), mean-impute PMI
# over donors within each cell type, as in the other covariate-adjusted analyses (n = 43)
for (ct_i in unique(long_all$cell_type)) {
  i_ct <- long_all$cell_type == ct_i
  long_all$PMI[i_ct & is.na(long_all$PMI)] <- mean(unique(long_all[i_ct, c("donor_id", "PMI")])$PMI, na.rm = TRUE)
}
long_all <- long_all %>% filter(!is.na(diagnosis) & !is.na(age) & !is.na(PMI) & !is.na(sex))
cat("  Rows after covariate filtering:", nrow(long_all), "(dropped", before - nrow(long_all), ")\n")


# =============================================================================
# SECTION 4: FIT lm(score ~ diagnosis + age + PMI + sex) PER MODULE × CELL TYPE
# =============================================================================

cat("\n========== FITTING LINEAR MODELS ==========\n")

fit_one <- function(df) {
  # Factorise
  df$diagnosis <- factor(df$diagnosis, levels = c("Control", "PD"))
  df$sex       <- factor(df$sex)

  fit <- tryCatch(
    lm(raw_score ~ diagnosis + age + PMI + sex, data = df),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  co <- summary(fit)$coefficients
  dx_row <- grep("^diagnosis", rownames(co))
  if (length(dx_row) == 0) return(NULL)

  # Adjusted score = raw - (age_contrib + PMI_contrib + sex_contrib) + grand_mean
  #                = residual + fitted_without_nuisance
  # Simpler: residual + intercept + beta_dx * (diagnosis == PD)
  # which keeps between-group separation on the original scale.
  b <- coef(fit)
  nuisance_contrib <- b["age"] * df$age + b["PMI"] * df$PMI
  sex_levels <- levels(df$sex)
  if (length(sex_levels) > 1) {
    for (sl in sex_levels[-1]) {
      term <- paste0("sex", sl)
      if (term %in% names(b)) {
        nuisance_contrib <- nuisance_contrib + b[term] * (df$sex == sl)
      }
    }
  }
  adjusted <- df$raw_score - nuisance_contrib + mean(nuisance_contrib, na.rm = TRUE)

  list(
    stats = data.frame(
      beta_diagnosis = co[dx_row, "Estimate"],
      SE             = co[dx_row, "Std. Error"],
      t_value        = co[dx_row, "t value"],
      p_raw          = co[dx_row, "Pr(>|t|)"],
      n_PD           = sum(df$diagnosis == "PD"),
      n_Ctrl         = sum(df$diagnosis == "Control"),
      stringsAsFactors = FALSE
    ),
    adjusted = adjusted
  )
}

lm_results_list <- list()
long_all$adjusted_score <- NA_real_

for (ct in unique(long_all$cell_type)) {
  for (mod in unique(long_all$module)) {
    idx <- which(long_all$cell_type == ct & long_all$module == mod)
    if (length(idx) < 5) next
    sub <- long_all[idx, ]
    res <- fit_one(sub)
    if (is.null(res)) next

    long_all$adjusted_score[idx] <- res$adjusted
    lm_results_list[[paste(ct, mod, sep = "_")]] <- data.frame(
      module    = mod,
      cell_type = ct,
      res$stats,
      stringsAsFactors = FALSE
    )
  }
}

lm_results <- do.call(rbind, lm_results_list)
rownames(lm_results) <- NULL

# BH correction across all tests (16 = 8 modules × 2 cell types)
lm_results$p_adj_BH <- p.adjust(lm_results$p_raw, method = "BH")
lm_results$sig <- ifelse(lm_results$p_adj_BH < 0.001, "***",
                  ifelse(lm_results$p_adj_BH < 0.01,  "**",
                  ifelse(lm_results$p_adj_BH < 0.05,  "*", "NS")))

cat("  Models fitted:", nrow(lm_results), "\n")
cat("  Significant (p_adj < 0.05):", sum(lm_results$p_adj_BH < 0.05, na.rm = TRUE), "\n")


# =============================================================================
# SECTION 5: WRITE CSVs
# =============================================================================

cat("\n========== WRITING CSVs ==========\n")

# Reorder columns for clarity
score_out <- long_all %>%
  dplyr::select(donor_id, cell_type, module, diagnosis,
                raw_score, adjusted_score, age, PMI, sex)

lm_out <- lm_results %>%
  dplyr::select(module, cell_type, beta_diagnosis, SE, t_value,
                p_raw, p_adj_BH, sig, n_PD, n_Ctrl)

write.csv(score_out,
          file.path(csv_dir, paste0("Fig4B_PDvsControl_ModuleScores_", REGION, ".csv")),
          row.names = FALSE)
write.csv(lm_out,
          file.path(csv_dir, paste0("Fig4B_LM_Results_", REGION, ".csv")),
          row.names = FALSE)

cat("  Fig4B CSVs saved to:\n   ", csv_dir, "\n")

cat("\n============================================================\n")
cat("=== FIG 4B SOURCE DATA EXPORT COMPLETE ===\n")
cat("============================================================\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
