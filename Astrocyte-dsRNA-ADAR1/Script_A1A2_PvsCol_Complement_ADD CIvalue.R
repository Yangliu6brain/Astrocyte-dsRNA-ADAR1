# ---------------------------------------------------------------------------------------------------------------
# Repository copy. Run from the repository root. Inputs: data/ (covariate-annotated Seurat objects from
# Full_merge_RDS_script_with_covariates.R) and results/<analysis output folder> (rename each timestamped output
# folder to the name used below, e.g. CX_TE_Control_<date> -> results/CX_TE_Control). R 4.4.3.
# ---------------------------------------------------------------------------------------------------------------
# ============================================================
# A1 & A2 REACTIVE ASTROCYTE MARKERS: PD vs CONTROL
# Complement vs Non-Complement Separation
# ============================================================
#
# PURPOSE:
#   Compare A1 (neurotoxic) and A2 (neuroprotective) reactive
#   astrocyte marker expression between PD and Control donors,
#   splitting A1 markers into complement and non-complement
#   subgroups.
#
# APPROACH:
#   1. Donor-level pseudobulk mean expression in astrocytes.
#   2. Covariate-adjust (age, PMI, sex) via linear regression.
#   3. Wilcoxon rank-sum test: PD vs Control for each marker.
#   4. Composite subgroup scores (A1_Complement, A1_NonComplement,
#      A2, Pan-reactive) compared PD vs Control.
#   5. Effect sizes (Cliff's delta) and fold-changes.
#   6. Correlation of A1/A2 subgroup scores with dsRNA/immune
#      modules, stratified by condition.
#
# OUTPUTS:
#   - Excel workbook with per-gene and subgroup-level results
#   - Box/violin plots per subgroup (PD vs Control)
#   - Heatmap of marker expression by donor (annotated by condition)
#   - Scatter plots: subgroup scores vs dsRNA modules by condition
#
# ============================================================

rm(list = ls()); gc()

# =============================================================================
# SECTION 0: Libraries
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
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
})

set.seed(42)

# =============================================================================
# SECTION 1: USER SETTINGS
# =============================================================================

REGION <- "CX"   # "CX" or "SN"

COVARIATES_ENABLED <- TRUE
COV_AGE_COL <- "age"
COV_PMI_COL <- "pmi"
COV_SEX_COL <- "sex"

# --- File paths (edit to match your system) ---
if (REGION == "CX") {
  astro_file <- astro_file  <- "data/CX_Astrocytes_withCov.rds"
} else {
  astro_file <- "data/SN_Astrocytes_withCov.rds"
}

output_base <- "results"

# --- Condition labels ---
PD_LABELS      <- c("PD", "Parkinson", "Parkinson's", "Case", "Disease")
CONTROL_LABELS <- c("Control", "Ctrl", "Normal", "Healthy")

# =============================================================================
# SECTION 2: GENE & MODULE DEFINITIONS
# =============================================================================

# --- A1 neurotoxic-reactive markers (Liddelow et al. 2017) ---
# Split into complement and non-complement subgroups
A1_COMPLEMENT     <- c("C3", "C1R", "C1S", "CFB", "CFI", "SERPING1")
A1_NON_COMPLEMENT <- c("GBP2", "FKBP5", "PSMB8", "SRGN", "AMIGO2")
A1_GENES          <- c(A1_COMPLEMENT, A1_NON_COMPLEMENT)

# --- A2 neuroprotective markers (Liddelow et al. 2017) ---
A2_GENES <- c(
  "EMP1", "S100A10", "SLC10A6", "TM4SF1", "CD109",
  "PTX3", "SPHK1", "CD14", "TGFB2", "B3GNT5"
)

# --- Pan-reactive markers ---
PAN_REACTIVE_GENES <- c("GFAP", "VIM", "SERPINA3", "LCN2")

ALL_MARKER_GENES <- unique(c(A1_GENES, A2_GENES, PAN_REACTIVE_GENES))

# --- dsRNA / immune modules (for correlation with subgroup scores) ---
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
    "DDX58", "IFIH1", "DHX58", "TLR3",
    "MAVS", "TICAM1", "TBK1", "IKBKE"
  ),
  Innate_Immune = c(
    "JAK1", "TYK2", "STAT1", "STAT2", "IRF9", "RELA",
    "ISG15", "USP18", "OAS1", "OAS2", "OAS3", "OASL",
    "RNASEL", "EIF2AK2", "PRKRA",
    "MX1", "IFIT1", "IFIT2", "IFIT3",
    "C3", "C4A"
  )
)

# =============================================================================
# SECTION 3: HELPER FUNCTIONS
# =============================================================================

find_col <- function(meta, candidates) {
  for (nm in candidates) if (nm %in% colnames(meta)) return(nm)
  NULL
}

# --- Covariate residualisation ---
residualise_donor <- function(y, df) {
  cov_list <- c("n_cells", "n_count")
  if (COVARIATES_ENABLED) {
    cov_list <- c(cov_list, COV_AGE_COL, COV_PMI_COL, COV_SEX_COL)
  }
  cov_list <- cov_list[cov_list %in% colnames(df)]
  cov_list <- cov_list[sapply(cov_list, function(cv) {
    v <- df[[cv]]
    if (all(is.na(v))) return(FALSE)
    if (is.numeric(v)) sd(v, na.rm = TRUE) > 0 else length(unique(na.omit(v))) > 1
  })]
  if (length(cov_list) == 0) return(y)

  fit_df <- data.frame(y = y)
  for (cv in cov_list) {
    v <- df[[cv]]
    if (is.numeric(v)) {
      v[is.na(v)] <- mean(v, na.rm = TRUE)
      fit_df[[cv]] <- v
    } else {
      fit_df[[cv]] <- as.factor(v)
    }
  }
  fit <- tryCatch(lm(y ~ ., data = fit_df), error = function(e) NULL)
  if (is.null(fit)) return(y)
  residuals(fit) + mean(y, na.rm = TRUE)
}

# --- Cliff's delta effect size ---
cliffs_delta <- function(x, y) {
  # x = group 1, y = group 2
  # Returns value in [-1, 1]; positive = group1 > group2
  n1 <- length(x); n2 <- length(y)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  comparisons <- outer(x, y, FUN = function(a, b) sign(a - b))
  sum(comparisons) / (n1 * n2)
}

# =============================================================================
# SECTION 4: LOAD DATA & BUILD PSEUDOBULK
# =============================================================================

cat("========== LOADING ASTROCYTE DATA ==========\n")
astro <- readRDS(astro_file)
cat("  Cells:", ncol(astro), "\n")

meta <- as.data.frame(astro@meta.data)
donor_col <- find_col(meta, c("donor_id","Donor","donor","subject_id","sample_id","orig.ident"))
cond_col  <- find_col(meta, c("condition","Condition","disease","Disease","group","Group","diagnosis"))
cat("  Donor column:", donor_col, "\n")
cat("  Condition column:", cond_col, "\n\n")

# --- Expression matrix ---
expr <- GetAssayData(astro, layer = "data")

# --- Check which marker genes are present ---
marker_genes_present <- ALL_MARKER_GENES[ALL_MARKER_GENES %in% rownames(expr)]
cat("  Marker genes found:", length(marker_genes_present), "/", length(ALL_MARKER_GENES), "\n")
cat("  Missing:", paste(setdiff(ALL_MARKER_GENES, rownames(expr)), collapse = ", "), "\n\n")

# --- Module scores (raw mean, within astrocytes only — no bg-sub needed) ---
cat("========== COMPUTING MODULE SCORES ==========\n")
module_scores_cell <- sapply(names(module_gene_sets), function(mod) {
  genes <- module_gene_sets[[mod]]
  genes <- genes[genes %in% rownames(expr)]
  if (length(genes) == 0) return(rep(NA_real_, ncol(astro)))
  if (length(genes) == 1) return(as.numeric(expr[genes, ]))
  colMeans(expr[genes, , drop = FALSE])
})
colnames(module_scores_cell) <- paste0(names(module_gene_sets), "_Score")

# --- Marker expression per cell ---
marker_expr_cell <- t(as.matrix(expr[marker_genes_present, , drop = FALSE]))

# --- Aggregate to donor level ---
cat("========== AGGREGATING TO DONOR LEVEL ==========\n")
donor_ids <- as.character(meta[[donor_col]])

cell_df <- data.frame(
  donor_id  = donor_ids,
  condition = if (!is.null(cond_col)) as.character(meta[[cond_col]]) else "All",
  n_count   = if ("nCount_RNA" %in% colnames(meta)) meta$nCount_RNA else NA_real_,
  stringsAsFactors = FALSE
)
cell_df <- cbind(cell_df, module_scores_cell, marker_expr_cell)

value_cols <- c(colnames(module_scores_cell), marker_genes_present, "n_count")

pb <- cell_df %>%
  group_by(donor_id, condition) %>%
  summarise(
    n_cells = n(),
    across(all_of(value_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  as.data.frame()

# Add demographic covariates
if (COVARIATES_ENABLED) {
  cov_candidates <- c(COV_AGE_COL, COV_PMI_COL, COV_SEX_COL)
  cov_cols <- cov_candidates[cov_candidates %in% colnames(meta)]
  if (length(cov_cols) > 0) {
    cov_df <- data.frame(
      donor_id = donor_ids,
      meta[, cov_cols, drop = FALSE],
      stringsAsFactors = FALSE
    ) %>%
      group_by(donor_id) %>%
      summarise(across(everything(), ~ first(na.omit(.))), .groups = "drop") %>%
      as.data.frame()
    pb <- merge(pb, cov_df, by = "donor_id", all.x = TRUE)
  }
}

# --- Standardise condition labels ---
pb$group <- case_when(
  toupper(pb$condition) %in% toupper(PD_LABELS)      ~ "PD",
  toupper(pb$condition) %in% toupper(CONTROL_LABELS)  ~ "Control",
  TRUE ~ pb$condition
)

n_PD   <- sum(pb$group == "PD")
n_Ctrl <- sum(pb$group == "Control")
cat("  Donors: PD =", n_PD, "| Control =", n_Ctrl, "| Total =", nrow(pb), "\n\n")

# =============================================================================
# SECTION 5: COMPUTE SUBGROUP COMPOSITE SCORES
# =============================================================================

cat("========== COMPUTING SUBGROUP COMPOSITE SCORES ==========\n")

# Subgroup definitions (only genes present in data)
subgroup_defs <- list(
  A1_Complement     = A1_COMPLEMENT[A1_COMPLEMENT %in% marker_genes_present],
  A1_NonComplement  = A1_NON_COMPLEMENT[A1_NON_COMPLEMENT %in% marker_genes_present],
  A1_All            = A1_GENES[A1_GENES %in% marker_genes_present],
  A2                = A2_GENES[A2_GENES %in% marker_genes_present],
  Pan_Reactive      = PAN_REACTIVE_GENES[PAN_REACTIVE_GENES %in% marker_genes_present]
)

for (sg_name in names(subgroup_defs)) {
  genes <- subgroup_defs[[sg_name]]
  if (length(genes) == 0) {
    pb[[paste0(sg_name, "_Score")]] <- NA_real_
  } else if (length(genes) == 1) {
    pb[[paste0(sg_name, "_Score")]] <- pb[[genes]]
  } else {
    pb[[paste0(sg_name, "_Score")]] <- rowMeans(pb[, genes, drop = FALSE], na.rm = TRUE)
  }
  cat("  ", sg_name, ": ", length(genes), " genes →",
      paste(genes, collapse = ", "), "\n")
}
cat("\n")

# =============================================================================
# SECTION 6: COVARIATE RESIDUALISATION
# =============================================================================

cat("========== RESIDUALISING FOR COVARIATES ==========\n")

# Residualise individual marker genes
for (gene in marker_genes_present) {
  pb[[paste0(gene, "_resid")]] <- residualise_donor(pb[[gene]], pb)
}

# Residualise subgroup scores
for (sg_name in names(subgroup_defs)) {
  score_col <- paste0(sg_name, "_Score")
  pb[[paste0(score_col, "_resid")]] <- residualise_donor(pb[[score_col]], pb)
}

# Residualise module scores
for (mod in names(module_gene_sets)) {
  score_col <- paste0(mod, "_Score")
  pb[[paste0(score_col, "_resid")]] <- residualise_donor(pb[[score_col]], pb)
}

cat("  Done.\n\n")

# =============================================================================
# SECTION 7: PD vs CONTROL — PER-GENE TESTS
# =============================================================================

cat("========== PD vs CONTROL: PER-GENE TESTS ==========\n")

gene_results <- list()

for (gene in marker_genes_present) {
  resid_col <- paste0(gene, "_resid")
  vals_PD   <- pb[[resid_col]][pb$group == "PD"]
  vals_Ctrl <- pb[[resid_col]][pb$group == "Control"]

  # Remove NAs
  vals_PD   <- vals_PD[!is.na(vals_PD)]
  vals_Ctrl <- vals_Ctrl[!is.na(vals_Ctrl)]

  if (length(vals_PD) < 3 || length(vals_Ctrl) < 3) next

  wt <- wilcox.test(vals_PD, vals_Ctrl, exact = FALSE)

  # Classify marker
  marker_class <- case_when(
    gene %in% A1_COMPLEMENT     ~ "A1_Complement",
    gene %in% A1_NON_COMPLEMENT ~ "A1_NonComplement",
    gene %in% A2_GENES          ~ "A2",
    gene %in% PAN_REACTIVE_GENES ~ "Pan_Reactive",
    TRUE                        ~ "Other"
  )

  gene_results[[length(gene_results) + 1]] <- data.frame(
    Gene          = gene,
    Marker_Class  = marker_class,
    mean_PD       = mean(vals_PD),
    mean_Control  = mean(vals_Ctrl),
    log2FC        = log2((mean(vals_PD) + 1e-6) / (mean(vals_Ctrl) + 1e-6)),
    median_PD     = median(vals_PD),
    median_Control = median(vals_Ctrl),
    Cliffs_delta  = cliffs_delta(vals_PD, vals_Ctrl),
    W_stat        = wt$statistic,
    p_value       = wt$p.value,
    n_PD          = length(vals_PD),
    n_Control     = length(vals_Ctrl),
    stringsAsFactors = FALSE
  )
}

gene_df <- bind_rows(gene_results)
gene_df$p_adj <- p.adjust(gene_df$p_value, method = "BH")
gene_df <- gene_df %>% arrange(p_adj, p_value)

cat("  Genes tested:", nrow(gene_df), "\n")
cat("  FDR < 0.05:", sum(gene_df$p_adj < 0.05), "\n")
cat("  FDR < 0.10:", sum(gene_df$p_adj < 0.10), "\n\n")

# =============================================================================
# SECTION 8: PD vs CONTROL — SUBGROUP COMPOSITE SCORES
# =============================================================================

cat("========== PD vs CONTROL: SUBGROUP COMPOSITE SCORES ==========\n")

subgroup_results <- list()

for (sg_name in names(subgroup_defs)) {
  resid_col <- paste0(sg_name, "_Score_resid")
  vals_PD   <- pb[[resid_col]][pb$group == "PD"]
  vals_Ctrl <- pb[[resid_col]][pb$group == "Control"]

  vals_PD   <- vals_PD[!is.na(vals_PD)]
  vals_Ctrl <- vals_Ctrl[!is.na(vals_Ctrl)]

  if (length(vals_PD) < 3 || length(vals_Ctrl) < 3) next

  wt <- wilcox.test(vals_PD, vals_Ctrl, exact = FALSE)

  subgroup_results[[length(subgroup_results) + 1]] <- data.frame(
    Subgroup       = sg_name,
    n_genes        = length(subgroup_defs[[sg_name]]),
    mean_PD        = mean(vals_PD),
    mean_Control   = mean(vals_Ctrl),
    log2FC         = log2((mean(vals_PD) + 1e-6) / (mean(vals_Ctrl) + 1e-6)),
    median_PD      = median(vals_PD),
    median_Control = median(vals_Ctrl),
    Cliffs_delta   = cliffs_delta(vals_PD, vals_Ctrl),
    W_stat         = wt$statistic,
    p_value        = wt$p.value,
    n_PD           = length(vals_PD),
    n_Control      = length(vals_Ctrl),
    stringsAsFactors = FALSE
  )
}

subgroup_df <- bind_rows(subgroup_results)
subgroup_df$p_adj <- p.adjust(subgroup_df$p_value, method = "BH")
subgroup_df <- subgroup_df %>% arrange(p_adj, p_value)

cat("  Subgroups tested:", nrow(subgroup_df), "\n")
print(subgroup_df %>% select(Subgroup, n_genes, log2FC, Cliffs_delta, p_value, p_adj) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))), row.names = FALSE)
cat("\n")

# =============================================================================
# SECTION 9: SUBGROUP SCORES × MODULE CORRELATIONS BY CONDITION
# =============================================================================

cat("========== SUBGROUP × MODULE CORRELATIONS BY CONDITION ==========\n")

corr_results <- list()

for (sg_name in names(subgroup_defs)) {
  sg_resid <- paste0(sg_name, "_Score_resid")
  for (mod_name in names(module_gene_sets)) {
    mod_resid <- paste0(mod_name, "_Score_resid")

    for (grp in c("ALL", "PD", "Control")) {
      if (grp == "ALL") {
        idx <- seq_len(nrow(pb))
      } else {
        idx <- which(pb$group == grp)
      }
      if (length(idx) < 5) next

      x <- pb[[sg_resid]][idx]
      y <- pb[[mod_resid]][idx]
      cc <- complete.cases(x, y)
      n <- sum(cc)
      if (n < 5) next

      ct <- cor.test(x[cc], y[cc], method = "spearman", exact = FALSE)

      corr_results[[length(corr_results) + 1]] <- data.frame(
        Subgroup  = sg_name,
        Module    = mod_name,
        Condition = grp,
        rho       = ct$estimate,
        p_value   = ct$p.value,
        n_donors  = n,
        stringsAsFactors = FALSE
      )
    }
  }
}

corr_df <- bind_rows(corr_results)
corr_df$p_adj <- p.adjust(corr_df$p_value, method = "BH")
corr_df <- corr_df %>% arrange(Condition, p_adj, p_value)

cat("  Total correlation tests:", nrow(corr_df), "\n")
cat("  FDR < 0.05:", sum(corr_df$p_adj < 0.05), "\n\n")


# =============================================================================
# ADD-ON to Script_A1_A2_PD_vs_Control_Complement_Split.R
# Bootstrap CIs for the A1-complement "collapse" (your strongest decoupling
# evidence, Figure 7c).
#
# PLACEMENT: paste AFTER Section 9 of that script, i.e. after `pb` exists with
# the *_Score_resid columns and `corr_df` has been built. It re-uses `pb`,
# `subgroup_defs`, `module_gene_sets`, and reproduces the script's OWN flow:
#   residualise on the FULL sample (already done -> the *_resid columns),
#   then split by group and Spearman. The bootstrap therefore resamples donors
#   and recomputes residualisation the same way, so point estimates match Fig7.
#
# WHY SEPARATE FROM THE Script4 add-on:
#   - A1 scores are raw colMeans (deterministic; no bg-sub, no seed issue).
#   - A1 residualisation adds mean back and includes n_cells + n_count as extra
#     covariates. This differs from Script4, so it needs its own engine to
#     reproduce the 0.78 collapse numbers faithfully.
#
# FOCUS PAIRS (Fig7c, A1_Complement arm, control -> PD):
#   the four couplings you report collapsing: EV transport, innate immune,
#   ADAR1, dsRNA sensing. Edit FOCUS_MODULES to match your module_gene_sets keys.
# =============================================================================

stopifnot(exists("pb"), exists("subgroup_defs"), exists("module_gene_sets"),
          exists("residualise_donor"))

N_BOOT   <- 10000
CI_LEVEL <- 0.95
set.seed(42)   # reseeded per analysis below

SUBGROUP_FOCUS <- "A1_Complement"
# module keys as they appear in names(module_gene_sets):
FOCUS_MODULES  <- c("dsRNA_EV_Transport", "Innate_Immune", "ADAR1", "dsRNA_Sensing")
FOCUS_MODULES  <- intersect(FOCUS_MODULES, names(module_gene_sets))

perc_ci <- function(v, level = CI_LEVEL) {
  a <- (1 - level)/2
  q <- quantile(v, c(a, 1 - a), na.rm = TRUE, names = FALSE)
  list(lo = q[1], hi = q[2], frac_na = mean(is.na(v)))
}
boot_p_two_sided <- function(v){ v<-v[!is.na(v)]; if(!length(v)) return(NA_real_)
min(1, 2*min(mean(v<=0), mean(v>=0))) }

# Reproduce the script's residualise-on-full-sample-then-correlate, on an
# arbitrary donor-row subset `d` (a bootstrap resample of pb rows).
# sg_raw / mod_raw are the RAW score column names (pre-resid).
rho_group_from_resample <- function(d, sg_raw, mod_raw, grp) {
  # residualise on the FULL resample d (matches original flow)
  sg_r  <- residualise_donor(d[[sg_raw]],  d)
  mod_r <- residualise_donor(d[[mod_raw]], d)
  keep  <- if (grp == "ALL") rep(TRUE, nrow(d)) else d$group == grp
  x <- sg_r[keep]; y <- mod_r[keep]
  cc <- is.finite(x) & is.finite(y)
  if (sum(cc) < 5) return(NA_real_)
  suppressWarnings(cor(x[cc], y[cc], method = "spearman"))
}

boot_A1 <- function(sg_name, mod_name, n_boot = N_BOOT, seed = 42) {
  sg_raw  <- paste0(sg_name,  "_Score")
  mod_raw <- paste0(mod_name, "_Score")
  stopifnot(sg_raw %in% colnames(pb), mod_raw %in% colnames(pb))
  
  # point estimates (reproduce Fig7 numbers on the real data)
  rC <- rho_group_from_resample(pb, sg_raw, mod_raw, "Control")
  rP <- rho_group_from_resample(pb, sg_raw, mod_raw, "PD")
  point_delta <- rC - rP   # positive = coupling higher in control (decoupling in PD)
  
  # stratified donor resample within each group
  idxC <- which(pb$group == "Control"); idxP <- which(pb$group == "PD")
  set.seed(seed)
  dvec <- numeric(n_boot); rCv <- numeric(n_boot); rPv <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    sC <- sample(idxC, length(idxC), replace = TRUE)
    sP <- sample(idxP, length(idxP), replace = TRUE)
    d  <- pb[c(sC, sP), , drop = FALSE]      # resampled full frame
    rc <- rho_group_from_resample(d, sg_raw, mod_raw, "Control")
    rp <- rho_group_from_resample(d, sg_raw, mod_raw, "PD")
    rCv[b] <- rc; rPv[b] <- rp; dvec[b] <- rc - rp
  }
  ciC <- perc_ci(rCv); ciP <- perc_ci(rPv); ciD <- perc_ci(dvec)
  data.frame(
    subgroup = sg_name, module = mod_name,
    rho_control = rC, ctrl_lo = ciC$lo, ctrl_hi = ciC$hi,
    rho_PD = rP,      pd_lo = ciP$lo,   pd_hi = ciP$hi,
    delta = point_delta, d_lo = ciD$lo, d_hi = ciD$hi,
    delta_excludes_zero = (ciD$lo > 0) | (ciD$hi < 0),
    boot_p = boot_p_two_sided(dvec),
    stringsAsFactors = FALSE
  )
}

cat("\n========== A1-COMPLEMENT BOOTSTRAP (Fig7c) ==========\n")
cat("Subgroup:", SUBGROUP_FOCUS, "| modules:", paste(FOCUS_MODULES, collapse=", "), "\n\n")

res_list <- lapply(seq_along(FOCUS_MODULES), function(i)
  boot_A1(SUBGROUP_FOCUS, FOCUS_MODULES[i], seed = 400 + i))
a1_boot <- do.call(rbind, res_list)

# BH across this A1-complement family (control-vs-PD deltas)
a1_boot$p_BH <- p.adjust(a1_boot$boot_p, method = "BH")

cat("--- A1_Complement: control vs PD coupling, with 95% CIs ---\n")
print(a1_boot[, c("module","rho_control","ctrl_lo","ctrl_hi",
                  "rho_PD","pd_lo","pd_hi",
                  "delta","d_lo","d_hi","delta_excludes_zero","boot_p","p_BH")],
      row.names = FALSE)

write.csv(a1_boot, file.path(getwd(), "boot_A1complement_covadj.csv"), row.names = FALSE)
cat("\nWrote boot_A1complement_covadj.csv\n")
# =============================================================================
# SECTION 10: OUTPUT DIRECTORY
# =============================================================================

timestamp   <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- file.path(output_base,
                         paste0(REGION, "_A1A2_PD_vs_Control_", timestamp))
dir.create(file.path(results_dir, "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SECTION 11: WRITE TABLES
# =============================================================================

cat("========== WRITING RESULTS ==========\n")

wb <- createWorkbook()

# Sheet 1: Per-gene PD vs Control
addWorksheet(wb, "PerGene_PD_vs_Control")
writeData(wb, "PerGene_PD_vs_Control", gene_df)

# Sheet 2: Subgroup composite PD vs Control
addWorksheet(wb, "Subgroup_PD_vs_Control")
writeData(wb, "Subgroup_PD_vs_Control", subgroup_df)

# Sheet 3: Subgroup × Module correlations (ALL)
addWorksheet(wb, "Corr_ALL")
writeData(wb, "Corr_ALL", corr_df %>% filter(Condition == "ALL"))

# Sheet 4: Subgroup × Module correlations (PD)
addWorksheet(wb, "Corr_PD")
writeData(wb, "Corr_PD", corr_df %>% filter(Condition == "PD"))

# Sheet 5: Subgroup × Module correlations (Control)
addWorksheet(wb, "Corr_Control")
writeData(wb, "Corr_Control", corr_df %>% filter(Condition == "Control"))

# Sheet 6: Correlation pivot (rho) — ALL
corr_pivot_ALL <- corr_df %>%
  filter(Condition == "ALL") %>%
  select(Subgroup, Module, rho) %>%
  pivot_wider(names_from = Module, values_from = rho)
addWorksheet(wb, "Corr_rho_pivot_ALL")
writeData(wb, "Corr_rho_pivot_ALL", corr_pivot_ALL)

# Sheet 7: Donor-level data for reproducibility
donor_export <- pb %>%
  select(donor_id, group, n_cells,
         all_of(paste0(names(subgroup_defs), "_Score_resid")),
         all_of(paste0(marker_genes_present, "_resid")))
addWorksheet(wb, "Donor_Level_Data")
writeData(wb, "Donor_Level_Data", donor_export)

saveWorkbook(wb, file.path(results_dir, "tables",
             paste0(REGION, "_A1A2_PD_vs_Control.xlsx")),
             overwrite = TRUE)

write.csv(gene_df,
          file.path(results_dir, "tables",
                    paste0(REGION, "_PerGene_PD_vs_Control.csv")),
          row.names = FALSE)

cat("  Tables saved.\n\n")

# =============================================================================
# SECTION 12: FIGURES
# =============================================================================

cat("========== GENERATING FIGURES ==========\n")

# Colour palette
COLOUR_PD   <- "#E64B35"
COLOUR_CTRL <- "#4DBBD5"

# ----- 12A: Box + jitter plots — subgroup composite scores PD vs Control -----

plot_list_box <- list()
for (sg_name in names(subgroup_defs)) {
  resid_col <- paste0(sg_name, "_Score_resid")
  plot_df <- pb %>%
    filter(group %in% c("PD", "Control")) %>%
    select(donor_id, group, value = all_of(resid_col)) %>%
    filter(!is.na(value))

  sg_stats <- subgroup_df %>% filter(Subgroup == sg_name)
  p_label <- if (nrow(sg_stats) > 0) {
    sprintf("p_adj = %.3g\nΔ = %.2f", sg_stats$p_adj[1], sg_stats$Cliffs_delta[1])
  } else { "" }

  p <- ggplot(plot_df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.7,
                aes(colour = group)) +
    scale_fill_manual(values = c("PD" = COLOUR_PD, "Control" = COLOUR_CTRL)) +
    scale_colour_manual(values = c("PD" = COLOUR_PD, "Control" = COLOUR_CTRL)) +
    labs(
      title = gsub("_", " ", sg_name),
      subtitle = p_label,
      y = "Score (residualised)", x = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")

  plot_list_box[[sg_name]] <- p
}

p_boxes <- wrap_plots(plot_list_box, nrow = 1)
ggsave(file.path(results_dir, "figures",
                 paste0(REGION, "_A1A2_Subgroup_BoxPlots.pdf")),
       p_boxes, width = 3.5 * length(plot_list_box), height = 4.5)

# ----- 12B: Per-gene box plots split by Marker_Class -----

for (mc in c("A1_Complement", "A1_NonComplement", "A2", "Pan_Reactive")) {
  genes_mc <- gene_df$Gene[gene_df$Marker_Class == mc]
  if (length(genes_mc) == 0) next

  plot_list_gene <- list()
  for (gene in genes_mc) {
    resid_col <- paste0(gene, "_resid")
    plot_df <- pb %>%
      filter(group %in% c("PD", "Control")) %>%
      select(donor_id, group, value = all_of(resid_col)) %>%
      filter(!is.na(value))

    g_stats <- gene_df %>% filter(Gene == gene)

    p <- ggplot(plot_df, aes(x = group, y = value, fill = group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.7,
                  aes(colour = group)) +
      scale_fill_manual(values = c("PD" = COLOUR_PD, "Control" = COLOUR_CTRL)) +
      scale_colour_manual(values = c("PD" = COLOUR_PD, "Control" = COLOUR_CTRL)) +
      labs(
        title = gene,
        subtitle = sprintf("p_adj=%.3g  Δ=%.2f",
                           g_stats$p_adj[1], g_stats$Cliffs_delta[1]),
        y = "Expression (resid.)", x = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(legend.position = "none")

    plot_list_gene[[gene]] <- p
  }

  p_gene <- wrap_plots(plot_list_gene, nrow = 1)
  ggsave(file.path(results_dir, "figures",
                   paste0(REGION, "_", mc, "_PerGene_BoxPlots.pdf")),
         p_gene,
         width = 3 * length(plot_list_gene),
         height = 4)
}

# ----- 12C: Heatmap — donor × marker genes, annotated by condition -----

heat_genes <- marker_genes_present[marker_genes_present %in%
                c(A1_COMPLEMENT, A1_NON_COMPLEMENT, A2_GENES)]
heat_cols   <- paste0(heat_genes, "_resid")
heat_cols   <- heat_cols[heat_cols %in% colnames(pb)]

if (length(heat_cols) > 2) {
  heat_mat <- pb %>%
    filter(group %in% c("PD", "Control")) %>%
    select(donor_id, group, all_of(heat_cols)) %>%
    arrange(group, donor_id) %>%
    column_to_rownames("donor_id")

  annot_row <- data.frame(
    Condition = heat_mat$group,
    row.names = rownames(heat_mat)
  )
  heat_mat$group <- NULL

  # Clean column names
  colnames(heat_mat) <- gsub("_resid$", "", colnames(heat_mat))

  # Gene annotation: subgroup
  gene_class <- case_when(
    colnames(heat_mat) %in% A1_COMPLEMENT     ~ "A1_Complement",
    colnames(heat_mat) %in% A1_NON_COMPLEMENT ~ "A1_NonComplement",
    colnames(heat_mat) %in% A2_GENES          ~ "A2",
    TRUE ~ "Other"
  )
  annot_col <- data.frame(
    Subgroup = gene_class,
    row.names = colnames(heat_mat)
  )

  ann_colors <- list(
    Condition = c("PD" = COLOUR_PD, "Control" = COLOUR_CTRL),
    Subgroup  = c("A1_Complement" = "#D62728", "A1_NonComplement" = "#FF7F0E",
                  "A2" = "#2CA02C", "Other" = "grey60")
  )

  p_heat <- pheatmap(
    as.matrix(heat_mat),
    scale = "column",
    annotation_row = annot_row,
    annotation_col = annot_col,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows = TRUE, cluster_cols = TRUE,
    show_rownames = FALSE,
    main = paste0(REGION, " Astrocytes: A1/A2 Marker Expression by Donor"),
    fontsize = 10,
    angle_col = 45,
    silent = TRUE
  )

  ggsave(file.path(results_dir, "figures",
                   paste0(REGION, "_A1A2_Donor_Heatmap.pdf")),
         p_heat, width = 10, height = 8)
}

# ----- 12D: Scatter — Subgroup scores vs dsRNA modules by condition -----

FOCUS_SUBGROUPS <- c("A1_Complement", "A1_NonComplement", "A2")
FOCUS_MODULES   <- c("dsRNA_Sensing", "Innate_Immune", "dsRNA_EV_Transport")

scatter_list <- list()
for (sg in FOCUS_SUBGROUPS) {
  for (mod in FOCUS_MODULES) {
    sg_col  <- paste0(sg, "_Score_resid")
    mod_col <- paste0(mod, "_Score_resid")

    plot_df <- pb %>%
      filter(group %in% c("PD", "Control")) %>%
      select(group, sg_val = all_of(sg_col), mod_val = all_of(mod_col)) %>%
      filter(complete.cases(.))

    if (nrow(plot_df) < 5) next

    # Get per-condition correlations
    corr_PD   <- corr_df %>% filter(Subgroup == sg, Module == mod, Condition == "PD")
    corr_Ctrl <- corr_df %>% filter(Subgroup == sg, Module == mod, Condition == "Control")

    sub_text <- paste0(
      "PD: ρ=", ifelse(nrow(corr_PD) > 0, round(corr_PD$rho[1], 2), "NA"),
      "  Ctrl: ρ=", ifelse(nrow(corr_Ctrl) > 0, round(corr_Ctrl$rho[1], 2), "NA")
    )

    p <- ggplot(plot_df, aes(x = mod_val, y = sg_val,
                              colour = group, shape = group)) +
      geom_point(size = 2.5, alpha = 0.75) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, alpha = 0.15) +
      scale_colour_manual(values = c("PD" = COLOUR_PD, "Control" = COLOUR_CTRL)) +
      scale_shape_manual(values = c("PD" = 16, "Control" = 17)) +
      labs(
        x = paste0(gsub("_", " ", mod), " (resid.)"),
        y = paste0(gsub("_", " ", sg), " (resid.)"),
        title = paste0(gsub("_", " ", sg), " vs ", gsub("_", " ", mod)),
        subtitle = sub_text
      ) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom", legend.title = element_blank())

    scatter_list[[paste0(sg, "_", mod)]] <- p
  }
}

if (length(scatter_list) > 0) {
  p_scatter <- wrap_plots(scatter_list, ncol = length(FOCUS_MODULES))
  ggsave(file.path(results_dir, "figures",
                   paste0(REGION, "_A1A2_vs_Modules_ByCondition_Scatters.pdf")),
         p_scatter,
         width = 5 * min(length(FOCUS_MODULES), 3),
         height = 5 * length(FOCUS_SUBGROUPS))
}

# ----- 12E: Grouped bar — log2FC by subgroup (overview) -----

bar_df <- gene_df %>%
  mutate(Marker_Class = factor(Marker_Class,
           levels = c("A1_Complement", "A1_NonComplement", "A2", "Pan_Reactive"))) %>%
  filter(!is.na(Marker_Class))

p_fc <- ggplot(bar_df, aes(x = reorder(Gene, log2FC), y = log2FC, fill = Marker_Class)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  coord_flip() +
  scale_fill_manual(values = c(
    "A1_Complement"    = "#D62728",
    "A1_NonComplement" = "#FF7F0E",
    "A2"               = "#2CA02C",
    "Pan_Reactive"     = "#9467BD"
  )) +
  labs(
    title = paste0(REGION, " Astrocytes: log2FC (PD / Control)"),
    subtitle = "Covariate-residualised donor-level means",
    x = NULL, y = "log2 Fold Change (PD / Control)", fill = "Marker Class"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

ggsave(file.path(results_dir, "figures",
                 paste0(REGION, "_A1A2_log2FC_barplot.pdf")),
       p_fc, width = 8, height = 6)

cat("  Figures saved.\n\n")

# =============================================================================
# SECTION 13: CONSOLE SUMMARY
# =============================================================================

cat("============================================================\n")
cat("=== RESULTS SUMMARY ===\n")
cat("============================================================\n\n")

cat("Region:", REGION, "\n")
cat("Donors: PD =", n_PD, "| Control =", n_Ctrl, "\n")
cat("Marker genes tested:", length(marker_genes_present), "\n\n")

cat("--- Subgroup Composite Scores: PD vs Control ---\n")
print(
  subgroup_df %>%
    select(Subgroup, n_genes, log2FC, Cliffs_delta, p_value, p_adj) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))),
  row.names = FALSE
)

cat("\n--- Per-Gene Results (top 15 by p_adj) ---\n")
print(
  gene_df %>%
    head(15) %>%
    select(Gene, Marker_Class, log2FC, Cliffs_delta, p_value, p_adj) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))),
  row.names = FALSE
)

cat("\n--- A1_Complement genes ---\n")
print(
  gene_df %>%
    filter(Marker_Class == "A1_Complement") %>%
    select(Gene, log2FC, Cliffs_delta, p_value, p_adj) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))),
  row.names = FALSE
)

cat("\n--- A1_NonComplement genes ---\n")
print(
  gene_df %>%
    filter(Marker_Class == "A1_NonComplement") %>%
    select(Gene, log2FC, Cliffs_delta, p_value, p_adj) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))),
  row.names = FALSE
)

cat("\n--- A2 genes ---\n")
print(
  gene_df %>%
    filter(Marker_Class == "A2") %>%
    select(Gene, log2FC, Cliffs_delta, p_value, p_adj) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))),
  row.names = FALSE
)

cat("\n--- Subgroup × Module Correlations: FDR < 0.10 (ALL donors) ---\n")
sig_corr <- corr_df %>%
  filter(Condition == "ALL", p_adj < 0.10) %>%
  select(Subgroup, Module, rho, p_adj, n_donors) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
if (nrow(sig_corr) > 0) print(sig_corr, row.names = FALSE) else cat("  None.\n")

cat("\n--- Condition-Stratified Highlights (PD vs Control) ---\n")
for (grp in c("PD", "Control")) {
  cat(sprintf("\n  %s:\n", grp))
  grp_hits <- corr_df %>%
    filter(Condition == grp, p_adj < 0.20,
           Subgroup %in% c("A1_Complement", "A1_NonComplement", "A2")) %>%
    select(Subgroup, Module, rho, p_adj) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
    arrange(p_adj)
  if (nrow(grp_hits) > 0) print(grp_hits, row.names = FALSE) else cat("    No hits at FDR < 0.20\n")
}

cat("\n\nDone. Output:", results_dir, "\n")
