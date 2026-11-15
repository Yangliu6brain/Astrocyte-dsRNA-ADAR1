# ---------------------------------------------------------------------------------------------------------------
# Repository copy. Run from the repository root. Inputs: data/ (covariate-annotated Seurat objects from
# Full_merge_RDS_script_with_covariates.R) and results/<analysis output folder> (rename each timestamped output
# folder to the name used below, e.g. CX_TE_Control_<date> -> results/CX_TE_Control). R 4.4.3.
# ---------------------------------------------------------------------------------------------------------------
# =============================================================================
# MS2 — FIGURE 9B + SUPPLEMENTARY TABLES SX/SY
# Cortex Braak stage analysis: per-cell-type AND astrocyte-neuron divergence
# =============================================================================
#
# Tests whether dsRNA / ADAR module activity changes with alpha-synuclein
# Braak pathology in postmortem cortex, at TWO complementary levels:
#
#   LAYER 1 (per-cell-type)  : does astrocyte / neuron module activity itself
#                              shift with Braak stage?
#   LAYER 2 (divergence)     : does the astrocyte-neuron gap shift with stage?
#
# Both null results together support a CONSTITUTIVE (rather than
# stage-escalating) framing of the cross-sectional astrocyte-specific
# dsRNA / innate-immune coupling identified in Figure 4.
#
# DATASETS:
#   GSE178146 (n=14 donors: 7 control, 7 PD Braak 5-6)
#   GSE253462 (n=17 donors: 8 control, 9 PD Braak 4-6)
#   Total available with Braak metadata: 31 donors
#   GSE202210 excluded (no Braak metadata)
#
# DESIGN:
#   Stage groups: Control (Braak 0, n=15)
#                 Braak 4-5  (pooled, n=5)
#                 Braak 6   (n=11)
#   8 modules   : ADAR1, ADAR2, SNCA, dsRNA_Processing, dsRNA_EV_Transport,
#                 dsRNA_Sensing, Innate_Immune, TE_Control
#                 (matches the 8-module framework used elsewhere in MS2)
#   PRIMARY    : partial Spearman | age + dataset (ppcor::pcor.test)
#   SECONDARY  : Jonckheere-Terpstra trend test (clinfun::jonckheere.test)
#   EXPLORATORY: Kruskal-Wallis + pairwise Wilcoxon
#   CI         : 10,000-resample percentile bootstrap on partial rho
#   FDR        : BH within each test family within each contrast
#
# CONTRASTS:
#   Per-cell-type : Astrocyte (separately), Total Neuron (separately),
#                   Excitatory Neuron, Inhibitory Neuron
#   Divergence    : Astro - Total Neuron (primary)
#                   Astro - Excitatory Neuron (sensitivity)
#                   Astro - Inhibitory Neuron (sensitivity)
#
# OUTPUT:
#   <OUT_DIR>/MS2_Figure9B_Braak_<timestamp>/
#     1_DonorScores/       per-donor wide tables, one per cell type
#     2_DivergenceScores/  per-donor divergence tables (3 contrasts)
#     3_Statistics/        stats CSVs per contrast
#     4_Figures/           Figure 9B + supplementary panels
#     Supplementary_Table_SX_PerCellType.xlsx
#     Supplementary_Table_SY_Divergence.xlsx
#     session_info.txt
# =============================================================================

rm(list = ls()); gc()

# =============================================================================
# 0. LIBRARIES
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
  library(ggbeeswarm)
})

for (pkg in c("ppcor", "clinfun")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, ask = FALSE)
}
library(ppcor)
library(clinfun)

set.seed(42)
N_BOOT <- 10000  # bootstrap resamples for partial-rho CI

# =============================================================================
# 1. USER SETTINGS - EDIT THESE PATHS BEFORE RUNNING
# =============================================================================
ASTRO_RDS  <- "data/CX_Astrocytes_withCov.rds"
NEURON_RDS <- "data/CX_Neurons_Total_withCov.rds"
OUT_DIR    <- "results"

DONOR_COL          <- "sample_id"
DATASET_COL        <- "dataset"
NEURON_SUBTYPE_COL <- "cell_type"

# =============================================================================
# 2. BRAAK STAGE METADATA
# =============================================================================
#
# Only Braak alpha-synuclein stage is hardcoded here. Donor age (and any
# other covariates such as PMI, sex) is ALWAYS taken from the Seurat object's
# own @meta.data — never from this hardcoded table. An earlier version of
# this script also hardcoded age, and the age vector for GSE178146 was
# mis-paired with sample_id, producing age errors of up to 27 years for some
# donors. The withCov RDS files are the canonical source for donor age, so
# we now defer to them and only carry Braak stage here.
#
# If your Seurat metadata does not contain age, add it before running this
# script (e.g. age <- match(sample_id, donor_metadata$sample_id)).
#
META_178146 <- data.frame(
  sample_id = c("C36","C48","PDC05","PDC22","PDC34","PDC87","PDC91",
                "PD732","PD747","PD413","PD416","PD523","PD666","PD683"),
  dataset = "GSE178146",
  braak_asyn_stage = c(0,0,0,0,0,0,0, 5,6,6,6,6,6,5),
  stringsAsFactors = FALSE
)

META_253462 <- data.frame(
  sample_id = c("T1619","T1952","T2017","T2465","T343","T4174","T4210","T551",
                "T1976","T233","T257","T290","T308","T362","T3860","T4162","T777"),
  dataset = "GSE253462",
  braak_asyn_stage = c(0,0,0,0,0,0,0,0, 6,5,6,5,6,4,6,6,6),
  stringsAsFactors = FALSE
)

STAGE_META <- rbind(META_178146, META_253462)
STAGE_META$stage_group <- factor(
  case_when(
    STAGE_META$braak_asyn_stage == 0          ~ "Control",
    STAGE_META$braak_asyn_stage %in% c(4, 5)  ~ "Braak_4-5",
    STAGE_META$braak_asyn_stage == 6          ~ "Braak_6"
  ),
  levels = c("Control","Braak_4-5","Braak_6")
)
STAGE_META$stage_numeric <- as.integer(STAGE_META$stage_group) - 1L

cat("Stage distribution by dataset:\n")
print(table(STAGE_META$stage_group, STAGE_META$dataset)); cat("\n")

# =============================================================================
# 3. MODULES AND KEY GENES
# =============================================================================
#
# 8-module framework matches the rest of MS2 (Script 5: TE_Control consolidated):
#   - 3 single-gene anchor modules : ADAR1, ADAR2, SNCA
#   - 5 multi-gene pathway modules : dsRNA_Processing, dsRNA_EV_Transport,
#                                    dsRNA_Sensing, Innate_Immune, TE_Control
#
# Note vs earlier versions:
#   - ADAR3 (ADARB2) is NOT in this framework.
#   - dsRNA_Processing replaces the older "dsRNA_Editing" label and excludes
#     DROSHA / DGCR8 (which belong to miRNA biogenesis rather than dsRNA
#     processing per se).
#   - dsRNA_Sensing is narrowed to sensors + adapters + kinases. Downstream
#     IFN signalling and ISGs are now their own module (Innate_Immune).
#
MODULE_GENES <- list(
  ADAR1              = "ADAR",
  ADAR2              = "ADARB1",
  SNCA               = "SNCA",
  dsRNA_Processing   = c("ILF2","ILF3","STRBP","ZFR","DHX9",
                         "DDX39B","SFPQ","NONO"),
  dsRNA_EV_Transport = c("SIDT1","SIDT2","LAMP2","RAB27A","RAB27B","RAB35",
                         "SMPD3","CD63","CD9","CD81","PDCD6IP","TSG101",
                         "SDCBP","MVP"),
  dsRNA_Sensing      = c("DDX58","IFIH1","DHX58",
                         "TLR3",
                         "MAVS","TICAM1",
                         "TBK1","IKBKE"),
  Innate_Immune      = c("JAK1","TYK2","STAT1","STAT2","IRF9","RELA",
                         "ISG15","USP18","OAS1","OAS2","OAS3","OASL",
                         "RNASEL","EIF2AK2","PRKRA",
                         "MX1","IFIT1","IFIT2","IFIT3",
                         "C3","C4A"),
  TE_Control         = c("MOV10","TRIM28","SETDB1","DNMT1",
                         "DICER1","AGO2","L1TD1","PIWIL1")
)
MODULE_NAMES <- names(MODULE_GENES)
stopifnot(length(MODULE_NAMES) == 8)

# =============================================================================
# 4. HELPER FUNCTIONS
# =============================================================================

detect_col <- function(meta, candidates, label) {
  for (nm in candidates) if (nm %in% colnames(meta)) {
    cat(sprintf("  [OK] %s column: '%s'\n", label, nm)); return(nm)
  }
  cat(sprintf("  [WARN] %s column not found.\n", label)); NULL
}

compute_module_score_per_donor <- function(seurat_obj, genes, donor_col) {
  genes_present <- genes[genes %in% rownames(seurat_obj)]
  if (length(genes_present) == 0) return(NULL)
  expr_mat <- GetAssayData(seurat_obj, layer = "data")[genes_present, , drop = FALSE]
  cell_scores <- Matrix::colMeans(expr_mat)
  meta <- seurat_obj@meta.data
  meta$donor <- as.character(meta[[donor_col]])
  meta$score <- cell_scores[rownames(meta)]
  donor_means <- meta %>%
    dplyr::group_by(donor) %>%
    dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")
  setNames(donor_means$mean_score, donor_means$donor)
}

add_stage_to_seurat <- function(sobj, donor_col, dataset_col, stage_meta) {
  meta <- sobj@meta.data
  
  # Verify the Seurat object already has age (canonical donor metadata).
  # If not, the user must add it before running this script — see SECTION 2.
  age_candidates <- c("age","Age","age_at_death","donor_age","AGE")
  age_col_match  <- intersect(age_candidates, colnames(meta))
  if (length(age_col_match) == 0) {
    stop("add_stage_to_seurat: no age column in Seurat metadata. ",
         "Tried: ", paste(age_candidates, collapse = ", "),
         ". Add donor age to the Seurat object before running this script.")
  }
  cat(sprintf("  [INFO] using age column '%s' from Seurat metadata.\n",
              age_col_match[1]))
  
  # Drop any pre-existing stage_* columns so the join overwrites cleanly
  # (relevant if this section is re-run on the same object).
  for (col in c("braak_asyn_stage","stage_group","stage_numeric")) {
    if (col %in% colnames(meta)) meta[[col]] <- NULL
  }
  
  meta$.donor_tmp   <- as.character(meta[[donor_col]])
  meta$.dataset_tmp <- as.character(meta[[dataset_col]])
  meta <- meta %>%
    dplyr::left_join(
      stage_meta %>%
        dplyr::select(sample_id, dataset,
                      braak_asyn_stage, stage_group, stage_numeric),
      by = c(".donor_tmp" = "sample_id", ".dataset_tmp" = "dataset"))
  rownames(meta)  <- rownames(sobj@meta.data)
  meta$.donor_tmp <- NULL; meta$.dataset_tmp <- NULL
  sobj@meta.data  <- meta
  sobj
}

# Compute partial Spearman rho with N_BOOT percentile bootstrap CI
partial_rho_with_boot <- function(y, x, covars, n_boot = N_BOOT) {
  idx <- complete.cases(y, x, covars)
  if (sum(idx) < 6) return(list(rho = NA, p = NA, ci_low = NA, ci_high = NA))
  y <- y[idx]; x <- x[idx]; covars <- covars[idx, , drop = FALSE]
  pc <- tryCatch(
    ppcor::pcor.test(y, x, covars, method = "spearman"),
    error = function(e) NULL)
  if (is.null(pc)) return(list(rho = NA, p = NA, ci_low = NA, ci_high = NA))
  rho0 <- pc$estimate; p0 <- pc$p.value
  n <- length(y)
  boot_rhos <- replicate(n_boot, {
    s <- sample.int(n, n, replace = TRUE)
    tryCatch(
      ppcor::pcor.test(y[s], x[s], covars[s, , drop = FALSE],
                       method = "spearman")$estimate,
      error = function(e) NA_real_)
  })
  ci <- quantile(boot_rhos, c(0.025, 0.975), na.rm = TRUE)
  list(rho = rho0, p = p0, ci_low = ci[1], ci_high = ci[2])
}

# =============================================================================
# 5. LOAD AND FILTER DATA
# =============================================================================
cat("========== LOADING DATA ==========\n")
astro <- readRDS(ASTRO_RDS)
neuron <- readRDS(NEURON_RDS)
cat(sprintf("  Astrocytes: %d cells | Neurons: %d cells\n", ncol(astro), ncol(neuron)))

donor_col_a   <- detect_col(astro@meta.data,
                            c(DONOR_COL,"donor_id","Donor_ID","donor","Donor","Sample_ID"),
                            "Astro donor")
donor_col_n   <- detect_col(neuron@meta.data,
                            c(DONOR_COL,"donor_id","Donor_ID","donor","Donor","Sample_ID"),
                            "Neuron donor")
dataset_col_a <- detect_col(astro@meta.data,
                            c(DATASET_COL,"Dataset","orig.ident"), "Astro dataset")
dataset_col_n <- detect_col(neuron@meta.data,
                            c(DATASET_COL,"Dataset","orig.ident"), "Neuron dataset")

astro  <- astro[,  astro@meta.data[[dataset_col_a]]  %in% c("GSE178146","GSE253462")]
neuron <- neuron[, neuron@meta.data[[dataset_col_n]] %in% c("GSE178146","GSE253462")]

astro  <- add_stage_to_seurat(astro,  donor_col_a, dataset_col_a, STAGE_META)
neuron <- add_stage_to_seurat(neuron, donor_col_n, dataset_col_n, STAGE_META)
astro  <- astro[,  !is.na(astro@meta.data$braak_asyn_stage)]
neuron <- neuron[, !is.na(neuron@meta.data$braak_asyn_stage)]
cat(sprintf("After stage filter - Astrocytes: %d | Neurons: %d\n\n",
            ncol(astro), ncol(neuron)))

# Split neurons into excitatory/inhibitory
subtype_col <- detect_col(neuron@meta.data,
                          c(NEURON_SUBTYPE_COL,"cell_type","subtype","Subtype"),
                          "Neuron subtype")
neuron_exc <- NULL; neuron_inh <- NULL
if (!is.null(subtype_col)) {
  subtypes_present <- unique(neuron@meta.data[[subtype_col]])
  exc_label <- subtypes_present[grepl("excit", subtypes_present, ignore.case = TRUE)][1]
  inh_label <- subtypes_present[grepl("inhib", subtypes_present, ignore.case = TRUE)][1]
  if (!is.na(exc_label)) neuron_exc <- neuron[, neuron@meta.data[[subtype_col]] == exc_label]
  if (!is.na(inh_label)) neuron_inh <- neuron[, neuron@meta.data[[subtype_col]] == inh_label]
}

# =============================================================================
# 6. COMPUTE DONOR-LEVEL PSEUDOBULK MODULE SCORES
# =============================================================================
cat("========== COMPUTING DONOR-LEVEL MODULE SCORES ==========\n")

build_donor_score_table <- function(sobj, donor_col, label) {
  if (is.null(sobj)) return(NULL)
  cat(sprintf("  %s\n", label))
  mod_scores <- lapply(MODULE_NAMES, function(m) {
    s <- compute_module_score_per_donor(sobj, MODULE_GENES[[m]], donor_col)
    if (is.null(s)) return(NULL)
    data.frame(donor_id = names(s), score = as.numeric(s),
               module = m, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), mod_scores))
  out$cell_type <- label
  out
}

scores_astro <- build_donor_score_table(astro,      donor_col_a, "Astrocyte")
scores_total <- build_donor_score_table(neuron,     donor_col_n, "Total_Neuron")
scores_exc   <- build_donor_score_table(neuron_exc, donor_col_n, "Excitatory_Neuron")
scores_inh   <- build_donor_score_table(neuron_inh, donor_col_n, "Inhibitory_Neuron")

to_wide <- function(df) {
  if (is.null(df)) return(NULL)
  tidyr::pivot_wider(df, names_from = module, values_from = score,
                     id_cols = c(donor_id, cell_type))
}

astro_wide <- to_wide(scores_astro)
total_wide <- to_wide(scores_total)
exc_wide   <- to_wide(scores_exc)
inh_wide   <- to_wide(scores_inh)

attach_stage <- function(wide_df) {
  if (is.null(wide_df)) return(NULL)
  
  # Detect the age column in the Seurat metadata. With the patched
  # add_stage_to_seurat() above, the user's pre-existing age column is
  # preserved as-is (could be 'age', 'Age', 'age_at_death', etc.).
  age_candidates <- c("age","Age","age_at_death","donor_age","AGE")
  age_col_match  <- intersect(age_candidates, colnames(astro@meta.data))
  if (length(age_col_match) == 0) {
    stop("attach_stage: no age column found in astro@meta.data. ",
         "Tried: ", paste(age_candidates, collapse = ", "),
         ". Existing columns: ", paste(colnames(astro@meta.data), collapse = ", "))
  }
  age_col <- age_col_match[1]
  
  meta_summary <- astro@meta.data %>%
    dplyr::mutate(donor_id = as.character(.data[[donor_col_a]]),
                  age      = .data[[age_col]]) %>%
    dplyr::group_by(donor_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(donor_id, dataset = !!dataset_col_a,
                  braak_asyn_stage, stage_group, stage_numeric, age)
  dplyr::left_join(wide_df, meta_summary, by = "donor_id")
}

astro_wide <- attach_stage(astro_wide)
total_wide <- attach_stage(total_wide)
exc_wide   <- attach_stage(exc_wide)
inh_wide   <- attach_stage(inh_wide)

# =============================================================================
# 7. STATISTICAL FRAMEWORK
# =============================================================================

run_braak_stats <- function(wide_df, contrast_label, n_boot = N_BOOT) {
  if (is.null(wide_df) || nrow(wide_df) == 0) return(NULL)
  cat(sprintf("--- %s ---\n", contrast_label))
  
  results <- lapply(MODULE_NAMES, function(mod) {
    if (!mod %in% colnames(wide_df)) return(NULL)
    dat <- wide_df[!is.na(wide_df[[mod]]) & !is.na(wide_df$stage_numeric), ]
    if (nrow(dat) < 6) return(NULL)
    
    n_total   <- nrow(dat)
    n_ctrl    <- sum(dat$stage_group == "Control",   na.rm = TRUE)
    n_b45     <- sum(dat$stage_group == "Braak_4-5", na.rm = TRUE)
    n_b6      <- sum(dat$stage_group == "Braak_6",   na.rm = TRUE)
    
    # ---- PRIMARY: partial Spearman + bootstrap CI ----
    dat$dataset_num <- as.numeric(factor(dat$dataset))
    covars <- dat[, c("age","dataset_num"), drop = FALSE]
    boot_res <- partial_rho_with_boot(dat[[mod]], dat$stage_numeric, covars, n_boot)
    partial_rho     <- boot_res$rho
    partial_p       <- boot_res$p
    partial_ci_low  <- boot_res$ci_low
    partial_ci_high <- boot_res$ci_high
    
    # ---- SECONDARY: Jonckheere-Terpstra ----
    jt_p <- NA; jt_stat <- NA
    tryCatch({
      jt <- clinfun::jonckheere.test(
        x = dat[[mod]],
        g = factor(dat$stage_numeric, ordered = TRUE),
        alternative = "two.sided")
      jt_stat <- as.numeric(jt$statistic)
      jt_p    <- jt$p.value
    }, error = function(e) {})
    
    # ---- EXPLORATORY: Kruskal-Wallis ----
    kw_p <- NA
    if (length(unique(dat$stage_group)) >= 2) {
      kw_p <- tryCatch(
        kruskal.test(as.formula(paste0("`", mod, "` ~ stage_group")), data = dat)$p.value,
        error = function(e) NA)
    }
    
    # ---- EXPLORATORY: pairwise Wilcoxon ----
    ctrl <- dat[[mod]][dat$stage_group == "Control"]
    b45  <- dat[[mod]][dat$stage_group == "Braak_4-5"]
    b6   <- dat[[mod]][dat$stage_group == "Braak_6"]
    w_ctrl_b45 <- if (length(ctrl) >= 2 && length(b45) >= 2)
      wilcox.test(ctrl, b45, exact = FALSE)$p.value else NA
    w_ctrl_b6  <- if (length(ctrl) >= 2 && length(b6) >= 2)
      wilcox.test(ctrl, b6, exact = FALSE)$p.value else NA
    w_b45_b6   <- if (length(b45) >= 2 && length(b6) >= 2)
      wilcox.test(b45, b6, exact = FALSE)$p.value else NA
    
    direction <- ifelse(mean(b6, na.rm = TRUE) > mean(ctrl, na.rm = TRUE),
                        "INCREASES_with_stage", "DECREASES_with_stage")
    
    cat(sprintf("  %-20s rho=%+.3f [%+.2f,%+.2f] p=%.3f JT_p=%.3f KW_p=%.3f %s\n",
                mod, partial_rho, partial_ci_low, partial_ci_high,
                partial_p, jt_p, kw_p, direction))
    
    data.frame(
      contrast            = contrast_label,
      module              = mod,
      n_total             = n_total,
      n_control           = n_ctrl,
      n_braak45           = n_b45,
      n_braak6            = n_b6,
      partial_rho         = partial_rho,
      partial_ci_low      = partial_ci_low,
      partial_ci_high     = partial_ci_high,
      partial_p           = partial_p,
      JT_statistic        = jt_stat,
      JT_p                = jt_p,
      kw_p_exploratory    = kw_p,
      wilcox_ctrl_vs_b45  = w_ctrl_b45,
      wilcox_ctrl_vs_b6   = w_ctrl_b6,
      wilcox_b45_vs_b6    = w_b45_b6,
      mean_ctrl           = mean(ctrl, na.rm = TRUE),
      mean_b45            = mean(b45,  na.rm = TRUE),
      mean_b6             = mean(b6,   na.rm = TRUE),
      direction           = direction,
      stringsAsFactors    = FALSE
    )
  })
  
  out <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(out) || nrow(out) == 0) return(out)
  
  # BH-FDR within each test family within this contrast
  out$partial_p_fdr        <- p.adjust(out$partial_p, method = "BH")
  out$JT_p_fdr             <- p.adjust(out$JT_p,      method = "BH")
  out$kw_p_fdr_exploratory <- p.adjust(out$kw_p_exploratory, method = "BH")
  out$wilcox_ctrl_vs_b45_fdr <- p.adjust(out$wilcox_ctrl_vs_b45, method = "BH")
  out$wilcox_ctrl_vs_b6_fdr  <- p.adjust(out$wilcox_ctrl_vs_b6,  method = "BH")
  out$wilcox_b45_vs_b6_fdr   <- p.adjust(out$wilcox_b45_vs_b6,   method = "BH")
  out
}

# ============================================================
# Per-cohort replication for per-cell-type analyses
# ============================================================
run_per_cohort_replication <- function(wide_df, contrast_label) {
  if (is.null(wide_df)) return(NULL)
  out <- lapply(MODULE_NAMES, function(mod) {
    if (!mod %in% colnames(wide_df)) return(NULL)
    dat <- wide_df[!is.na(wide_df[[mod]]) & !is.na(wide_df$stage_numeric), ]
    if (nrow(dat) < 6) return(NULL)
    per_ds <- lapply(unique(dat$dataset), function(ds) {
      d <- dat[dat$dataset == ds, ]
      if (length(unique(d$stage_numeric)) < 2) return(NULL)
      rho_p <- tryCatch(
        cor.test(d[[mod]], d$stage_numeric, method = "spearman", exact = FALSE),
        error = function(e) NULL)
      if (is.null(rho_p)) return(NULL)
      data.frame(dataset = ds, n = nrow(d),
                 rho = rho_p$estimate, p = rho_p$p.value,
                 stringsAsFactors = FALSE)
    })
    per_ds <- do.call(rbind, Filter(Negate(is.null), per_ds))
    if (is.null(per_ds) || nrow(per_ds) < 2) return(NULL)
    same_dir <- length(unique(sign(per_ds$rho))) == 1
    both_nom <- all(per_ds$p < 0.05, na.rm = TRUE)
    rep_status <- if (both_nom && same_dir) "REPLICATED" else
      if (same_dir) "SAME_DIRECTION_NS" else "OPPOSITE_DIRECTION"
    data.frame(
      contrast     = contrast_label,
      module       = mod,
      rho_GSE178146 = per_ds$rho[per_ds$dataset == "GSE178146"][1],
      p_GSE178146   = per_ds$p[per_ds$dataset == "GSE178146"][1],
      rho_GSE253462 = per_ds$rho[per_ds$dataset == "GSE253462"][1],
      p_GSE253462   = per_ds$p[per_ds$dataset == "GSE253462"][1],
      replication_status = rep_status,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), out))
}

# =============================================================================
# 8. RUN LAYER 1 - PER-CELL-TYPE
# =============================================================================
cat("\n========== LAYER 1: PER-CELL-TYPE STAGE ANALYSIS ==========\n\n")
stats_astro <- run_braak_stats(astro_wide, "Astrocyte")
stats_total <- run_braak_stats(total_wide, "Total_Neuron")
stats_exc   <- run_braak_stats(exc_wide,   "Excitatory_Neuron")
stats_inh   <- run_braak_stats(inh_wide,   "Inhibitory_Neuron")
cat("\n")

cat("========== Per-cohort replication checks ==========\n")
rep_astro <- run_per_cohort_replication(astro_wide, "Astrocyte")
rep_total <- run_per_cohort_replication(total_wide, "Total_Neuron")
rep_exc   <- run_per_cohort_replication(exc_wide,   "Excitatory_Neuron")
rep_inh   <- run_per_cohort_replication(inh_wide,   "Inhibitory_Neuron")
cat("\n")

# =============================================================================
# 9. RUN LAYER 2 - DIVERGENCE
# =============================================================================
cat("========== LAYER 2: ASTROCYTE-NEURON DIVERGENCE ==========\n\n")

compute_divergence <- function(astro_w, neuron_w, neuron_label) {
  if (is.null(astro_w) || is.null(neuron_w)) return(NULL)
  common <- intersect(astro_w$donor_id, neuron_w$donor_id)
  if (length(common) == 0) return(NULL)
  a <- astro_w[match(common, astro_w$donor_id), ]
  n <- neuron_w[match(common, neuron_w$donor_id), ]
  shared_mods <- intersect(MODULE_NAMES, intersect(colnames(a), colnames(n)))
  div_df <- data.frame(donor_id = a$donor_id, stringsAsFactors = FALSE)
  for (m in shared_mods) div_df[[m]] <- a[[m]] - n[[m]]
  for (mc in c("braak_asyn_stage","stage_group","stage_numeric","dataset","age"))
    if (mc %in% colnames(a)) div_df[[mc]] <- a[[mc]]
  div_df
}

div_total <- compute_divergence(astro_wide, total_wide, "Total_Neuron")
div_exc   <- compute_divergence(astro_wide, exc_wide,   "Excitatory_Neuron")
div_inh   <- compute_divergence(astro_wide, inh_wide,   "Inhibitory_Neuron")

stats_div_total <- run_braak_stats(div_total, "Astro_vs_Total_Neuron")
stats_div_exc   <- run_braak_stats(div_exc,   "Astro_vs_Excitatory_Neuron")
stats_div_inh   <- run_braak_stats(div_inh,   "Astro_vs_Inhibitory_Neuron")
cat("\n")

# =============================================================================
# 10. OUTPUT TREE
# =============================================================================
timestamp   <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- file.path(OUT_DIR, paste0("MS2_Figure9B_Braak_", timestamp))
for (sd in c("1_DonorScores","2_DivergenceScores","3_Statistics","4_Figures"))
  dir.create(file.path(results_dir, sd), recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(results_dir, "4_Figures")

# Save donor scores
if (!is.null(astro_wide)) write.csv(astro_wide,
                                    file.path(results_dir, "1_DonorScores", "Donor_Scores_Astrocyte.csv"), row.names = FALSE)
if (!is.null(total_wide)) write.csv(total_wide,
                                    file.path(results_dir, "1_DonorScores", "Donor_Scores_TotalNeuron.csv"), row.names = FALSE)
if (!is.null(exc_wide))   write.csv(exc_wide,
                                    file.path(results_dir, "1_DonorScores", "Donor_Scores_ExcitatoryNeuron.csv"), row.names = FALSE)
if (!is.null(inh_wide))   write.csv(inh_wide,
                                    file.path(results_dir, "1_DonorScores", "Donor_Scores_InhibitoryNeuron.csv"), row.names = FALSE)

# Save divergence tables
if (!is.null(div_total)) write.csv(div_total,
                                   file.path(results_dir, "2_DivergenceScores", "Divergence_Astro_vs_Total.csv"), row.names = FALSE)
if (!is.null(div_exc)) write.csv(div_exc,
                                 file.path(results_dir, "2_DivergenceScores", "Divergence_Astro_vs_Exc.csv"), row.names = FALSE)
if (!is.null(div_inh)) write.csv(div_inh,
                                 file.path(results_dir, "2_DivergenceScores", "Divergence_Astro_vs_Inh.csv"), row.names = FALSE)

# Save stats CSVs
for (nm in c("stats_astro","stats_total","stats_exc","stats_inh",
             "stats_div_total","stats_div_exc","stats_div_inh",
             "rep_astro","rep_total","rep_exc","rep_inh")) {
  v <- get(nm)
  if (!is.null(v)) write.csv(v,
                             file.path(results_dir, "3_Statistics", paste0(nm, ".csv")), row.names = FALSE)
}

# =============================================================================
# 11. FIGURE 9B - MAIN-TEXT NULL-RESULT PANEL
# =============================================================================
# Design rationale: Results sentence makes two claims:
#   (1) per-cell-type null in astrocytes and total neurons
#   (2) divergence null in astro vs total neuron
# Figure 9B encodes both with a single visual idiom (forest plot of
# partial Spearman rho with bootstrap 95% CI) across two stacked sub-panels.
# Modules are ordered consistently across both sub-panels for direct comparison.
# A point that straddles rho = 0 (vertical reference line) is null.

# Color palette (colorblind-safe)
COL_NULL <- "grey45"
COL_AST  <- "#2C7A4A"  # astrocyte green
COL_NEU  <- "#5D478B"  # neuron purple
COL_DIV  <- "#B85C00"  # divergence orange

# ---- Sub-panel A: per-cell-type (astrocyte and total neuron stacked) ----

prep_forest <- function(stats_df, cell_label, color_hex) {
  if (is.null(stats_df)) return(NULL)
  stats_df %>%
    dplyr::filter(module %in% MODULE_NAMES) %>%
    dplyr::mutate(
      module     = factor(module, levels = rev(MODULE_NAMES)),
      cell_type  = cell_label,
      ci_low     = pmax(-1, partial_ci_low),
      ci_high    = pmin( 1, partial_ci_high),
      label_col  = color_hex,
      sig_label  = ifelse(!is.na(partial_p_fdr) & partial_p_fdr < 0.05, "*", "")
    )
}

forest_astro <- prep_forest(stats_astro, "Astrocyte",     COL_AST)
forest_total <- prep_forest(stats_total, "Total Neuron",  COL_NEU)
forest_div   <- prep_forest(stats_div_total, "Astrocyte − Total Neuron\n(divergence)", COL_DIV)

build_forest_panel <- function(df, panel_title, color_hex) {
  if (is.null(df)) return(NULL)
  ggplot(df, aes(x = partial_rho, y = module)) +
    annotate("rect", xmin = -0.05, xmax = 0.05, ymin = -Inf, ymax = Inf,
             alpha = 0.06, fill = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                   height = 0.25, colour = color_hex, linewidth = 0.7) +
    geom_point(size = 3, colour = color_hex) +
    geom_text(aes(x = ci_high + 0.04, label = sig_label),
              hjust = 0, size = 4, fontface = "bold") +
    coord_cartesian(xlim = c(-0.8, 0.8)) +
    labs(title = panel_title,
         x = "Partial Spearman ρ (95% bootstrap CI)\nstage_numeric ~ score | age + dataset",
         y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 11,
                                   colour = color_hex, margin = margin(b = 4)),
      axis.text.y   = element_text(size = 10),
      axis.title.x  = element_text(size = 9, colour = "grey25"),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3)
    )
}

p_astro <- build_forest_panel(forest_astro,
                              "(i) Astrocyte module score vs Braak stage", COL_AST)
p_total <- build_forest_panel(forest_total,
                              "(ii) Total Neuron module score vs Braak stage", COL_NEU)
p_div   <- build_forest_panel(forest_div,
                              "(iii) Astrocyte − Total Neuron divergence vs Braak stage", COL_DIV)

# Compose Figure 9B
figure_9B <- (p_astro | p_total | p_div) +
  patchwork::plot_annotation(
    title = "Figure 9B. dsRNA / ADAR module activity does not scale with α-synuclein Braak stage in cortex",
    subtitle = sprintf("n = 31 donors (15 Control / 5 Braak 4–5 / 11 Braak 6). Two cohorts (GSE178146 + GSE253462). Asterisk: BH-FDR < 0.05."),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, colour = "grey30",
                                   margin = margin(b = 8))))

ggsave(file.path(fig_dir, "Figure_9B.pdf"), figure_9B, width = 14, height = 6)
ggsave(file.path(fig_dir, "Figure_9B.png"), figure_9B, width = 14, height = 6, dpi = 300)

# =============================================================================
# 12. SUPPLEMENTARY FIGURE - DIVERGENCE BOXPLOTS (primary contrast)
# =============================================================================
stage_colours <- c("Control"="#4E9BB9","Braak_4-5"="#E8A838","Braak_6"="#C0392B")

plot_div_boxbee <- function(div_df, col, title_label) {
  if (is.null(div_df) || !col %in% colnames(div_df)) return(NULL)
  d <- div_df[!is.na(div_df[[col]]) & !is.na(div_df$stage_group), ]
  d$stage_group <- factor(d$stage_group, levels = c("Control","Braak_4-5","Braak_6"))
  ggplot(d, aes(x = stage_group, y = .data[[col]],
                fill = stage_group, colour = stage_group)) +
    geom_boxplot(alpha = 0.25, outlier.shape = NA, width = 0.45) +
    ggbeeswarm::geom_beeswarm(size = 2.2, cex = 2.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_fill_manual(values = stage_colours) +
    scale_colour_manual(values = stage_colours) +
    labs(title = title_label, x = NULL,
         y = "Astrocyte − Total Neuron") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 10),
          axis.text.x = element_text(size = 8))
}

if (!is.null(div_total)) {
  panels <- lapply(MODULE_NAMES, function(m) plot_div_boxbee(div_total, m, m))
  panels <- Filter(Negate(is.null), panels)
  if (length(panels) > 0) {
    pp <- patchwork::wrap_plots(panels, ncol = 4) +
      patchwork::plot_annotation(
        title = "Supplementary Figure. Astrocyte – Total Neuron divergence by Braak stage",
        subtitle = "Primary contrast — eight modules. All comparisons null after FDR.",
        theme = theme(plot.title = element_text(face = "bold", size = 13)))
    ggsave(file.path(fig_dir, "SuppFig_Divergence_Boxplots.pdf"),
           pp, width = 14, height = 7)
  }
}

# =============================================================================
# 13. SUPPLEMENTARY TABLE SX - PER-CELL-TYPE
# =============================================================================
cat("Building Supplementary Table SX (per-cell-type)...\n")
wb_sx <- createWorkbook()

addWorksheet(wb_sx, "README")
readme_sx <- data.frame(
  Item = c("Title","Purpose","Donor n","Stage encoding","Primary test",
           "CI method","Secondary tests","Exploratory tests","FDR correction",
           "Sheet 'Stats_Astrocyte'","Sheet 'Stats_TotalNeuron'",
           "Sheet 'Stats_ExcitatoryNeuron'","Sheet 'Stats_InhibitoryNeuron'",
           "Sheet 'Replication_*'","Modules","Reference figure"),
  Description = c(
    "Supplementary Table SX - per-cell-type Braak stage statistics for MS2 Figure 9B",
    "Tests whether each module score changes with alpha-synuclein Braak stage WITHIN each cell type, separately",
    "31 cortical donors total: 15 Control (Braak 0), 5 PD Braak 4-5, 11 PD Braak 6 (GSE178146 + GSE253462)",
    "Ordinal numeric: 0 = Control, 1 = Braak 4-5, 2 = Braak 6",
    "Partial Spearman correlation between stage_numeric and module score, controlling for age and dataset (ppcor::pcor.test)",
    "10,000-resample percentile bootstrap on the partial Spearman rho",
    "Jonckheere-Terpstra ordered trend test (clinfun::jonckheere.test, two-sided)",
    "Kruskal-Wallis across the three stage groups; pairwise Wilcoxon (Control vs B4-5; Control vs B6; B4-5 vs B6)",
    "BH within each test family within each contrast (8 modules)",
    "Astrocyte module scores vs Braak stage",
    "Total Neuron module scores vs Braak stage",
    "Excitatory Neuron module scores vs Braak stage (sensitivity, given prior cell-type heterogeneity)",
    "Inhibitory Neuron module scores vs Braak stage (sensitivity)",
    "Per-cohort Spearman correlations + replication status (REPLICATED requires both cohorts p<0.05 and same direction)",
    "Eight modules: ADAR1 (ADAR), ADAR2 (ADARB1), SNCA, dsRNA_Processing, dsRNA_EV_Transport, dsRNA_Sensing, Innate_Immune, TE_Control",
    "Cited in main-text Figure 9B"
  ),
  stringsAsFactors = FALSE
)
writeDataTable(wb_sx, "README", readme_sx, tableStyle = "TableStyleLight2")

add_stats_sheet <- function(wb, name, df) {
  if (is.null(df) || nrow(df) == 0) return()
  addWorksheet(wb, name)
  writeDataTable(wb, name, df, tableStyle = "TableStyleLight9")
  pidx <- which(colnames(df) == "partial_p_fdr")
  if (length(pidx) > 0)
    conditionalFormatting(wb, name, cols = pidx, rows = 2:(nrow(df) + 1),
                          rule = "<0.05",
                          style = createStyle(fontColour = "#006100", bgFill = "#C6EFCE"))
}

add_stats_sheet(wb_sx, "Stats_Astrocyte",         stats_astro)
add_stats_sheet(wb_sx, "Stats_TotalNeuron",       stats_total)
add_stats_sheet(wb_sx, "Stats_ExcitatoryNeuron",  stats_exc)
add_stats_sheet(wb_sx, "Stats_InhibitoryNeuron",  stats_inh)

add_simple_sheet <- function(wb, name, df) {
  if (is.null(df) || nrow(df) == 0) return()
  addWorksheet(wb, name)
  writeDataTable(wb, name, df, tableStyle = "TableStyleLight2")
}
add_simple_sheet(wb_sx, "Replication_Astrocyte",        rep_astro)
add_simple_sheet(wb_sx, "Replication_TotalNeuron",      rep_total)
add_simple_sheet(wb_sx, "Replication_ExcitatoryNeuron", rep_exc)
add_simple_sheet(wb_sx, "Replication_InhibitoryNeuron", rep_inh)

# Per-donor scores (one sheet per cell type)
add_donor_sheet <- function(wb, name, df) {
  if (is.null(df) || nrow(df) == 0) return()
  out_cols <- intersect(c("donor_id","dataset","braak_asyn_stage","stage_group","age", MODULE_NAMES),
                        colnames(df))
  addWorksheet(wb, name)
  writeDataTable(wb, name,
                 df[order(df$stage_group, df$donor_id), out_cols],
                 tableStyle = "TableStyleLight2")
}
add_donor_sheet(wb_sx, "DonorScores_Astrocyte",        astro_wide)
add_donor_sheet(wb_sx, "DonorScores_TotalNeuron",      total_wide)
add_donor_sheet(wb_sx, "DonorScores_ExcitatoryNeuron", exc_wide)
add_donor_sheet(wb_sx, "DonorScores_InhibitoryNeuron", inh_wide)

sx_path <- file.path(results_dir, "Supplementary_Table_SX_PerCellType.xlsx")
saveWorkbook(wb_sx, sx_path, overwrite = TRUE)
cat("Saved:", sx_path, "\n")

# =============================================================================
# 14. SUPPLEMENTARY TABLE SY - DIVERGENCE
# =============================================================================
cat("Building Supplementary Table SY (divergence)...\n")
wb_sy <- createWorkbook()

addWorksheet(wb_sy, "README")
readme_sy <- data.frame(
  Item = c("Title","Purpose","Donor n","Stage encoding","Divergence definition",
           "Primary test","CI method","Secondary tests","Exploratory tests",
           "FDR correction","Primary contrast","Sensitivity contrasts",
           "Sheet 'Stats_Astro_vs_TotalNeuron'",
           "Sheet 'Stats_Astro_vs_ExcitatoryNeuron'",
           "Sheet 'Stats_Astro_vs_InhibitoryNeuron'",
           "Modules","Reference figure"),
  Description = c(
    "Supplementary Table SY - astrocyte-neuron divergence Braak stage statistics for MS2 Figure 9B",
    "Tests whether the GAP between astrocyte and neuron module score changes with Braak stage",
    "31 cortical donors total: 15 Control, 5 PD Braak 4-5, 11 PD Braak 6 (GSE178146 + GSE253462)",
    "Ordinal numeric: 0 = Control, 1 = Braak 4-5, 2 = Braak 6",
    "Within-donor astrocyte score MINUS matched neuron score, per module",
    "Partial Spearman correlation between stage_numeric and divergence, controlling for age and dataset (ppcor::pcor.test)",
    "10,000-resample percentile bootstrap on the partial Spearman rho",
    "Jonckheere-Terpstra ordered trend test (clinfun::jonckheere.test, two-sided)",
    "Kruskal-Wallis across stage groups; pairwise Wilcoxon",
    "BH within each test family within each contrast (8 modules)",
    "Astrocyte vs Total Neuron (matches Figure 4 and other MS2 cell-type comparisons)",
    "Astrocyte vs Excitatory Neuron AND Astrocyte vs Inhibitory Neuron (sensitivity contrasts; excitatory neurons are the primary cortical population affected by alpha-synuclein pathology, motivating the Exc-specific contrast)",
    "Astrocyte minus Total Neuron divergence vs Braak stage (primary)",
    "Astrocyte minus Excitatory Neuron divergence vs Braak stage (sensitivity)",
    "Astrocyte minus Inhibitory Neuron divergence vs Braak stage (sensitivity)",
    "Eight modules: ADAR1, ADAR2, SNCA, dsRNA_Processing, dsRNA_EV_Transport, dsRNA_Sensing, Innate_Immune, TE_Control",
    "Cited in main-text Figure 9B"
  ),
  stringsAsFactors = FALSE
)
writeDataTable(wb_sy, "README", readme_sy, tableStyle = "TableStyleLight2")

add_stats_sheet(wb_sy, "Stats_Astro_vs_TotalNeuron",      stats_div_total)
add_stats_sheet(wb_sy, "Stats_Astro_vs_ExcitatoryNeuron", stats_div_exc)
add_stats_sheet(wb_sy, "Stats_Astro_vs_InhibitoryNeuron", stats_div_inh)

add_donor_div_sheet <- function(wb, name, df) {
  if (is.null(df) || nrow(df) == 0) return()
  out_cols <- intersect(c("donor_id","dataset","braak_asyn_stage","stage_group","age", MODULE_NAMES),
                        colnames(df))
  # Excel hard limit: sheet names ≤ 31 chars. Truncate with warning rather than crash.
  if (nchar(name) > 31) {
    warning(sprintf("Sheet name '%s' (%d chars) truncated to 31.", name, nchar(name)))
    name <- substr(name, 1, 31)
  }
  addWorksheet(wb, name)
  writeDataTable(wb, name,
                 df[order(df$stage_group, df$donor_id), out_cols],
                 tableStyle = "TableStyleLight2")
}
add_donor_div_sheet(wb_sy, "DivScores_Astro_vs_TotalNeuron", div_total)
add_donor_div_sheet(wb_sy, "DivScores_Astro_vs_Exc",         div_exc)
add_donor_div_sheet(wb_sy, "DivScores_Astro_vs_Inh",         div_inh)

sy_path <- file.path(results_dir, "Supplementary_Table_SY_Divergence.xlsx")
saveWorkbook(wb_sy, sy_path, overwrite = TRUE)
cat("Saved:", sy_path, "\n")

# =============================================================================
# 15. SESSION INFO
# =============================================================================
writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "session_info.txt"))

cat("\n============================================================\n")
cat("MS2 FIGURE 9B + SUPPLEMENTARY TABLES SX/SY COMPLETE\n")
cat("============================================================\n")
cat("Output directory:", results_dir, "\n\n")
cat("Files for the manuscript:\n")
cat("  - 4_Figures/Figure_9B.pdf (and .png)\n")
cat("  - Supplementary_Table_SX_PerCellType.xlsx\n")
cat("  - Supplementary_Table_SY_Divergence.xlsx\n")
cat("\nQuick summary (primary contrast partial rho with bootstrap CI):\n")
if (!is.null(stats_astro)) {
  cat("\nAstrocyte:\n")
  print(stats_astro[, c("module","partial_rho","partial_ci_low","partial_ci_high",
                        "partial_p","partial_p_fdr")])
}
if (!is.null(stats_total)) {
  cat("\nTotal Neuron:\n")
  print(stats_total[, c("module","partial_rho","partial_ci_low","partial_ci_high",
                        "partial_p","partial_p_fdr")])
}
if (!is.null(stats_div_total)) {
  cat("\nAstrocyte - Total Neuron divergence:\n")
  print(stats_div_total[, c("module","partial_rho","partial_ci_low","partial_ci_high",
                            "partial_p","partial_p_fdr")])
}
cat("\n")
