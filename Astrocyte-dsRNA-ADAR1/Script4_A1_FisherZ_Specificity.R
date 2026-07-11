# ============================================================
# A1 REACTIVE MARKERS × dsRNA/IMMUNE MODULES
# Astrocyte-Specificity Test via Fisher Z Transformation
# ============================================================
#
# PURPOSE:
#   For each A1 marker gene × dsRNA/immune module pair, compute
#   donor-level Spearman correlations SEPARATELY in astrocytes
#   and neurons, then test whether the astrocyte correlation is
#   significantly stronger using Fisher Z transformation.
#
# METHODOLOGY:
#   1. Background-subtracted module scores (removes global
#      expression artefact for cross-cell-type comparison).
#   2. Donor-level pseudobulk (mean per donor).
#   3. QC residualisation (n_cells, mean_nCount, age, PMI, sex)
#      to remove technical confounds.
#   4. Spearman ρ per cell type; Fisher Z for ρ_astro vs ρ_neuron.
#   5. BH correction across all tests.
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
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
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

library(here)
astro_file  <- here::here("data", "CX_Astrocytes_withCov.rds")
neuron_file <- here::here("data", "CX_Neurons_Total_withCov.rds")
results_rds <- here::here("data", "CX_8Module_AllResults.rds")
OUT_DIR     <- here::here("results")

# =============================================================================
# SECTION 2: GENE & MODULE DEFINITIONS
# =============================================================================

# --- A1 neurotoxic-reactive markers (Liddelow et al. 2017) ---
A1_GENES <- c(
  "C3", "SERPING1", "GBP2", "FKBP5", "PSMB8",
  "SRGN", "AMIGO2", "C1R", "C1S", "CFB", "CFI"
)

# --- A2 neuroprotective markers (contrast set) ---
A2_GENES <- c(
  "EMP1", "S100A10", "SLC10A6", "TM4SF1", "CD109",
  "PTX3", "SPHK1", "CD14", "TGFB2", "B3GNT5"
)

# --- Pan-reactive ---
PAN_REACTIVE_GENES <- c("GFAP", "VIM", "SERPINA3", "LCN2")

ALL_MARKER_GENES <- unique(c(A1_GENES, A2_GENES, PAN_REACTIVE_GENES))

# --- Modules ---
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

# --- Background-subtracted module scoring (removes library-size artefact) ---
calculate_module_scores_bgsub <- function(seurat_obj, module_list, n_ctrl = 100) {
  expr <- GetAssayData(seurat_obj, layer = "data")
  all_genes <- rownames(expr)
  gene_means <- Matrix::rowMeans(expr)
  gene_bins  <- cut(gene_means, breaks = 25, labels = FALSE)
  names(gene_bins) <- names(gene_means)

  scores <- sapply(names(module_list), function(mod) {
    genes <- module_list[[mod]]
    genes <- genes[genes %in% all_genes]
    if (length(genes) == 0) return(rep(NA_real_, ncol(seurat_obj)))

    # Foreground
    fg <- if (length(genes) == 1) {
      as.numeric(expr[genes, ])
    } else {
      colMeans(expr[genes, , drop = FALSE])
    }

    # Background: expression-matched control genes
    module_bins <- unique(gene_bins[genes])
    candidates  <- names(gene_bins)[gene_bins %in% module_bins]
    candidates  <- setdiff(candidates, genes)
    if (length(candidates) > n_ctrl) candidates <- sample(candidates, n_ctrl)
    bg <- if (length(candidates) > 0) colMeans(expr[candidates, , drop = FALSE]) else 0

    fg - bg
  })
  colnames(scores) <- paste0(names(module_list), "_Score")
  as.data.frame(scores)
}

# --- Donor-level pseudobulk aggregation ---
aggregate_to_donor <- function(seurat_obj, module_list, marker_genes, donor_col, cond_col) {
  expr <- GetAssayData(seurat_obj, layer = "data")
  meta <- as.data.frame(seurat_obj@meta.data)

  # Module scores (background-subtracted)
  mod_scores <- calculate_module_scores_bgsub(seurat_obj, module_list)

  # Marker gene expression (raw, no bg-sub needed for individual genes
  # because we residualise at donor level below)
  marker_present <- marker_genes[marker_genes %in% rownames(expr)]
  marker_expr    <- t(as.matrix(expr[marker_present, , drop = FALSE]))

  # Combine
  cell_df <- data.frame(
    donor_id  = as.character(meta[[donor_col]]),
    condition = if (!is.null(cond_col)) as.character(meta[[cond_col]]) else "All",
    n_count   = if ("nCount_RNA" %in% colnames(meta)) meta$nCount_RNA else NA_real_,
    stringsAsFactors = FALSE
  )
  cell_df <- cbind(cell_df, mod_scores, marker_expr)

  # Aggregate to donor
  value_cols <- c(colnames(mod_scores), marker_present, "n_count")
  donor_df <- cell_df %>%
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
        donor_id = as.character(meta[[donor_col]]),
        meta[, cov_cols, drop = FALSE],
        stringsAsFactors = FALSE
      ) %>%
        group_by(donor_id) %>%
        summarise(across(everything(), ~ first(na.omit(.))), .groups = "drop") %>%
        as.data.frame()
      donor_df <- merge(donor_df, cov_df, by = "donor_id", all.x = TRUE)
    }
  }

  list(donor_df = donor_df, marker_present = marker_present)
}

# --- QC + covariate residualisation ---
residualise_donor <- function(y, df) {
  # Regress out: n_cells, mean_nCount, age, PMI, sex
  cov_list <- c("n_cells", "n_count")
  if (COVARIATES_ENABLED) {
    cov_list <- c(cov_list, COV_AGE_COL, COV_PMI_COL, COV_SEX_COL)
  }
  cov_list <- cov_list[cov_list %in% colnames(df)]

  # Filter to usable covariates
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

# --- Fisher Z test ---
fisher_z_test <- function(rho1, n1, rho2, n2) {
  # Two-tailed test: H0: rho1 == rho2
  z1 <- atanh(rho1)
  z2 <- atanh(rho2)
  se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  z_diff <- (z1 - z2) / se
  p <- 2 * pnorm(-abs(z_diff))
  list(z_diff = z_diff, p_value = p)
}

# =============================================================================
# SECTION 4: LOAD DATA & COMPUTE DONOR PSEUDOBULK
# =============================================================================

cat("========== LOADING DATA ==========\n")

astro  <- readRDS(astro_file)
neuron <- readRDS(neuron_file)
cat("  Astrocytes:", ncol(astro), "cells\n")
cat("  Neurons:   ", ncol(neuron), "cells\n")

# Detect columns
find_col <- function(meta, candidates) {
  for (nm in candidates) if (nm %in% colnames(meta)) return(nm)
  NULL
}

donor_col_a <- find_col(astro@meta.data,  c("donor_id","Donor","donor","subject_id","sample_id","orig.ident"))
donor_col_n <- find_col(neuron@meta.data, c("donor_id","Donor","donor","subject_id","sample_id","orig.ident"))
cond_col_a  <- find_col(astro@meta.data,  c("condition","Condition","disease","Disease","group","Group","diagnosis"))
cond_col_n  <- find_col(neuron@meta.data, c("condition","Condition","disease","Disease","group","Group","diagnosis"))

cat("  Donor cols:  astro =", donor_col_a, "| neuron =", donor_col_n, "\n")
cat("  Cond cols:   astro =", cond_col_a,  "| neuron =", cond_col_n, "\n\n")

cat("========== COMPUTING DONOR-LEVEL PSEUDOBULK ==========\n")

res_a <- aggregate_to_donor(astro,  module_gene_sets, ALL_MARKER_GENES, donor_col_a, cond_col_a)
res_n <- aggregate_to_donor(neuron, module_gene_sets, ALL_MARKER_GENES, donor_col_n, cond_col_n)

pb_astro  <- res_a$donor_df
pb_neuron <- res_n$donor_df

# Use intersection of marker genes present in both cell types
markers_both <- intersect(res_a$marker_present, res_n$marker_present)
cat("  Marker genes in both cell types:", length(markers_both), "\n")
cat("  Donors: astro =", nrow(pb_astro), "| neuron =", nrow(pb_neuron), "\n\n")

# =============================================================================
# SECTION 5: RESIDUALISE
# =============================================================================

cat("========== RESIDUALISING ==========\n")

mod_score_cols <- paste0(names(module_gene_sets), "_Score")

for (col in c(mod_score_cols, markers_both)) {
  if (col %in% colnames(pb_astro))
    pb_astro[[paste0(col, "_resid")]]  <- residualise_donor(pb_astro[[col]],  pb_astro)
  if (col %in% colnames(pb_neuron))
    pb_neuron[[paste0(col, "_resid")]] <- residualise_donor(pb_neuron[[col]], pb_neuron)
}

cat("  Done.\n\n")

# =============================================================================
# SECTION 6: SPEARMAN CORRELATIONS + FISHER Z
# =============================================================================

cat("========== SPEARMAN + FISHER Z ==========\n")

results_list <- list()

for (mod in names(module_gene_sets)) {
  mod_col <- paste0(mod, "_Score_resid")

  for (gene in markers_both) {
    gene_col <- paste0(gene, "_resid")

    # --- Astrocyte ---
    xa <- pb_astro[[mod_col]];  ya <- pb_astro[[gene_col]]
    ok_a <- complete.cases(xa, ya); n_a <- sum(ok_a)

    # --- Neuron ---
    xn <- pb_neuron[[mod_col]]; yn <- pb_neuron[[gene_col]]
    ok_n <- complete.cases(xn, yn); n_n <- sum(ok_n)

    if (n_a < 5 || n_n < 5) next

    ct_a <- cor.test(xa[ok_a], ya[ok_a], method = "spearman", exact = FALSE)
    ct_n <- cor.test(xn[ok_n], yn[ok_n], method = "spearman", exact = FALSE)

    fz <- fisher_z_test(ct_a$estimate, n_a, ct_n$estimate, n_n)

    marker_class <- case_when(
      gene %in% A1_GENES         ~ "A1",
      gene %in% A2_GENES         ~ "A2",
      gene %in% PAN_REACTIVE_GENES ~ "Pan-reactive",
      TRUE                       ~ "Other"
    )

    results_list[[length(results_list) + 1]] <- data.frame(
      Module         = mod,
      Marker_Gene    = gene,
      Marker_Class   = marker_class,
      rho_Astro      = ct_a$estimate,
      p_Astro        = ct_a$p.value,
      n_Astro        = n_a,
      rho_Neuron     = ct_n$estimate,
      p_Neuron       = ct_n$p.value,
      n_Neuron       = n_n,
      Fisher_Z       = fz$z_diff,
      Fisher_p       = fz$p_value,
      rho_diff       = ct_a$estimate - ct_n$estimate,
      stringsAsFactors = FALSE
    )
  }
}

results_df <- bind_rows(results_list)

# BH adjustment — separate for within-cell-type and Fisher Z
results_df$p_adj_Astro  <- p.adjust(results_df$p_Astro,  method = "BH")
results_df$p_adj_Neuron <- p.adjust(results_df$p_Neuron, method = "BH")
results_df$Fisher_p_adj <- p.adjust(results_df$Fisher_p, method = "BH")

results_df <- results_df %>% arrange(Fisher_p_adj, Fisher_p)

cat("  Total tests:", nrow(results_df), "\n")
cat("  Fisher Z FDR < 0.05:", sum(results_df$Fisher_p_adj < 0.05), "\n")
cat("  Fisher Z FDR < 0.10:", sum(results_df$Fisher_p_adj < 0.10), "\n\n")

# =============================================================================
# SECTION 7: OUTPUT TABLES
# =============================================================================

cat("========== WRITING RESULTS ==========\n")

timestamp   <- format(Sys.time(), "%Y%m%d_%H%M")
results_dir <- file.path(output_base, paste0(REGION, "_A1_FisherZ_", timestamp))
dir.create(file.path(results_dir, "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

wb <- createWorkbook()

# Sheet 1: Full results
addWorksheet(wb, "All_Results")
writeData(wb, "All_Results", results_df)

# Sheet 2: A1 genes only, sorted by Fisher Z
a1_df <- results_df %>% filter(Marker_Class == "A1") %>% arrange(Fisher_p_adj)
addWorksheet(wb, "A1_FisherZ")
writeData(wb, "A1_FisherZ", a1_df)

# Sheet 3: Focus — A1 × dsRNA_Sensing + Innate_Immune
focus_df <- a1_df %>% filter(Module %in% c("dsRNA_Sensing", "Innate_Immune"))
addWorksheet(wb, "A1_FocusModules")
writeData(wb, "A1_FocusModules", focus_df)

# Sheet 4: Significance summary pivot
sig_summary <- results_df %>%
  filter(Marker_Class == "A1") %>%
  mutate(
    astro_sig  = ifelse(p_adj_Astro  < 0.05, sprintf("%.2f*", rho_Astro),  sprintf("%.2f",  rho_Astro)),
    neuron_sig = ifelse(p_adj_Neuron < 0.05, sprintf("%.2f*", rho_Neuron), sprintf("%.2f",  rho_Neuron)),
    fisher_sig = ifelse(Fisher_p_adj < 0.05, "***",
                 ifelse(Fisher_p_adj < 0.10, "**",
                 ifelse(Fisher_p_adj < 0.20, "*", "")))
  ) %>%
  select(Module, Marker_Gene, rho_Astro, rho_Neuron, rho_diff, Fisher_Z, Fisher_p_adj, fisher_sig)

addWorksheet(wb, "A1_Summary")
writeData(wb, "A1_Summary", sig_summary)

saveWorkbook(wb, file.path(results_dir, "tables",
             paste0(REGION, "_A1_FisherZ_Results.xlsx")),
             overwrite = TRUE)

write.csv(results_df,
          file.path(results_dir, "tables",
                    paste0(REGION, "_A1_FisherZ_ALL.csv")),
          row.names = FALSE)

cat("  Tables saved.\n\n")

# =============================================================================
# SECTION 8: FIGURES
# =============================================================================

cat("========== GENERATING FIGURES ==========\n")

# --- 8A: Paired bar chart (rho_Astro vs rho_Neuron) for A1 × focus modules ---
bar_df <- results_df %>%
  filter(Marker_Class == "A1", Module %in% c("dsRNA_Sensing", "Innate_Immune")) %>%
  select(Module, Marker_Gene, rho_Astro, rho_Neuron, Fisher_p_adj) %>%
  pivot_longer(cols = c(rho_Astro, rho_Neuron),
               names_to = "Cell_Type", values_to = "rho") %>%
  mutate(
    Cell_Type = gsub("rho_", "", Cell_Type),
    star = ifelse(Fisher_p_adj < 0.05, "***",
           ifelse(Fisher_p_adj < 0.10, "**",
           ifelse(Fisher_p_adj < 0.20, "*", "")))
  )

p_bar <- ggplot(bar_df, aes(x = Marker_Gene, y = rho, fill = Cell_Type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_wrap(~ Module, scales = "free_x") +
  scale_fill_manual(values = c("Astro" = "#E64B35", "Neuron" = "#4DBBD5")) +
  labs(
    title = paste0(REGION, ": A1 Marker ~ dsRNA/Immune Module Correlation by Cell Type"),
    subtitle = "Donor-level Spearman ρ (BG-subtracted, QC-residualised) | Fisher Z significance shown",
    y = "Spearman ρ (residualised)", x = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, "figures",
                 paste0(REGION, "_A1_FisherZ_paired_bars.pdf")),
       p_bar, width = 12, height = 5)

# --- 8B: Heatmap of rho_diff (astro − neuron) with Fisher Z stars ---
rho_diff_mat <- a1_df %>%
  select(Module, Marker_Gene, rho_diff) %>%
  pivot_wider(names_from = Module, values_from = rho_diff) %>%
  column_to_rownames("Marker_Gene") %>%
  as.matrix()

fisher_star_mat <- a1_df %>%
  mutate(star = ifelse(Fisher_p_adj < 0.001, "***",
                ifelse(Fisher_p_adj < 0.01,  "**",
                ifelse(Fisher_p_adj < 0.05,  "*",
                ifelse(Fisher_p_adj < 0.10,  "†", ""))))) %>%
  select(Module, Marker_Gene, star) %>%
  pivot_wider(names_from = Module, values_from = star) %>%
  column_to_rownames("Marker_Gene") %>%
  as.matrix()

if (nrow(rho_diff_mat) > 1 && ncol(rho_diff_mat) > 1) {
  p_heat <- pheatmap::pheatmap(
    rho_diff_mat,
    display_numbers = fisher_star_mat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    cluster_rows = TRUE, cluster_cols = TRUE,
    main = paste0(REGION, ": Δρ (Astro − Neuron) | A1 Markers × Modules"),
    fontsize = 10,
    angle_col = 45,
    silent = TRUE
  )
  ggsave(file.path(results_dir, "figures",
                   paste0(REGION, "_A1_FisherZ_delta_heatmap.pdf")),
         p_heat, width = 9, height = 6)
}

# --- 8C: Scatter panels for top Fisher Z hits ---
top_fz <- results_df %>%
  filter(Marker_Class == "A1", Fisher_p_adj < 0.20) %>%
  arrange(Fisher_p_adj) %>%
  head(9)

if (nrow(top_fz) == 0) {
  cat("  No Fisher Z hits at FDR<0.20; plotting top 6 by raw Fisher p.\n")
  top_fz <- results_df %>%
    filter(Marker_Class == "A1") %>%
    arrange(Fisher_p) %>%
    head(6)
}

scatter_list <- list()
for (i in seq_len(nrow(top_fz))) {
  row <- top_fz[i, ]
  mod_col  <- paste0(row$Module, "_Score_resid")
  gene_col <- paste0(row$Marker_Gene, "_resid")

  df_a <- data.frame(
    module = pb_astro[[mod_col]], marker = pb_astro[[gene_col]], cell_type = "Astrocyte"
  )
  df_n <- data.frame(
    module = pb_neuron[[mod_col]], marker = pb_neuron[[gene_col]], cell_type = "Neuron"
  )
  plot_df <- rbind(df_a, df_n) %>% filter(complete.cases(.))

  p <- ggplot(plot_df, aes(x = module, y = marker, colour = cell_type)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
    scale_colour_manual(values = c("Astrocyte" = "#E64B35", "Neuron" = "#4DBBD5")) +
    labs(
      x = paste0(row$Module, " score (resid.)"),
      y = paste0(row$Marker_Gene, " (resid.)"),
      title = paste0(row$Marker_Gene, " ~ ", row$Module),
      subtitle = sprintf("ρ_A=%.2f  ρ_N=%.2f  FisherZ FDR=%.3g",
                         row$rho_Astro, row$rho_Neuron, row$Fisher_p_adj)
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", legend.title = element_blank())

  scatter_list[[i]] <- p
}

if (length(scatter_list) > 0) {
  combined <- wrap_plots(scatter_list, ncol = min(3, length(scatter_list)))
  ggsave(file.path(results_dir, "figures",
                   paste0(REGION, "_A1_FisherZ_scatters.pdf")),
         combined,
         width  = 4.5 * min(3, length(scatter_list)),
         height = 4.5 * ceiling(length(scatter_list) / 3))
}

cat("  Figures saved.\n\n")

# =============================================================================
# SECTION 9: CONSOLE SUMMARY
# =============================================================================

cat("============================================================\n")
cat("=== RESULTS SUMMARY ===\n")
cat("============================================================\n\n")

cat("Region:", REGION, "\n")
cat("Donors: astro =", nrow(pb_astro), "| neuron =", nrow(pb_neuron), "\n")
cat("Marker genes tested:", length(markers_both), "\n")
cat("Modules:", length(module_gene_sets), "\n")
cat("Total tests:", nrow(results_df), "\n\n")

cat("--- A1 × Focus Modules (dsRNA_Sensing + Innate_Immune) ---\n")
print(
  focus_df %>%
    select(Module, Marker_Gene, rho_Astro, rho_Neuron, rho_diff,
           Fisher_Z, Fisher_p_adj) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))),
  row.names = FALSE
)

cat("\n--- All Fisher Z FDR < 0.10 (any marker class) ---\n")
sig_fz <- results_df %>%
  filter(Fisher_p_adj < 0.10) %>%
  select(Module, Marker_Gene, Marker_Class, rho_Astro, rho_Neuron,
         rho_diff, Fisher_Z, Fisher_p_adj) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
if (nrow(sig_fz) > 0) print(sig_fz, row.names = FALSE) else cat("  None.\n")

cat("\n\nDone. Output:", results_dir, "\n")
