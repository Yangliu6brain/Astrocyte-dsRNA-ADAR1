# ============================================================
# CX MERGE PIPELINE + COVARIATE INJECTION
# COMBINED SCRIPT
# UPDATED: 04/06/2026
# ============================================================
# 
# PURPOSE: 
#   1. Merge 3 CX datasets → Harmony integration → annotate
#      4 cell types → subcluster with Harmony
#   2. Save clean RDS objects for 4 cell types
#   3. Inject age, sex, PMI covariates into the RDS files
#   4. Save final annotated + covariate-enriched RDS objects
#
# WORKFLOW:
#   STEP 1: Run Sections 0-10 (merge pipeline from original script)
#   STEP 2: Run Section 11 DIAGNOSTIC to see sample_ids in your RDS
#   STEP 3: Fill in Section 12 lookup table with age/sex/pmi
#   STEP 4: Run Section 13 to inject covariates
#
# FINAL OUTPUTS: RDS files with "_withCov" suffix
#   - CX_Astrocytes_withCov.rds
#   - CX_Neurons_Total_withCov.rds
#   - CX_Excitatory_Neurons_withCov.rds
#   - CX_Inhibitory_Neurons_withCov.rds
#
# ============================================================

# =============================================================================
# SECTION 0: Setup
# =============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(harmony)
})

set.seed(42)

# =============================================================================
# SECTION 1: User Settings (MERGE PIPELINE)
# =============================================================================

library(here)

DATA_DIR <- here::here("data")   
OUT_DIR  <- here::here("data")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

IN_DIR  <- here::here("data")  # Read from data/
dataset_rds <- c(
  file.path(IN_DIR, "GSE202210_merged_step1.rds"),
  file.path(IN_DIR, "GSE178146_merged_step1.rds"),
  file.path(IN_DIR, "GSE253462_merged_step1.rds")
)

infer_condition <- function(sample_id) {
  sample_id <- as.character(sample_id)
  cond <- rep(NA_character_, length(sample_id))
  cond[grepl("^PD[0-9]|_PD|Parkinson", sample_id, ignore.case = TRUE)] <- "PD"
  cond[grepl("^C[0-9]+$|^PDC[0-9]|^HC|^NO|Control|^Ctrl|Normal", sample_id, ignore.case = TRUE)] <- "Control"
  return(cond)
}

NFEATURES <- 2000
N_PCS <- 30
HARMONY_DIMS <- 1:20
CLUSTER_RESOLUTION <- 0.5

cat("=== CX Merge Pipeline with Covariate Injection ===\n")
cat("Output:", results_dir, "\n")

# =============================================================================
# SECTION 2: Load and Prepare Datasets
# =============================================================================

cat("\n========== SECTION 2: Loading Datasets ==========\n")

for (f in dataset_rds) {
  cat("  ", basename(f), ":", ifelse(file.exists(f), "FOUND", "MISSING"), "\n")
}

objs <- lapply(dataset_rds, function(f) {
  if (!file.exists(f)) {
    cat("  WARNING: Skipping missing file:", basename(f), "\n")
    return(NULL)
  }
  cat("  Loading:", basename(f), "\n")
  readRDS(f)
})

objs <- objs[!sapply(objs, is.null)]
stopifnot(length(objs) >= 2)

for (i in seq_along(objs)) {
  o <- objs[[i]]
  md <- o@meta.data
  
  if (!"dataset" %in% colnames(md))
    md$dataset <- sub("_merged_step1\\.rds$", "", basename(dataset_rds[i]))
  
  if (!"sample_id" %in% colnames(md)) {
    md$sample_id <- if ("orig.ident" %in% colnames(md)) as.character(md$orig.ident) else "unknown"
  }
  
  if (!"condition" %in% colnames(md) || any(is.na(md$condition)) || any(md$condition == "Unknown")) {
    inferred <- infer_condition(md$sample_id)
    if ("condition" %in% colnames(md)) {
      md$condition[is.na(md$condition) | md$condition == "Unknown"] <- 
        inferred[is.na(md$condition) | md$condition == "Unknown"]
    } else {
      md$condition <- inferred
    }
  }
  
  cat("  Dataset:", unique(md$dataset), "| Cells:", nrow(md),
      "| Samples:", length(unique(md$sample_id)), "\n")
  cat("    Conditions:", paste(names(table(md$condition)), table(md$condition), collapse = ", "), "\n")
  
  o@meta.data <- md
  objs[[i]] <- o
}

# =============================================================================
# SECTION 3: Common Gene Intersection + Merge
# =============================================================================

cat("\n========== SECTION 3: Gene Intersection & Merge ==========\n")

gene_lists <- lapply(objs, rownames)
common_genes <- Reduce(intersect, gene_lists)

cat("Genes per dataset:\n")
for (i in seq_along(objs)) {
  cat("  ", unique(objs[[i]]$dataset), ":", length(rownames(objs[[i]])), "\n")
}
cat("Common genes:", length(common_genes), "\n")

# Verify key genes
key_genes <- c("ADAR", "ADARB1", "ADARB2", "SNCA", "IFIH1", "DDX58",
               "EIF2AK2", "GFAP", "AQP4", "SLC17A7", "GAD1", "MBP")
key_missing <- key_genes[!key_genes %in% common_genes]
if (length(key_missing) > 0) cat("  WARNING — key genes missing:", paste(key_missing, collapse = ", "), "\n")

# Remove Ensembl ID contamination
ensembl_ids <- grep("^ENSG[0-9]", common_genes, value = TRUE)
if (length(ensembl_ids) > 0) {
  common_genes <- common_genes[!grepl("^ENSG[0-9]", common_genes)]
  cat("  Removed", length(ensembl_ids), "Ensembl IDs ->", length(common_genes), "genes\n")
}

objs <- lapply(objs, function(o) o[common_genes, ])

combined <- merge(objs[[1]], y = objs[-1],
                  add.cell.ids = sapply(objs, function(x) unique(x$dataset)[1]))
rm(objs); gc()

cat("Total cells:", ncol(combined), "\n")
print(table(combined$dataset, combined$condition, useNA = "ifany"))

# Remove NA conditions
na_cond <- sum(is.na(combined$condition))
if (na_cond > 0) {
  cat("Removing", na_cond, "cells with NA condition\n")
  combined <- subset(combined, subset = !is.na(condition))
}

saveRDS(combined, file.path(results_dir, "00_merged_raw.rds"))

# =============================================================================
# SECTION 4: Post-Merge QC Verification (NO additional filtering)
# Step 1 already applied unified QC (500/500/10% + scDblFinder).
# This section only generates diagnostic plots.
# =============================================================================

cat("\n========== SECTION 4: QC Verification ==========\n")

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^MT-")

cat("Total cells (no additional filtering):", ncol(combined), "\n")

qc_summary <- combined@meta.data %>%
  group_by(dataset) %>%
  summarise(n_cells = n(), med_feat = median(nFeature_RNA),
            med_count = median(nCount_RNA), med_mt = median(percent.mt))
print(as.data.frame(qc_summary))

pdf(file.path(results_dir, "QC_post_merge_verification.pdf"), width = 15, height = 10)
p1 <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "dataset", pt.size = 0, ncol = 3)
p2 <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              group.by = "condition", pt.size = 0, ncol = 3)
print(p1 / p2)
dev.off()

# =============================================================================
# SECTION 5: Data Integration with Harmony
# =============================================================================

cat("\n========== SECTION 5: Harmony Integration ==========\n")

combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, nfeatures = NFEATURES, verbose = FALSE)
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = N_PCS, verbose = FALSE)
combined <- RunHarmony(combined, group.by.vars = "dataset", dims.use = HARMONY_DIMS)

cat("Harmony integration complete\n")

combined <- RunUMAP(combined, reduction = "harmony", dims = HARMONY_DIMS, verbose = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = HARMONY_DIMS, verbose = FALSE)
combined <- FindClusters(combined, resolution = CLUSTER_RESOLUTION, verbose = FALSE)

pdf(file.path(results_dir, "UMAP_integration.pdf"), width = 15, height = 5)
p1 <- DimPlot(combined, reduction = "umap", group.by = "dataset", label = FALSE)
p2 <- DimPlot(combined, reduction = "umap", group.by = "condition", label = FALSE)
p3 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
print(p1 + p2 + p3)
dev.off()

# Batch correction assessment
batch_summary <- combined@meta.data %>%
  group_by(seurat_clusters, dataset) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from = dataset, values_from = n_cells, values_fill = 0)

write.csv(as.data.frame(batch_summary),
          file.path(results_dir, "Batch_mixing_per_cluster.csv"), row.names = FALSE)

saveRDS(combined, file.path(results_dir, "01_integrated.rds"))

# =============================================================================
# SECTION 6: Cell Type Annotation
# =============================================================================

cat("\n========== SECTION 6: Cell Type Annotation ==========\n")

# Define cell type markers
cell_markers <- list(
  Astrocyte = c("GFAP", "AQP4", "SLC1A2", "SOX9", "ALDH1L1"),
  Excitatory_Neuron = c("SLC17A7", "CAMK2A", "VGLUT1", "CaMKII"),
  Inhibitory_Neuron = c("GAD1", "GAD2", "SLC32A1", "GABA"),
  Oligodendrocyte = c("MBP", "MOBP", "MOG", "PLP1"),
  OPC = c("PDGFRA", "CSPG4", "NG2"),
  Microglia = c("CX3CR1", "IBA1", "TMEM119", "P2RY12"),
  Endothelial = c("PECAM1", "VWF", "CDH5")
)

# Calculate marker-based scores (mean log-normalised expression)
for (marker_set in names(cell_markers)) {
  markers <- cell_markers[[marker_set]]
  markers_present <- markers[markers %in% rownames(combined)]
  if (length(markers_present) > 0) {
    combined[[paste0(marker_set, "_Score1")]] <- 
      colMeans(combined@assays$RNA@data[markers_present, , drop = FALSE])
  }
}

score_cols <- colnames(combined@meta.data)[grepl("_Score1$", colnames(combined@meta.data))]
all_genes <- rownames(combined)

scores <- as.matrix(combined@meta.data[, score_cols])
best_idx <- apply(scores, 1, which.max)
best_score <- apply(scores, 1, function(x) sort(x, decreasing = TRUE)[1])
second_score <- apply(scores, 1, function(x) sort(x, decreasing = TRUE)[2])

best_type <- sub("_Score1$", "", score_cols[best_idx])
confidence_pass <- (best_score > 0.1) & (best_score > 1.5 * second_score)

combined$cell_type <- best_type
combined$cell_type_confidence <- ifelse(confidence_pass, "High", "Low")

# Broad cell type labels
combined$cell_type_broad <- case_when(
  combined$cell_type %in% c("Excitatory_Neuron", "Inhibitory_Neuron") ~ "Neuron",
  combined$cell_type == "Astrocyte" ~ "Astrocyte",
  combined$cell_type == "Oligodendrocyte" ~ "Oligodendrocyte",
  combined$cell_type == "OPC" ~ "OPC",
  combined$cell_type == "Microglia" ~ "Microglia",
  combined$cell_type == "Endothelial" ~ "Endothelial",
  TRUE ~ "Other"
)

cat("\nAnnotation confidence:\n")
print(table(combined$cell_type_confidence))
cat("\nCell type (detailed):\n")
print(table(combined$cell_type))
cat("\nCell type (broad):\n")
print(table(combined$cell_type_broad))
cat("\nCell type by condition:\n")
print(table(combined$cell_type_broad, combined$condition))

# Supplementary table
cell_counts <- table(combined$dataset, combined$cell_type_broad)
write.csv(as.data.frame.matrix(cell_counts),
          file.path(results_dir, "Cell_counts_dataset_x_celltype.csv"))

saveRDS(combined, file.path(results_dir, "02_annotated.rds"))

# Annotation plots
pdf(file.path(results_dir, "UMAP_celltypes.pdf"), width = 14, height = 12)
p1 <- DimPlot(combined, group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend()
p2 <- DimPlot(combined, group.by = "cell_type_broad", label = TRUE)
p3 <- DimPlot(combined, group.by = "cell_type_confidence") + ggtitle("Annotation Confidence")
print(p1 + p2 + p3)
dev.off()

# Marker DotPlot validation
marker_genes <- unique(unlist(cell_markers))
marker_genes <- marker_genes[marker_genes %in% all_genes]

pdf(file.path(results_dir, "DotPlot_markers_validation.pdf"), width = 16, height = 8)
print(DotPlot(combined, features = marker_genes[1:min(35, length(marker_genes))],
              group.by = "cell_type_broad") + RotatedAxis())
dev.off()

# =============================================================================
# SECTION 7: Subset Cell Types
# =============================================================================

cat("\n========== SECTION 7: Subsetting Cell Types ==========\n")

astro <- subset(combined, subset = cell_type_broad == "Astrocyte")
cat("Astrocytes:", ncol(astro),
    "| PD:", sum(astro$condition == "PD"),
    "| Control:", sum(astro$condition == "Control"), "\n")

neurons <- subset(combined, subset = cell_type_broad == "Neuron")
cat("Neurons (total):", ncol(neurons),
    "| PD:", sum(neurons$condition == "PD"),
    "| Control:", sum(neurons$condition == "Control"), "\n")

exc_neurons <- subset(combined, subset = cell_type == "Excitatory_Neuron")
cat("Excitatory Neurons:", ncol(exc_neurons),
    "| PD:", sum(exc_neurons$condition == "PD"),
    "| Control:", sum(exc_neurons$condition == "Control"), "\n")

inh_neurons <- subset(combined, subset = cell_type == "Inhibitory_Neuron")
cat("Inhibitory Neurons:", ncol(inh_neurons),
    "| PD:", sum(inh_neurons$condition == "PD"),
    "| Control:", sum(inh_neurons$condition == "Control"), "\n")

# =============================================================================
# SECTION 8: Astrocyte Subclustering (with Harmony)
# =============================================================================

cat("\n========== SECTION 8: Astrocyte Subclustering ==========\n")

astro <- NormalizeData(astro, verbose = FALSE)
astro <- FindVariableFeatures(astro, nfeatures = 2000, verbose = FALSE)
astro <- ScaleData(astro, features = VariableFeatures(astro), verbose = FALSE)
astro <- RunPCA(astro, npcs = 20, verbose = FALSE)
astro <- RunHarmony(astro, group.by.vars = "dataset", dims.use = 1:15)
astro <- RunUMAP(astro, reduction = "harmony", dims = 1:15, verbose = FALSE)
astro <- FindNeighbors(astro, reduction = "harmony", dims = 1:15, verbose = FALSE)
astro <- FindClusters(astro, resolution = 0.3, verbose = FALSE)

cat("Astrocyte clusters:\n")
print(table(astro$seurat_clusters, astro$condition))

# Purity check
purity_markers <- c("RBFOX3", "SYT1", "SNAP25", "SLC17A7", "GAD1",
                    "AQP4", "GFAP", "SLC1A2", "ALDH1L1")
purity_markers <- purity_markers[purity_markers %in% rownames(astro)]

pdf(file.path(results_dir, "Astro_purity_check.pdf"), width = 14, height = 8)
print(VlnPlot(astro, features = purity_markers, group.by = "seurat_clusters",
              pt.size = 0, ncol = 3))
dev.off()

pdf(file.path(results_dir, "Astro_subclusters.pdf"), width = 15, height = 5)
p1 <- DimPlot(astro, group.by = "seurat_clusters", label = TRUE) + ggtitle("Subclusters")
p2 <- DimPlot(astro, group.by = "condition") + ggtitle("Condition")
p3 <- DimPlot(astro, group.by = "dataset") + ggtitle("Dataset")
print(p1 + p2 + p3)
dev.off()

cat("Check Astro_purity_check.pdf — clusters with high RBFOX3/SYT1 = neuronal contamination\n")

# =============================================================================
# SECTION 9: Neuron Subclustering (with Harmony)
# =============================================================================

cat("\n========== SECTION 9: Neuron Subclustering ==========\n")

neurons <- NormalizeData(neurons, verbose = FALSE)
neurons <- FindVariableFeatures(neurons, nfeatures = 2000, verbose = FALSE)
neurons <- ScaleData(neurons, features = VariableFeatures(neurons), verbose = FALSE)
neurons <- RunPCA(neurons, npcs = 20, verbose = FALSE)
neurons <- RunHarmony(neurons, group.by.vars = "dataset", dims.use = 1:15)
neurons <- RunUMAP(neurons, reduction = "harmony", dims = 1:15, verbose = FALSE)
neurons <- FindNeighbors(neurons, reduction = "harmony", dims = 1:15, verbose = FALSE)
neurons <- FindClusters(neurons, resolution = 0.5, verbose = FALSE)

cat("Neuron clusters:\n")
print(table(neurons$seurat_clusters, neurons$condition))

# Purity check
neuron_purity <- c("RBFOX3", "SYT1", "SNAP25", "MAP2",
                   "GFAP", "AQP4", "MBP", "PLP1", "PDGFRA", "CX3CR1")
neuron_purity <- neuron_purity[neuron_purity %in% rownames(neurons)]

pdf(file.path(results_dir, "Neuron_purity_check.pdf"), width = 14, height = 8)
print(VlnPlot(neurons, features = neuron_purity, group.by = "seurat_clusters",
              pt.size = 0, ncol = 3))
dev.off()

pdf(file.path(results_dir, "Neuron_subclusters.pdf"), width = 15, height = 10)
p1 <- DimPlot(neurons, group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
p2 <- DimPlot(neurons, group.by = "condition") + ggtitle("Condition")
p3 <- DimPlot(neurons, group.by = "cell_type") + ggtitle("Exc vs Inh")
p4 <- DimPlot(neurons, group.by = "dataset") + ggtitle("Dataset")
print((p1 + p2) / (p3 + p4))
dev.off()

cat("Check Neuron_purity_check.pdf — clusters with high GFAP/MBP = glial contamination\n")

# =============================================================================
# SECTION 10: Save All Objects (BEFORE covariate injection)
# =============================================================================

cat("\n========== SECTION 10: Saving Objects ==========\n")

saveRDS(astrocytes_obj, file.path(OUT_DIR, "CX_Astrocytes_withCov.rds"))
saveRDS(neurons_obj,    file.path(OUT_DIR, "CX_Neurons_Total_withCov.rds"))

saveRDS(exc_neurons, file.path(results_dir, "CX_Excitatory_Neurons.rds"))
cat("Saved: CX_Excitatory_Neurons.rds (", ncol(exc_neurons), "cells)\n")

saveRDS(inh_neurons, file.path(results_dir, "CX_Inhibitory_Neurons.rds"))
cat("Saved: CX_Inhibitory_Neurons.rds (", ncol(inh_neurons), "cells)\n")

saveRDS(combined, file.path(results_dir, "CX_All_Celltypes_Annotated.rds"))
cat("Saved: CX_All_Celltypes_Annotated.rds (", ncol(combined), "cells)\n")

# Session info
writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info.txt"))

cat("\n============================================================\n")
cat("=== MERGE PIPELINE COMPLETE ===\n")
cat("============================================================\n\n")

# =============================================================================
# SECTION 11: DIAGNOSTIC — Extract sample_ids from saved RDS
# =============================================================================
# Run this section to see the exact sample_id values in your merged RDS.
# Copy the output and fill in Section 12 lookup table below.

cat("\n========== SECTION 11: DIAGNOSTIC (Covariate Setup) ==========\n")

cat(">>> DIAGNOSTIC: Reading sample_id values from CX_Astrocytes.rds\n\n")

astro_diag <- readRDS(file.path(results_dir, "CX_Astrocytes.rds"))

cat("--- Metadata columns available ---\n")
cat(paste(colnames(astro_diag@meta.data), collapse = ", "), "\n\n")

# Find the donor ID column
donor_col_found <- NULL
for (nm in c("sample_id", "donor_id", "Donor_ID", "orig.ident", "Individual")) {
  if (nm %in% colnames(astro_diag@meta.data)) {
    donor_col_found <- nm
    break
  }
}
cat("Donor ID column detected:", donor_col_found, "\n\n")

# Print all unique sample_ids grouped by condition and dataset
diag_summary <- astro_diag@meta.data %>%
  group_by(dataset, condition, .data[[donor_col_found]]) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(dataset, condition, .data[[donor_col_found]])

cat("=== ALL DONOR IDs IN YOUR MERGED RDS (copy these for Section 12) ===\n")
print(as.data.frame(diag_summary), row.names = FALSE)

cat("\n=== EXACT sample_id values to paste into Section 12 ===\n")
cat('sample_id = c(\n')
by_dataset <- split(diag_summary, diag_summary$dataset)
for (ds in names(by_dataset)) {
  ids <- by_dataset[[ds]][[donor_col_found]]
  cat(sprintf('  # --- %s ---\n', ds))
  cat(sprintf('  "%s",\n', paste(ids, collapse = '", "')))
}
cat(')\n\n')

cat(">>> STOP HERE TEMPORARILY.\n")
cat("    Fill in Section 12 lookup table using the sample_id values above,\n")
cat("    then run Section 13 to inject covariates.\n\n")

rm(astro_diag); gc()

# =============================================================================
# SECTION 12: DONOR LOOKUP TABLE — FILL IN AFTER RUNNING SECTION 11
# =============================================================================
# Instructions:
#   - Copy the sample_id values printed above into the sample_id column
#   - Fill in age, sex, pmi from your GSE_list_Demographic.xlsx file
#   - Use exact same strings as printed (case-sensitive!)
#   - sex: use "Male" or "Female"
#   - pmi: in HOURS. If not available for a dataset, use NA
#   - age: in YEARS. For age ranges (GSE253462), use midpoint
#           e.g. "70-75" → 72.5, "<65" → 62, "80-85" → 82.5

cat("============================================================\n")
cat("SECTION 12: Donor lookup table\n")
cat("============================================================\n\n")

# ============================================================
# >>> EDIT THIS TABLE — paste your sample_ids from Section 11 output
# ============================================================

donor_metadata <- data.frame(
  
  sample_id = c(
    # ---- PASTE YOUR SAMPLE IDs HERE (from Section 11 output) ----
    # Example format — replace with your actual values:
    
    # --- GSE178146 ---
    "C36", "C48", "PDC05", "PDC22", "PDC34", "PDC87", "PDC91", "PD413", "PD416", "PD523", "PD666", "PD683", "PD732", "PD747",
    # --- GSE202210 ---
    "HSDG07HC", "HSDG10HC", "hsDG101HC", "hsDG13HC", "hsDG30HC", "hsDG99HC", "HSDG148PD", "HSDG199PD", "hsDG142PD", "hsDG151PD", "hsDG197PD", "hsDG208PD",
    # --- GSE253462 ---
    "T1619", "T1952", "T2017", "T2465", "T343", "T4174", "T4210", "T551", "T1976", "T233", "T257", "T290", "T308", "T362", "T3860", "T4162", "T777"
    
    # Note: The exact IDs above are PLACEHOLDERS.
    # Run Section 11 first to get the real values from your RDS.
  ),
  
  age = c(
    # GSE202210 Controls (from your demographic sheet)
    81, 70, 72, 65, 63, 69,
    # GSE202210 PD
    96, 86, 79, 81, 72, 77, 
    
    # GSE178146 Controls
    68, 68, 58, 65, 90, 92, 85,
    # GSE178146 PD
    72, 85, 86, 93, 83, 81, 83,
    
    # GSE253462 Controls (midpoints: 75-80→77.5, 80-85→82.5, <65→62, 70-75→72.5)
    82.5, 82.5, 72.5, 82.5, 62, 82.5, 72.5, 72.5,
    # GSE253462 PD
    62, 72.5, 67.5, 72.5, 77.5, 72.5, 72.5, 72.5, 77.5
  ),
  
  sex = c(
    # GSE202210 Controls
    "Female", "Male", "Female", "Female", "Male", "Male",
    # GSE202210 PD
    "Male", "Male", "Male", "Female", "Female","Female", 
    
    # GSE178146 Controls
    "Male", "Male", "Male", "Male", "Male", "Female", "Male",
    # GSE178146 PD
    "Female", "Female", "Female", "Female", "Female", "Male", "Male",
    
    # GSE253462 Controls
    "Female", "Male", "Male", "Female", "Male", "Male", "Male", "Female",
    # GSE253462 PD
    "Male", "Male", "Female", "Male", "Male", "Female", "Male", "Male", "Female"
  ),
  
  pmi = c(
    # GSE202210 Controls (hours)
    8, 10.5, 7, 19, 5, 11,
    # GSE202210 PD
    6, 13, 5.5, 8, 6, 7,
    
    # GSE178146 Controls
    30, 10, 9, 12, 12, 24, 29,
    # GSE178146 PD
    14, 17, 9, 24, 24, 14, 22,
    
    # GSE253462: PMI NOT AVAILABLE — all NA
    NA, NA, NA, NA, NA, NA, NA, NA,   # Controls
    NA, NA, NA, NA, NA, NA, NA, NA, NA # PD
  ),
  
  stringsAsFactors = FALSE
)

# Validate table
cat(sprintf("Lookup table: %d donors defined\n", nrow(donor_metadata)))
if (any(duplicated(donor_metadata$sample_id))) {
  stop("DUPLICATE sample_ids found: ",
       paste(donor_metadata$sample_id[duplicated(donor_metadata$sample_id)], collapse = ", "))
}
cat("No duplicate sample_ids — OK\n\n")

# =============================================================================
# SECTION 13: INJECT COVARIATES INTO SAVED RDS FILES
# =============================================================================

cat("============================================================\n")
cat("SECTION 13: Injecting covariates into RDS files\n")
cat("============================================================\n\n")

# Helper function -------------------------------------------------------
inject_covariates <- function(rds_path, donor_meta, output_dir) {
  
  fname <- basename(rds_path)
  cat(sprintf("--- %s ---\n", fname))
  
  obj <- readRDS(rds_path)
  meta <- obj@meta.data
  
  # Find donor column
  donor_col <- NULL
  for (nm in c("sample_id", "donor_id", "Donor_ID", "orig.ident", "Individual")) {
    if (nm %in% colnames(meta)) { donor_col <- nm; break }
  }
  if (is.null(donor_col)) stop("No donor column found in: ", fname)
  cat(sprintf("  Donor column: '%s' (%d unique donors)\n",
              donor_col, length(unique(meta[[donor_col]]))))
  
  # Check match rate
  donors_in_rds   <- unique(as.character(meta[[donor_col]]))
  donors_in_table <- donor_meta$sample_id
  matched   <- intersect(donors_in_rds, donors_in_table)
  unmatched <- setdiff(donors_in_rds, donors_in_table)
  
  cat(sprintf("  Matched: %d / %d donors\n", length(matched), length(donors_in_rds)))
  
  if (length(unmatched) > 0) {
    cat(sprintf("  WARNING: %d donors NOT in lookup table:\n", length(unmatched)))
    cat(sprintf("    %s\n", paste(unmatched, collapse = ", ")))
    cat("  These cells will get NA. Fix sample_id values in Section 12.\n")
  }
  
  if (length(matched) == 0) {
    cat("  ERROR: Zero donors matched. Skipping this file.\n")
    cat(sprintf("  RDS has: %s\n", paste(head(donors_in_rds, 5), collapse=", ")))
    cat(sprintf("  Table has: %s\n", paste(head(donors_in_table, 5), collapse=", ")))
    return(invisible(NULL))
  }
  
  # Build named lookup vectors
  age_lut <- setNames(donor_meta$age, donor_meta$sample_id)
  sex_lut <- setNames(donor_meta$sex, donor_meta$sample_id)
  pmi_lut <- setNames(donor_meta$pmi, donor_meta$sample_id)
  
  cell_ids <- as.character(meta[[donor_col]])
  
  obj$age <- as.numeric(age_lut[cell_ids])
  obj$sex <- as.character(sex_lut[cell_ids])
  obj$pmi <- as.numeric(pmi_lut[cell_ids])
  
  # Report coverage
  n <- ncol(obj)
  cat(sprintf("  age: %d/%d cells assigned (%.0f%%)\n",
              sum(!is.na(obj$age)), n, 100*mean(!is.na(obj$age))))
  cat(sprintf("  sex: %d/%d cells assigned (%.0f%%)\n",
              sum(!is.na(obj$sex)), n, 100*mean(!is.na(obj$sex))))
  cat(sprintf("  pmi: %d/%d cells assigned (%.0f%%)\n",
              sum(!is.na(obj$pmi)), n, 100*mean(!is.na(obj$pmi))))
  
  # Show donor-level summary to verify
  cat("\n  Donor-level verification (first 10 rows):\n")
  donor_check <- obj@meta.data %>%
    group_by(.data[[donor_col]], condition) %>%
    summarise(age = first(age), sex = first(sex), pmi = first(pmi),
              n_cells = n(), .groups = "drop") %>%
    arrange(condition, .data[[donor_col]])
  print(head(as.data.frame(donor_check), 10))
  
  # Save with "_withCov" suffix
  out_name <- sub("\\.rds$", "_withCov.rds", fname)
  out_path <- file.path(output_dir, out_name)
  saveRDS(obj, out_path)
  cat(sprintf("\n  Saved: %s\n\n", out_name))
  
  rm(obj); gc()
  return(out_path)
}

# -----------------------------------------------------------------------
# Run on all CX cell type RDS files
# -----------------------------------------------------------------------

# Inject covariates into the 4 cell-type files
cat("Injecting covariates into CX cell-type files...\n\n")

inject_covariates(file.path(results_dir, "CX_Astrocytes.rds"), 
                  donor_metadata, results_dir)
inject_covariates(file.path(results_dir, "CX_Neurons_Total.rds"), 
                  donor_metadata, results_dir)
inject_covariates(file.path(results_dir, "CX_Excitatory_Neurons.rds"), 
                  donor_metadata, results_dir)
inject_covariates(file.path(results_dir, "CX_Inhibitory_Neurons.rds"), 
                  donor_metadata, results_dir)

# Also save the full annotated object with covariates
inject_covariates(file.path(results_dir, "CX_All_Celltypes_Annotated.rds"), 
                  donor_metadata, results_dir)

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n============================================================\n")
cat("=== FULL PIPELINE COMPLETE ===\n")
cat("=== Merge + Covariate Injection ===\n")
cat("============================================================\n\n")

cat("Output directory:", results_dir, "\n\n")

cat("FINAL RDS FILES (with covariates injected):\n")
cat("  CX_Astrocytes_withCov.rds\n")
cat("  CX_Neurons_Total_withCov.rds\n")
cat("  CX_Excitatory_Neurons_withCov.rds\n")
cat("  CX_Inhibitory_Neurons_withCov.rds\n")
cat("  CX_All_Celltypes_Annotated_withCov.rds\n\n")

cat("All files now contain metadata columns:\n")
cat("  - age (years)\n")
cat("  - sex (Male/Female)\n")
cat("  - pmi (hours)\n")
cat("  - (plus original columns: sample_id, condition, dataset, cell_type, etc.)\n\n")

cat("Ready for downstream analysis with covariates!\n")
cat("Done:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
