# ============================================================
# STEP 1 — GSE178146 raw 10x (Cingulate Cortex - folder structure) -> Seurat
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scDblFinder)
  library(SingleCellExperiment)
})

set.seed(42)

# ---------------------------
# 1) USER SETTINGS
# ---------------------------
DATASET_ID <- "GSE178146"

# Folder that contains the sample subfolders like:
# GSM5380921_C36_raw_feature_bc_matrix/
# GSM5380930_PD413_raw_feature_bc_matrix/
# Each subfolder contains: barcodes.tsv(.gz), features.tsv(.gz), matrix.mtx(.gz)
library(here)
DATA_DIR <- here::here("data", "raw", "GSE178146")   # GEO files go here
OUT_DIR  <- here::here("data")          # QC'd output
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
# <-- CHANGE THIS


# Clinical file (optional)
CLINICAL_FILE <- NULL

# --- Unified QC thresholds (consistent across all 6 datasets) ---
# CHANGED: Was MIN_FEATURES=100, MIN_COUNTS=200, UPPER_Q=0.995
# This was the most permissive dataset, creating systematic bias at merge.
# If GSE178146 uses 10x v2 chemistry with genuinely lower sensitivity,
# document this justification in your methods section.
# However, unified thresholds are preferred for cross-dataset integration.
MIN_FEATURES <- 500    # Unified (was 100)
MIN_COUNTS   <- 500    # Unified (was 200)
MAX_MT       <- 10
MITO_PATTERN <- "^MT-"

# Upper cutoff strategy
UPPER_Q <- 0.99        # Unified (was 0.995)

# --- Gene removal settings ---
REMOVE_MALAT1 <- TRUE
REMOVE_RIBO   <- TRUE

# --- Doublet detection ---
RUN_DOUBLET_DETECTION <- TRUE

# Dim reduction for QC visualization only
N_PCS_QC     <- 20
UMAP_DIMS_QC <- 1:20

SAVE_PLOTS <- TRUE


# ---------------------------
# 2) HELPER FUNCTIONS
# ---------------------------

find_first_existing <- function(paths) {
  p <- paths[file.exists(paths)]
  if (length(p) == 0) return(NA_character_)
  p[[1]]
}

# Read 10x mtx from a folder containing barcodes.tsv, features.tsv, matrix.mtx
# Handles nested folder structure: GSM###_XXX_raw_feature_bc_matrix/XXX_raw_feature_bc_matrix/
read_10x_folder <- function(folder_path) {
  
  subfolders <- list.dirs(folder_path, full.names = TRUE, recursive = FALSE)
  nested_folder <- subfolders[grepl("_raw_feature_bc_matrix$", basename(subfolders))]
  
  if (length(nested_folder) > 0) {
    data_folder <- nested_folder[1]
    message("  Found nested folder: ", basename(data_folder))
  } else {
    data_folder <- folder_path
  }
  
  mtx_path <- find_first_existing(c(
    file.path(data_folder, "matrix.mtx.gz"),
    file.path(data_folder, "matrix.mtx")
  ))
  features_path <- find_first_existing(c(
    file.path(data_folder, "features.tsv.gz"),
    file.path(data_folder, "features.tsv"),
    file.path(data_folder, "genes.tsv.gz"),
    file.path(data_folder, "genes.tsv")
  ))
  barcodes_path <- find_first_existing(c(
    file.path(data_folder, "barcodes.tsv.gz"),
    file.path(data_folder, "barcodes.tsv")
  ))
  
  if (is.na(mtx_path) || is.na(features_path) || is.na(barcodes_path)) {
    stop("Missing required files in folder: ", data_folder, 
         "\n  Looked for matrix.mtx(.gz), features.tsv(.gz), barcodes.tsv(.gz)")
  }
  
  m <- Matrix::readMM(mtx_path)
  feats <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)
  bars  <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
  
  gene_names <- if (ncol(feats) >= 2) feats[[2]] else feats[[1]]
  gene_names <- make.unique(as.character(gene_names))
  
  colnames(m) <- as.character(bars[[1]])
  rownames(m) <- gene_names
  
  m <- as(m, "dgCMatrix")
  return(m)
}

# Extract sample_id and GSM from folder names
extract_ids_from_folder <- function(folder_name) {
  m <- stringr::str_match(folder_name, "^(GSM[0-9]+)_([A-Za-z0-9]+)_raw_feature_bc_matrix$")
  if (is.na(m[1, 2])) return(list(gsm = NA_character_, sample_id = NA_character_))
  return(list(gsm = m[1, 2], sample_id = m[1, 3]))
}

# Infer condition from sample ID
# PD### -> "PD", C## -> "Control", PDC## -> "Control" (VERIFY!)
infer_condition <- function(sample_id) {
  if (is.na(sample_id)) return(NA_character_)
  
  if (grepl("^PDC[0-9]+$", sample_id)) {
    # PDC samples — assumed Control based on Al-Dalahmah et al. naming
    # IMPORTANT: Verify with Supplementary Table 1!
    return("Control")
  } else if (grepl("^PD[0-9]+$", sample_id)) {
    return("PD")
  } else if (grepl("^C[0-9]+$", sample_id)) {
    return("Control")
  } else {
    return(NA_character_)
  }
}


# ---------------------------
# 3) LOAD CLINICAL METADATA (if available)
# ---------------------------
condition_lookup <- NULL

if (!is.null(CLINICAL_FILE) && file.exists(CLINICAL_FILE)) {
  message("Loading clinical metadata from: ", CLINICAL_FILE)
  
  clinical_df <- read.csv(CLINICAL_FILE)
  
  possible_id_cols <- c("Sample.ID", "sample_id", "SampleID", "Sample", "GSM")
  possible_cond_cols <- c("Clinical.diagnosis", "diagnosis", "Condition", "condition", "Disease")
  
  id_col <- intersect(possible_id_cols, colnames(clinical_df))[1]
  cond_col <- intersect(possible_cond_cols, colnames(clinical_df))[1]
  
  if (!is.na(id_col) && !is.na(cond_col)) {
    clinical_df$sample_id <- trimws(as.character(clinical_df[[id_col]]))
    clinical_df$diagnosis <- trimws(as.character(clinical_df[[cond_col]]))
    
    clinical_df$condition <- case_when(
      grepl("parkinson", clinical_df$diagnosis, ignore.case = TRUE) ~ "PD",
      grepl("control", clinical_df$diagnosis, ignore.case = TRUE) ~ "Control",
      grepl("^PD$", clinical_df$diagnosis) ~ "PD",
      TRUE ~ clinical_df$diagnosis
    )
    
    condition_lookup <- setNames(clinical_df$condition, clinical_df$sample_id)
    message("Loaded condition mapping for ", length(condition_lookup), " samples")
  } else {
    message("Could not find expected columns in clinical file. Will infer conditions from sample names.")
  }
} else {
  message("No clinical file provided. Conditions will be inferred from sample naming (PD/C/PDC prefixes).")
}


# ---------------------------
# 4) DISCOVER SAMPLE FOLDERS
# ---------------------------
all_folders <- list.dirs(DATA_DIR, full.names = TRUE, recursive = FALSE)
sample_folders <- all_folders[grepl("_raw_feature_bc_matrix$", basename(all_folders))]

if (length(sample_folders) == 0) {
  stop("No sample folders found matching '*_raw_feature_bc_matrix' pattern in DATA_DIR: ", DATA_DIR)
}

message("Found ", length(sample_folders), " sample folders.")


# ---------------------------
# 5) PROCESS EACH SAMPLE
# ---------------------------
message("\n========== PROCESSING SAMPLES ==========")

seurat_list <- list()
qc_summaries <- list()
load_errors <- character()

for (folder in sample_folders) {
  folder_name <- basename(folder)
  ids <- extract_ids_from_folder(folder_name)
  sid <- ids$sample_id
  gsm <- ids$gsm
  
  if (is.na(sid)) {
    message("Skipping (could not parse sample id): ", folder_name)
    load_errors <- c(load_errors, folder_name)
    next
  }
  
  # Get condition
  if (!is.null(condition_lookup) && sid %in% names(condition_lookup)) {
    condition <- unname(condition_lookup[sid])
  } else {
    condition <- infer_condition(sid)
  }
  
  message("\n========== ", sid, " (", gsm, ", ", condition, ") ==========")
  
  # -------------------------------------------------------
  # 5a) Load count matrix
  # -------------------------------------------------------
  counts <- tryCatch({
    read_10x_folder(folder)
  }, error = function(e) {
    message("  ERROR loading sample: ", e$message)
    load_errors <<- c(load_errors, folder_name)
    return(NULL)
  })
  
  if (is.null(counts)) next
  
  message("  Raw matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")
  
  # -------------------------------------------------------
  # 5b) Remove MALAT1 and ribosomal genes BEFORE Seurat object creation
  # -------------------------------------------------------
  genes_to_remove <- character(0)
  
  if (REMOVE_MALAT1) {
    malat1 <- grep("^MALAT1$", rownames(counts), value = TRUE)
    genes_to_remove <- c(genes_to_remove, malat1)
    if (length(malat1) > 0) message("  Removing MALAT1")
  }
  
  if (REMOVE_RIBO) {
    ribo_genes <- grep("^RP[SL][0-9]|^MT-RP", rownames(counts), value = TRUE)
    genes_to_remove <- c(genes_to_remove, ribo_genes)
    if (length(ribo_genes) > 0) message("  Removing ", length(ribo_genes), " ribosomal genes")
  }
  
  if (length(genes_to_remove) > 0) {
    counts <- counts[!rownames(counts) %in% genes_to_remove, ]
    message("  After gene removal: ", nrow(counts), " genes")
  }
  
  # -------------------------------------------------------
  # 5c) Create Seurat object
  # -------------------------------------------------------
  obj <- CreateSeuratObject(
    counts = counts,
    project = paste0(DATASET_ID, "_", sid),
    min.cells = 3,
    min.features = 1  # Drop completely empty barcodes early (raw 10x has 737K barcodes)
  )
  
  obj$dataset   <- DATASET_ID
  obj$gsm       <- gsm
  obj$sample_id <- sid
  obj$condition <- condition
  
  obj <- RenameCells(obj, add.cell.id = paste0(DATASET_ID, "_", sid))
  
  # Calculate percent.mt
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = MITO_PATTERN)
  
  # Handle NA values in percent.mt (can occur with empty/near-empty droplets)
  if (any(is.na(obj$percent.mt))) {
    n_na <- sum(is.na(obj$percent.mt))
    message("  NOTE: ", n_na, " cells have NA percent.mt — setting to 100 (will be filtered)")
    obj$percent.mt[is.na(obj$percent.mt)] <- 100
  }
  
  if (max(obj$percent.mt, na.rm = TRUE) == 0) {
    message("  NOTE: No MT genes found with '^MT-', trying alternatives...")
    for (pattern in c("^mt-", "^Mt-")) {
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = pattern)
      if (max(obj$percent.mt, na.rm = TRUE) > 0) {
        message("  Found MT genes with pattern: ", pattern)
        break
      }
    }
  }
  
  # -------------------------------------------------------
  # 5d) QC filtering with unified thresholds
  # -------------------------------------------------------
  MAX_COUNTS_Q <- as.numeric(quantile(obj$nCount_RNA, probs = UPPER_Q))
  MAX_FEATURES_Q <- as.numeric(quantile(obj$nFeature_RNA, probs = UPPER_Q))
  
  message("  QC cutoffs:")
  message("    MIN_FEATURES = ", MIN_FEATURES, " (unified; was 100)")
  message("    MIN_COUNTS = ", MIN_COUNTS, " (unified; was 200)")
  message("    MAX_MT = ", MAX_MT, "%")
  message("    MAX_COUNTS (", UPPER_Q*100, "th %ile) = ", round(MAX_COUNTS_Q))
  message("    MAX_FEATURES (", UPPER_Q*100, "th %ile) = ", round(MAX_FEATURES_Q))
  
  pre_qc_n <- ncol(obj)
  
  if (MAX_FEATURES_Q < MIN_FEATURES || MAX_COUNTS_Q < MIN_COUNTS) {
    message("  WARNING: Upper quantile cutoffs are below minimum thresholds!")
    message("  This may indicate very low-quality data for this sample.")
    message("  Applying only lower thresholds...")
    obj <- subset(
      obj,
      subset = nFeature_RNA >= MIN_FEATURES &
        nCount_RNA >= MIN_COUNTS &
        percent.mt <= MAX_MT
    )
  } else {
    obj <- subset(
      obj,
      subset = nFeature_RNA >= MIN_FEATURES &
        nCount_RNA >= MIN_COUNTS &
        percent.mt <= MAX_MT &
        nCount_RNA <= MAX_COUNTS_Q &
        nFeature_RNA <= MAX_FEATURES_Q
    )
  }
  
  post_qc_n <- ncol(obj)
  
  if (post_qc_n == 0) {
    message("  WARNING: No cells passed QC! Skipping this sample.")
    message("  Consider checking if this dataset has lower sensitivity (10x v2?)")
    next
  }
  
  message("  After QC: ", pre_qc_n, " -> ", post_qc_n, " (", round(post_qc_n/pre_qc_n*100, 1), "% retained)")
  
  # -------------------------------------------------------
  # 5e) Doublet detection with scDblFinder
  # -------------------------------------------------------
  n_doublets <- 0
  
  if (RUN_DOUBLET_DETECTION && post_qc_n >= 100) {
    message("  Running scDblFinder...")
    
    tryCatch({
      obj_tmp <- NormalizeData(obj, verbose = FALSE)
      obj_tmp <- FindVariableFeatures(obj_tmp, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      obj_tmp <- ScaleData(obj_tmp, verbose = FALSE)
      obj_tmp <- RunPCA(obj_tmp, npcs = N_PCS_QC, verbose = FALSE)
      
      sce <- as.SingleCellExperiment(obj_tmp)
      sce <- scDblFinder(sce)
      
      obj$scDblFinder.class <- sce$scDblFinder.class
      obj$scDblFinder.score <- sce$scDblFinder.score
      
      n_doublets <- sum(obj$scDblFinder.class == "doublet")
      message("  scDblFinder: ", n_doublets, " doublets detected (",
              round(n_doublets/ncol(obj)*100, 1), "%)")
      
      obj <- subset(obj, subset = scDblFinder.class == "singlet")
      message("  After doublet removal: ", ncol(obj), " cells")
      
      rm(obj_tmp, sce)
      gc(verbose = FALSE)
      
    }, error = function(e) {
      message("  WARNING: scDblFinder failed for ", sid, ": ", e$message)
      message("  Proceeding without doublet removal for this sample.")
    })
  } else if (post_qc_n < 100) {
    message("  Skipping doublet detection: too few cells (", post_qc_n, ")")
  }
  
  post_dbl_n <- ncol(obj)
  
  if (post_dbl_n == 0) {
    message("  WARNING: No cells remain after doublet removal! Skipping.")
    next
  }
  
  # -------------------------------------------------------
  # 5f) Save QC summary
  # -------------------------------------------------------
  qc_summaries[[sid]] <- data.frame(
    sample_id = sid,
    gsm = gsm,
    condition = condition,
    pre_qc_cells = pre_qc_n,
    post_qc_cells = post_qc_n,
    qc_kept_frac = round(post_qc_n / pre_qc_n, 4),
    n_doublets = n_doublets,
    post_doublet_cells = post_dbl_n,
    final_kept_frac = round(post_dbl_n / pre_qc_n, 4),
    max_counts_q = round(MAX_COUNTS_Q),
    max_features_q = round(MAX_FEATURES_Q),
    stringsAsFactors = FALSE
  )
  
  seurat_list[[sid]] <- obj
}

# Check results
if (length(seurat_list) == 0) {
  stop("No samples loaded successfully. Check folder structure and file naming.")
}

# Report loading summary
message("\n========== LOADING SUMMARY ==========")
message("Successfully processed: ", length(seurat_list), " samples")
if (length(load_errors) > 0) {
  message("Failed to load: ", length(load_errors), " folders")
  print(load_errors)
}

# Save QC summary
qc_df <- bind_rows(qc_summaries)
message("\n========== QC SUMMARY ==========")
print(qc_df)
write.csv(qc_df, file.path(OUT_DIR, "GSE178146_QC_summary_per_sample.csv"), row.names = FALSE)

# Check condition mapping
cond_check <- sapply(seurat_list, function(o) unique(o$condition))
message("\nCondition mapping per sample:")
print(cond_check)
message("\nCondition distribution:")
print(table(unlist(cond_check), useNA = "ifany"))

unmapped <- names(cond_check)[is.na(cond_check)]
if (length(unmapped) > 0) {
  message("\nWARNING: Samples with unknown condition (NA):")
  print(unmapped)
  message("These samples will be excluded unless you provide clinical metadata.")
}


# ---------------------------
# 6) MERGE SAMPLES
# ---------------------------
message("\nMerging ", length(seurat_list), " samples...")

if (length(seurat_list) == 1) {
  obj_merged <- seurat_list[[1]]
} else {
  obj_merged <- Reduce(function(x, y) merge(x, y), seurat_list)
}

message("Merged cells: ", ncol(obj_merged))
message("\nCells per condition:")
print(table(obj_merged$condition, useNA = "ifany"))

# Drop NA condition cells
na_cells <- sum(is.na(obj_merged$condition))
if (na_cells > 0) {
  message("\nWARNING: ", na_cells, " cells have NA condition. Removing.")
  obj_merged <- subset(obj_merged, subset = !is.na(condition))
}


# ---------------------------
# 7) MINIMAL PROCESSING FOR QC VISUALIZATION
# Note: PCA/UMAP/Clustering here is ONLY for QC plots.
# The merge script will re-normalize, re-integrate (Harmony),
# and re-cluster from the raw counts.
# ---------------------------
message("\nRunning minimal Seurat workflow for QC visualization...")

obj_merged <- NormalizeData(obj_merged, verbose = FALSE)
obj_merged <- FindVariableFeatures(obj_merged, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj_merged <- ScaleData(obj_merged, features = VariableFeatures(obj_merged), verbose = FALSE)
obj_merged <- RunPCA(obj_merged, features = VariableFeatures(obj_merged), npcs = N_PCS_QC, verbose = FALSE)
obj_merged <- RunUMAP(obj_merged, dims = UMAP_DIMS_QC, verbose = FALSE)


# ---------------------------
# 8) QC PLOTS
# ---------------------------
if (SAVE_PLOTS) {
  message("\nSaving QC plots...")
  
  p1 <- VlnPlot(obj_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = "condition", pt.size = 0, ncol = 3)
  ggsave(file.path(OUT_DIR, "QC_violin_by_condition.png"), p1, width = 12, height = 4, dpi = 300)
  
  p2 <- VlnPlot(obj_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = "sample_id", pt.size = 0, ncol = 3)
  ggsave(file.path(OUT_DIR, "QC_violin_by_sample.png"), p2, width = 14, height = 6, dpi = 300)
  
  p3 <- DimPlot(obj_merged, group.by = "condition") + 
    ggtitle("GSE178146 - UMAP by condition (pre-Harmony)")
  ggsave(file.path(OUT_DIR, "UMAP_by_condition_preHarmony.png"), p3, width = 7, height = 6, dpi = 300)
  
  p4 <- DimPlot(obj_merged, group.by = "sample_id", label = FALSE) + 
    ggtitle("GSE178146 - UMAP by sample (pre-Harmony)")
  ggsave(file.path(OUT_DIR, "UMAP_by_sample_preHarmony.png"), p4, width = 9, height = 7, dpi = 300)
  
  sample_counts <- as.data.frame(table(obj_merged$sample_id, obj_merged$condition))
  colnames(sample_counts) <- c("sample_id", "condition", "cells")
  
  p5 <- ggplot(sample_counts, aes(x = sample_id, y = cells, fill = condition)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "GSE178146 - Cells per sample (post-QC + doublet removal)",
         x = "Sample", y = "Number of cells")
  ggsave(file.path(OUT_DIR, "Cells_per_sample_barplot.png"), p5, width = 12, height = 5, dpi = 300)
  
  if ("scDblFinder.score" %in% colnames(obj_merged@meta.data)) {
    p6 <- ggplot(obj_merged@meta.data, aes(x = scDblFinder.score)) +
      geom_histogram(bins = 50, fill = "steelblue", color = "white") +
      theme_minimal() +
      labs(title = "scDblFinder score distribution (retained singlets)",
           x = "Doublet score", y = "Count")
    ggsave(file.path(OUT_DIR, "scDblFinder_score_distribution.png"), p6, width = 7, height = 5, dpi = 300)
  }
}


# ---------------------------
# 9) SAVE OUTPUT
# ---------------------------
output_file <- file.path(OUT_DIR, paste0(DATASET_ID, "_merged_step1.rds"))
saveRDS(obj_merged, file = output_file)
message("\nSaved: ", output_file)


# ---------------------------
# 10) SESSION INFO
# ---------------------------
session_file <- file.path(OUT_DIR, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_file)
message("Session info saved: ", session_file)


# ---------------------------
# 11) FINAL SUMMARY
# ---------------------------
message("\n========== FINAL SUMMARY ==========")
print(obj_merged)
message("\nCells per sample_id:")
print(table(obj_merged$sample_id))
message("\nCells per condition:")
print(table(obj_merged$condition))

message("\nDONE.")
