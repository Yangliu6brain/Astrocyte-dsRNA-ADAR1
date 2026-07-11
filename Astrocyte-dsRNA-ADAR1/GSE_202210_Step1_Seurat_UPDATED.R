# ============================================================
# STEP 1 — GSE202210 snRNA (Cortex - flat 10x files) -> Seurat
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
DATASET_ID <- "GSE202210"

# Folder containing the flat 10x files like:
# GSM6106340_HSDG07HC_barcodes.tsv.gz
# GSM6106340_HSDG07HC_genes.tsv.gz
# GSM6106340_HSDG07HC_matrix.mtx.gz

library(here)
DATA_DIR <- here::here("data", "raw", "GSE202210")   # GEO files go here
OUT_DIR  <- here::here("data")            # QC'd output
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)  
# <-- CHANGE THIS


# --- Unified QC thresholds (consistent across all 6 datasets) ---
# Rationale: Prevents systematic bias when merging datasets with
# different QC stringency. See MS2_Script_Review Section 1.2.
MIN_FEATURES <- 500    # Unified (was 200)
MIN_COUNTS   <- 500    # Unified (was 500 — already matched)
MAX_MT       <- 10     # Consistent with original papers (Lee 10%, Smajić 10%)
MITO_PATTERN <- "^MT-" # human

# Upper cutoff strategy (99th percentile)
UPPER_Q <- 0.99

# --- Gene removal settings ---
REMOVE_MALAT1 <- TRUE   # Nuclear-enriched lncRNA that dominates snRNAseq clustering
REMOVE_RIBO   <- TRUE   # Ribosomal genes (RPL/RPS) — not informative in snRNAseq

# --- Doublet detection ---
RUN_DOUBLET_DETECTION <- TRUE  # scDblFinder — current gold standard

# Dim reduction for QC visualization only (NOT used downstream — merge script re-does this)
N_PCS_QC     <- 20
UMAP_DIMS_QC <- 1:20

SAVE_PLOTS <- TRUE


# ---------------------------
# 2) HELPER FUNCTIONS
# ---------------------------

# Find files even if gzipped
find_first_existing <- function(paths) {
  p <- paths[file.exists(paths)]
  if (length(p) == 0) return(NA_character_)
  p[[1]]
}

# Read 10x triplet files (matrix, genes/features, barcodes)
read_10x_triplet <- function(mtx_path, genes_path, barcodes_path) {
  m <- Matrix::readMM(mtx_path)
  genes <- read.delim(genes_path, header = FALSE, stringsAsFactors = FALSE)
  bars  <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)
  
  # genes.tsv typically has 2 columns: gene_id, gene_name
  # Use gene_name if available (column 2), otherwise gene_id (column 1)
  gene_names <- if (ncol(genes) >= 2) genes[[2]] else genes[[1]]
  gene_names <- make.unique(as.character(gene_names))
  
  colnames(m) <- as.character(bars[[1]])
  rownames(m) <- gene_names
  
  # Convert to sparse dgCMatrix
  m <- as(m, "dgCMatrix")
  return(m)
}

# Extract sample ID and GSM from filename
extract_ids_from_filename <- function(filename) {
  m <- str_match(basename(filename), "^(GSM[0-9]+)_([A-Za-z0-9]+)_(matrix|barcodes|genes|features)")
  if (is.na(m[1, 2])) return(list(gsm = NA_character_, sample_id = NA_character_))
  return(list(gsm = m[1, 2], sample_id = m[1, 3]))
}

# Infer condition from sample ID
infer_condition <- function(sample_id) {
  if (is.na(sample_id)) return(NA_character_)
  
  sid_upper <- toupper(sample_id)
  
  if (grepl("HC$", sid_upper)) {
    return("Control")  # HC = Healthy Control
  } else if (grepl("PD", sid_upper)) {
    return("PD")
  } else if (grepl("AD", sid_upper)) {
    return("AD")
  } else if (grepl("CTRL|CTL|CON", sid_upper)) {
    return("Control")
  } else {
    return("Unknown")
  }
}


# ---------------------------
# 3) DISCOVER SAMPLE FILES
# ---------------------------
message("Scanning for files in: ", DATA_DIR)

mtx_files <- list.files(DATA_DIR, pattern = "_matrix\\.mtx(\\.gz)?$", full.names = TRUE)

if (length(mtx_files) == 0) {
  stop("No matrix.mtx files found in DATA_DIR: ", DATA_DIR)
}

message("Found ", length(mtx_files), " matrix files")

samples_found <- data.frame(
  mtx_file = mtx_files,
  stringsAsFactors = FALSE
)

samples_found$gsm <- sapply(samples_found$mtx_file, function(f) extract_ids_from_filename(f)$gsm)
samples_found$sample_id <- sapply(samples_found$mtx_file, function(f) extract_ids_from_filename(f)$sample_id)
samples_found$condition <- sapply(samples_found$sample_id, infer_condition)

message("\nSamples discovered:")
print(samples_found[, c("gsm", "sample_id", "condition")])

message("\nCondition distribution (inferred):")
print(table(samples_found$condition, useNA = "ifany"))


# ---------------------------
# 4) PROCESS EACH SAMPLE
# ---------------------------
message("\n========== PROCESSING SAMPLES ==========")

seurat_list <- list()
qc_summaries <- list()

for (i in seq_len(nrow(samples_found))) {
  mtx_file <- samples_found$mtx_file[i]
  gsm <- samples_found$gsm[i]
  sid <- samples_found$sample_id[i]
  condition <- samples_found$condition[i]
  
  message("\n========== ", sid, " (", gsm, ", ", condition, ") ==========")
  
  # Determine file prefix
  prefix <- sub("_matrix\\.mtx(\\.gz)?$", "", mtx_file)
  
  # Find corresponding genes and barcodes files
  genes_file <- find_first_existing(c(
    paste0(prefix, "_genes.tsv.gz"),
    paste0(prefix, "_genes.tsv"),
    paste0(prefix, "_features.tsv.gz"),
    paste0(prefix, "_features.tsv")
  ))
  
  barcodes_file <- find_first_existing(c(
    paste0(prefix, "_barcodes.tsv.gz"),
    paste0(prefix, "_barcodes.tsv")
  ))
  
  if (is.na(genes_file)) {
    message("  SKIPPING - genes/features file not found")
    next
  }
  if (is.na(barcodes_file)) {
    message("  SKIPPING - barcodes file not found")
    next
  }
  
  message("  Matrix: ", basename(mtx_file))
  message("  Genes: ", basename(genes_file))
  message("  Barcodes: ", basename(barcodes_file))
  
  # Load count matrix
  counts <- tryCatch({
    read_10x_triplet(mtx_file, genes_file, barcodes_file)
  }, error = function(e) {
    message("  ERROR reading files: ", e$message)
    return(NULL)
  })
  
  if (is.null(counts)) next
  
  message("  Raw matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")
  
  # -------------------------------------------------------
  # 4a) Remove MALAT1 and ribosomal genes BEFORE Seurat object creation
  # Rationale: MALAT1 is nuclear-enriched and dominates snRNAseq clustering.
  # Ribosomal genes (RPL/RPS) are not informative in snRNAseq.
  # Martirosyan et al. (GSE243639) excluded MALAT1;
  # Smajić et al. (GSE157783) removed ribosomal genes.
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
  # 4b) Create Seurat object
  # -------------------------------------------------------
  obj <- CreateSeuratObject(
    counts = counts,
    project = paste0(DATASET_ID, "_", sid),
    min.cells = 3,
    min.features = 0
  )
  
  # Add metadata
  obj$dataset   <- DATASET_ID
  obj$gsm       <- gsm
  obj$sample_id <- sid
  obj$condition <- condition
  
  # Rename cells to avoid collisions when merging
  obj <- RenameCells(obj, add.cell.id = paste0(DATASET_ID, "_", sid))
  
  # Calculate percent.mt (MT genes kept for QC metric, just not in HVGs)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = MITO_PATTERN)
  
  # Check for MT genes with alternative patterns if needed
  if (max(obj$percent.mt) == 0) {
    message("  NOTE: No MT genes found with '^MT-', trying alternatives...")
    for (pattern in c("^mt-", "^Mt-")) {
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = pattern)
      if (max(obj$percent.mt) > 0) {
        message("  Found MT genes with pattern: ", pattern)
        break
      }
    }
  }
  
  # -------------------------------------------------------
  # 4c) QC filtering with unified thresholds
  # -------------------------------------------------------
  MAX_COUNTS_Q <- as.numeric(quantile(obj$nCount_RNA, probs = UPPER_Q))
  MAX_FEATURES_Q <- as.numeric(quantile(obj$nFeature_RNA, probs = UPPER_Q))
  
  message("  QC cutoffs:")
  message("    MIN_FEATURES = ", MIN_FEATURES, " (unified)")
  message("    MIN_COUNTS = ", MIN_COUNTS, " (unified)")
  message("    MAX_MT = ", MAX_MT, "%")
  message("    MAX_COUNTS (", UPPER_Q*100, "th %ile) = ", round(MAX_COUNTS_Q))
  message("    MAX_FEATURES (", UPPER_Q*100, "th %ile) = ", round(MAX_FEATURES_Q))
  
  pre_qc_n <- ncol(obj)
  
  if (MAX_FEATURES_Q < MIN_FEATURES || MAX_COUNTS_Q < MIN_COUNTS) {
    message("  WARNING: Upper quantile cutoffs are below minimum thresholds!")
    message("  Adjusting to use only lower thresholds for this sample...")
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
    next
  }
  
  message("  After QC: ", pre_qc_n, " -> ", post_qc_n, " (", round(post_qc_n/pre_qc_n*100, 1), "% retained)")
  
  # -------------------------------------------------------
  # 4d) Doublet detection with scDblFinder
  # Rationale: All original papers used Scrublet or DoubletFinder.
  # scDblFinder is the current gold standard in R.
  # The 99th percentile upper cutoff is NOT a substitute for
  # dedicated doublet detection (MS2_Script_Review Section 1.2).
  # -------------------------------------------------------
  pre_dbl_n <- ncol(obj)
  n_doublets <- 0
  
  if (RUN_DOUBLET_DETECTION && post_qc_n >= 100) {
    message("  Running scDblFinder...")
    
    tryCatch({
      # Minimal processing for doublet detection
      obj_tmp <- NormalizeData(obj, verbose = FALSE)
      obj_tmp <- FindVariableFeatures(obj_tmp, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      obj_tmp <- ScaleData(obj_tmp, verbose = FALSE)
      obj_tmp <- RunPCA(obj_tmp, npcs = N_PCS_QC, verbose = FALSE)
      
      # Convert to SingleCellExperiment for scDblFinder
      sce <- as.SingleCellExperiment(obj_tmp)
      sce <- scDblFinder(sce)
      
      # Store doublet classifications back in Seurat object
      obj$scDblFinder.class <- sce$scDblFinder.class
      obj$scDblFinder.score <- sce$scDblFinder.score
      
      n_doublets <- sum(obj$scDblFinder.class == "doublet")
      message("  scDblFinder: ", n_doublets, " doublets detected (",
              round(n_doublets/ncol(obj)*100, 1), "%)")
      
      # Remove doublets
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
    message("  WARNING: No cells remain after doublet removal! Skipping this sample.")
    next
  }
  
  # -------------------------------------------------------
  # 4e) Save QC summary
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

# Check if we have any samples
if (length(seurat_list) == 0) {
  stop("No samples were successfully processed! Check file paths and formats.")
}

# Save QC summary
qc_df <- bind_rows(qc_summaries)
message("\n========== QC SUMMARY ==========")
print(qc_df)
write.csv(qc_df, file.path(OUT_DIR, "GSE202210_QC_summary_per_sample.csv"), row.names = FALSE)


# ---------------------------
# 5) MERGE ALL SAMPLES
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

# Warn about unknown conditions
if ("Unknown" %in% obj_merged$condition) {
  message("\n*** WARNING ***")
  message("Some samples have 'Unknown' condition!")
  message("Please check GEO metadata and update the infer_condition() function")
  message("***************\n")
}


# ---------------------------
# 6) MINIMAL PROCESSING FOR QC VISUALIZATION
# Note: PCA/UMAP/Clustering is done here ONLY for QC plots.
# The merge script will re-normalize, re-integrate (Harmony),
# and re-cluster from the raw counts in this object.
# ---------------------------
message("\nRunning minimal Seurat workflow for QC visualization...")

obj_merged <- NormalizeData(obj_merged, verbose = FALSE)
obj_merged <- FindVariableFeatures(obj_merged, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj_merged <- ScaleData(obj_merged, features = VariableFeatures(obj_merged), verbose = FALSE)
obj_merged <- RunPCA(obj_merged, features = VariableFeatures(obj_merged), npcs = N_PCS_QC, verbose = FALSE)
obj_merged <- RunUMAP(obj_merged, dims = UMAP_DIMS_QC, verbose = FALSE)


# ---------------------------
# 7) QC PLOTS
# ---------------------------
if (SAVE_PLOTS) {
  message("\nSaving QC plots...")
  
  # Per-condition violin plots
  p1 <- VlnPlot(obj_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = "condition", pt.size = 0, ncol = 3)
  ggsave(file.path(OUT_DIR, "QC_violin_by_condition.png"), p1, width = 12, height = 4, dpi = 300)
  
  # Per-sample violin plots — important for checking distribution consistency
  p2 <- VlnPlot(obj_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = "sample_id", pt.size = 0, ncol = 3)
  ggsave(file.path(OUT_DIR, "QC_violin_by_sample.png"), p2, width = 16, height = 6, dpi = 300)
  
  # UMAP colored by condition and sample (pre-Harmony, for QC only)
  p3 <- DimPlot(obj_merged, group.by = "condition") + 
    ggtitle("GSE202210 - UMAP by condition (pre-Harmony)")
  ggsave(file.path(OUT_DIR, "UMAP_by_condition_preHarmony.png"), p3, width = 7, height = 6, dpi = 300)
  
  p4 <- DimPlot(obj_merged, group.by = "sample_id") + 
    ggtitle("GSE202210 - UMAP by sample (pre-Harmony)")
  ggsave(file.path(OUT_DIR, "UMAP_by_sample_preHarmony.png"), p4, width = 10, height = 7, dpi = 300)
  
  # Sample composition bar plot
  sample_counts <- as.data.frame(table(obj_merged$sample_id, obj_merged$condition))
  colnames(sample_counts) <- c("sample_id", "condition", "cells")
  
  p5 <- ggplot(sample_counts, aes(x = sample_id, y = cells, fill = condition)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "GSE202210 - Cells per sample (post-QC + doublet removal)", 
         x = "Sample", y = "Number of cells")
  ggsave(file.path(OUT_DIR, "Cells_per_sample_barplot.png"), p5, width = 12, height = 5, dpi = 300)
  
  # Doublet score distribution (if available)
  if ("scDblFinder.score" %in% colnames(obj_merged@meta.data)) {
    # Note: after filtering, all remaining cells are singlets
    # This shows the score distribution of retained cells
    p6 <- ggplot(obj_merged@meta.data, aes(x = scDblFinder.score)) +
      geom_histogram(bins = 50, fill = "steelblue", color = "white") +
      theme_minimal() +
      labs(title = "scDblFinder score distribution (retained singlets)",
           x = "Doublet score", y = "Count")
    ggsave(file.path(OUT_DIR, "scDblFinder_score_distribution.png"), p6, width = 7, height = 5, dpi = 300)
  }
}


# ---------------------------
# 8) SAVE OUTPUT
# ---------------------------
output_file <- file.path(OUT_DIR, paste0(DATASET_ID, "_merged_step1.rds"))
saveRDS(obj_merged, file = output_file)
message("\nSaved: ", output_file)


# ---------------------------
# 9) SESSION INFO (for reproducibility)
# ---------------------------
session_file <- file.path(OUT_DIR, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_file)
message("Session info saved: ", session_file)


# ---------------------------
# 10) FINAL SUMMARY
# ---------------------------
message("\n========== FINAL SUMMARY ==========")
print(obj_merged)
message("\nCells per sample:")
print(table(obj_merged$sample_id))
message("\nCells per condition:")
print(table(obj_merged$condition))

message("\nDONE.")
