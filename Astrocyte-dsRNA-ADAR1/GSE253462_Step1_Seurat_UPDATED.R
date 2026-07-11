# ============================================================
# STEP 1 — GSE253462 snRNA (Cingulate Cortex - text matrices) -> Seurat
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(data.table)
  library(scDblFinder)
  library(SingleCellExperiment)
})

set.seed(42)

# ---------------------------
# 1) USER SETTINGS
# ---------------------------
DATASET_ID <- "GSE253462"

library(here)
DATA_DIR <- here::here("data", "raw", "GSE253462")   # GEO files go here
OUT_DIR  <- here::here("data")           # QC'd output
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
# <-- CHANGE THIS


# Clinical/sample metadata

# Check: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253462
SAMPLES <- data.frame(
  gsm = c("GSM8021001", "GSM8021002", "GSM8021003", "GSM8021004", 
          "GSM8021005", "GSM8021006", "GSM8021007", "GSM8021008",
          "GSM8021009", "GSM8021010", "GSM8021011", "GSM8021012",
          "GSM8021013", "GSM8021014", "GSM8021015", "GSM8021040", "GSM8021041"),
  sample_id = c("T362", "T343", "T257", "T551", 
                "T290", "T1619", "T777", "T1952",
                "T308", "T2017", "T1976", "T233",
                "T4162", "T3860", "T4174", "T2465", "T4210"),

  condition = c("PD", "Control", "PD", "Control",
                "PD", "Control", "PD", "Control",
                "PD", "Control", "PD", "PD",
                "PD", "PD", "Control", "Control", "Control"),
  stringsAsFactors = FALSE
)

# --- Unified QC thresholds (consistent across all 6 datasets) ---
MIN_FEATURES <- 500    # Unified (was 200)
MIN_COUNTS   <- 500    # Unified (already 500)
MAX_MT       <- 10
MITO_PATTERN <- "^MT-"

# Upper cutoff strategy
UPPER_Q <- 0.99

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

# Read a count matrix from .txt or .txt.gz file
# Handles GSE253462 format: quoted cell barcodes in header, genes as rows
read_count_matrix <- function(file_path) {
  message("  Reading: ", basename(file_path))
  
  if (grepl("\\.gz$", file_path)) {
    con <- gzfile(file_path, "rt")
  } else {
    con <- file(file_path, "rt")
  }
  
  first_line <- readLines(con, n = 1)
  close(con)
  
  has_quotes <- grepl('^"', first_line)
  message("  Header has quotes: ", has_quotes)
  
  if (grepl("\t", first_line)) {
    sep_char <- "\t"
    message("  Separator: tab")
  } else {
    sep_char <- " "
    message("  Separator: space")
  }
  
  dt <- NULL
  
  # Strategy 1: fread
  tryCatch({
    dt <- fread(file_path, header = TRUE, sep = sep_char, data.table = FALSE, 
                fill = TRUE, strip.white = TRUE)
    message("  Read with fread, dimensions: ", nrow(dt), " x ", ncol(dt))
  }, error = function(e) {
    message("  fread failed: ", e$message)
  })
  
  # If only 1 column, try space separator
  if (!is.null(dt) && ncol(dt) <= 1) {
    message("  Only 1 column detected, trying alternative parsing...")
    dt <- NULL
    tryCatch({
      dt <- fread(file_path, header = TRUE, sep = " ", data.table = FALSE,
                  fill = TRUE, strip.white = TRUE)
      message("  Retry with space separator: ", nrow(dt), " x ", ncol(dt))
    }, error = function(e) {
      message("  Space separator failed: ", e$message)
    })
  }
  
  # Strategy 2: base R read.table
  if (is.null(dt) || ncol(dt) <= 1) {
    message("  Trying base R read.table...")
    tryCatch({
      dt <- read.table(file_path, header = TRUE, sep = "", 
                       row.names = 1, check.names = FALSE,
                       stringsAsFactors = FALSE, fill = TRUE,
                       strip.white = TRUE)
      message("  read.table dimensions: ", nrow(dt), " x ", ncol(dt))
      
      gene_names <- rownames(dt)
      colnames(dt) <- gsub('^"|"$', '', colnames(dt))
      colnames(dt) <- gsub("^X_", "", colnames(dt))
      
      mat <- as.matrix(dt)
      storage.mode(mat) <- "numeric"
      rownames(mat) <- make.unique(gene_names)
      
      mat <- as(mat, "dgCMatrix")
      message("  Final matrix: ", nrow(mat), " genes x ", ncol(mat), " cells")
      return(mat)
      
    }, error = function(e) {
      message("  read.table failed: ", e$message)
    })
  }
  
  # Strategy 3: manual line-by-line parsing
  if (is.null(dt) || ncol(dt) <= 1) {
    message("  Trying manual line-by-line parsing...")
    
    if (grepl("\\.gz$", file_path)) {
      con <- gzfile(file_path, "rt")
    } else {
      con <- file(file_path, "rt")
    }
    
    all_lines <- readLines(con)
    close(con)
    
    header_parts <- strsplit(all_lines[1], "\\s+")[[1]]
    header_parts <- gsub('^"|"$', '', header_parts)
    header_parts <- gsub("^X_", "", header_parts)
    
    message("  Header has ", length(header_parts), " cell barcodes")
    
    data_lines <- all_lines[-1]
    
    parsed_data <- lapply(data_lines, function(line) {
      parts <- strsplit(line, "\\s+")[[1]]
      parts <- gsub('^"|"$', '', parts)
      return(parts)
    })
    
    gene_names <- sapply(parsed_data, function(x) x[1])
    
    count_data <- lapply(parsed_data, function(x) {
      as.numeric(x[-1])
    })
    
    mat <- do.call(rbind, count_data)
    rownames(mat) <- make.unique(gene_names)
    colnames(mat) <- header_parts
    
    mat <- as(mat, "dgCMatrix")
    message("  Final matrix: ", nrow(mat), " genes x ", ncol(mat), " cells")
    return(mat)
  }
  
  # Process the data.table/data.frame result
  message("  Raw dimensions: ", nrow(dt), " rows x ", ncol(dt), " cols")
  
  first_col <- dt[, 1]
  
  if (is.character(first_col)) {
    gene_names <- as.character(first_col)
    count_cols <- dt[, -1, drop = FALSE]
    
    mat <- as.matrix(count_cols)
    storage.mode(mat) <- "numeric"
    rownames(mat) <- make.unique(gene_names)
    
    colnames(mat) <- gsub('^"|"$', '', colnames(mat))
    colnames(mat) <- gsub("^X_", "", colnames(mat))
    
    message("  Format: genes as rows (", nrow(mat), "), cells as columns (", ncol(mat), ")")
  } else {
    message("  Format: cells as rows, genes as columns - transposing")
    mat <- t(as.matrix(dt))
  }
  
  mat <- as(mat, "dgCMatrix")
  message("  Final matrix: ", nrow(mat), " genes x ", ncol(mat), " cells")
  return(mat)
}

# Extract sample ID from filename
extract_sample_id <- function(filename) {
  m <- str_match(basename(filename), "^GSM[0-9]+_([A-Za-z0-9]+)_cingulate")
  if (is.na(m[1, 2])) return(NA_character_)
  return(m[1, 2])
}

# Extract GSM ID from filename
extract_gsm <- function(filename) {
  m <- str_match(basename(filename), "^(GSM[0-9]+)_")
  if (is.na(m[1, 2])) return(NA_character_)
  return(m[1, 2])
}


# ---------------------------
# 3) DISCOVER AND VERIFY FILES
# ---------------------------
message("Checking files in: ", DATA_DIR)

all_files <- list.files(DATA_DIR, pattern = "_cingulate\\.txt(\\.gz)?$", full.names = TRUE)

if (length(all_files) == 0) {
  all_files <- list.files(DATA_DIR, pattern = "\\.txt\\.gz$", full.names = TRUE)
  all_files <- all_files[!grepl("series.matrix|series_matrix", all_files, ignore.case = TRUE)]
}

message("Found ", length(all_files), " sample files:")
for (f in all_files) {
  message("  ", basename(f))
}

if (any(SAMPLES$condition == "Unknown")) {
  message("\n*** WARNING ***")
  message("T2465 and T4210 have 'Unknown' condition!")
  message("MUST verify from GEO supplementary tables before manuscript submission.")
  message("***************\n")
}


# ---------------------------
# 4) PROCESS EACH SAMPLE
# ---------------------------
message("\n========== PROCESSING SAMPLES ==========")

seurat_list <- list()
qc_summaries <- list()

for (file_path in all_files) {
  gsm <- extract_gsm(file_path)
  sid <- extract_sample_id(file_path)
  
  if (is.na(sid)) {
    message("\nSkipping (could not parse sample ID): ", basename(file_path))
    next
  }
  
  # Look up condition from SAMPLES table
  idx <- which(SAMPLES$sample_id == sid | SAMPLES$gsm == gsm)
  if (length(idx) > 0) {
    condition <- SAMPLES$condition[idx[1]]
  } else {
    condition <- "Unknown"
  }
  
  message("\n========== ", sid, " (", gsm, ", ", condition, ") ==========")
  
  # -------------------------------------------------------
  # 4a) Load count matrix
  # -------------------------------------------------------
  counts <- tryCatch({
    read_count_matrix(file_path)
  }, error = function(e) {
    message("  ERROR reading file: ", e$message)
    return(NULL)
  })
  
  if (is.null(counts)) next
  
  if (ncol(counts) < 10 || nrow(counts) < 100) {
    message("  Skipping - matrix too small: ", nrow(counts), " x ", ncol(counts))
    next
  }
  
  # -------------------------------------------------------
  # 4b) Remove MALAT1 and ribosomal genes BEFORE Seurat object creation
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
  # 4c) Create Seurat object
  # -------------------------------------------------------
  obj <- CreateSeuratObject(
    counts = counts,
    project = paste0(DATASET_ID, "_", sid),
    min.cells = 3,
    min.features = 1
  )
  
  obj$dataset   <- DATASET_ID
  obj$gsm       <- gsm
  obj$sample_id <- sid
  obj$condition <- condition
  
  obj <- RenameCells(obj, add.cell.id = paste0(DATASET_ID, "_", sid))
  
  # -------------------------------------------------------
  # 4d) Calculate percent.mt (with NA handling)
  # -------------------------------------------------------
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = MITO_PATTERN)
  
  # Handle NA values in percent.mt
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
  # 4e) QC filtering with unified thresholds
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
    next
  }
  
  message("  After QC: ", pre_qc_n, " -> ", post_qc_n, " (", round(post_qc_n/pre_qc_n*100, 1), "% retained)")
  
  # -------------------------------------------------------
  # 4f) Doublet detection with scDblFinder
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
  # 4g) Save QC summary
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
write.csv(qc_df, file.path(OUT_DIR, "GSE253462_QC_summary_per_sample.csv"), row.names = FALSE)


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

if ("Unknown" %in% obj_merged$condition) {
  message("\n*** WARNING ***")
  message("Samples with 'Unknown' condition: T2465, T4210")
  message("These will introduce noise at merge unless resolved.")
  message("Please update the SAMPLES data frame with correct labels.")
  message("***************\n")
}


# ---------------------------
# 6) MINIMAL PROCESSING FOR QC VISUALIZATION
# Note: PCA/UMAP here is ONLY for QC plots.
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
# 7) QC PLOTS
# ---------------------------
if (SAVE_PLOTS) {
  message("\nSaving QC plots...")
  
  p1 <- VlnPlot(obj_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = "condition", pt.size = 0, ncol = 3)
  ggsave(file.path(OUT_DIR, "QC_violin_by_condition.png"), p1, width = 12, height = 4, dpi = 300)
  
  p2 <- VlnPlot(obj_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                group.by = "sample_id", pt.size = 0, ncol = 3)
  ggsave(file.path(OUT_DIR, "QC_violin_by_sample.png"), p2, width = 16, height = 6, dpi = 300)
  
  p3 <- DimPlot(obj_merged, group.by = "condition") + 
    ggtitle("GSE253462 - UMAP by condition (pre-Harmony)")
  ggsave(file.path(OUT_DIR, "UMAP_by_condition_preHarmony.png"), p3, width = 7, height = 6, dpi = 300)
  
  p4 <- DimPlot(obj_merged, group.by = "sample_id") + 
    ggtitle("GSE253462 - UMAP by sample (pre-Harmony)")
  ggsave(file.path(OUT_DIR, "UMAP_by_sample_preHarmony.png"), p4, width = 10, height = 7, dpi = 300)
  
  sample_counts <- as.data.frame(table(obj_merged$sample_id, obj_merged$condition))
  colnames(sample_counts) <- c("sample_id", "condition", "cells")
  
  p5 <- ggplot(sample_counts, aes(x = sample_id, y = cells, fill = condition)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "GSE253462 - Cells per sample (post-QC + doublet removal)",
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
# 8) SAVE OUTPUT
# ---------------------------
output_file <- file.path(OUT_DIR, paste0(DATASET_ID, "_merged_step1.rds"))
saveRDS(obj_merged, file = output_file)
message("\nSaved: ", output_file)


# ---------------------------
# 9) SESSION INFO
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
