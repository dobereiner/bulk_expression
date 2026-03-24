suppressPackageStartupMessages({
  ensure_cran <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }

  ensure_bioc <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
  }

  for (pkg in c("data.table", "dplyr", "stringr", "readxl")) {
    ensure_cran(pkg)
  }

  for (pkg in c("DESeq2", "AnnotationDbi", "org.Hs.eg.db")) {
    ensure_bioc(pkg)
  }

  library(data.table)
  library(dplyr)
  library(stringr)
  library(readxl)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = FALSE) else normalizePath(".", winslash = "/", mustWork = TRUE)

counts_path <- file.path(project_dir, "data", "raw", "gene_counts.txt")
metadata_path <- file.path(project_dir, "data", "raw", "Organoids_RNAseq_May2025.xlsx")
output_dir <- file.path(project_dir, "data", "processed")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(counts_path)) {
  stop(sprintf("Missing input counts file: %s", counts_path))
}

if (!file.exists(metadata_path)) {
  stop(sprintf("Missing metadata workbook: %s", metadata_path))
}

extract_sample_metadata <- function(xlsx_path) {
  raw_meta <- read_excel(xlsx_path, sheet = "исход образцы", col_names = FALSE)
  colnames(raw_meta) <- paste0("col", seq_len(ncol(raw_meta)))

  meta <- raw_meta %>%
    dplyr::transmute(
      sample_id_raw = col1,
      sample_name = col2,
      rna_conc_nanodrop_ng_ul = suppressWarnings(as.numeric(col3)),
      rna_conc_qubit_ng_ul = suppressWarnings(as.numeric(col4)),
      rin = suppressWarnings(as.numeric(col5)),
      magicpure_sample_ul = suppressWarnings(as.numeric(col6)),
      magicpure_water_ul = suppressWarnings(as.numeric(col7)),
      dna_conc_qubit_ng_ul = suppressWarnings(as.numeric(col9)),
      fragment_length_bp = suppressWarnings(as.numeric(col10))
    ) %>%
    dplyr::mutate(
      sample_id = str_extract(as.character(sample_id_raw), "\\d+"),
      sample_id = suppressWarnings(as.integer(sample_id)),
      sample_name = str_squish(as.character(sample_name))
    ) %>%
    dplyr::filter(!is.na(sample_id), sample_id >= 1, sample_id <= 48, !is.na(sample_name), sample_name != "") %>%
    dplyr::distinct(sample_id, .keep_all = TRUE) %>%
    dplyr::arrange(sample_id) %>%
    dplyr::mutate(
      sample_key = as.character(sample_id),
      display_name = sprintf("%02d | %s", sample_id, sample_name),
      passage = str_match(sample_name, "\\(p(\\d+)[^\\)]*\\)")[, 2],
      passage = suppressWarnings(as.integer(passage)),
      sample_name_no_passage = str_squish(str_remove(sample_name, "\\s*\\(p\\d+[^\\)]*\\)")),
      condition = case_when(
        str_detect(sample_name_no_passage, " full$") ~ "full",
        str_detect(sample_name_no_passage, " -pp$") ~ "-pp",
        str_detect(sample_name_no_passage, " -nog -rspo$") ~ "-nog -rspo",
        str_detect(sample_name_no_passage, " -nog-rspo$") ~ "-nog -rspo",
        str_detect(sample_name_no_passage, " -nog$") ~ "-nog",
        str_detect(sample_name_no_passage, " -rspo$") ~ "-rspo",
        TRUE ~ "baseline"
      ),
      line = case_when(
        condition == "baseline" ~ sample_name_no_passage,
        condition == "full" ~ str_remove(sample_name_no_passage, " full$"),
        condition == "-pp" ~ str_remove(sample_name_no_passage, " -pp$"),
        condition == "-nog -rspo" ~ str_remove(sample_name_no_passage, " -(nog -rspo|nog-rspo)$"),
        condition == "-nog" ~ str_remove(sample_name_no_passage, " -nog$"),
        condition == "-rspo" ~ str_remove(sample_name_no_passage, " -rspo$"),
        TRUE ~ sample_name_no_passage
      ),
      project_group = case_when(
        str_starts(line, "PRC_") ~ "PRC",
        str_starts(line, "CC_") ~ "CC",
        TRUE ~ "other"
      )
    ) %>%
    dplyr::select(
      sample_id,
      sample_key,
      display_name,
      sample_name,
      line,
      condition,
      project_group,
      passage,
      rna_conc_nanodrop_ng_ul,
      rna_conc_qubit_ng_ul,
      rin,
      magicpure_sample_ul,
      magicpure_water_ul,
      dna_conc_qubit_ng_ul,
      fragment_length_bp
    )

  if (nrow(meta) != 48) {
    stop(sprintf("Expected 48 samples in workbook, found %d", nrow(meta)))
  }

  meta
}

calc_tpm <- function(count_matrix, gene_length_bp) {
  valid <- is.finite(gene_length_bp) & gene_length_bp > 0
  if (!any(valid)) {
    stop("No genes with positive length available for TPM calculation")
  }

  count_matrix <- count_matrix[valid, , drop = FALSE]
  gene_length_kb <- gene_length_bp[valid] / 1000
  rate <- sweep(count_matrix, 1, gene_length_kb, "/")
  rate[!is.finite(rate)] <- 0

  scaling_factors <- colSums(rate)
  if (any(!is.finite(scaling_factors) | scaling_factors <= 0)) {
    bad_samples <- colnames(count_matrix)[!is.finite(scaling_factors) | scaling_factors <= 0]
    stop(sprintf("TPM scaling factor is invalid for sample(s): %s", paste(bad_samples, collapse = ", ")))
  }

  sweep(rate, 2, scaling_factors, "/") * 1e6
}

sample_meta <- extract_sample_metadata(metadata_path)
sample_ids <- sample_meta$sample_key

raw_counts <- fread(counts_path)
count_sample_cols <- colnames(raw_counts)[7:ncol(raw_counts)]
missing_ids <- setdiff(sample_ids, count_sample_cols)

if (length(missing_ids) > 0) {
  stop(sprintf("Missing sample columns in gene_counts.txt: %s", paste(missing_ids, collapse = ", ")))
}

counts_dt <- raw_counts[, c("Geneid", "Length", sample_ids), with = FALSE]
counts_dt[, gene_id := gsub("\\.\\d+$", "", Geneid)]
counts_dt[, Length := suppressWarnings(as.numeric(Length))]

collapsed_counts <- counts_dt[
  ,
  lapply(.SD, sum, na.rm = TRUE),
  by = gene_id,
  .SDcols = sample_ids
]

collapsed_lengths <- counts_dt[
  ,
  .(
    gene_length_bp = {
      valid_lengths <- unique(Length[is.finite(Length) & Length > 0])
      if (length(valid_lengths) == 0) {
        NA_real_
      } else {
        max(valid_lengths)
      }
    }
  ),
  by = gene_id
]

counts_merged <- merge(
  collapsed_lengths,
  collapsed_counts,
  by = "gene_id",
  all.y = TRUE,
  sort = FALSE
)
counts_merged <- counts_merged[match(unique(collapsed_counts$gene_id), counts_merged$gene_id)]

gene_symbol_map <- withCallingHandlers(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = counts_merged$gene_id,
    columns = "SYMBOL",
    keytype = "ENSEMBL"
  ),
  warning = function(w) {
    if (grepl("1:many mapping", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  }
) %>%
  as.data.table() %>%
  .[!is.na(SYMBOL) & SYMBOL != ""] %>%
  .[, .SD[1], by = ENSEMBL]

gene_symbol_lookup <- setNames(gene_symbol_map$SYMBOL, gene_symbol_map$ENSEMBL)
counts_merged$gene_symbol <- unname(gene_symbol_lookup[counts_merged$gene_id])
counts_merged <- counts_merged[!is.na(gene_symbol) & gene_symbol != ""]

if (nrow(counts_merged) == 0) {
  stop("No Ensembl IDs were mapped to gene symbols")
}

count_mat <- as.matrix(counts_merged[, ..sample_ids])
mode(count_mat) <- "numeric"
rownames(count_mat) <- counts_merged$gene_id

symbol_count_mat <- rowsum(
  count_mat,
  group = counts_merged$gene_symbol[match(rownames(count_mat), counts_merged$gene_id)],
  reorder = FALSE
)

tpm_gene_mat <- calc_tpm(count_mat, counts_merged$gene_length_bp)
tpm_symbol_mat <- rowsum(
  tpm_gene_mat,
  group = counts_merged$gene_symbol[match(rownames(tpm_gene_mat), counts_merged$gene_id)],
  reorder = FALSE
)

dds <- DESeqDataSetFromMatrix(
  countData = round(symbol_count_mat),
  colData = data.frame(row.names = sample_ids),
  design = ~ 1
)
dds <- dds[rowSums(counts(dds)) > 0, ]
vst_symbol_mat <- assay(vst(dds, blind = TRUE))

keep_symbols <- rownames(tpm_symbol_mat)[rowSums(tpm_symbol_mat > 0) > 0]
symbol_count_mat <- symbol_count_mat[keep_symbols, , drop = FALSE]
tpm_symbol_mat <- tpm_symbol_mat[keep_symbols, , drop = FALSE]
vst_symbol_mat <- vst_symbol_mat[keep_symbols, , drop = FALSE]

symbol_count_mat <- symbol_count_mat[order(rownames(symbol_count_mat)), sample_ids, drop = FALSE]
tpm_symbol_mat <- tpm_symbol_mat[order(rownames(tpm_symbol_mat)), sample_ids, drop = FALSE]
vst_symbol_mat <- vst_symbol_mat[order(rownames(vst_symbol_mat)), sample_ids, drop = FALSE]

export_matrix <- function(mat, out_path) {
  export_df <- data.frame(
    gene_symbol = rownames(mat),
    mat,
    check.names = FALSE
  )
  fwrite(export_df, out_path)
}

export_matrix(symbol_count_mat, file.path(output_dir, "expression_counts.csv.gz"))
export_matrix(tpm_symbol_mat, file.path(output_dir, "expression_tpm.csv.gz"))
export_matrix(vst_symbol_mat, file.path(output_dir, "expression_vst.csv.gz"))
fwrite(sample_meta, file.path(output_dir, "sample_metadata.csv"))
writeLines(sort(rownames(tpm_symbol_mat)), file.path(output_dir, "gene_symbols.txt"))

cat("Saved processed data to:", output_dir, "\n")
cat("Genes:", nrow(tpm_symbol_mat), "\n")
cat("Samples:", ncol(tpm_symbol_mat), "\n")
