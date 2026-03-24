suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(openxlsx))

counts_path <- "gene_counts.txt"
output_xlsx <- "PRC_TPM_VST_symbols.xlsx"

if (!file.exists(counts_path)) {
  stop("Missing input: gene_counts.txt")
}

sample_map <- c(
  "29" = "PRC_U_N1 (p2)",
  "30" = "PRC_U_N1 (p3)",
  "31" = "PRC_U_N1 (p4)",
  "32" = "PRC_U_N1 (p5)",
  "33" = "PRC_U_N1 (p6)",
  "34" = "PRC_U_N2 (p2)",
  "35" = "PRC_U_N2 (p3)",
  "36" = "PRC_U_N2 (p4)",
  "37" = "PRC_U_N2 (p5)",
  "38" = "PRC_U_N2 (p6 or p7?)",
  "39" = "PRC_U_N4 (p2)",
  "40" = "PRC_U_N4 (p3)",
  "41" = "PRC_U_N4 (p4)",
  "42" = "PRC_U_N4 (p5)",
  "43" = "PRC_U_N4 (p6)",
  "44" = "PRC_P2 (p2)",
  "45" = "PRC_P2 (p3)",
  "46" = "PRC_P2 (p4)",
  "47" = "PRC_P2 (p5)",
  "48" = "PRC_P2 (p6)"
)
prc_ids <- names(sample_map)

sample_meta <- data.frame(
  sample_id = prc_ids,
  sample_name = unname(sample_map[prc_ids]),
  line = str_extract(unname(sample_map[prc_ids]), "^PRC_[^ ]+"),
  passage = str_match(unname(sample_map[prc_ids]), "\\(p(\\d+)")[, 2],
  stringsAsFactors = FALSE
) %>%
  mutate(
    passage_num = suppressWarnings(as.integer(passage))
  )

raw <- fread(counts_path)
sample_cols <- colnames(raw)[7:ncol(raw)]
missing_ids <- setdiff(prc_ids, sample_cols)

if (length(missing_ids) > 0) {
  stop(sprintf(
    "Missing PRC sample columns in gene_counts.txt: %s",
    paste(missing_ids, collapse = ", ")
  ))
}

counts_dt <- raw[, c("Geneid", "Length", prc_ids), with = FALSE]
counts_dt[, gene_id := gsub("\\.\\d+$", "", Geneid)]
counts_dt[, Length := suppressWarnings(as.numeric(Length))]

# Collapse rows if multiple Ensembl versioned IDs map to the same stable gene ID.
collapsed_counts <- counts_dt[
  ,
  lapply(.SD, sum, na.rm = TRUE),
  by = gene_id,
  .SDcols = prc_ids
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
count_mat <- as.matrix(counts_merged[, ..prc_ids])
mode(count_mat) <- "integer"
rownames(count_mat) <- counts_merged$gene_id

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

count_mat <- as.matrix(counts_merged[, ..prc_ids])
mode(count_mat) <- "integer"
rownames(count_mat) <- counts_merged$gene_id

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
    stop(sprintf(
      "TPM scaling factor is invalid for sample(s): %s",
      paste(bad_samples, collapse = ", ")
    ))
  }

  tpm <- sweep(rate, 2, scaling_factors, "/") * 1e6
  tpm
}

tpm_gene_mat <- calc_tpm(count_mat, counts_merged$gene_length_bp)
tpm_symbol_mat <- rowsum(
  tpm_gene_mat,
  group = counts_merged$gene_symbol[match(rownames(tpm_gene_mat), counts_merged$gene_id)],
  reorder = FALSE
)

symbol_count_mat <- rowsum(
  count_mat,
  group = counts_merged$gene_symbol[match(rownames(count_mat), counts_merged$gene_id)],
  reorder = FALSE
)

dds <- DESeqDataSetFromMatrix(
  countData = symbol_count_mat,
  colData = data.frame(row.names = prc_ids),
  design = ~ 1
)
dds <- dds[rowSums(counts(dds)) > 0, ]
vst_symbol_mat <- assay(vst(dds, blind = TRUE))

# Keep only genes with non-zero TPM in at least one sample.
keep_symbols <- rownames(tpm_symbol_mat)[rowSums(tpm_symbol_mat != 0) > 0]
tpm_symbol_mat <- tpm_symbol_mat[keep_symbols, , drop = FALSE]
vst_symbol_mat <- vst_symbol_mat[rownames(vst_symbol_mat) %in% keep_symbols, , drop = FALSE]

sample_name_lookup <- setNames(sample_meta$sample_name, sample_meta$sample_id)
ordered_sample_names <- unname(sample_name_lookup[prc_ids])

colnames(tpm_symbol_mat) <- ordered_sample_names
colnames(vst_symbol_mat) <- ordered_sample_names

tpm_symbol_mat <- tpm_symbol_mat[order(rownames(tpm_symbol_mat)), , drop = FALSE]
vst_symbol_mat <- vst_symbol_mat[order(rownames(vst_symbol_mat)), , drop = FALSE]

tpm_export <- data.frame(
  gene_symbol = rownames(tpm_symbol_mat),
  tpm_symbol_mat,
  check.names = FALSE
)

vst_export <- data.frame(
  gene_symbol = rownames(vst_symbol_mat),
  vst_symbol_mat,
  check.names = FALSE
)

wb <- createWorkbook()
modifyBaseFont(wb, fontName = "Arial", fontSize = 11)
addWorksheet(wb, "TPM")
addWorksheet(wb, "VST")

writeData(wb, sheet = "TPM", x = tpm_export, withFilter = TRUE)
writeData(wb, sheet = "VST", x = vst_export, withFilter = TRUE)

header_style <- createStyle(
  textDecoration = "bold",
  fontName = "Arial",
  halign = "center"
)

addStyle(
  wb, "TPM", header_style,
  rows = 1, cols = seq_len(ncol(tpm_export)),
  gridExpand = TRUE, stack = TRUE
)
addStyle(
  wb, "VST", header_style,
  rows = 1, cols = seq_len(ncol(vst_export)),
  gridExpand = TRUE, stack = TRUE
)

freezePane(wb, "TPM", firstRow = TRUE, firstCol = TRUE)
freezePane(wb, "VST", firstRow = TRUE, firstCol = TRUE)
setColWidths(wb, "TPM", cols = 1, widths = 18)
setColWidths(wb, "TPM", cols = 2:ncol(tpm_export), widths = 16)
setColWidths(wb, "VST", cols = 1, widths = 18)
setColWidths(wb, "VST", cols = 2:ncol(vst_export), widths = 16)
saveWorkbook(wb, output_xlsx, overwrite = TRUE)

cat("Saved: ", output_xlsx, "\n", sep = "")
cat("Sheets: TPM, VST\n")
