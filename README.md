# Bulk Expression Viewer

Professional repository layout for exploring organoid bulk RNA-seq expression in Google Colab.

The project separates two concerns:

1. local preprocessing of raw featureCounts-style data into compact expression matrices
2. interactive visualization of one gene or multiple genes in Colab

## What is included

- `data/raw/`
  Raw inputs used to build the project dataset.
- `data/processed/`
  Ready-to-load matrices and metadata for visualization.
- `scripts/preprocess_expression.R`
  Reproducible preprocessing pipeline for all 48 samples.
- `scripts/legacy/`
  Older analysis scripts kept as reference.
- `src/bulk_expression_viewer/`
  Python helpers for loading data and drawing plots.
- `notebooks/bulk_expression_viewer_colab.ipynb`
  Technical Colab notebook with inline controls.
- `notebooks/bulk_expression_colab_app.ipynb`
  Colleague-facing Colab app notebook with a cleaner launch flow.
- `requirements.txt`
  Python dependencies for Colab or local notebook runs.

## Processed outputs

After running preprocessing, the following files are created in `data/processed/`:

- `expression_counts.csv.gz`
- `expression_tpm.csv.gz`
- `expression_vst.csv.gz`
- `sample_metadata.csv`
- `gene_symbols.txt`

Matrices are stored as:

- rows: `gene_symbol`
- columns: sample numbers as strings (`"1"` ... `"48"`)

Metadata contains:

- `sample_id`, `sample_key`, `display_name`, `sample_name`
- `line`, `condition`, `project_group`, `passage`
- QC fields from the Excel workbook

## Preprocessing

Run locally from the repository root:

```bash
Rscript scripts/preprocess_expression.R .
```

The script:

- reads raw counts from `data/raw/gene_counts.txt`
- reads sample metadata from `data/raw/Organoids_RNAseq_May2025.xlsx`
- maps Ensembl IDs to HGNC gene symbols with `org.Hs.eg.db`
- collapses duplicated Ensembl versions
- exports symbol-level `counts`, `TPM`, and `VST`

## Google Colab workflow

1. Open the repository in Colab.
2. Prefer opening `notebooks/bulk_expression_colab_app.ipynb` for colleagues.
3. Run the setup form cell.
4. Run the launch form cell.

The notebook supports:

- one gene or multiple genes
- all samples or custom sample number lists
- `counts`, `TPM`, `log2(TPM+1)`, and `VST`
- sample-level barplots
- passage trend plots for single genes
- heatmaps for selected genes
- aggregated heatmaps by line, condition, project group, or passage
- export of the current selection as CSV tables and PNG figures into `data/exports/`

## Recommended repo usage

For a private GitHub repository:

- keep `data/processed/` committed if the matrices are the canonical dataset for colleagues
- keep `data/raw/` only if you want full reproducibility in the same repo
- if raw inputs are sensitive or too large, move them to private storage and keep only `data/processed/`

## Notes

- The viewer is optimized for speed by loading processed matrices instead of the original raw counts file.
- Sample selection is currently designed around sample numbers because this is the most stable and compact UI for colleagues.
- The notebook can be extended later with differential-expression overlays, export buttons, or preset marker panels.
