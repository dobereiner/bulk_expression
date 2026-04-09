from __future__ import annotations

from dataclasses import dataclass
from difflib import get_close_matches
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


@dataclass
class ProjectData:
    counts: pd.DataFrame
    tpm: pd.DataFrame
    vst: pd.DataFrame
    metadata: pd.DataFrame


def _read_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["gene_symbol"] = df["gene_symbol"].astype(str)
    return df.set_index("gene_symbol")


def load_project_data(data_dir: str | Path = "data/processed") -> ProjectData:
    data_path = Path(data_dir)
    metadata = pd.read_csv(data_path / "sample_metadata.csv", dtype={"sample_key": str})
    metadata["sample_key"] = metadata["sample_key"].astype(str)
    metadata = metadata.sort_values("sample_id").reset_index(drop=True)

    counts = _read_matrix(data_path / "expression_counts.csv.gz")
    tpm = _read_matrix(data_path / "expression_tpm.csv.gz")
    vst = _read_matrix(data_path / "expression_vst.csv.gz")

    ordered_cols = metadata["sample_key"].tolist()
    for frame in (counts, tpm, vst):
        missing = set(ordered_cols) - set(frame.columns)
        if missing:
            raise ValueError(f"Missing columns in expression matrix: {sorted(missing)}")
        frame = frame.loc[:, ordered_cols]

    counts = counts.loc[:, ordered_cols]
    tpm = tpm.loc[:, ordered_cols]
    vst = vst.loc[:, ordered_cols]

    return ProjectData(counts=counts, tpm=tpm, vst=vst, metadata=metadata)


def parse_gene_input(raw_text: str) -> list[str]:
    if not raw_text:
        return []

    tokens = (
        pd.Series(raw_text.replace("\n", " ").split(" "))
        .astype(str)
        .str.strip()
        .replace("", pd.NA)
        .dropna()
        .tolist()
    )
    return list(dict.fromkeys(tokens))


def parse_sample_input(raw_text: str, available_sample_ids: Iterable[int]) -> list[str]:
    available = {int(sample_id) for sample_id in available_sample_ids}
    if not raw_text or not raw_text.strip():
        return [str(sample_id) for sample_id in sorted(available)]

    tokens = (
        pd.Series(raw_text.replace("\n", " ").split(" "))
        .astype(str)
        .str.strip()
        .replace("", pd.NA)
        .dropna()
        .tolist()
    )

    sample_ids: list[str] = []
    invalid: list[str] = []
    for token in tokens:
        if token.isdigit() and int(token) in available:
            sample_ids.append(str(int(token)))
        else:
            invalid.append(token)

    if invalid:
        raise ValueError(f"Unknown sample numbers: {', '.join(invalid)}")

    return list(dict.fromkeys(sample_ids))


def resolve_genes(requested_genes: list[str], matrix_index: Iterable[str]) -> tuple[list[str], dict[str, list[str]]]:
    available = list(matrix_index)
    uppercase_lookup = {gene.upper(): gene for gene in available}
    resolved: list[str] = []
    suggestions: dict[str, list[str]] = {}

    for gene in requested_genes:
        canonical = uppercase_lookup.get(gene.upper())
        if canonical is not None:
            resolved.append(canonical)
            continue
        suggestions[gene] = get_close_matches(gene.upper(), uppercase_lookup.keys(), n=5, cutoff=0.5)

    resolved = list(dict.fromkeys(resolved))
    pretty_suggestions = {
        key: [uppercase_lookup[item] for item in values]
        for key, values in suggestions.items()
        if values
    }
    return resolved, pretty_suggestions


def get_scale_matrix(scale: str, project_data: ProjectData) -> pd.DataFrame:
    normalized_scale = scale.strip().lower()
    if normalized_scale == "counts":
        return project_data.counts
    if normalized_scale == "tpm":
        return project_data.tpm
    if normalized_scale in {"log2(tpm+1)", "log2_tpm", "log2 tpm"}:
        return np.log2(project_data.tpm + 1.0)
    if normalized_scale == "vst":
        return project_data.vst
    raise ValueError(f"Unsupported scale: {scale}")


def matrix_to_long(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
    sample_keys: list[str],
) -> pd.DataFrame:
    subset = matrix.loc[genes, sample_keys]
    long_df = (
        subset.rename_axis(index="gene_symbol", columns="sample_key")
        .stack()
        .rename("expression")
        .reset_index()
        .merge(metadata, on="sample_key", how="left")
    )
    return long_df


def summarize_selection(long_df: pd.DataFrame, genes: list[str], sample_keys: list[str]) -> pd.DataFrame:
    subset = long_df.loc[long_df["gene_symbol"].isin(genes) & long_df["sample_key"].isin(sample_keys)].copy()
    subset["expression"] = subset["expression"].astype(float)
    return subset
