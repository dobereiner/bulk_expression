from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd


def _safe_slug(parts: Iterable[str]) -> str:
    raw = "_".join(parts)
    keep = [char if char.isalnum() or char in {"_", "-"} else "_" for char in raw]
    slug = "".join(keep)
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug.strip("_") or "selection"


def export_render_bundle(
    long_df: pd.DataFrame,
    matrix_subset: pd.DataFrame,
    metadata_subset: pd.DataFrame,
    figures: dict[str, object],
    output_dir: str | Path = "exports",
) -> Path:
    export_root = Path(output_dir)
    export_root.mkdir(parents=True, exist_ok=True)

    gene_slug = _safe_slug(long_df["gene_symbol"].astype(str).unique().tolist()[:5])
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    bundle_dir = export_root / f"{timestamp}_{gene_slug}"
    bundle_dir.mkdir(parents=True, exist_ok=False)

    long_df.to_csv(bundle_dir / "expression_long.csv", index=False)
    matrix_subset.to_csv(bundle_dir / "expression_matrix.csv")
    metadata_subset.to_csv(bundle_dir / "sample_metadata.csv", index=False)

    for figure_name, figure in figures.items():
        if figure is None:
            continue
        figure.savefig(bundle_dir / f"{_safe_slug([figure_name])}.png", dpi=200, bbox_inches="tight")

    return bundle_dir
