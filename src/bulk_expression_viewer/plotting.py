from __future__ import annotations

from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def setup_plot_style() -> None:
    sns.set_theme(
        style="whitegrid",
        context="talk",
        rc={
            "axes.facecolor": "#FAF7F2",
            "figure.facecolor": "#FFFDF8",
            "figure.dpi": 160,
            "savefig.dpi": 300,
            "grid.color": "#D8D2C8",
            "axes.edgecolor": "#B5ADA1",
            "axes.labelcolor": "#2E2A26",
            "xtick.color": "#2E2A26",
            "ytick.color": "#2E2A26",
            "text.color": "#2E2A26",
        },
    )


def _palette(categories: Iterable[str]) -> list[str]:
    unique_count = len(pd.Index(categories).unique())
    palette = sns.color_palette("Spectral", n_colors=max(unique_count, 1))
    return palette


def _sample_labels(metadata: pd.DataFrame) -> list[str]:
    duplicate_name_mask = metadata["sample_name"].duplicated(keep=False)
    labels = metadata["sample_name"].where(~duplicate_name_mask, metadata["display_name"])
    return labels.astype(str).tolist()


def _clip_bounds(values: np.ndarray, clip_quantile: float) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return -1.0, 1.0

    if not 0.5 < clip_quantile < 1.0:
        raise ValueError("clip_quantile must be between 0.5 and 1.0")

    lower = float(np.nanquantile(finite, 1.0 - clip_quantile))
    upper = float(np.nanquantile(finite, clip_quantile))
    if lower == upper:
        lower = float(np.nanmin(finite))
        upper = float(np.nanmax(finite))
        if lower == upper:
            lower -= 1.0
            upper += 1.0
    return lower, upper


def _prepare_heatmap_frame(
    frame: pd.DataFrame,
    scale_label: str,
    transform: str = "quantile_clip",
    clip_quantile: float = 0.99,
) -> tuple[pd.DataFrame, str, dict[str, float]]:
    numeric_frame = frame.astype(float).copy()
    normalized_transform = transform.strip().lower()

    if normalized_transform == "raw":
        return numeric_frame, scale_label, {}

    if normalized_transform == "quantile_clip":
        lower, upper = _clip_bounds(numeric_frame.to_numpy(), clip_quantile)
        clipped = numeric_frame.clip(lower=lower, upper=upper, axis=None)
        label = f"{scale_label} (q{int(round(clip_quantile * 100))} clipped)"
        return clipped, label, {"vmin": lower, "vmax": upper}

    if normalized_transform == "row_zscore":
        row_means = numeric_frame.mean(axis=1)
        row_stds = numeric_frame.std(axis=1).replace(0, np.nan)
        z_frame = numeric_frame.sub(row_means, axis=0).div(row_stds, axis=0).fillna(0.0)
        bound = float(np.nanquantile(np.abs(z_frame.to_numpy()), clip_quantile))
        if not np.isfinite(bound) or bound == 0:
            bound = 1.0
        z_frame = z_frame.clip(lower=-bound, upper=bound, axis=None)
        label = f"{scale_label} row z-score (q{int(round(clip_quantile * 100))})"
        return z_frame, label, {"center": 0.0, "vmin": -bound, "vmax": bound}

    raise ValueError(f"Unsupported heatmap transform: {transform}")


def plot_gene_barplot(long_df: pd.DataFrame, scale_label: str) -> plt.Figure:
    genes = long_df["gene_symbol"].unique().tolist()
    if len(genes) != 1:
        raise ValueError("Barplot requires exactly one gene")

    plot_df = long_df.sort_values(["sample_id"]).copy()
    plot_df["label"] = plot_df["display_name"]

    fig_width = max(12, min(24, len(plot_df) * 0.45))
    fig, ax = plt.subplots(figsize=(fig_width, 6.5))
    sns.barplot(
        data=plot_df,
        x="label",
        y="expression",
        hue="project_group",
        dodge=False,
        palette=_palette(plot_df["project_group"]),
        edgecolor="#5C554A",
        ax=ax,
    )
    ax.set_title(f"{genes[0]} expression across samples", loc="left", pad=16, fontsize=20, fontweight="bold")
    ax.set_xlabel("Sample")
    ax.set_ylabel(scale_label)
    ax.tick_params(axis="x", rotation=75, labelsize=10)
    ax.legend(title="Project group", frameon=False)
    fig.tight_layout()
    return fig


def plot_passage_trend(long_df: pd.DataFrame, scale_label: str) -> plt.Figure:
    genes = long_df["gene_symbol"].unique().tolist()
    if len(genes) != 1:
        raise ValueError("Passage plot requires exactly one gene")

    plot_df = long_df.dropna(subset=["passage"]).sort_values(["line", "passage", "sample_id"]).copy()
    fig, ax = plt.subplots(figsize=(11, 6))
    sns.lineplot(
        data=plot_df,
        x="passage",
        y="expression",
        hue="line",
        style="line",
        markers=True,
        dashes=False,
        linewidth=2.6,
        markersize=8,
        palette=_palette(plot_df["line"]),
        ax=ax,
    )
    ax.set_title(f"{genes[0]} expression by passage", loc="left", pad=16, fontsize=20, fontweight="bold")
    ax.set_xlabel("Passage")
    ax.set_ylabel(scale_label)
    ax.legend(title="Line", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    return fig


def plot_multi_gene_heatmap(
    matrix_subset: pd.DataFrame,
    metadata: pd.DataFrame,
    scale_label: str,
    transform: str = "quantile_clip",
    clip_quantile: float = 0.99,
) -> plt.Figure:
    ordered_meta = metadata.sort_values(["project_group", "line", "passage", "sample_id"]).copy()
    ordered_cols = ordered_meta["sample_key"].tolist()
    ordered_matrix = matrix_subset.loc[:, ordered_cols].copy()
    sample_labels = _sample_labels(ordered_meta)
    display_matrix, cbar_label, heatmap_kwargs = _prepare_heatmap_frame(
        ordered_matrix,
        scale_label=scale_label,
        transform=transform,
        clip_quantile=clip_quantile,
    )

    fig_width = max(12, min(24, len(ordered_cols) * 0.35))
    fig_height = max(4, min(14, display_matrix.shape[0] * 0.55))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        display_matrix,
        cmap="bwr",
        linewidths=0.4,
        linecolor="#EFEAE0",
        cbar_kws={"label": cbar_label},
        **heatmap_kwargs,
        ax=ax,
    )
    ax.set_title("Expression heatmap", loc="left", pad=16, fontsize=20, fontweight="bold")
    ax.set_xlabel("Samples")
    ax.set_ylabel("Genes")
    ax.set_xticks(np.arange(display_matrix.shape[1]) + 0.5)
    ax.set_xticklabels(sample_labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(np.arange(display_matrix.shape[0]) + 0.5)
    ax.set_yticklabels(display_matrix.index.tolist(), rotation=0, va="center", fontsize=10)
    fig.tight_layout()
    return fig


def plot_aggregated_heatmap(
    long_df: pd.DataFrame,
    group_by: str,
    agg_func: str,
    scale_label: str,
    transform: str = "quantile_clip",
    clip_quantile: float = 0.99,
) -> plt.Figure:
    if group_by not in long_df.columns:
        raise ValueError(f"Unknown aggregation column: {group_by}")

    grouped = (
        long_df.groupby(["gene_symbol", group_by], dropna=False)["expression"]
        .agg(agg_func)
        .unstack(group_by)
        .sort_index(axis=1)
    )
    display_grouped, cbar_label, heatmap_kwargs = _prepare_heatmap_frame(
        grouped,
        scale_label=f"{agg_func} {scale_label}",
        transform=transform,
        clip_quantile=clip_quantile,
    )

    fig_width = max(8, min(18, display_grouped.shape[1] * 1.1))
    fig_height = max(4, min(14, display_grouped.shape[0] * 0.7))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        display_grouped,
        cmap="bwr",
        linewidths=0.5,
        linecolor="#EFEAE0",
        cbar_kws={"label": cbar_label},
        **heatmap_kwargs,
        ax=ax,
    )
    ax.set_title(f"Aggregated expression by {group_by}", loc="left", pad=16, fontsize=20, fontweight="bold")
    ax.set_xlabel(group_by)
    ax.set_ylabel("Genes")
    ax.set_xticks(np.arange(display_grouped.shape[1]) + 0.5)
    ax.set_xticklabels(display_grouped.columns.astype(str).tolist(), rotation=45, ha="right")
    ax.set_yticks(np.arange(display_grouped.shape[0]) + 0.5)
    ax.set_yticklabels(display_grouped.index.tolist(), rotation=0, va="center")
    fig.tight_layout()
    return fig
