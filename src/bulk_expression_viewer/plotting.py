from __future__ import annotations

from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def setup_plot_style() -> None:
    sns.set_theme(
        style="whitegrid",
        context="talk",
        rc={
            "axes.facecolor": "#FAF7F2",
            "figure.facecolor": "#FFFDF8",
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


def plot_multi_gene_heatmap(matrix_subset: pd.DataFrame, metadata: pd.DataFrame, scale_label: str) -> plt.Figure:
    ordered_meta = metadata.sort_values(["project_group", "line", "passage", "sample_id"]).copy()
    ordered_cols = ordered_meta["sample_key"].tolist()
    ordered_matrix = matrix_subset.loc[:, ordered_cols].copy()

    fig_width = max(12, min(24, len(ordered_cols) * 0.35))
    fig_height = max(4, min(14, ordered_matrix.shape[0] * 0.55))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        ordered_matrix,
        cmap="mako",
        linewidths=0.4,
        linecolor="#EFEAE0",
        cbar_kws={"label": scale_label},
        ax=ax,
    )
    ax.set_title("Expression heatmap", loc="left", pad=16, fontsize=20, fontweight="bold")
    ax.set_xlabel("Samples")
    ax.set_ylabel("Genes")
    ax.set_xticklabels(ordered_meta["display_name"], rotation=75, ha="right", fontsize=9)
    fig.tight_layout()
    return fig


def plot_aggregated_heatmap(long_df: pd.DataFrame, group_by: str, agg_func: str, scale_label: str) -> plt.Figure:
    if group_by not in long_df.columns:
        raise ValueError(f"Unknown aggregation column: {group_by}")

    grouped = (
        long_df.groupby(["gene_symbol", group_by], dropna=False)["expression"]
        .agg(agg_func)
        .unstack(group_by)
        .sort_index(axis=1)
    )

    fig_width = max(8, min(18, grouped.shape[1] * 1.1))
    fig_height = max(4, min(14, grouped.shape[0] * 0.7))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        grouped,
        cmap="rocket",
        linewidths=0.5,
        linecolor="#EFEAE0",
        cbar_kws={"label": f"{agg_func} {scale_label}"},
        ax=ax,
    )
    ax.set_title(f"Aggregated expression by {group_by}", loc="left", pad=16, fontsize=20, fontweight="bold")
    ax.set_xlabel(group_by)
    ax.set_ylabel("Genes")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    return fig
