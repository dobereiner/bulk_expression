from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import ipywidgets as widgets
from IPython.display import HTML, Markdown, clear_output, display

from .data import (
    get_scale_matrix,
    load_project_data,
    matrix_to_long,
    parse_gene_input,
    parse_sample_input,
    resolve_genes,
)
from .exporting import export_render_bundle
from .plotting import (
    plot_aggregated_heatmap,
    plot_gene_barplot,
    plot_multi_gene_heatmap,
    plot_passage_trend,
    setup_plot_style,
)


def launch_colab_viewer(data_dir: str | Path = "data/processed", export_dir: str | Path = "data/exports") -> None:
    setup_plot_style()
    project = load_project_data(data_dir)
    export_dir = Path(export_dir)

    metadata = project.metadata.copy()
    project_groups = sorted(metadata["project_group"].dropna().astype(str).unique().tolist())
    lines = sorted(metadata["line"].dropna().astype(str).unique().tolist())
    conditions = sorted(metadata["condition"].dropna().astype(str).unique().tolist())

    title = widgets.HTML(
        """
        <div style="padding: 16px 18px; border-radius: 16px; background: linear-gradient(135deg, #FFF3D6 0%, #F6E7C9 100%);
                    border: 1px solid #E7D2A7; margin-bottom: 12px;">
          <div style="font-size: 28px; font-weight: 700; color: #2F2618;">Bulk Expression Viewer</div>
          <div style="font-size: 14px; color: #5F523D; margin-top: 6px;">
            Query one gene or multiple genes, filter samples, switch expression scales, and export figures for colleagues.
          </div>
        </div>
        """
    )

    genes_widget = widgets.Textarea(
        value="EPCAM, VIM, MKI67",
        description="Genes",
        placeholder="Enter one gene or a comma-separated list",
        layout=widgets.Layout(width="100%", height="90px"),
        style={"description_width": "initial"},
    )
    scale_widget = widgets.Dropdown(
        options=["counts", "TPM", "log2(TPM+1)", "VST"],
        value="TPM",
        description="Scale",
        style={"description_width": "initial"},
    )
    sample_mode_widget = widgets.ToggleButtons(
        options=["All samples", "Filter widgets", "Custom numbers"],
        value="All samples",
        description="Samples",
        style={"description_width": "initial"},
        layout=widgets.Layout(width="100%"),
    )
    project_group_widget = widgets.SelectMultiple(
        options=project_groups,
        value=tuple(project_groups),
        description="Project group",
        rows=min(5, max(len(project_groups), 2)),
        style={"description_width": "initial"},
        layout=widgets.Layout(width="100%"),
    )
    line_widget = widgets.SelectMultiple(
        options=lines,
        value=tuple(lines[: min(8, len(lines))]),
        description="Line",
        rows=10,
        style={"description_width": "initial"},
        layout=widgets.Layout(width="100%"),
    )
    condition_widget = widgets.SelectMultiple(
        options=conditions,
        value=tuple(conditions),
        description="Condition",
        rows=min(6, max(len(conditions), 2)),
        style={"description_width": "initial"},
        layout=widgets.Layout(width="100%"),
    )
    sample_numbers_widget = widgets.Textarea(
        value="29, 30, 31, 32, 33, 34, 35, 36, 37, 38",
        description="Sample numbers",
        placeholder="Example: 1, 2, 3, 29, 30",
        layout=widgets.Layout(width="100%", height="80px"),
        style={"description_width": "initial"},
        disabled=True,
    )
    aggregate_widget = widgets.Checkbox(
        value=False,
        description="Add aggregated heatmap",
        indent=False,
    )
    groupby_widget = widgets.Dropdown(
        options=[("Line", "line"), ("Condition", "condition"), ("Project group", "project_group"), ("Passage", "passage")],
        value="line",
        description="Aggregate by",
        style={"description_width": "initial"},
    )
    aggfunc_widget = widgets.Dropdown(
        options=["mean", "median"],
        value="mean",
        description="Reducer",
        style={"description_width": "initial"},
    )
    table_rows_widget = widgets.IntSlider(
        value=25,
        min=5,
        max=100,
        step=5,
        description="Preview rows",
        style={"description_width": "initial"},
    )
    render_button = widgets.Button(
        description="Render",
        button_style="warning",
        icon="bar-chart",
        layout=widgets.Layout(width="180px"),
    )
    export_button = widgets.Button(
        description="Export bundle",
        button_style="success",
        icon="download",
        layout=widgets.Layout(width="180px"),
    )
    reset_button = widgets.Button(
        description="Reset filters",
        button_style="",
        icon="refresh",
        layout=widgets.Layout(width="180px"),
    )

    status_output = widgets.Output()
    plot_output = widgets.Output()
    table_output = widgets.Output()
    last_render: dict[str, object] = {}

    help_box = widgets.HTML(
        """
        <div style="padding: 12px 14px; border-radius: 14px; background: #F9F4EA; border: 1px solid #E6DCCB;">
          <b>Usage</b><br>
          1. Enter one or more gene symbols.<br>
          2. Choose the expression scale.<br>
          3. Show all samples, filter with widgets, or enter sample numbers manually.<br>
          4. Click <b>Render</b>.<br>
          5. Click <b>Export bundle</b> to save CSV tables and PNG figures.
        </div>
        """
    )

    filter_box = widgets.VBox(
        [
            widgets.HTML("<b>Metadata filters</b>"),
            project_group_widget,
            condition_widget,
            line_widget,
        ],
        layout=widgets.Layout(width="100%"),
    )

    settings_box = widgets.VBox(
        [
            genes_widget,
            scale_widget,
            sample_mode_widget,
            sample_numbers_widget,
            aggregate_widget,
            widgets.HBox([groupby_widget, aggfunc_widget]),
            table_rows_widget,
            widgets.HBox([render_button, export_button, reset_button]),
        ],
        layout=widgets.Layout(width="100%"),
    )

    layout = widgets.VBox(
        [
            title,
            help_box,
            widgets.HBox(
                [
                    widgets.VBox([settings_box], layout=widgets.Layout(width="55%")),
                    widgets.VBox([filter_box], layout=widgets.Layout(width="45%")),
                ]
            ),
            status_output,
            plot_output,
            table_output,
        ]
    )

    def _toggle_sample_mode(change: dict) -> None:
        sample_numbers_widget.disabled = change["new"] != "Custom numbers"
        filters_disabled = change["new"] != "Filter widgets"
        project_group_widget.disabled = filters_disabled
        line_widget.disabled = filters_disabled
        condition_widget.disabled = filters_disabled

    def _selected_sample_keys() -> list[str]:
        if sample_mode_widget.value == "All samples":
            return metadata["sample_key"].tolist()

        if sample_mode_widget.value == "Custom numbers":
            return parse_sample_input(sample_numbers_widget.value, metadata["sample_id"].tolist())

        filtered = metadata.copy()
        project_groups_selected = list(project_group_widget.value)
        lines_selected = list(line_widget.value)
        conditions_selected = list(condition_widget.value)

        if project_groups_selected:
            filtered = filtered[filtered["project_group"].isin(project_groups_selected)]
        if lines_selected:
            filtered = filtered[filtered["line"].isin(lines_selected)]
        if conditions_selected:
            filtered = filtered[filtered["condition"].isin(conditions_selected)]

        if filtered.empty:
            raise ValueError("Filter widgets returned zero samples.")

        return filtered.sort_values("sample_id")["sample_key"].tolist()

    def _render(_=None) -> None:
        nonlocal last_render

        with status_output:
            clear_output(wait=True)
            try:
                requested_genes = parse_gene_input(genes_widget.value)
                if not requested_genes:
                    raise ValueError("Please provide at least one gene symbol.")

                sample_keys = _selected_sample_keys()
                selected_meta = metadata[metadata["sample_key"].isin(sample_keys)].copy().sort_values("sample_id")
                sample_keys = selected_meta["sample_key"].tolist()

                matrix = get_scale_matrix(scale_widget.value, project)
                resolved_genes, suggestions = resolve_genes(requested_genes, matrix.index)
                if not resolved_genes:
                    suggestion_lines = []
                    for gene, matches in suggestions.items():
                        if matches:
                            suggestion_lines.append(f"- `{gene}`: try {', '.join(matches)}")
                    suggestion_text = "\n".join(suggestion_lines) if suggestion_lines else "No close matches found."
                    raise ValueError(f"No gene symbols were resolved.\n\n{suggestion_text}")

                long_df = matrix_to_long(matrix, selected_meta, resolved_genes, sample_keys)
                matrix_subset = matrix.loc[resolved_genes, sample_keys]

                summary_md = (
                    f"**Resolved genes:** {', '.join(resolved_genes)}  \n"
                    f"**Samples shown:** {len(sample_keys)}  \n"
                    f"**Scale:** `{scale_widget.value}`  \n"
                    f"**Export folder:** `{export_dir}`"
                )
                display(Markdown(summary_md))

                unresolved = [gene for gene in requested_genes if gene not in resolved_genes]
                if unresolved:
                    unresolved_lines = []
                    for gene in unresolved:
                        matches = suggestions.get(gene, [])
                        unresolved_lines.append(f"- `{gene}`: {', '.join(matches) if matches else 'no close match'}")
                    display(Markdown("**Unresolved genes**\n" + "\n".join(unresolved_lines)))

                figures: dict[str, object] = {}
                with plot_output:
                    clear_output(wait=True)

                    if len(resolved_genes) == 1:
                        fig = plot_gene_barplot(long_df, scale_widget.value)
                        figures["sample_barplot"] = fig
                        display(fig)
                        plt.close(fig)

                        if long_df["passage"].notna().sum() >= 2:
                            fig = plot_passage_trend(long_df, scale_widget.value)
                            figures["passage_trend"] = fig
                            display(fig)
                            plt.close(fig)

                    fig = plot_multi_gene_heatmap(matrix_subset, selected_meta, scale_widget.value)
                    figures["heatmap"] = fig
                    display(fig)
                    plt.close(fig)

                    if aggregate_widget.value:
                        fig = plot_aggregated_heatmap(long_df, groupby_widget.value, aggfunc_widget.value, scale_widget.value)
                        figures["aggregated_heatmap"] = fig
                        display(fig)
                        plt.close(fig)

                preview = long_df[
                    ["gene_symbol", "sample_id", "sample_name", "line", "condition", "project_group", "passage", "expression"]
                ].sort_values(["gene_symbol", "sample_id"])

                with table_output:
                    clear_output(wait=True)
                    display(preview.head(table_rows_widget.value).reset_index(drop=True))

                last_render = {
                    "long_df": preview.copy(),
                    "matrix_subset": matrix_subset.copy(),
                    "metadata_subset": selected_meta.copy(),
                    "figures": figures,
                }

            except Exception as exc:
                with plot_output:
                    clear_output(wait=True)
                with table_output:
                    clear_output(wait=True)
                display(HTML(f"<div style='color:#8A1C1C; font-weight:600;'>Error: {exc}</div>"))

    def _export(_=None) -> None:
        with status_output:
            if not last_render:
                display(Markdown("**Error:** render plots first, then export."))
                return

            bundle_dir = export_render_bundle(
                long_df=last_render["long_df"],
                matrix_subset=last_render["matrix_subset"],
                metadata_subset=last_render["metadata_subset"],
                figures=last_render["figures"],
                output_dir=export_dir,
            )
            display(Markdown(f"**Exported to:** `{bundle_dir}`"))

    def _reset(_=None) -> None:
        genes_widget.value = "EPCAM, VIM, MKI67"
        scale_widget.value = "TPM"
        sample_mode_widget.value = "All samples"
        project_group_widget.value = tuple(project_groups)
        condition_widget.value = tuple(conditions)
        line_widget.value = tuple(lines[: min(8, len(lines))])
        sample_numbers_widget.value = "29, 30, 31, 32, 33, 34, 35, 36, 37, 38"
        aggregate_widget.value = False
        groupby_widget.value = "line"
        aggfunc_widget.value = "mean"
        table_rows_widget.value = 25

    sample_mode_widget.observe(_toggle_sample_mode, names="value")
    render_button.on_click(_render)
    export_button.on_click(_export)
    reset_button.on_click(_reset)

    _toggle_sample_mode({"new": sample_mode_widget.value})
    display(layout)
    _render()
