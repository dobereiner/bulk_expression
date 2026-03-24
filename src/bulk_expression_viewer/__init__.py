from .data import (
    get_scale_matrix,
    load_project_data,
    matrix_to_long,
    parse_gene_input,
    parse_sample_input,
    resolve_genes,
    summarize_selection,
)
from .exporting import export_render_bundle
from .plotting import (
    plot_aggregated_heatmap,
    plot_gene_barplot,
    plot_multi_gene_heatmap,
    plot_passage_trend,
    setup_plot_style,
)

__all__ = [
    "get_scale_matrix",
    "launch_colab_viewer",
    "load_project_data",
    "matrix_to_long",
    "parse_gene_input",
    "parse_sample_input",
    "export_render_bundle",
    "plot_aggregated_heatmap",
    "plot_gene_barplot",
    "plot_multi_gene_heatmap",
    "plot_passage_trend",
    "resolve_genes",
    "setup_plot_style",
    "summarize_selection",
]


def launch_colab_viewer(*args, **kwargs):
    from .colab_app import launch_colab_viewer as _launch_colab_viewer

    return _launch_colab_viewer(*args, **kwargs)
