import logging

import numpy as np
import pandas as pd
import plotly.express as px
from scipy.sparse import issparse

from dash import Input, Output, State, html, dash_table, dcc, no_update, ctx
import dash_bootstrap_components as dbc

from viral_platform.state.dataset_store import (
    cache_results,
    get_state_store,
    get_working_dataset,
    reset_state_store,
    update_state_store,
)
from viral_platform.analysis.pseudobulk import subset_cells
from viral_platform.analysis.differential_expression import run_differential_expression
from viral_platform.scmovir.parsers.DEGS_parser import DEGsParser
from viral_platform.scmovir.parsers.geneheatmap_parser import GeneHeatmapParser
from viral_platform.utils.sample_column_utils import resolve_sample_column
from viral_platform.utils.reference_file_utils import (
    read_downloaded_reference_file_map,
    read_downloaded_reference_filenames,
)
from viral_platform.utils.runtime_tracker import get_tracker

logger = logging.getLogger(__name__)

_de_tracker = get_tracker("Differential Expression")

MAX_DROPDOWN_CATEGORY_VALUES = 500
DE_REFERENCE_SUFFIXES = ("_degs.json", "_degs_top.json", "_gene_heatmap.json")
DEGS_SUFFIX = "_degs.json"
DEGS_TOP_SUFFIX = "_degs_top.json"
GENE_HEATMAP_SUFFIX = "_gene_heatmap.json"


def _read_downloaded_de_reference_filenames():
    return read_downloaded_reference_filenames(
        lambda filename: filename.lower().endswith(DE_REFERENCE_SUFFIXES)
    )


def _build_options(values):
    """Normalize iterable values into unique Dash dropdown option dictionaries."""
    seen = set()
    options = []
    for value in values:
        normalized = str(value)
        if normalized in seen:
            continue
        seen.add(normalized)
        options.append({"label": normalized, "value": normalized})
    return options


def _guess_cell_type_columns(groupable_columns):
    """Find likely cell-type metadata columns from available groupable columns."""
    matches = []
    for col in groupable_columns:
        lowered = str(col).lower()
        if "cell_type" in lowered or "celltype" in lowered or "cell type" in lowered:
            matches.append(col)
    return matches


def _normalize_de_result_gene_column(df, celltype_label):
    """Normalize DE result gene identifier column names to a single 'gene' field."""
    if "gene" in df.columns:
        return df

    rename_candidates = ("var_names", "index", "gene_name", "names")
    for candidate in rename_candidates:
        if candidate in df.columns:
            print(
                f"[DE heatmap debug] celltype={celltype_label} using '{candidate}' as gene column"
            )
            return df.rename(columns={candidate: "gene"})

    print(
        f"[DE heatmap debug] celltype={celltype_label} no recognizable gene column; "
        f"available={list(df.columns)}"
    )
    return df


def _normalize_reference_de_columns(df):
    """Normalize scMOVIR reference DE columns to the app's expected schema."""
    normalized = df.copy()
    if "log2FoldChange" not in normalized.columns and "logFC" in normalized.columns:
        normalized["log2FoldChange"] = normalized["logFC"]
    if "padj" not in normalized.columns:
        for candidate in ("pvals_adj", "adj_pval", "p_adj"):
            if candidate in normalized.columns:
                normalized["padj"] = normalized[candidate]
                break
    if "gene" not in normalized.columns and "names" in normalized.columns:
        normalized = normalized.rename(columns={"names": "gene"})
    return normalized


def _build_de_table_component(results_df, table_id):
    return dash_table.DataTable(
        id=table_id,
        columns=[{"name": col, "id": col} for col in results_df.columns],
        data=results_df.to_dict("records"),
        sort_action="native",
        sort_mode="single",
        page_size=20,
        style_table={"overflowX": "auto"},
    )


def _deduplicate_comparison_directions(df):
    """When both (A vs B) and (B vs A) are present, keep only the first direction."""
    if "group" not in df.columns or "reference" not in df.columns:
        return df
    seen = set()
    keep_pairs = []
    for g, r in zip(df["group"].astype(str), df["reference"].astype(str)):
        pair = (g, r)
        reverse = (r, g)
        if pair not in seen and reverse not in seen:
            seen.add(pair)
            keep_pairs.append(pair)
    keep_set = set(keep_pairs)
    mask = [
        (g, r) in keep_set
        for g, r in zip(df["group"].astype(str), df["reference"].astype(str))
    ]
    return df[mask].copy()


def _build_volcano_component(results_df, title, graph_id=None):
    volcano_df = _normalize_de_result_gene_column(results_df.copy(), title)
    volcano_df = _normalize_reference_de_columns(volcano_df)
    # Remove mirror comparison directions before plotting.
    volcano_df = _deduplicate_comparison_directions(volcano_df)
    volcano_df["log2FoldChange"] = pd.to_numeric(
        volcano_df.get("log2FoldChange"), errors="coerce"
    )
    volcano_df["padj"] = pd.to_numeric(volcano_df.get("padj"), errors="coerce")
    volcano_df = volcano_df[
        volcano_df["log2FoldChange"].notna() & volcano_df["padj"].notna()
    ].copy()

    if volcano_df.empty:
        return html.Div(f"Volcano plot for {title} is not available.")

    """
    Compute -log10(padj) for non-zero rows, then cap zero-padj rows at 5% 
    above the highest non-zero value instead of using an arbitrary floor.
    This avoids a flat line of points floating at y=300 
    """
    nonzero_mask = volcano_df["padj"] > 0
    volcano_df.loc[nonzero_mask, "neg_log10_padj"] = -np.log10(
        volcano_df.loc[nonzero_mask, "padj"]
    )
    if nonzero_mask.any():
        y_cap = float(volcano_df.loc[nonzero_mask, "neg_log10_padj"].max()) * 1.05
    else:
        y_cap = 10.0
    volcano_df["neg_log10_padj"] = volcano_df["neg_log10_padj"].fillna(y_cap)

    volcano_df["regulation"] = "not significant"
    sig_mask = (volcano_df["padj"] < 0.05) | (volcano_df["padj"] == 0)
    volcano_df.loc[sig_mask & (volcano_df["log2FoldChange"] > 0), "regulation"] = (
        "upregulated"
    )
    volcano_df.loc[sig_mask & (volcano_df["log2FoldChange"] < 0), "regulation"] = (
        "downregulated"
    )

    hover_fields = {
        "padj": ":.3e",
        "log2FoldChange": ":.3f",
        "regulation": True,
    }
    if "gene" in volcano_df.columns:
        hover_fields["gene"] = True

    fig = px.scatter(
        volcano_df,
        x="log2FoldChange",
        y="neg_log10_padj",
        color="regulation",
        color_discrete_map={
            "not significant": "grey",
            "downregulated": "blue",
            "upregulated": "red",
        },
        title=f"Volcano plot for {title}",
        hover_data=hover_fields,
    )
    fig.update_layout(template="plotly_white", legend_title_text="", yaxis_title="-log10(padj)")
    return dcc.Graph(id=graph_id, figure=fig)


def _build_top_genes_heatmap_from_results(
    adata,
    results_df,
    grouping,
    celltype_label,
    graph_id=None,
):
    heatmap_df = _normalize_de_result_gene_column(results_df.copy(), celltype_label)
    heatmap_df = _normalize_reference_de_columns(heatmap_df)

    if "gene" not in heatmap_df.columns or "padj" not in heatmap_df.columns:
        return html.Div(f"Heatmap for {celltype_label} is not available.")

    heatmap_df["padj"] = pd.to_numeric(heatmap_df["padj"], errors="coerce")
    if "log2FoldChange" in heatmap_df.columns:
        heatmap_df["log2FoldChange"] = pd.to_numeric(
            heatmap_df["log2FoldChange"], errors="coerce"
        )
        heatmap_df["abs_log2fc"] = heatmap_df["log2FoldChange"].abs()
    else:
        heatmap_df["abs_log2fc"] = 0.0

    top_genes_df = heatmap_df[
        heatmap_df["gene"].notna() & heatmap_df["padj"].notna()
    ].copy()
    top_genes_df = top_genes_df.sort_values(
        ["padj", "abs_log2fc"], ascending=[True, False]
    )
    top_genes = top_genes_df["gene"].astype(str).tolist()[:20]
    available_top_genes = [gene for gene in top_genes if gene in adata.var_names]

    if not available_top_genes or grouping not in adata.obs.columns:
        return html.Div(f"Heatmap for {celltype_label} is not available.")

    heatmap_subset = adata[:, available_top_genes].copy()

    if issparse(heatmap_subset.X):
        heatmap_values = heatmap_subset.X.toarray()
    else:
        heatmap_values = np.asarray(heatmap_subset.X)

    heatmap_values = np.log1p(heatmap_values.astype(float))
    sample_order = np.argsort(
        heatmap_subset.obs[grouping].astype(str).to_numpy(),
        kind="stable",
    )
    heatmap_values = heatmap_values[sample_order, :]
    ordered_obs = heatmap_subset.obs.iloc[sample_order]

    gene_means = heatmap_values.mean(axis=0, keepdims=True)
    gene_stds = heatmap_values.std(axis=0, keepdims=True)
    gene_stds[gene_stds == 0] = 1.0
    heatmap_z = (heatmap_values - gene_means) / gene_stds

    sample_labels = [
        f"{grp}_{i+1}"
        for i, grp in enumerate(
            ordered_obs[grouping].astype(str)
        )
    ]

    heatmap_fig = px.imshow(
        heatmap_z.T,
        x=sample_labels,
        y=available_top_genes,
        color_continuous_scale="RdBu_r",
        aspect="auto",
        labels={
            "x": "Pseudobulk Samples",
            "y": "Top DE Genes",
            "color": "Z-score",
        },
        title=f"Top 20 DE genes heatmap for {celltype_label}",
    )
    heatmap_fig.update_layout(
        template="plotly_white",
        height=max(500, 40 * len(available_top_genes) + 180),
        margin=dict(
            l=150,
            r=30,
            t=70,
            b=120,
        ),
        xaxis_tickangle=45,
    )
    heatmap_fig.update_yaxes(
        tickmode="linear",
        dtick=1,
    )
    sample_names = ordered_obs.index.astype(str).to_numpy()

    customdata = np.tile(
        sample_names,
        (len(available_top_genes), 1)
    )

    heatmap_fig.update_traces(
        customdata=customdata,
        hovertemplate=(
            "Sample: %{customdata}<br>"
            "Gene: %{y}<br>"
            "Z-score: %{z:.2f}<extra></extra>"
        )
    )
    return dcc.Graph(id=graph_id, figure=heatmap_fig)


def _build_cross_population_heatmap(gene_heatmap_df, title, graph_id=None):
    if gene_heatmap_df is None or gene_heatmap_df.empty:
        return html.Div("Cross-cell-population heatmap is not available.")

    required = {"cluster_name", "gene_name", "expression"}
    if not required.issubset(set(gene_heatmap_df.columns)):
        return html.Div("Cross-cell-population heatmap is not available.")

    cross_df = gene_heatmap_df.copy()
    cross_df["expression"] = pd.to_numeric(cross_df["expression"], errors="coerce")
    cross_df = cross_df.dropna(subset=["cluster_name", "gene_name", "expression"])

    if cross_df.empty:
        return html.Div("Cross-cell-population heatmap is not available.")

    if "gene_order" in cross_df.columns:
        cross_df["gene_order"] = pd.to_numeric(cross_df["gene_order"], errors="coerce")
        gene_order_df = (
            cross_df.dropna(subset=["gene_order"])
            .sort_values("gene_order")
            .drop_duplicates(subset=["gene_name"], keep="first")
        )
        gene_order = gene_order_df["gene_name"].astype(str).tolist()
    else:
        gene_order = pd.unique(cross_df["gene_name"].astype(str)).tolist()

    """
    Keep the cross-population heatmap consistent with the runtime DE heatmap
    and prevent unusually large gene lists from exploding the plot.
    """
    gene_order = gene_order[:20]

    if "cluster_order" in cross_df.columns:
        cross_df["cluster_order"] = pd.to_numeric(
            cross_df["cluster_order"], errors="coerce"
        )
        cluster_order_df = (
            cross_df.dropna(subset=["cluster_order"])
            .sort_values("cluster_order")
            .drop_duplicates(subset=["cluster_name"], keep="first")
        )
        cluster_order = cluster_order_df["cluster_name"].astype(str).tolist()
    else:
        cluster_order = pd.unique(cross_df["cluster_name"].astype(str)).tolist()

    matrix_df = cross_df.pivot_table(
        index="gene_name",
        columns="cluster_name",
        values="expression",
        aggfunc="mean",
    ).fillna(0.0)

    matrix_df = matrix_df.reindex(index=gene_order, columns=cluster_order)

    if matrix_df.empty:
        return html.Div("Cross-cell-population heatmap is not available.")

    matrix_values = matrix_df.to_numpy(dtype=float)
    vmax = float(np.nanpercentile(matrix_values, 99)) if matrix_values.size else 0.0
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0

    heatmap_height = max(600, min(1400, int(30 * len(matrix_df.index)) + 180))

    fig = px.imshow(
        matrix_values,
        x=matrix_df.columns.tolist(),
        y=matrix_df.index.tolist(),
        color_continuous_scale="YlOrRd",
        zmin=0,
        zmax=vmax,
        aspect="auto",
        labels={
            "x": "Cell populations",
            "y": "Top DE genes",
            "color": "Expression",
        },
        title=(
            "Heatmap of Top 20 Genes across Cell Populations"
            "<br><sup>Values represent averaged log10(CP10K + 1) normalized expression</sup>"
        ),
    )
    fig.update_layout(
        template="plotly_white",
        height=heatmap_height,
        xaxis_tickangle=45,
        margin=dict(l=150, r=30, t=90, b=120),
    )
    fig.update_yaxes(
        tickmode="linear",
        dtick=1,
    )
    return dcc.Graph(id=graph_id, figure=fig)


def _infer_reference_prefix(filename):
    lowered = filename.lower()
    for suffix in DE_REFERENCE_SUFFIXES:
        if lowered.endswith(suffix):
            return filename[: -len(suffix)]
    return None


def _find_case_insensitive_file(file_map, expected_filename):
    expected_lower = expected_filename.lower()
    for name, path in file_map.items():
        if name.lower() == expected_lower:
            return name, path
    return None, None


def _collect_reference_de_files(selected_filename):
    file_map = read_downloaded_reference_file_map(
        lambda filename: filename.lower().endswith(DE_REFERENCE_SUFFIXES)
    )
    if not file_map:
        return None

    selected_name, selected_path = _find_case_insensitive_file(
        file_map, selected_filename
    )
    if selected_name is None:
        return None

    prefix = _infer_reference_prefix(selected_name)
    if prefix is None:
        return None

    degs_name, degs_path = _find_case_insensitive_file(
        file_map, f"{prefix}{DEGS_SUFFIX}"
    )
    degs_top_name, degs_top_path = _find_case_insensitive_file(
        file_map, f"{prefix}{DEGS_TOP_SUFFIX}"
    )
    gene_heatmap_name, gene_heatmap_path = _find_case_insensitive_file(
        file_map, f"{prefix}{GENE_HEATMAP_SUFFIX}"
    )

    # If selected is a DE file but sibling was not found by convention, still use selected file.
    if degs_path is None and selected_name.lower().endswith(DEGS_SUFFIX):
        degs_name, degs_path = selected_name, selected_path
    if degs_top_path is None and selected_name.lower().endswith(DEGS_TOP_SUFFIX):
        degs_top_name, degs_top_path = selected_name, selected_path
    if gene_heatmap_path is None and selected_name.lower().endswith(
        GENE_HEATMAP_SUFFIX
    ):
        gene_heatmap_name, gene_heatmap_path = selected_name, selected_path

    return {
        "degs": (degs_name, degs_path),
        "degs_top": (degs_top_name, degs_top_path),
        "gene_heatmap": (gene_heatmap_name, gene_heatmap_path),
        "selected": (selected_name, selected_path),
        "prefix": prefix,
    }


def _build_cross_population_heatmap_from_runtime_de(
    de_results_by_celltype,
    adata,
    celltype_column,
    graph_id=None,
):
    if not de_results_by_celltype:
        return html.Div("Cross-cell-population heatmap is not available.")
    if adata is None or celltype_column not in adata.obs.columns:
        return html.Div("Cross-cell-population heatmap is not available.")

    selected_genes = []
    for _, results_df in de_results_by_celltype.items():
        normalized_df = _normalize_reference_de_columns(
            _normalize_de_result_gene_column(results_df.copy(), "All Cells")
        )
        if "gene" not in normalized_df.columns or "padj" not in normalized_df.columns:
            continue
        normalized_df["padj"] = pd.to_numeric(normalized_df["padj"], errors="coerce")
        if "log2FoldChange" in normalized_df.columns:
            normalized_df["log2FoldChange"] = pd.to_numeric(
                normalized_df["log2FoldChange"], errors="coerce"
            )
            normalized_df["abs_log2fc"] = normalized_df["log2FoldChange"].abs()
        else:
            normalized_df["abs_log2fc"] = 0.0
        ranked = normalized_df.dropna(subset=["gene", "padj"]).sort_values(
            ["padj", "abs_log2fc"],
            ascending=[True, False],
        )
        selected_genes.extend(ranked["gene"].astype(str).tolist()[:4])

    selected_genes = list(dict.fromkeys(selected_genes))[:20]
    selected_genes = [gene for gene in selected_genes if gene in adata.var_names]

    if not selected_genes:
        return html.Div("Cross-cell-population heatmap is not available.")

    gene_indices = [int(adata.var_names.get_loc(gene)) for gene in selected_genes]
    source_labels = adata.obs[celltype_column].astype(str)
    cell_populations = [
        value for value in source_labels.unique().tolist() if value.strip()
    ]
    if not cell_populations:
        return html.Div("Cross-cell-population heatmap is not available.")

    layer_matrix = (
        adata.layers["log_normalized"] if "log_normalized" in adata.layers else adata.X
    )
    if issparse(layer_matrix):
        expression_matrix = layer_matrix[:, gene_indices].toarray().astype(float)
    else:
        expression_matrix = np.asarray(layer_matrix[:, gene_indices], dtype=float)

    if "log_normalized" not in adata.layers:
        expression_matrix = np.log1p(expression_matrix)

    averaged_rows = []
    for population in cell_populations:
        mask = source_labels.to_numpy() == population
        if not np.any(mask):
            continue
        averaged_rows.append(expression_matrix[mask].mean(axis=0))

    if not averaged_rows:
        return html.Div("Cross-cell-population heatmap is not available.")

    averaged_matrix = np.vstack(averaged_rows)
    fig = px.imshow(
        averaged_matrix.T,
        x=cell_populations,
        y=selected_genes,
        color_continuous_scale="YlOrRd",
        aspect="auto",
        labels={
            "x": "Cell populations",
            "y": "Top DE genes",
            "color": "Average expression",
        },
        title="Heatmap of top genes across cell populations",
    )
    fig.update_layout(template="plotly_white", xaxis_tickangle=45)
    return dcc.Graph(id=graph_id, figure=fig)


def register_differential_expression_callbacks(app):
    """Register callbacks that keep DE dropdown options synchronized with dataset metadata/state."""

    @app.callback(
        Output("differential-expression-advanced-options-collapse", "is_open"),
        Input("differential-expression-advanced-options-button", "n_clicks"),
        State("differential-expression-advanced-options-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_de_advanced_options(n_clicks, is_open):
        if not n_clicks:
            return is_open
        return not is_open

    @app.callback(
        Output("de-min-psbulk-cells-slider", "value"),
        Output("de-min-psbulk-counts-slider", "value"),
        Output("de-min-gene-count-slider", "value"),
        Output("de-min-samples-per-gene-slider", "value"),
        Input("de-advanced-reset-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_de_advanced_defaults(n_clicks):
        if not n_clicks:
            return no_update, no_update, no_update, no_update
        return 10, 1000, 10, 2

    @app.callback(
        Output("de-reference-file-radio", "options"),
        Output("de-reference-file-radio", "value"),
        Input("scmovir-refresh-token", "data"),
        prevent_initial_call=False,
    )
    def populate_de_reference_files(_refresh_token):
        print("Refreshing DE reference files")
        filenames = _read_downloaded_de_reference_filenames()
        if not filenames:
            return (
                [{"label": "No downloaded DE reference files found.", "value": ""}],
                "",
            )

        options = [{"label": name, "value": name} for name in filenames]
        return options, options[0]["value"]

    @app.callback(
        Output("de-reference-results-container", "children"),
        Input("de-reference-display-button", "n_clicks"),
        Input("de-reference-file-radio", "value"),
        State("de-reference-display-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def display_reference_de_results(_button_click, selected_filename, n_clicks):
        trigger_id = ctx.triggered_id
        if trigger_id == "de-reference-file-radio" and not n_clicks:
            return no_update
        if not n_clicks:
            return no_update

        if not selected_filename:
            return dbc.Alert("Select a reference file first.", color="warning")

        reference_files = _collect_reference_de_files(selected_filename)
        if reference_files is None:
            return dbc.Alert(
                "Selected reference file is not available locally.", color="danger"
            )

        degs_name, degs_path = reference_files["degs"]
        degs_top_name, degs_top_path = reference_files["degs_top"]
        gene_heatmap_name, gene_heatmap_path = reference_files["gene_heatmap"]
        selected_name, selected_path = reference_files["selected"]

        reference_results_df = None
        selected_lower = selected_name.lower()

        try:
            if selected_lower.endswith(DEGS_TOP_SUFFIX) and selected_path is not None:
                reference_results_df = _normalize_reference_de_columns(
                    DEGsParser.parse(selected_path)
                )
            elif selected_lower.endswith(DEGS_SUFFIX) and selected_path is not None:
                reference_results_df = _normalize_reference_de_columns(
                    DEGsParser.parse(selected_path)
                )
        except Exception as exc:
            logger.exception(
                "Failed to parse DE reference file '%s': %s", selected_name, exc
            )
            reference_results_df = None

        if reference_results_df is None or reference_results_df.empty:
            de_table_output = html.Div(
                "No DE table data available for the selected reference file."
            )
        else:
            # Deduplicate mirror comparison directions before rendering.
            display_df = _deduplicate_comparison_directions(reference_results_df)
            de_table_output = _build_de_table_component(
                display_df.reset_index(drop=True),
                "de-reference-results-table",
            )

        gene_heatmap_df = None
        if selected_lower.endswith(GENE_HEATMAP_SUFFIX) and selected_path is not None:
            try:
                gene_heatmap_df = GeneHeatmapParser.parse(selected_path)
                gene_heatmap_name = selected_name
            except Exception as exc:
                logger.exception(
                    "Failed to parse gene heatmap reference file '%s': %s",
                    selected_name,
                    exc,
                )
        elif gene_heatmap_path is not None:
            try:
                gene_heatmap_df = GeneHeatmapParser.parse(gene_heatmap_path)
            except Exception as exc:
                logger.exception(
                    "Failed to parse gene heatmap reference file '%s': %s",
                    gene_heatmap_name,
                    exc,
                )

        cross_population_heatmap = _build_cross_population_heatmap(
            gene_heatmap_df,
            "Heatmap of top genes across cell populations",
            graph_id={"type": "de-reference-cross-population-heatmap"},
        )

        heatmap_output = dbc.Accordion(
            [
                dbc.AccordionItem(
                    cross_population_heatmap,
                    title="Heatmap of top genes across cell populations",
                    item_id="reference-cross-cell-population-heatmap",
                )
            ],
            start_collapsed=True,
            always_open=True,
            flush=True,
            id="de-reference-heatmap-accordion",
        )

        table_source = (
            selected_name if reference_results_df is not None else "Not available"
        )
        heatmap_source = (
            gene_heatmap_name if gene_heatmap_df is not None else "Not available"
        )

        result = html.Div(
            [
                dbc.Alert(
                    f"Loaded reference outputs for {selected_name}.",
                    color="success",
                    className="mb-3",
                ),
                html.P(
                    f"DE table source: {table_source}",
                    style={"marginBottom": "4px"},
                ),
                html.P(
                    f"Cross-population heatmap source: {heatmap_source}",
                    style={"color": "#6c757d", "marginBottom": "12px"},
                ),
                html.H5("DE Table"),
                de_table_output,
                html.Hr(),
                html.H5("Heatmap"),
                heatmap_output,
            ]
        )
        cache_results(**{"de-reference-results-container": result})
        return result

    @app.callback(
        Output("grouping-variable-dropdown", "options"),
        Output("grouping-variable-dropdown", "value"),
        Output("celltype-variable-dropdown", "options"),
        Output("celltype-variable-dropdown", "value"),
        Output("celltype-dropdown", "options"),
        Output("celltype-dropdown", "value"),
        Input("active-dataset-version", "data"),
    )
    def populate_metadata_driven_dropdowns(_dataset_version):
        """Populate grouping-variable and celltype dropdowns from current metadata."""
        if _dataset_version is None:
            # Browser was refreshed or is a new session with no upload — clear any
            # stale server-side state left over from a previous upload in this process.
            reset_state_store()

        state = get_state_store()
        metadata_info = state.get("metadata_info", {})

        adata = get_working_dataset()
        groupable_columns = list(metadata_info.get("groupable_columns", []))
        if adata is not None:
            groupable_columns = [
                column for column in groupable_columns if column in adata.obs.columns
            ]
        grouping_options = _build_options(groupable_columns)
        grouping_value = grouping_options[0]["value"] if grouping_options else None

        # Allow any groupable metadata column to be used as the cell-type field.
        # Some datasets use non-standard names for cell type annotations.
        celltype_variable_options = _build_options(groupable_columns)
        likely_celltype_columns = _guess_cell_type_columns(groupable_columns)
        if likely_celltype_columns:
            celltype_variable_value = str(likely_celltype_columns[0])
        else:
            celltype_variable_value = (
                celltype_variable_options[0]["value"]
                if celltype_variable_options
                else None
            )

        # Values for this dropdown are populated by the celltype-column callback.
        celltype_options = [
            {"label": "Select a cell type column to load cell types", "value": ""}
        ]
        celltype_value = ""

        if not grouping_options:
            grouping_options = [
                {"label": "Upload a dataset to select grouping variable", "value": ""}
            ]
            grouping_value = ""

        if not celltype_variable_options:
            celltype_variable_options = [
                {"label": "Upload a dataset to select cell type column", "value": ""}
            ]
            celltype_variable_value = ""

        logger.info(
            "Populated DE metadata dropdowns with %s grouping columns and %s cell-type columns.",
            len(grouping_options),
            len(celltype_variable_options),
        )

        return (
            grouping_options,
            grouping_value,
            celltype_variable_options,
            celltype_variable_value,
            celltype_options,
            celltype_value,
        )

    @app.callback(
        Output("ccc-grouping-dropdown", "options"),
        Output("ccc-grouping-dropdown", "value"),
        Input("active-dataset-version", "data"),
    )
    def populate_ccc_grouping_dropdown(_dataset_version):
        """Populate the CCC grouping dropdown from current metadata state."""
        if _dataset_version is None:
            reset_state_store()

        state = get_state_store()
        metadata_info = state.get("metadata_info", {})
        adata = get_working_dataset()
        groupable_columns = list(metadata_info.get("groupable_columns", []))
        if adata is not None:
            groupable_columns = [
                column for column in groupable_columns if column in adata.obs.columns
            ]
        grouping_options = _build_options(groupable_columns)

        if not grouping_options:
            return [
                {"label": "Upload a dataset to select grouping variable", "value": ""}
            ], ""

        return grouping_options, grouping_options[0]["value"]

    @app.callback(
        Output("celltype-dropdown", "options", allow_duplicate=True),
        Output("celltype-dropdown", "value", allow_duplicate=True),
        Input("celltype-variable-dropdown", "value"),
        Input("active-dataset-version", "data"),
        prevent_initial_call=True,
    )
    def populate_celltype_values(celltype_column, _dataset_version):
        """Populate cell-type values based on the selected cell-type metadata column."""
        adata = get_working_dataset()
        if adata is None:
            return (
                [{"label": "Upload a dataset to select cell type", "value": ""}],
                "",
            )

        if not celltype_column or celltype_column not in adata.obs.columns:
            return (
                [
                    {
                        "label": "Select a cell type column to load cell types",
                        "value": "",
                    }
                ],
                "",
            )

        unique_values = (
            adata.obs[celltype_column].dropna().astype(str).unique().tolist()
        )
        unique_values = [v for v in unique_values if v.strip()]

        if len(unique_values) > MAX_DROPDOWN_CATEGORY_VALUES:
            logger.warning(
                "Skipping DE cell-type option expansion for '%s': %d values exceed max %d.",
                celltype_column,
                len(unique_values),
                MAX_DROPDOWN_CATEGORY_VALUES,
            )
            return (
                [
                    {
                        "label": (
                            f"Selected column has {len(unique_values)} distinct values. "
                            "Choose a lower-cardinality annotation column."
                        ),
                        "value": "",
                    }
                ],
                "",
            )

        ordered_values = ["All Cells"] + [
            v for v in unique_values if v.lower() != "all cells"
        ]
        options = _build_options(ordered_values)

        if not options:
            return (
                [{"label": "No cell types found for selected column", "value": ""}],
                "",
            )

        logger.info(
            "Populated DE cell-type dropdown for '%s' with %s values.",
            celltype_column,
            len(options),
        )
        return options, options[0]["value"]

    @app.callback(
        Output("group1-dropdown", "options"),
        Output("group1-dropdown", "value"),
        Output("group2-dropdown", "options"),
        Output("group2-dropdown", "value"),
        Input("grouping-variable-dropdown", "value"),
        Input("active-dataset-version", "data"),
    )
    def populate_group_comparison_dropdowns(grouping_column, _dataset_version):
        """Populate group1/group2 dropdowns with unique values from the selected grouping column."""
        adata = get_working_dataset()
        if adata is None:
            return (
                [
                    {
                        "label": "Select a grouping variable to select group 1",
                        "value": "",
                    }
                ],
                "",
                [
                    {
                        "label": "Select a grouping variable to select group 2",
                        "value": "",
                    }
                ],
                "",
            )

        if not grouping_column or grouping_column not in adata.obs.columns:
            return (
                [
                    {
                        "label": "Select a grouping variable to select group 1",
                        "value": "",
                    }
                ],
                "",
                [
                    {
                        "label": "Select a grouping variable to select group 2",
                        "value": "",
                    }
                ],
                "",
            )

        unique_values = (
            adata.obs[grouping_column].dropna().astype(str).unique().tolist()
        )
        unique_values = [v for v in unique_values if v.strip()]

        if len(unique_values) > MAX_DROPDOWN_CATEGORY_VALUES:
            logger.warning(
                "Skipping DE group option expansion for '%s': %d values exceed max %d.",
                grouping_column,
                len(unique_values),
                MAX_DROPDOWN_CATEGORY_VALUES,
            )
            warning_option = [
                {
                    "label": (
                        f"Selected column has {len(unique_values)} distinct values. "
                        "Choose a lower-cardinality grouping column."
                    ),
                    "value": "",
                }
            ]
            return warning_option, "", warning_option, ""

        options = _build_options(unique_values)

        if not options:
            return (
                [{"label": "No groups found for selected variable", "value": ""}],
                "",
                [{"label": "No groups found for selected variable", "value": ""}],
                "",
            )

        group1_value = options[0]["value"]
        group2_value = options[1]["value"] if len(options) > 1 else options[0]["value"]

        logger.info(
            "Populated DE group dropdowns for '%s' with %s group values.",
            grouping_column,
            len(options),
        )
        return options, group1_value, options, group2_value

    @app.callback(
        Output("differential-expression-loading-signal", "children"),
        Output("pseudobulk-container", "children"),
        Output("de-table-container", "children"),
        Output("volcano-plot-container", "children"),
        Output("de-heatmap-container", "children"),
        Input("run-differential-expression-analysis-button", "n_clicks"),
        State("grouping-variable-dropdown", "value"),
        State("group1-dropdown", "value"),
        State("group2-dropdown", "value"),
        State("celltype-variable-dropdown", "value"),
        State("celltype-dropdown", "value"),
        State("de-min-psbulk-cells-slider", "value"),
        State("de-min-psbulk-counts-slider", "value"),
        State("de-min-gene-count-slider", "value"),
        State("de-min-samples-per-gene-slider", "value"),
    )
    def run_DE_analysis(
        n_clicks,
        grouping,
        group1,
        group2,
        celltype_column,
        celltype,
        min_psbulk_cells,
        min_psbulk_counts,
        min_gene_count,
        min_samples_per_gene,
    ):
        """Run differential expression analysis when the button is clicked."""
        if n_clicks == 0:
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
            )

        _de_tracker.begin()

        adata_for_values = get_working_dataset()
        if adata_for_values is None:
            _de_tracker.finish()
            return (
                "Upload a dataset to run differential expression analysis.",
                "Upload a dataset to run differential expression analysis.",
                "Upload a dataset to run differential expression analysis.",
                "Upload a dataset to run differential expression analysis.",
                "Upload a dataset to run differential expression analysis.",
            )

        if not celltype_column or celltype_column not in adata_for_values.obs.columns:
            _de_tracker.finish()
            return (
                "Select a valid cell type column before running differential expression.",
                "Select a valid cell type column before running differential expression.",
                "Select a valid cell type column before running differential expression.",
                "Select a valid cell type column before running differential expression.",
                "Select a valid cell type column before running differential expression.",
            )

        # Build the analysis target list: a single selected cell type, or every
        # value in the selected cell-type column when user selects "All Cells".
        if celltype != "All Cells":
            target_celltypes = [celltype]
        else:
            target_celltypes = [
                ct
                for ct in adata_for_values.obs[celltype_column]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
                if str(ct).strip().lower() != "all cells" and str(ct).strip() != ""
            ]

        if not target_celltypes:
            _de_tracker.finish()
            return (
                "No cell types available for differential expression analysis.",
                "No cell types available for differential expression analysis.",
                "No cell types available for differential expression analysis.",
                "No cell types available for differential expression analysis.",
                "No cell types available for differential expression analysis.",
            )

        # Keep unique order in case metadata_info has duplicates.
        target_celltypes = list(dict.fromkeys(target_celltypes))

        state = get_state_store()
        metadata_info = state.get("metadata_info", {})
        metadata_sample_columns = metadata_info.get("sample_columns", [])

        pseudobulk_items = []
        de_items = []
        volcano_items = []
        heatmap_items = []
        de_results_by_celltype = {}
        completed = 0

        # Run the full DE pipeline independently for each target cell type so
        # each section can be rendered as its own collapsible output.
        for ct in target_celltypes:
            print(f"Running DE for {ct}")
            with _de_tracker.phase("analysis"):
                adata = subset_cells(grouping, group1, group2, ct, celltype_column)
                if adata is not None:
                    adata, results = run_differential_expression(
                        adata,
                        grouping,
                        group1,
                        group2,
                        ct,
                        min_psbulk_cells=min_psbulk_cells if min_psbulk_cells is not None else 10,
                        min_psbulk_counts=min_psbulk_counts if min_psbulk_counts is not None else 1000,
                        min_gene_count=min_gene_count if min_gene_count is not None else 10,
                        min_samples_per_gene=min_samples_per_gene if min_samples_per_gene is not None else 2,
                    )
            if adata is None:
                logger.warning(
                    "Skipping DE for cell type '%s' due to missing/invalid subset.", ct
                )
                continue
            if adata is None or not hasattr(adata, "obs"):
                logger.warning(
                    "Skipping DE for cell type '%s' due to failed DE analysis.", ct
                )
                continue

            completed += 1
            if hasattr(results, "reset_index"):
                de_results_by_celltype[ct] = results.reset_index().copy()

            if grouping in adata.obs.columns:
                sample_id_col = resolve_sample_column(
                    adata.obs.columns,
                    metadata_sample_columns=metadata_sample_columns,
                    obs_df=adata.obs,
                )
                if sample_id_col is not None:
                    group1_n = adata.obs[adata.obs[grouping] == group1][
                        sample_id_col
                    ].nunique()
                    group2_n = adata.obs[adata.obs[grouping] == group2][
                        sample_id_col
                    ].nunique()
                else:
                    group1_n = "N/A"
                    group2_n = "N/A"
            else:
                group1_n = "N/A"
                group2_n = "N/A"

            median_cells = (
                int(adata.obs["psbulk_cells"].median())
                if "psbulk_cells" in adata.obs.columns
                else "N/A"
            )
            median_counts = (
                int(adata.obs["psbulk_counts"].median())
                if "psbulk_counts" in adata.obs.columns
                else "N/A"
            )
            min_cells = (
                int(adata.obs["psbulk_cells"].min())
                if "psbulk_cells" in adata.obs.columns
                else "N/A"
            )
            min_counts = (
                int(adata.obs["psbulk_counts"].min())
                if "psbulk_counts" in adata.obs.columns
                else "N/A"
            )

            # Per-celltype pseudobulk summary table.
            pseudobulk_summary = dbc.Table(
                [
                    html.Tbody(
                        [
                            html.Tr([html.Td("Grouping Variable"), html.Td(grouping)]),
                            html.Tr([html.Td("Group 1"), html.Td(group1)]),
                            html.Tr([html.Td("Group 2"), html.Td(group2)]),
                            html.Tr([html.Td("Cell Type Filter"), html.Td(ct)]),
                            html.Tr(
                                [html.Td("Biological Samples"), html.Td(adata.n_obs)]
                            ),
                            html.Tr([html.Td("Samples (Group 1)"), html.Td(group1_n)]),
                            html.Tr([html.Td("Samples (Group 2)"), html.Td(group2_n)]),
                            html.Tr(
                                [
                                    html.Td("Median Cells / Sample"),
                                    html.Td(median_cells),
                                ]
                            ),
                            html.Tr(
                                [
                                    html.Td("Median Counts / Sample"),
                                    html.Td(median_counts),
                                ]
                            ),
                            html.Tr(
                                [html.Td("Minimum Cells / Sample"), html.Td(min_cells)]
                            ),
                            html.Tr(
                                [
                                    html.Td("Minimum Counts / Sample"),
                                    html.Td(min_counts),
                                ]
                            ),
                        ]
                    )
                ]
            )

            pseudobulk_items.append(
                dbc.AccordionItem(
                    pseudobulk_summary,
                    title=f"Pseudobulk summary for {ct}",
                    item_id=f"pseudobulk-{completed}",
                )
            )

            # Per-celltype DE result table with native sortable columns.
            if hasattr(results, "reset_index"):
                results_table_df = results.reset_index()
                with _de_tracker.phase("processing"):
                    de_table_content = _build_de_table_component(
                        results_table_df,
                        {"type": "de-results-table", "celltype": str(ct)},
                    )
            else:
                de_table_content = html.Div(
                    f"DE results for {ct} are not available in table format."
                )

            de_items.append(
                dbc.AccordionItem(
                    de_table_content,
                    title=f"DE results for {ct}",
                    item_id=f"de-{completed}",
                )
            )

            # Convert the DE results into volcano-plot inputs and label each point
            # by significance / direction so the plot can be color-coded clearly.
            # Build a volcano plot for this cell type using DE outputs.
            volcano_content = html.Div(f"Volcano plot for {ct} is not available.")
            if hasattr(results, "reset_index"):
                with _de_tracker.phase("processing"):
                    volcano_content = _build_volcano_component(
                        results.reset_index().copy(),
                        ct,
                        graph_id={"type": "de-volcano-plot", "celltype": str(ct)},
                    )

            volcano_items.append(
                dbc.AccordionItem(
                    volcano_content,
                    title=f"Volcano plot for {ct}",
                    item_id=f"volcano-{completed}",
                )
            )

            # Build a top-20 DE gene heatmap for this cell type.
            heatmap_content = html.Div(f"Heatmap for {ct} is not available.")
            if hasattr(results, "reset_index"):
                with _de_tracker.phase("processing"):
                    heatmap_content = _build_top_genes_heatmap_from_results(
                        adata,
                        results.reset_index().copy(),
                        grouping,
                        ct,
                        graph_id={"type": "de-heatmap-plot", "celltype": str(ct)},
                    )

            heatmap_items.append(
                dbc.AccordionItem(
                    heatmap_content,
                    title=f"Heatmap for {ct}",
                    item_id=f"heatmap-{completed}",
                )
            )

        if celltype == "All Cells":
            cross_population_heatmap = _build_cross_population_heatmap_from_runtime_de(
                de_results_by_celltype,
                adata_for_values,
                celltype_column,
                graph_id={
                    "type": "de-cross-population-heatmap",
                    "celltype": "All_Cells",
                },
            )
            heatmap_items.append(
                dbc.AccordionItem(
                    cross_population_heatmap,
                    title="Heatmap of top genes across cell populations",
                    item_id="cross-cell-population-heatmap",
                )
            )

        if completed == 0:
            _de_tracker.finish()
            return "No valid cell type analyses could be completed.", "No valid cell type analyses could be completed.", "No valid cell type analyses could be completed.", "No valid cell type analyses could be completed.", "No valid cell type analyses could be completed."

        logger.info(
            "Running differential expression analysis for grouping '%s', comparing '%s' vs '%s', filtered by cell type '%s'.",
            grouping,
            group1,
            group2,
            celltype,
        )

        de_results = f"Differential expression analysis completed successfully for {completed} cell type(s)."
        # Wrap each result type in its own accordion so every cell type remains
        # collapsed by default and the output stays compact.
        with _de_tracker.phase("visualization"):
            pseudobulk_output = dbc.Accordion(
                pseudobulk_items,
                start_collapsed=True,
                always_open=True,
                flush=True,
                id="pseudobulk-accordion",
            )
            de_output = dbc.Accordion(
                de_items,
                start_collapsed=True,
                always_open=True,
                flush=True,
                id="de-results-accordion",
            )
            volcano_output = dbc.Accordion(
                volcano_items,
                start_collapsed=True,
                always_open=True,
                flush=True,
                id="volcano-results-accordion",
            )
            heatmap_output = dbc.Accordion(
                heatmap_items,
                start_collapsed=True,
                always_open=True,
                flush=True,
                id="de-heatmap-accordion",
            )

            cache_results(**{
                "pseudobulk-container": pseudobulk_output,
                "de-table-container": de_output,
                "volcano-plot-container": volcano_output,
                "de-heatmap-container": heatmap_output,
            })


        # For each cell type dataframe, normalize the gene column and replace the updated frame in the dictionary (for pathway enrichment later)
        for celltype, results_df in de_results_by_celltype.items():
            normalized_df = _normalize_reference_de_columns(
                _normalize_de_result_gene_column(results_df.copy(), celltype)
            )
            de_results_by_celltype[celltype] = normalized_df
        
        update_state_store(**{"DE_results": {"results_by_celltype": de_results_by_celltype}})
        _de_tracker.finish()
        return de_results, pseudobulk_output, de_output, volcano_output, heatmap_output
