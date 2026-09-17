from dash import Input, Output, State, ctx, dash_table, html, no_update, dcc
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
from viral_platform.state.dataset_store import (
    cache_results,
    get_working_dataset,
    get_state_store,
    update_state_store,
)
from viral_platform.analysis.CellCellLIANA import (
    run_liana,
    filter_liana_results,
    liana_output_table,
    summarise_celltype_interactions,
)
from viral_platform.scmovir.parsers.cellphonedb_parser import (
    CellPhoneDBParser,
    CellPhoneDBDCCParser,
)
from viral_platform.utils.reference_file_utils import (
    read_downloaded_reference_file_map,
    read_downloaded_reference_filenames,
)
import logging
import math
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import networkx as nx
from viral_platform.utils.runtime_tracker import get_tracker

logger = logging.getLogger(__name__)

_ccc_tracker = get_tracker("Cell-Cell Communication")

MAX_DROPDOWN_CATEGORY_VALUES = 500  # Maximum number of unique values allowed in the source/target dropdowns to prevent performance issues.
DEFAULT_CCC_RENDER_ROWS = (
    1000  # Default number of rows to render in the CCC results table and bubble plot.
)
MAX_CCC_RENDER_ROWS = 5000  # Maximum number of rows to render in the CCC results table and bubble plot to prevent performance issues.


def _read_downloaded_cellphonedb_filenames():
    return read_downloaded_reference_filenames(
        lambda filename: "cellphonedb" in filename.lower()
    )


# CellPhoneDB reference rendering helpers
CCC_REFERENCE_SUFFIXES = (
    "_cellphonedb.json",
    "_cellphonedb_top.json",
    "_cellphonedb_dcc.csv",
)


def _load_ccc_reference_file(filename, path):
    """Parse a CellPhoneDB reference file and return (df, file_type)."""
    lowered = filename.lower()
    if lowered.endswith(".csv"):
        return CellPhoneDBDCCParser.parse(path), "dcc"
    return CellPhoneDBParser.parse(path), "json"


def _build_ccc_reference_table(df, file_type):
    if file_type == "json":
        cols = [
            c
            for c in [
                "source",
                "target",
                "ligand_complex",
                "receptor_complex",
                "lr_means",
                "cellphone_pvals",
                "condition",
                "top_status",
            ]
            if c in df.columns
        ]
    else:
        cols = list(df.columns)
    display_df = df[cols].copy()
    return dash_table.DataTable(
        id="ccc-reference-table",
        columns=[{"name": c, "id": c} for c in cols],
        data=display_df.to_dict("records"),
        sort_action="native",
        sort_mode="single",
        page_size=20,
        style_table={"overflowX": "auto"},
        style_cell={"textAlign": "left", "padding": "6px", "fontSize": "13px"},
        style_header={"fontWeight": "600"},
    )


def _build_ccc_reference_summary(df, file_type):
    """Return (summary_df, mean_magnitude_col, mean_specificity_col) for network."""
    if file_type == "json":
        if not {"source", "target", "ligand_complex", "cellphone_pvals"}.issubset(
            df.columns
        ):
            return None
        pvals = pd.to_numeric(df["cellphone_pvals"], errors="coerce")
        pval_max = pvals[pvals > 0].max() if (pvals > 0).any() else 1.0
        df = df.copy()
        df["_norm_pval"] = pvals.clip(lower=0).fillna(pval_max) / pval_max
        lr_means = pd.to_numeric(
            df.get("lr_means", pd.Series(dtype=float)), errors="coerce"
        ).fillna(0)
        lr_max = lr_means.max() if lr_means.max() > 0 else 1.0
        df["_norm_lr"] = lr_means / lr_max
        summary = (
            df.groupby(["source", "target"])
            .agg(
                interaction_count=("ligand_complex", "count"),
                mean_magnitude=("_norm_pval", "mean"),
                mean_specificity=("_norm_lr", "mean"),
            )
            .reset_index()
        )
        return summary
    else:  # dcc
        if not {"source", "target"}.issubset(df.columns):
            return None
        df = df.copy()
        lfc = pd.to_numeric(
            df.get("abs_log2FC", pd.Series(dtype=float)), errors="coerce"
        ).fillna(0)
        lfc_max = lfc.max() if lfc.max() > 0 else 1.0
        df["_norm_lfc"] = lfc / lfc_max
        summary = (
            df.groupby(["source", "target"])
            .agg(
                interaction_count=("source", "count"),
                mean_magnitude=("_norm_lfc", "mean"),
                mean_specificity=("_norm_lfc", "mean"),
            )
            .reset_index()
        )
        return summary


def _build_ccc_reference_bubble(df, file_type, title):
    if file_type == "json":
        required = {
            "source",
            "target",
            "ligand_complex",
            "receptor_complex",
            "lr_means",
            "cellphone_pvals",
        }
        if not required.issubset(df.columns):
            return html.Div("Bubble plot is not available for this reference file.")
        plot_df = df.copy()
        plot_df["interaction"] = (
            plot_df["ligand_complex"].astype(str)
            + " -> "
            + plot_df["receptor_complex"].astype(str)
        )
        plot_df["lr_means"] = pd.to_numeric(
            plot_df["lr_means"], errors="coerce"
        ).fillna(0)
        plot_df["cellphone_pvals"] = pd.to_numeric(
            plot_df["cellphone_pvals"], errors="coerce"
        ).clip(lower=1e-300)
        plot_df["neg_log10_pval"] = -np.log10(plot_df["cellphone_pvals"])
        # Cap display rows for performance
        plot_df = plot_df.sort_values("cellphone_pvals", ascending=True).head(1000)
        fig = px.scatter(
            plot_df,
            x="target",
            y="interaction",
            size="lr_means",
            color="neg_log10_pval",
            color_continuous_scale="RdYlBu_r",
            hover_data=[
                "source",
                "ligand_complex",
                "receptor_complex",
                "lr_means",
                "cellphone_pvals",
            ],
            title=title,
            labels={
                "neg_log10_pval": "-log10(p-value)",
                "target": "Target cell type",
                "interaction": "Interaction (Ligand → Receptor)",
            },
        )
        fig.update_layout(template="plotly_white", xaxis_tickangle=45)
        return dcc.Graph(figure=fig)
    else:  # dcc
        required = {"source", "target", "lr_pair", "abs_log2FC"}
        if not required.issubset(df.columns):
            return html.Div("Bubble plot is not available for this reference file.")
        plot_df = df.copy()
        plot_df["cell_pair"] = (
            plot_df["source"].astype(str) + " → " + plot_df["target"].astype(str)
        )
        plot_df["abs_log2FC"] = pd.to_numeric(
            plot_df["abs_log2FC"], errors="coerce"
        ).fillna(0)
        plot_df["log2FC"] = pd.to_numeric(plot_df["log2FC"], errors="coerce").fillna(0)
        fig = px.scatter(
            plot_df,
            x="cell_pair",
            y="lr_pair",
            size="abs_log2FC",
            color="log2FC",
            color_continuous_scale="RdBu_r",
            color_continuous_midpoint=0,
            hover_data=["source", "target", "ligand", "receptor", "abs_log2FC"],
            title=title,
            labels={
                "cell_pair": "Cell pair (Source → Target)",
                "lr_pair": "Interaction (Ligand → Receptor)",
                "log2FC": "log2 Fold Change",
            },
        )
        fig.update_layout(template="plotly_white", xaxis_tickangle=45)
        return dcc.Graph(figure=fig)


def _build_ccc_reference_filter_options(df):
    """Build source/target filter options from reference data if available."""
    if df is None or df.empty:
        fallback = [{"label": "All", "value": "_all_"}]
        return fallback, "_all_", fallback, "_all_"

    source_values = []
    target_values = []
    if "source" in df.columns:
        source_values = [
            str(v) for v in df["source"].dropna().astype(str).tolist() if v.strip()
        ]
    if "target" in df.columns:
        target_values = [
            str(v) for v in df["target"].dropna().astype(str).tolist() if v.strip()
        ]

    source_options = _build_options(list(dict.fromkeys(source_values)))
    target_options = _build_options(list(dict.fromkeys(target_values)))
    return source_options, "_all_", target_options, "_all_"


def _filter_ccc_reference_df(
    df, file_type, source_filter, target_filter, interaction_filter
):
    """Filter and optionally cap reference CCC rows for display components."""
    filtered = df.copy()

    if "source" in filtered.columns and source_filter and source_filter != "_all_":
        filtered = filtered[filtered["source"].astype(str) == str(source_filter)]
    if "target" in filtered.columns and target_filter and target_filter != "_all_":
        filtered = filtered[filtered["target"].astype(str) == str(target_filter)]

    render_limit = _resolve_render_limit(interaction_filter)
    display_df = filtered.copy()
    if file_type == "json":
        if "cellphone_pvals" in display_df.columns:
            display_df = display_df.assign(
                _sort_pval=pd.to_numeric(display_df["cellphone_pvals"], errors="coerce")
            ).sort_values("_sort_pval", ascending=True)
            display_df = display_df.drop(columns=["_sort_pval"])
        elif "lr_means" in display_df.columns:
            display_df = display_df.assign(
                _sort_lr=pd.to_numeric(display_df["lr_means"], errors="coerce")
            ).sort_values("_sort_lr", ascending=False)
            display_df = display_df.drop(columns=["_sort_lr"])
    else:
        if "abs_log2FC" in display_df.columns:
            display_df = display_df.assign(
                _sort_fc=pd.to_numeric(display_df["abs_log2FC"], errors="coerce")
            ).sort_values("_sort_fc", ascending=False)
            display_df = display_df.drop(columns=["_sort_fc"])

    return filtered, display_df.head(render_limit)


def _render_ccc_reference_outputs(df, file_type, selected_filename):
    """Render summary table, bubble plot, and network plot for reference CCC data."""
    table_output = _build_ccc_reference_table(df, file_type)
    bubble_output = _build_ccc_reference_bubble(
        df, file_type, f"Cell-Cell Communication - {selected_filename}"
    )

    summary_df = _build_ccc_reference_summary(df, file_type)
    if summary_df is not None and not summary_df.empty:
        network_fig = create_network_plot(summary_df)
        network_output = dcc.Graph(figure=network_fig)
    else:
        network_output = html.Div(
            "Network plot is not available for this reference file."
        )

    return html.Div(
        [
            dbc.Alert(
                f"Loaded reference CCC outputs for {selected_filename}.",
                color="success",
                className="mb-3",
            ),
            html.H5("Summary Table"),
            table_output,
            html.Hr(),
            html.H5("Bubble Plot"),
            bubble_output,
            html.Hr(),
            html.H5("Network Plot"),
            network_output,
        ]
    )


def make_sortable_table(df, table_id):
    """Create a sortable Dash DataTable from a DataFrame."""
    return dash_table.DataTable(
        id=table_id,
        columns=[{"name": col, "id": col} for col in df.columns],
        data=df.to_dict("records"),
        sort_action="native",
        page_action="native",
        page_size=20,
        style_table={"overflowX": "auto"},
        style_cell={"textAlign": "left", "padding": "6px", "fontSize": "13px"},
        style_header={"fontWeight": "600"},
    )


def _build_options(values):
    """Normalize iterable values into unique Dash dropdown option dictionaries."""
    seen = set()
    options = [
        {"label": "All", "value": "_all_"},
    ]
    for value in values:
        normalized = str(value)
        if normalized in seen:
            continue
        seen.add(normalized)
        options.append({"label": normalized, "value": normalized})
    return options


def _resolve_render_limit(requested_limit):
    """Normalize requested row limit and cap to a safe upper bound."""
    if requested_limit is None:
        return DEFAULT_CCC_RENDER_ROWS
    try:
        limit = int(requested_limit)
    except (TypeError, ValueError):
        return DEFAULT_CCC_RENDER_ROWS
    if limit <= 0:
        return DEFAULT_CCC_RENDER_ROWS
    return min(limit, MAX_CCC_RENDER_ROWS)


def _prepare_liana_for_display(liana_results):
    """Precompute CCC display columns once to avoid repeated per-click work."""
    prepared = liana_results.copy()
    prepared["interaction"] = (
        prepared["ligand_complex"].astype(str)
        + " -> "
        + prepared["receptor_complex"].astype(str)
    )
    prepared["bubble_size"] = 1 - prepared["magnitude_rank"]
    prepared["bubble_color"] = 1 - prepared["specificity_rank"]
    return prepared.sort_values(by="magnitude_rank", ascending=True)


def _display_rows(results, requested_limit=None):
    """Return only the strongest rows for components that render each interaction."""
    return results.head(_resolve_render_limit(requested_limit))


def _regular_polygon_layout(nodes):
    """Place nodes at equally spaced points around a regular polygon."""
    ordered_nodes = sorted(nodes, key=lambda value: str(value))
    count = len(ordered_nodes)

    if count == 1:
        return {ordered_nodes[0]: (0.0, 0.0)}, ordered_nodes

    radius = 1.0
    step = (2 * math.pi) / count
    start_angle = math.pi / 2

    positions = {}
    for index, node in enumerate(ordered_nodes):
        angle = start_angle - (index * step)
        positions[node] = (radius * math.cos(angle), radius * math.sin(angle))

    return positions, ordered_nodes


def _label_rotation(x, y, counter_clockwise_offset=8):
    """Point a label away from the centre, with a small counter-clockwise turn."""
    radial_angle = math.degrees(math.atan2(y, x))
    # Plotly rotates annotation text in screen coordinates, whose vertical
    # direction is the reverse of the graph's y-axis.
    return ((-radial_angle - counter_clockwise_offset + 180) % 360) - 180


def _label_position(x, y, node_size, label):
    """Place a label on an outer ring and return its position."""
    distance_from_origin = math.hypot(x, y)
    if distance_from_origin == 0:
        return 0.0, 1.58

    # Keep labels outside the graph. Longer labels receive slightly more room
    # so adjacent labels on the ring have less chance of colliding.
    radius = 1.56 + (node_size / 240) + (0.012 * len(str(label)))
    if x < 0:
        # Left-side text begins at the anchor position. Push the anchor farther
        # out so left-side labels have room and avoid touching the network.
        radius += 0.025 * len(str(label))
    return (
        radius * x / distance_from_origin,
        radius * y / distance_from_origin,
    )


def create_network_plot(summary, show_labels=True):
    """
    Create a cell-cell communication network.

    Parameters
    ----------
    summary : pd.DataFrame

    Required columns:
        source
        target
        interaction_count
        mean_magnitude
        mean_specificity
    """

    if summary is None or summary.empty:
        fig = go.Figure()
        fig.update_layout(
            title="Cell-Cell Communication Network",
            template="plotly_white",
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            annotations=[
                dict(
                    text="No interactions available for selected filters.",
                    x=0.5,
                    y=0.5,
                    xref="paper",
                    yref="paper",
                    showarrow=False,
                )
            ],
        )
        return fig

    # Build graph
    G = nx.DiGraph()

    for _, row in summary.iterrows():
        G.add_edge(
            row["source"],
            row["target"],
            interactions=row["interaction_count"],
            magnitude=1 - row["mean_magnitude"],
            specificity=1 - row["mean_specificity"],
        )

    # Layout
    pos, ordered_nodes = _regular_polygon_layout(G.nodes())

    # Node statistics
    degree = dict(G.degree())

    weighted_degree = dict(G.degree(weight="interactions"))

    # Normalise node sizes
    min_node = 12
    max_node = 30

    dmin = min(degree.values())
    dmax = max(degree.values())

    # Normalise edge widths
    min_edge = 1
    max_edge = 8

    wmin = summary["interaction_count"].min()
    wmax = summary["interaction_count"].max()

    # Figure
    fig = go.Figure()

    colours = px.colors.sequential.Turbo

    # Normalise edge colours
    strength = 1 - summary["mean_magnitude"]

    smin = strength.min()
    smax = strength.max()

    # Draw edges
    for u, v, data in G.edges(data=True):
        x0, y0 = pos[u]
        x1, y1 = pos[v]

        # Normalise colour across this dataset
        if smax == smin:
            scaled = 0.5

        else:
            scaled = (data["magnitude"] - smin) / (smax - smin)

        scaled = max(0, min(1, scaled))

        colour_index = int(scaled * (len(colours) - 1))

        hex_colour = colours[colour_index]

        r = int(hex_colour[1:3], 16)
        g = int(hex_colour[3:5], 16)
        b = int(hex_colour[5:7], 16)

        alpha = 0.35 + 0.45 * scaled  # 0.25–0.75

        colour = f"rgba({r},{g},{b},{alpha:.2f})"

        # Edge width
        if wmax == wmin:
            edge_width = (min_edge + max_edge) / 2

        else:
            edge_width = min_edge + (data["interactions"] - wmin) * (
                max_edge - min_edge
            ) / (wmax - wmin)

        # Edge trace
        fig.add_trace(
            go.Scatter(
                x=[x0, x1],
                y=[y0, y1],
                mode="lines",
                line=dict(
                    width=edge_width,
                    color=colour,
                ),
                hoverinfo="text",
                text=(
                    f"<b>{u}</b> → <b>{v}</b><br>"
                    f"Interactions: {data['interactions']}<br>"
                    f"Mean Magnitude: {data['magnitude']:.3f}<br>"
                    f"Mean Specificity: {data['specificity']:.3f}"
                ),
                showlegend=False,
            )
        )

        # Arrow
        fig.add_annotation(
            x=x1,
            y=y1,
            ax=x0,
            ay=y0,
            xref="x",
            yref="y",
            axref="x",
            ayref="y",
            showarrow=True,
            arrowhead=2,
            arrowsize=0.6,
            arrowwidth=max(1.5, edge_width * 0.6),
            arrowcolor=colour,
            text="",
        )

    # Nodes
    node_x = []
    node_y = []

    node_size = []
    node_sizes = {}

    node_colour = []

    node_text = []

    for node in ordered_nodes:
        x, y = pos[node]

        node_x.append(x)
        node_y.append(y)

        # Node size
        if dmax == dmin:
            size = (min_node + max_node) / 2

        else:
            size = min_node + (degree[node] - dmin) * (max_node - min_node) / (
                dmax - dmin
            )

        node_size.append(size)
        node_sizes[node] = size

        # Colour
        node_colour.append(weighted_degree[node])

        node_text.append(
            f"<b>{node}</b><br>"
            f"Connected Cell Types: {degree[node]}<br>"
            f"Total Interactions: {weighted_degree[node]}"
        )

    # Draw nodes
    fig.add_trace(
        go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers",
            hoverinfo="text",
            hovertext=node_text,
            marker=dict(
                size=node_size,
                color=node_colour,
                colorscale="Viridis",
                showscale=True,
                colorbar=dict(
                    title="Node<br>Connectivity",
                    x=1.02,
                ),
                line=dict(
                    width=2,
                    color="black",
                ),
            ),
            showlegend=False,
        )
    )

    # Node Labels
    if show_labels:
        for node in ordered_nodes:
            x, y = pos[node]
            rotation = _label_rotation(x, y)
            label_x, label_y = _label_position(x, y, node_sizes[node], node)

            fig.add_annotation(
                x=label_x,
                y=label_y,
                ax=x,
                ay=y,
                xref="x",
                yref="y",
                axref="x",
                ayref="y",
                text=node,
                # Keep annotation anchor behavior for stable label placement,
                # but hide the leader line itself.
                showarrow=True,
                arrowhead=0,
                arrowwidth=1,
                arrowcolor="rgba(0, 0, 0, 0)",
                standoff=14,
                font=dict(
                    size=12,
                    color="black",
                ),
                xanchor="left",
                yanchor="middle",
                align="left",
                textangle=rotation,
            )

    # Dummy trace for edge colourbar
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                colorscale="Turbo",
                cmin=smin,
                cmax=smax,
                color=[smin],
                size=0,
                showscale=True,
                colorbar=dict(
                    title="Edge<br>Strength",
                    tickvals=[
                        smin,
                        (smin + smax) / 2,
                        smax,
                    ],
                    ticktext=[
                        "Weak",
                        "Moderate",
                        "Strong",
                    ],
                    thickness=18,
                    len=0.75,
                    x=1.12,
                ),
            ),
            hoverinfo="skip",
            showlegend=False,
        )
    )
    # Layout
    fig.update_layout(
        title="Cell-Cell Communication Network",
        template="plotly_white",
        height=780,
        hovermode="closest",
        xaxis=dict(
            visible=False,
            zeroline=False,
            showgrid=False,
            range=[-2.45, 2.45],
        ),
        yaxis=dict(
            visible=False,
            zeroline=False,
            showgrid=False,
            range=[-2.45, 2.45],
        ),
        margin=dict(
            l=20,
            r=20,
            t=60,
            b=20,
        ),
    )

    return fig


def register_ccc_callbacks(app):
    @app.callback(
        Output("ccc-advanced-options-collapse", "is_open"),
        Input("ccc-advanced-options-button", "n_clicks"),
        State("ccc-advanced-options-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_ccc_advanced_options(n_clicks, is_open):
        if not n_clicks:
            return is_open
        return not is_open

    @app.callback(
        Output("ccc-max-magnitude-rank-input", "value"),
        Output("ccc-max-specificity-rank-input", "value"),
        Output("ccc-default-render-rows-input", "value"),
        Output("ccc-interaction-filter-input", "value"),
        Input("ccc-advanced-reset-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_ccc_advanced_defaults(n_clicks):
        if not n_clicks:
            return no_update, no_update, no_update, no_update
        return None, None, DEFAULT_CCC_RENDER_ROWS, DEFAULT_CCC_RENDER_ROWS

    @app.callback(
        Output("ccc-reference-file-radio", "options"),
        Output("ccc-reference-file-radio", "value"),
        Input("scmovir-refresh-token", "data"),
        prevent_initial_call=False,
    )
    def populate_ccc_reference_files(_refresh_token):
        print("Refreshing CCC reference files")
        filenames = _read_downloaded_cellphonedb_filenames()
        if not filenames:
            return (
                [{"label": "No downloaded CellPhoneDB files found.", "value": ""}],
                "",
            )

        options = [{"label": name, "value": name} for name in filenames]
        return options, options[0]["value"]

    @app.callback(
        Output("ccc-reference-results-container", "children"),
        Output("ccc-source-filter-dropdown", "options", allow_duplicate=True),
        Output("ccc-source-filter-dropdown", "value", allow_duplicate=True),
        Output("ccc-target-filter-dropdown", "options", allow_duplicate=True),
        Output("ccc-target-filter-dropdown", "value", allow_duplicate=True),
        Input("ccc-reference-select-button", "n_clicks"),
        Input("ccc-reference-file-radio", "value"),
        State("ccc-reference-select-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def display_reference_ccc_results(_button_click, selected_filename, n_clicks):
        trigger_id = ctx.triggered_id
        if trigger_id == "ccc-reference-file-radio" and not n_clicks:
            return no_update, no_update, no_update, no_update, no_update
        if not n_clicks:
            return no_update, no_update, no_update, no_update, no_update

        if not selected_filename:
            return (
                dbc.Alert("Select a reference file first.", color="warning"),
                no_update,
                no_update,
                no_update,
                no_update,
            )

        file_map = read_downloaded_reference_file_map(
            lambda fn: "cellphonedb" in fn.lower()
        )
        selected_lower = selected_filename.lower()
        file_path = None
        for name, path in file_map.items():
            if name.lower() == selected_lower:
                file_path = path
                break

        if file_path is None:
            return (
                dbc.Alert(
                    "Selected reference file is not available locally.", color="danger"
                ),
                no_update,
                no_update,
                no_update,
                no_update,
            )

        try:
            df, file_type = _load_ccc_reference_file(selected_filename, file_path)
        except Exception as exc:
            logger.exception(
                "Failed to parse CCC reference file '%s': %s", selected_filename, exc
            )
            return (
                dbc.Alert(f"Failed to parse file: {exc}", color="danger"),
                no_update,
                no_update,
                no_update,
                no_update,
            )

        update_state_store(
            CCC_reference_results={
                "results": df,
                "file_type": file_type,
                "selected_filename": selected_filename,
            },
            CCC_active_context="reference",
        )

        source_options, source_value, target_options, target_value = (
            _build_ccc_reference_filter_options(df)
        )

        return (
            _render_ccc_reference_outputs(df, file_type, selected_filename),
            source_options,
            source_value,
            target_options,
            target_value,
        )

    @app.callback(
        Output("ccc-source-filter-dropdown", "options"),
        Output("ccc-source-filter-dropdown", "value"),
        Output("ccc-target-filter-dropdown", "options"),
        Output("ccc-target-filter-dropdown", "value"),
        Input("ccc-grouping-dropdown", "value"),
        Input("active-dataset-version", "data"),
    )
    def populate_group_comparison_dropdowns(grouping_column, _dataset_version):
        """Populate source/target dropdowns with unique values from the selected grouping column."""
        adata = get_working_dataset()
        if adata is None:
            return (
                [{"label": "Select source cell type", "value": ""}],
                "",
                [{"label": "Select target cell type", "value": ""}],
                "",
            )

        if not grouping_column or grouping_column not in adata.obs.columns:
            return (
                [{"label": "Select source cell type", "value": ""}],
                "",
                [{"label": "Select target cell type", "value": ""}],
                "",
            )

        unique_values = (
            adata.obs[grouping_column].dropna().astype(str).unique().tolist()
        )
        unique_values = [v for v in unique_values if v.strip()]

        if len(unique_values) > MAX_DROPDOWN_CATEGORY_VALUES:
            logger.warning(
                "Skipping CCC source/target option expansion for '%s': %d values exceed max %d.",
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
        group2_value = options[0]["value"]

        logger.info(
            "Populated DE group dropdowns for '%s' with %s group values.",
            grouping_column,
            len(options),
        )
        return options, group1_value, options, group2_value

    @app.callback(
        Output("ccc-loading-signal", "children"),
        Output("ccc-summary-container", "children", allow_duplicate=True),
        Output("ccc-bubble-plot-container", "children", allow_duplicate=True),
        Output("ccc-network-plot-container", "children", allow_duplicate=True),
        Input("run-ccc-button", "n_clicks"),
        State("active-dataset-version", "data"),
        State("ccc-grouping-dropdown", "value"),
        State("ccc-method-dropdown", "value"),
        State("ccc-resource-dropdown", "value"),
        State("ccc-max-magnitude-rank-input", "value"),
        State("ccc-max-specificity-rank-input", "value"),
        State("ccc-default-render-rows-input", "value"),
        prevent_initial_call=True,
    )
    @_ccc_tracker.track
    def run_ccc_analysis(
        n_clicks,
        dataset_version,
        grouping_variable,
        method,
        resource,
        max_magnitude_rank,
        max_specificity_rank,
        default_render_rows,
    ):
        """
        Runs the Cell-Cell Communication (CCC) analysis using the selected parameters.

        Parameters
        ----------
        n_clicks : int
            Number of times the "Run Cell-Cell Communication" button has been clicked.
        grouping_variable : str
            The selected grouping variable for the CCC analysis.
        method : str
            The selected method for the CCC analysis.
        resource : str
            The selected resource for the CCC analysis.

        Returns
        -------
        tuple
            A tuple containing:
            - A string indicating the loading status.
            - A Dash DataTable displaying the CCC analysis results.
        """
        if not n_clicks or n_clicks == 0:
            return (
                no_update,
                html.P(
                    "No cell-cell communication results yet. Run the analysis to see results here.",
                    style={"color": "#6c757d", "fontSize": "14px", "marginTop": "10px"},
                ),
                "",
                "",
            )

        adata = get_working_dataset()
        if adata is None:
            return (
                "done",
                "No dataset available for cell-cell communication analysis.",
                "",
                "",
            )

        cache_key = {
            "dataset_version": str(dataset_version),
            "grouping_variable": str(grouping_variable),
            "method": str(method),
            "resource": str(resource),
        }

        history = get_state_store()
        cached_ccc = history.get("CCC_results", {})
        if (
            cached_ccc.get("cache_key") == cache_key
            and cached_ccc.get("results") is not None
        ):
            prepared_results = cached_ccc["results"]
            logger.info("Reusing cached CCC results for key=%s", cache_key)
        else:
            # Run the CCC analysis using the provided parameters.
            with _ccc_tracker.phase("analysis"):
                liana_results = run_liana(adata, grouping_variable, method, resource)
                prepared_results = _prepare_liana_for_display(liana_results)
            update_state_store(
                CCC_results={"results": prepared_results, "cache_key": cache_key}
            )

        update_state_store(CCC_active_context="uploaded")

        with _ccc_tracker.phase("processing"):
            filtered_results = filter_liana_results(
                prepared_results,
                max_magnitude_rank=max_magnitude_rank,
                max_specificity_rank=max_specificity_rank,
            )

            display_results = _display_rows(filtered_results, default_render_rows)

            results_table = liana_output_table(display_results)
            results_table = make_sortable_table(results_table, "ccc-results-table")
            # The summary represents every LIANA result.  Only the interaction-level
            # table is capped, because rendering thousands of individual marks is
            # what slows the browser down.
            summary = summarise_celltype_interactions(filtered_results)

            # Network plot generation
            network_fig = create_network_plot(summary)

            fig = px.scatter(
                summary,
                x="target",
                y="source",
                size="bubble_size",
                color="bubble_color",
                labels={"bubble_color": "Interaction strength"},
                hover_data=[
                    "interaction_count",
                    "mean_magnitude",
                    "mean_specificity",
                ],
            )
            fig.update_coloraxes(colorbar_title_text="Interaction<br>strength")

            fig.update_yaxes(categoryorder="category ascending")
            fig.update_layout(
                template="plotly_white",
                title="Cell-Cell Communication Bubble Plot",
            )

        with _ccc_tracker.phase("visualization"):
            bubble_graph = dcc.Graph(figure=fig)
            network_graph = dcc.Graph(figure=network_fig)
            cache_results(**{
                "ccc-summary-container": results_table,
                "ccc-bubble-plot-container": bubble_graph,
                "ccc-network-plot-container": network_graph,
            })
        return (
            "done",
            results_table,
            bubble_graph,
            network_graph,
        )

    @app.callback(
        Output("ccc-loading-signal", "children", allow_duplicate=True),
        Output("ccc-summary-container", "children", allow_duplicate=True),
        Output("ccc-bubble-plot-container", "children", allow_duplicate=True),
        Output("ccc-network-plot-container", "children", allow_duplicate=True),
        Output("ccc-reference-results-container", "children", allow_duplicate=True),
        Input("ccc-apply-filters-button", "n_clicks"),
        State("ccc-source-filter-dropdown", "value"),
        State("ccc-target-filter-dropdown", "value"),
        State("ccc-interaction-filter-input", "value"),
        State("ccc-show-network-labels", "value"),
        State("ccc-max-magnitude-rank-input", "value"),
        State("ccc-max-specificity-rank-input", "value"),
        State("ccc-default-render-rows-input", "value"),
        prevent_initial_call=True,
    )
    def apply_liana_filters(
        n_clicks,
        source_filter,
        target_filter,
        interaction_filter,
        show_labels_value,
        max_magnitude_rank,
        max_specificity_rank,
        default_render_rows,
    ):
        """
        Filters the LIANA results based on user-selected criteria.

        Parameters
        ----------
        n_clicks : int
            Number of times the "Apply Filters" button has been clicked.
        source_filter : str
            The selected source cell type for filtering.
        target_filter : str
            The selected target cell type for filtering.
        interaction_filter : int
            The input number for filtering interactions.

        Returns
        -------
        tuple
            A tuple containing:
            - A Dash DataTable displaying the filtered LIANA results.
            - A Dash Graph displaying the filtered bubble plot.
        """
        if not n_clicks or n_clicks == 0:
            return no_update, no_update, no_update, no_update, no_update

        history = get_state_store()
        active_context = history.get("CCC_active_context", "uploaded")

        if active_context == "reference":
            ref_payload = history.get("CCC_reference_results", {})
            reference_results = ref_payload.get("results")
            file_type = ref_payload.get("file_type")
            selected_filename = ref_payload.get("selected_filename", "reference file")

            if reference_results is None or file_type is None:
                return (
                    "done",
                    no_update,
                    no_update,
                    no_update,
                    dbc.Alert(
                        "No reference CCC results available to filter.",
                        color="warning",
                    ),
                )

            filtered_results, display_results = _filter_ccc_reference_df(
                reference_results,
                file_type,
                source_filter,
                target_filter,
                interaction_filter,
            )

            reference_output = _render_ccc_reference_outputs(
                display_results,
                file_type,
                selected_filename,
            )
            return "done", no_update, no_update, no_update, reference_output

        liana_results = history.get("CCC_results", {}).get("results")

        if liana_results is None:
            logger.warning("No LIANA results found in state store for filtering.")
            return (
                "done",
                html.P(
                    "No cell-cell communication results available to filter.",
                    style={"color": "#6c757d", "fontSize": "14px", "marginTop": "10px"},
                ),
                "",
                "",
                no_update,
            )

        filtered_results = filter_liana_results(
            liana_results,
            source=source_filter if source_filter != "_all_" else None,
            target=target_filter if target_filter != "_all_" else None,
            max_magnitude_rank=max_magnitude_rank,
            max_specificity_rank=max_specificity_rank,
        )

        requested_limit = interaction_filter
        if requested_limit is None:
            requested_limit = default_render_rows
        render_limit = _resolve_render_limit(requested_limit)
        display_results = _display_rows(filtered_results, render_limit)

        logger.info(
            "CCC filter render rows: requested=%s resolved=%s remaining=%s",
            interaction_filter,
            render_limit,
            len(display_results),
        )

        show_labels = True in (show_labels_value or [])
        network_fig = create_network_plot(
            summarise_celltype_interactions(filtered_results), show_labels=show_labels
        )

        results_table = liana_output_table(display_results)
        results_table = make_sortable_table(results_table, "ccc-results-table")

        if target_filter == "_all_":
            fig = px.scatter(
                display_results,
                x="target",
                y="interaction",
                size="bubble_size",
                color="bubble_color",
                labels={"bubble_color": "Interaction strength"},
                hover_data=[
                    "ligand_complex",
                    "receptor_complex",
                    "magnitude_rank",
                    "specificity_rank",
                ],
            )
            fig.update_coloraxes(colorbar_title_text="Interaction<br>strength")

            fig.update_traces(marker=dict(sizemode="area"))
            fig.update_layout(
                template="plotly_white",
                title="Cell-Cell Communication Bubble Plot",
                xaxis_title="Target Cell Type",
                yaxis_title="Interaction (Ligand → Receptor)",
            )
        else:
            fig = px.scatter(
                display_results,
                x="magnitude_rank",
                y="interaction",
                size="bubble_size",
                color="bubble_color",
                labels={"bubble_color": "Interaction strength"},
                hover_data=[
                    "ligand_complex",
                    "receptor_complex",
                    "magnitude_rank",
                    "specificity_rank",
                ],
            )
            fig.update_coloraxes(colorbar_title_text="Interaction<br>strength")

            fig.update_xaxes(autorange="reversed")

            fig.update_traces(marker=dict(sizemode="area"))
            fig.update_layout(
                template="plotly_white",
                title="Cell-Cell Communication Bubble Plot",
                xaxis_title="Magnitude Rank (lower is stronger)",
                yaxis_title="Interaction (Ligand → Receptor)",
            )

        return (
            "done",
            results_table,
            dcc.Graph(figure=fig),
            dcc.Graph(figure=network_fig),
            no_update,
        )

    @app.callback(
        Output("ccc-loading-signal", "children", allow_duplicate=True),
        Output("ccc-summary-container", "children", allow_duplicate=True),
        Output("ccc-bubble-plot-container", "children", allow_duplicate=True),
        Output("ccc-network-plot-container", "children", allow_duplicate=True),
        Output("ccc-reference-results-container", "children", allow_duplicate=True),
        Output("ccc-show-network-labels", "value"),
        Input("ccc-reset-filters-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_liana_filters(n_clicks):
        """ "
        Reset the LIANA results filters and restore the original unfiltered results.

        Parameters
        ----------
        n_clicks : int
            Number of times the "Reset Filters" button has been clicked.

        Returns
        -------

        """
        if not n_clicks or n_clicks == 0:
            return no_update, no_update, no_update, no_update, no_update, no_update

        history = get_state_store()
        active_context = history.get("CCC_active_context", "uploaded")

        if active_context == "reference":
            ref_payload = history.get("CCC_reference_results", {})
            reference_results = ref_payload.get("results")
            file_type = ref_payload.get("file_type")
            selected_filename = ref_payload.get("selected_filename", "reference file")

            if reference_results is None or file_type is None:
                return (
                    "done",
                    no_update,
                    no_update,
                    no_update,
                    dbc.Alert(
                        "No reference CCC results available to reset.",
                        color="warning",
                    ),
                    [True],
                )

            reference_output = _render_ccc_reference_outputs(
                reference_results,
                file_type,
                selected_filename,
            )

            return (
                "done",
                no_update,
                no_update,
                no_update,
                reference_output,
                [True],
            )

        liana_results = history.get("CCC_results", {}).get("results")

        if liana_results is None:
            logger.warning("No LIANA results found in state store for filtering.")
            return (
                "done",
                html.P(
                    "No cell-cell communication results available to filter.",
                    style={"color": "#6c757d", "fontSize": "14px", "marginTop": "10px"},
                ),
                "",
                "",
                no_update,
                [True],
            )

        display_results = _display_rows(liana_results)

        results_table = liana_output_table(display_results)
        results_table = make_sortable_table(results_table, "ccc-results-table")
        summary = summarise_celltype_interactions(liana_results)

        network_fig = create_network_plot(summary, show_labels=True)

        fig = px.scatter(
            summary,
            x="target",
            y="source",
            size="bubble_size",
            color="bubble_color",
            labels={"bubble_color": "Interaction strength"},
            hover_data=[
                "interaction_count",
                "mean_magnitude",
                "mean_specificity",
            ],
        )
        fig.update_coloraxes(colorbar_title_text="Interaction<br>strength")

        fig.update_yaxes(categoryorder="category ascending")
        fig.update_layout(
            template="plotly_white",
            title="Cell-Cell Communication Bubble Plot",
        )

        return (
            "done",
            results_table,
            dcc.Graph(figure=fig),
            dcc.Graph(figure=network_fig),
            no_update,
            [True],
        )

    @app.callback(
        Output("ccc-export-loading-signal", "children"),
        Output("ccc-network-html-download", "data"),
        Input("ccc-export-network-html-button", "n_clicks"),
        State("ccc-source-filter-dropdown", "value"),
        State("ccc-target-filter-dropdown", "value"),
        State("ccc-show-network-labels", "value"),
        prevent_initial_call=True,
    )
    def export_network_fullscreen_html(
        n_clicks, source_filter, target_filter, show_labels_value
    ):
        """Export a standalone fullscreen HTML for the current CCC network view."""
        if not n_clicks or n_clicks == 0:
            return no_update, no_update

        history = get_state_store()
        liana_results = history.get("CCC_results", {}).get("results")
        if liana_results is None:
            logger.warning("No LIANA results found while exporting CCC network HTML.")
            return "No CCC results available. Run analysis first.", no_update

        filtered_results = filter_liana_results(
            liana_results,
            source=source_filter
            if source_filter and source_filter != "_all_"
            else None,
            target=target_filter
            if target_filter and target_filter != "_all_"
            else None,
        )

        show_labels = True in (show_labels_value or [])
        network_fig = create_network_plot(
            summarise_celltype_interactions(filtered_results),
            show_labels=show_labels,
        )

        network_fig.update_layout(template="plotly_white", autosize=True, margin=dict(l=20, r=20, t=60, b=20))

        figure_html = pio.to_html(
            network_fig,
            include_plotlyjs="cdn",
            full_html=False,
            config={"responsive": True, "displaylogo": False},
        )

        full_html = f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>CCC Network Fullscreen</title>
  <style>
    html, body {{
      width: 100%;
      height: 100%;
      margin: 0;
      background: #ffffff;
      overflow: hidden;
    }}
    .plotly-graph-div {{
      width: 100vw !important;
      height: 100vh !important;
    }}
  </style>
</head>
<body>
  {figure_html}
</body>
</html>
"""

        return "Export ready. Your download should begin automatically.", {
            "content": full_html,
            "filename": "ccc_network_fullscreen.html",
            "type": "text/html",
        }
