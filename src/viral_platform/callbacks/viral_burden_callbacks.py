from dash import Input, Output, State, dash_table, html, no_update, dcc
import dash_bootstrap_components as dbc
from viral_platform.analysis.viral_burden_associations import calculate_viral_burden_associations, identify_significant_associations
from viral_platform.state.dataset_store import cache_results, get_working_dataset, get_state_store, set_working_dataset, sync_state_with_dataset
import numpy as np
import pandas as pd
import plotly.express as px
import logging
from viral_platform.utils.sample_column_utils import normalize_column_name, resolve_sample_column
from viral_platform.utils.runtime_tracker import get_tracker

logger = logging.getLogger(__name__)

_viral_burden_tracker = get_tracker("Viral Burden: Analysis")
_viral_burden_assoc_tracker = get_tracker("Viral Burden: Associations")


def _build_options(values):
    """Build unique dropdown options from an iterable of values."""
    seen = set()
    options = []
    for value in values:
        normalized = str(value)
        if normalized in seen:
            continue
        seen.add(normalized)
        options.append({"label": normalized, "value": normalized})
    return options


def _resolve_column(columns, priority_names, contains_names=None):
    """Resolve a metadata column using exact normalized names then contains checks."""
    contains_names = contains_names or []
    normalized_to_original = {normalize_column_name(col): col for col in list(columns)}

    for name in priority_names:
        normalized = normalize_column_name(name)
        if normalized in normalized_to_original:
            return normalized_to_original[normalized]

    for normalized_col, original_col in normalized_to_original.items():
        for token in contains_names:
            normalized_token = normalize_column_name(token)
            if normalized_token and normalized_token in normalized_col:
                return original_col

    return None


def _resolve_celltype_column(columns):
    """Resolve the most likely cell-type annotation column."""
    return _resolve_column(
        columns,
        priority_names=["cell_type", "celltype", "CellType"],
        contains_names=["celltype", "cell_type", "cell type", "annotation", "cluster"],
    )


def _resolve_condition_column(columns):
    """Resolve the most likely condition/group column."""
    return _resolve_column(
        columns,
        priority_names=["condition", "group", "status", "disease", "treatment"],
        contains_names=["condition", "group", "status", "disease", "treat"],
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

def viral_burden_results(adata):
    """Render summary metrics for viral burden analysis as a card/table."""
    infected_cells = (adata.obs['viral_counts'] > 0).sum()
    return dbc.Card(
        dbc.CardBody(
            [
                html.H5("Viral Burden Analysis Results"),
                dbc.Table(
                    [
                        html.Tbody([
                            html.Tr([html.Td("Cells with Viral Reads"), html.Td(f"{infected_cells}/{adata.n_obs}")]),
                            html.Tr([html.Td("Maximum Viral Counts"), html.Td(f"{adata.obs['viral_counts'].max()}")]),
                            html.Tr([html.Td("Average Viral Counts"), html.Td(f"{adata.obs['viral_counts'].mean():.2f}")]),
                            html.Tr([html.Td("Maximum Viral Burden"), html.Td(f"{adata.obs['viral_burden'].max():.4f}")]),
                            html.Tr([html.Td("Average Viral Burden"), html.Td(f"{adata.obs['viral_burden'].mean():.4f}")]),
                            html.Tr([html.Td("Maximum Viral Burden (%)"), html.Td(f"{adata.obs['viral_burden_percent'].max():.2f}%")]),
                            html.Tr([html.Td("Average Viral Burden (%)"), html.Td(f"{adata.obs['viral_burden_percent'].mean():.2f}%")]),
                        ])
                    ]
                )
            ]
        )
    )


def _build_infection_umap(adata):
    """Create UMAP colored by binary infection status (infected vs non-infected)."""
    if "X_umap" not in adata.obsm:
        return html.Div("UMAP coordinates are not available. Run preprocessing/clustering first.")

    umap_df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    umap_df["infection_status"] = adata.obs["infection_status"].astype(str).values

    # Ensure all non-infected labels map to grey while infected remains green.
    umap_df["infection_status_plot"] = np.where(
        umap_df["infection_status"].str.lower().isin(["infected"]),
        "Infected",
        "Non-infected",
    )

    fig = px.scatter(
        umap_df,
        x="UMAP1",
        y="UMAP2",
        color="infection_status_plot",
        color_discrete_map={"Infected": "green", "Non-infected": "grey"},
        title="Infection UMAP",
        opacity=0.85,
    )
    fig.update_traces(marker={"size": 4})
    fig.update_layout(template="plotly_white", legend_title_text="")
    return dcc.Graph(figure=fig)


def _build_viral_burden_umap(adata):
    """Create UMAP colored by continuous viral burden values."""
    if "X_umap" not in adata.obsm:
        return html.Div("UMAP coordinates are not available. Run preprocessing/clustering first.")

    umap_df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    umap_df["viral_burden"] = pd.to_numeric(adata.obs["viral_burden"], errors="coerce").fillna(0.0).values

    fig = px.scatter(
        umap_df,
        x="UMAP1",
        y="UMAP2",
        color="viral_burden",
        color_continuous_scale="YlOrRd",
        title="Viral Burden UMAP",
        opacity=0.85,
    )
    fig.update_layout(template="plotly_white")
    fig.update_traces(marker={"size": 4})
    fig.update_coloraxes(colorbar_title="Viral burden")
    return dcc.Graph(figure=fig)


def _build_violin_plots(
    adata,
    celltype_col=None,
    condition_col=None,
    sample_col=None,
    metadata_sample_columns=None,
):
    """Create violin plots of viral burden percent by cell type, condition, and sample."""
    obs = adata.obs.copy()
    obs["viral_burden_percent"] = pd.to_numeric(obs["viral_burden_percent"], errors="coerce").fillna(0.0)

    if celltype_col is None:
        celltype_col = _resolve_celltype_column(obs.columns)
    if condition_col is None:
        condition_col = _resolve_condition_column(obs.columns)
    if sample_col is None:
        sample_col = resolve_sample_column(
            obs.columns,
            metadata_sample_columns=metadata_sample_columns,
            obs_df=obs,
        )

    sections = []

    if celltype_col is not None:
        fig_celltype = px.violin(
            obs,
            x=celltype_col,
            y="viral_burden_percent",
            box=True,
            points=False,
            title=f"Viral burden (%) by cell type ({celltype_col})",
        )
        fig_celltype.update_layout(template="plotly_white", xaxis_tickangle=45)
        sections.append(dcc.Graph(figure=fig_celltype))
    else:
        sections.append(html.Div("No cell type column found for cell type violin plot."))

    if condition_col is not None:
        fig_condition = px.violin(
            obs,
            x=condition_col,
            y="viral_burden_percent",
            box=True,
            points=False,
            title=f"Viral burden (%) by condition ({condition_col})",
        )
        fig_condition.update_layout(template="plotly_white", xaxis_tickangle=45)
        sections.append(dcc.Graph(figure=fig_condition))
    else:
        sections.append(html.Div("No condition/group column found for condition violin plot."))

    if sample_col is not None:
        fig_sample = px.violin(
            obs,
            x=sample_col,
            y="viral_burden_percent",
            box=True,
            points=False,
            title=f"Viral burden (%) by sample ({sample_col})",
        )
        fig_sample.update_layout(template="plotly_white", xaxis_tickangle=45)
        sections.append(dcc.Graph(figure=fig_sample))
    else:
        sections.append(html.Div("No sample column found for sample violin plot."))

    return html.Div(sections)


def _build_celltype_infection_fraction_plot(adata, celltype_col=None):
    """Create per-cell-type infection fraction bar chart."""
    obs = adata.obs.copy()
    if celltype_col is None:
        celltype_col = _resolve_celltype_column(obs.columns)
    if celltype_col is None:
        return html.Div("No cell type column found for infection fraction plot.")

    summary = (
        obs.assign(
            infected=obs["infection_status"].astype(str).str.lower().eq("infected").astype(int)
        )
        .groupby(celltype_col, observed=False)
        .agg(total_cells=("infected", "size"), infected_cells=("infected", "sum"))
        .reset_index()
    )

    if summary.empty:
        return html.Div("No cells available for infection fraction plot.")

    summary["infection_fraction"] = summary["infected_cells"] / summary["total_cells"]
    summary = summary.sort_values("infection_fraction", ascending=False)

    fig = px.bar(
        summary,
        x=celltype_col,
        y="infection_fraction",
        color="infection_fraction",
        color_continuous_scale="Greens",
        title="Cell type infection fraction",
        hover_data={
            "infected_cells": True,
            "total_cells": True,
            "infection_fraction": ":.2%",
        },
    )
    fig.update_layout(template="plotly_white",xaxis_tickangle=45, yaxis_tickformat=".0%")
    return dcc.Graph(figure=fig)

def register_viral_burden_callbacks(app):
    """Register callbacks for viral burden analysis, plotting, and associations."""
    @app.callback(
        Output("viral-burden-advanced-options-collapse", "is_open"),
        Input("viral-burden-advanced-options-button", "n_clicks"),
        State("viral-burden-advanced-options-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_viral_burden_advanced_options(n_clicks, is_open):
        """Toggle visibility of advanced viral burden column selection controls."""
        if not n_clicks:
            return is_open
        return not is_open

    @app.callback(
        Output("viral-burden-infected-threshold-input", "value"),
        Output("viral-burden-association-min-cells-input", "value"),
        Output("viral-burden-association-corr-threshold-input", "value"),
        Output("viral-burden-association-fdr-threshold-input", "value"),
        Input("viral-burden-advanced-reset-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_viral_burden_advanced_defaults(n_clicks):
        if not n_clicks:
            return no_update, no_update, no_update, no_update
        return 0, 10, 0.3, 0.05

    @app.callback(
        Output("viral-burden-celltype-column-dropdown", "options"),
        Output("viral-burden-celltype-column-dropdown", "value"),
        Output("viral-burden-condition-column-dropdown", "options"),
        Output("viral-burden-condition-column-dropdown", "value"),
        Output("viral-burden-sample-column-dropdown", "options"),
        Output("viral-burden-sample-column-dropdown", "value"),
        Input("active-dataset-version", "data"),
    )
    def populate_viral_burden_column_dropdowns(_dataset_version):
        """Populate advanced column dropdowns from the current AnnData metadata."""
        # Always use the live working dataset.  ``metadata_info`` is refreshed
        # whenever a module updates it, but deriving from ``obs`` here prevents
        # a stale state snapshot from hiding newly-created annotation columns.
        adata = get_working_dataset()
        state = get_state_store()
        metadata_info = state.get("metadata_info", {})
        groupable_columns = list(metadata_info.get("groupable_columns", []))
        if adata is not None:
            groupable_columns = [
                column for column in groupable_columns if column in adata.obs.columns
            ]
        metadata_sample_columns = metadata_info.get("sample_columns", [])

        celltype_options = _build_options(groupable_columns)
        condition_options = _build_options(groupable_columns)

        sample_candidates = list(groupable_columns)
        sample_candidates.extend([c for c in metadata_sample_columns if c not in sample_candidates])
        sample_options = _build_options(sample_candidates)

        default_celltype = _resolve_celltype_column(groupable_columns)
        default_condition = _resolve_condition_column(groupable_columns)
        obs_df = adata.obs if adata is not None else None
        default_sample = resolve_sample_column(
            sample_candidates,
            metadata_sample_columns=metadata_sample_columns,
            obs_df=obs_df,
        )

        if not celltype_options:
            celltype_options = [{"label": "Upload a dataset to select cell type column", "value": ""}]
            default_celltype = ""
        if not condition_options:
            condition_options = [{"label": "Upload a dataset to select condition column", "value": ""}]
            default_condition = ""
        if not sample_options:
            sample_options = [{"label": "Upload a dataset to select sample column", "value": ""}]
            default_sample = ""

        if default_celltype not in {opt["value"] for opt in celltype_options}:
            default_celltype = celltype_options[0]["value"] if celltype_options else ""
        if default_condition not in {opt["value"] for opt in condition_options}:
            default_condition = condition_options[0]["value"] if condition_options else ""
        if default_sample not in {opt["value"] for opt in sample_options}:
            default_sample = ""

        if sample_options and default_sample == "":
            sample_options = [{
                "label": "No confident sample ID detected automatically. Please choose a sample column.",
                "value": "",
            }] + [opt for opt in sample_options if opt.get("value") != ""]

        return (
            celltype_options,
            default_celltype,
            condition_options,
            default_condition,
            sample_options,
            default_sample,
        )

    @app.callback(
        Output("viral-burden-loading-signal", "children"),
        Output("viral-burden-results-container", "children"),
        Output("viral-burden-infection-umap-container", "children"),
        Output("viral-burden-umap-container", "children"),
        Output("viral-burden-violin-container", "children"),
        Output("viral-burden-celltype-fraction-container", "children"),
        Output("viral-burden-associations-container", "hidden"),
        Input("run-viral-burden-analysis-button", "n_clicks"),
        State("viral-burden-celltype-column-dropdown", "value"),
        State("viral-burden-condition-column-dropdown", "value"),
        State("viral-burden-sample-column-dropdown", "value"),
        State("viral-burden-infected-threshold-input", "value"),
        prevent_initial_call=True,
    )
    @_viral_burden_tracker.track
    def run_viral_burden_analysis(
        n_clicks,
        selected_celltype_column,
        selected_condition_column,
        selected_sample_column,
        infected_threshold,
    ):
        """
        Calculate viral burden for each cell and render tabbed outputs.

        Parameters
        ----------
        adata : AnnData
            Dataset containing expression matrix.

        detected_features : list
            List of viral feature names returned by the viral
            detection module.

        Returns
        -------
        tuple
            Loading signal, rendered results/plots, and association-section visibility.
        """
        if n_clicks is None or n_clicks == 0:
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
            )
        
        history = get_state_store()
        adata = get_working_dataset()
        if adata is None:
            message = "No dataset available for viral burden analysis."
            return "done", message, message, message, message, message, True
        
        features = history.get("viral_detection", {}).get("viral_features", "")
        metadata_info = history.get("metadata_info", {})
        metadata_sample_columns = metadata_info.get("sample_columns", [])
        history = None # remove local history object to free memory
        if not features:
            message = "No viral features detected. Please run viral gene detection first."
            return "done", message, message, message, message, message, True

        if isinstance(features, str):
            features = [f.strip() for f in features.split(",") if f.strip()]
        else:
            features = [f for f in features if f]

        if not features:
            message = "No valid viral features detected. Please run viral gene detection first."
            return "done", message, message, message, message, message, True
        
        # grabbing raw counts matrix from adata.layers["counts"]
        with _viral_burden_tracker.phase("analysis"):
            matrix = (adata.layers["counts"])

            # extracting only viral features
            viral_matrix = matrix[:, adata.var_names.isin(features)].copy()

            # Sum viral counts per cell
            viral_counts = viral_matrix.sum(axis=1)
            if hasattr(viral_counts, "A1"):  # Check if it's a sparse matrix
                viral_counts = viral_counts.A1  # Convert to 1D array
            else:
                viral_counts = np.array(viral_counts).flatten()  # Ensure it's a 1D array

            
            # store raw viral counts in adata.obs
            adata.obs["viral_counts"] = viral_counts


            # calculate viral burden
            adata.obs["viral_burden"] = (
                adata.obs["viral_counts"] / adata.obs["total_counts"]
            )

            # percentage burden
            adata.obs["viral_burden_percent"] = adata.obs["viral_burden"] * 100

            #infection status
            resolved_infected_threshold = (
                float(infected_threshold) if infected_threshold is not None else 0.0
            )
            adata.obs["infection_status"] = np.where(
                adata.obs["viral_counts"] > resolved_infected_threshold,
                "Infected",
                "Bystander",
            )

            # log1p transformation of viral counts
            adata.obs["log1p_viral_counts"] = np.log1p(adata.obs["viral_counts"])

            # update the working dataset in the state store
            set_working_dataset(adata)
            sync_state_with_dataset(adata)


        logger.info("Viral burden analysis completed successfully.")

        celltype_col = (
            selected_celltype_column
            if selected_celltype_column and selected_celltype_column in adata.obs.columns
            else _resolve_celltype_column(adata.obs.columns)
        )
        condition_col = (
            selected_condition_column
            if selected_condition_column and selected_condition_column in adata.obs.columns
            else _resolve_condition_column(adata.obs.columns)
        )
        sample_col = (
            selected_sample_column
            if selected_sample_column and selected_sample_column in adata.obs.columns
            else resolve_sample_column(
                adata.obs.columns,
                metadata_sample_columns=metadata_sample_columns,
                obs_df=adata.obs,
            )
        )

        with _viral_burden_tracker.phase("processing"):
            infection_umap = _build_infection_umap(adata)
            viral_burden_umap = _build_viral_burden_umap(adata)
            violin_plots = _build_violin_plots(
                adata,
                celltype_col=celltype_col,
                condition_col=condition_col,
                sample_col=sample_col,
                metadata_sample_columns=metadata_sample_columns,
            )
            celltype_fraction_plot = _build_celltype_infection_fraction_plot(
                adata,
                celltype_col=celltype_col,
            )

        with _viral_burden_tracker.phase("visualization"):
            result_components = {
                "viral-burden-results-container": viral_burden_results(adata),
                "viral-burden-infection-umap-container": infection_umap,
                "viral-burden-umap-container": viral_burden_umap,
                "viral-burden-violin-container": violin_plots,
                "viral-burden-celltype-fraction-container": celltype_fraction_plot,
            }
            cache_results(**result_components)
        return (
            "done",
            result_components["viral-burden-results-container"],
            result_components["viral-burden-infection-umap-container"],
            result_components["viral-burden-umap-container"],
            result_components["viral-burden-violin-container"],
            result_components["viral-burden-celltype-fraction-container"],
            False,
        )
    
    @app.callback(
        Output("viral-burden-associations-loading-signal", "children"),
        Output("viral-burden-associations-results-container", "children"),
        Output("viral-burden-associations-significant-results-container", "children"),
        Input("run-viral-burden-association-button", "n_clicks"),
        State("viral-burden-association-min-cells-input", "value"),
        State("viral-burden-association-corr-threshold-input", "value"),
        State("viral-burden-association-fdr-threshold-input", "value"),
        prevent_initial_call=True,
    )
    @_viral_burden_assoc_tracker.track
    def run_viral_burden_associations(n_clicks, min_cells, corr_threshold, fdr_threshold):
        """Run viral burden association analysis and render full/significant result tables."""
        if n_clicks is None or n_clicks == 0:
            return no_update, no_update, no_update
        
        adata = get_working_dataset()
        if adata is None:
            return no_update, "No dataset available for viral burden association analysis.", ""
        
        features = get_state_store().get("viral_detection", {}).get("viral_features", "")
        if not features:
            return no_update, "No viral features detected. Please run viral gene detection first.", ""

        if isinstance(features, str):
            features = [f.strip() for f in features.split(",") if f.strip()]
        else:
            features = [f for f in features if f]

        if not features:
            return no_update, "No valid viral features detected. Please run viral gene detection first.", ""
        
        try:
            resolved_min_cells = int(min_cells) if min_cells is not None else 10
            resolved_corr_threshold = float(corr_threshold) if corr_threshold is not None else 0.3
            resolved_fdr_threshold = float(fdr_threshold) if fdr_threshold is not None else 0.05

            with _viral_burden_assoc_tracker.phase("analysis"):
                associations_df = calculate_viral_burden_associations(
                    features,
                    min_cells=resolved_min_cells,
                )
                logger.info("Viral burden association analysis completed successfully.")
                significant_associations_df = identify_significant_associations(
                    associations_df,
                    corr_threshold=resolved_corr_threshold,
                    fdr_threshold=resolved_fdr_threshold,
                )
            logger.info("Significant viral burden associations identified successfully.")
            with _viral_burden_assoc_tracker.phase("processing"):
                full_table = make_sortable_table(associations_df, "viral-burden-associations-table")
                significant_table = make_sortable_table(significant_associations_df, "viral-burden-significant-associations-table")

                full_table_section = html.Div(
                    [
                        html.H5("Full Associations Table", style={"marginBottom": "8px"}),
                        full_table,
                    ],
                    style={"marginBottom": "20px"},
                )
                significant_table_section = html.Div(
                    [
                        html.H5("Significant Associations Table", style={"marginBottom": "8px"}),
                        significant_table,
                    ]
                )

            with _viral_burden_assoc_tracker.phase("visualization"):
                cache_results(**{
                    "viral-burden-associations-results-container": full_table_section,
                    "viral-burden-associations-significant-results-container": significant_table_section,
                })
            return (
                "done",
                full_table_section,
                significant_table_section,
            )
        except Exception as e:
            logger.exception("Failed to calculate viral burden associations: %s", str(e))
            return no_update, f"Failed to calculate viral burden associations: {str(e)}", ""
