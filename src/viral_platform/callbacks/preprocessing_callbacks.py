from dash import Input, Output, State, no_update
from viral_platform.state.dataset_store import cache_results, get_working_dataset
from viral_platform.analysis.preprocessing import preprocess_data, run_clustering
from viral_platform.plotting.elbow_plot import create_elbow_plot
from viral_platform.plotting.clustering import create_umap_plot
from viral_platform.utils.runtime_tracker import get_tracker
import plotly.express as px
import plotly.graph_objects as go
import logging


logger = logging.getLogger(__name__)

_preprocess_tracker = get_tracker("Preprocessing: Run")
_pca_tracker = get_tracker("Preprocessing: PCA Selection")


def register_preprocessing_callbacks(app):
    @app.callback(
        Output("preprocess-advanced-options-collapse", "is_open"),
        Input("advanced-options-button", "n_clicks"),
        State("preprocess-advanced-options-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_preprocess_advanced_options(n_clicks, is_open):
        if not n_clicks:
            return is_open
        return not is_open

    @app.callback(
        Output("preprocess-hvg-slider", "value"),
        Output("preprocess-scale-max-slider", "value"),
        Output("preprocess-neighbors-slider", "value"),
        Input("preprocess-advanced-reset-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_preprocess_advanced_defaults(n_clicks):
        if not n_clicks:
            return no_update, no_update, no_update
        return 2000, 10, 10

    @app.callback(
        Output("preprocess-loading-signal", "children"),
        Output("preprocess-temp-container", "children"),
        Input("run-preprocess-button", "n_clicks"),
        State("preprocess-hvg-slider", "value"),
        State("preprocess-scale-max-slider", "value"),
        prevent_initial_call=True,
    )
    def run_preprocessing(n_clicks, n_top_genes, scale_max_value):
        if not n_clicks:
            return no_update, no_update
        _preprocess_tracker.begin()
        try:
            with _preprocess_tracker.phase("analysis"):
                preprocess_data(
                    n_top_genes=n_top_genes if n_top_genes is not None else 2000,
                    scale_max_value=scale_max_value if scale_max_value is not None else 10,
                )
                adata = get_working_dataset()
            if adata is None:
                logger.warning("Preprocessing completed but no dataset found in state store.")
                _preprocess_tracker.finish()
                return no_update, "Preprocessing completed, but no dataset found."
            logger.info("Preprocessing completed successfully.")
            with _preprocess_tracker.phase("processing"):
                result = create_elbow_plot(adata)
            with _preprocess_tracker.phase("visualization"):
                cache_results(**{"preprocess-temp-container": result})
            _preprocess_tracker.finish()
            return no_update, result
        except Exception as exc:
            logger.exception("Preprocessing failed: %s", str(exc))
            _preprocess_tracker.finish()
            return no_update, f"Preprocessing failed: {str(exc)}"
        
    @app.callback(
        Output("elbow-plot", "figure"),
        Input("pc-slider", "value")
    )
    def update_elbow_plot(n_pcs):
        adata = get_working_dataset()
        if adata is None:
            return px.bar(title="Elbow Plot")
        try:
            elbow_layout = create_elbow_plot(adata)
            elbow_fig = go.Figure(elbow_layout.children[0].figure)
            if n_pcs is not None:
                elbow_fig.add_vline(
                    x=n_pcs-1,
                    line_dash="dash",
                    line_color="red",
                    annotation_text=f"Selected PCs: {n_pcs}",
                    annotation_position="top right",
                )
            return elbow_fig
        except Exception as exc:
            logger.exception("Failed to update elbow plot: %s", str(exc))
            return px.bar(title="Elbow Plot")
    
    @app.callback(
        Output("clustering-loading", "children"),
        Output("selected-pcs-output", "children"),
        Input("select-pcs-button", "n_clicks"),
        State("pc-slider", "value"),
        State("preprocess-neighbors-slider", "value"),
        prevent_initial_call=True,
    )
    def apply_pca_selection(n_clicks, n_pcs, n_neighbors):
        if not n_clicks or n_pcs is None:
            return no_update, no_update
        _pca_tracker.begin()
        try:
            logger.info("Apply PCA selection clicked. n_clicks=%s, n_pcs=%s", n_clicks, n_pcs)
            with _pca_tracker.phase("analysis"):
                run_clustering(
                    n_dims=n_pcs,
                    n_neighbors=n_neighbors if n_neighbors is not None else 10,
                )
            logger.info("PCA selection applied successfully with %d PCs.", n_pcs)
            with _pca_tracker.phase("processing"):
                graph = create_umap_plot()
            with _pca_tracker.phase("visualization"):
                cache_results(**{"selected-pcs-output": graph})
            _pca_tracker.finish()
            return f"PCA selection applied with {n_pcs} PCs.", graph
        except Exception as exc:
            logger.exception("Failed to apply PCA selection: %s", str(exc))
            return no_update, f"Failed to apply PCA selection: {str(exc)}"
