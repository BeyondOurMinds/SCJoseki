import logging

from dash import Input, Output, State, html, no_update
import plotly.express as px

from viral_platform.plotting.QC_plots import create_qc_plots
from viral_platform.state.dataset_store import cache_results, get_state_store, get_working_dataset, set_working_dataset, update_state_store
from viral_platform.utils.runtime_tracker import get_tracker

logger = logging.getLogger(__name__)

_qc_tracker = get_tracker("QC Plots")


def register_qc_callbacks(app):
    @app.callback(
        Output("qc-generate-loading-signal", "children"),
        Output("qc-temp-container", "children"),
        Input("generate-qc-plots-button", "n_clicks")
    )
    def render_qc_plots(n_clicks):
        if n_clicks == 0:
            # Page remounts invoke this callback with the initial button value.
            # Preserve the component restored from the shared results cache.
            return no_update, no_update
        _qc_tracker.begin()
        adata = get_working_dataset()
        if adata is None:
            logger.warning("QC plot generation requested without an active dataset.")
            _qc_tracker.finish()
            return "done", "Upload a dataset to view QC plots."
        try:
            with _qc_tracker.phase("processing"):
                result = create_qc_plots(adata)
            with _qc_tracker.phase("visualization"):
                cache_results(**{"qc-temp-container": result})
            _qc_tracker.finish()
            return "done", result
        except Exception:
            logger.exception("Failed to render QC plots.")
            _qc_tracker.finish()
            return "done", "An error occurred while generating QC plots."
    
    @app.callback(
        Output("ncount-violin", "figure"),
        Input("min-counts-slider", "value")
    )
    def update_ncount_violin(min_counts):
        adata = get_working_dataset()
        state = get_state_store()
        if adata is None:
            logger.warning("nCount violin update requested without dataset.")
            return px.violin(title="nCount_RNA")

        if not min_counts or len(min_counts) != 2:
            stored = state.get("nCount_RNA", {})
            stored_min = stored.get("min")
            stored_max = stored.get("max")
            if stored_min is not None and stored_max is not None:
                min_counts = [stored_min, stored_max]
            else:
                min_counts = [adata.obs["nCount_RNA"].min(), adata.obs["nCount_RNA"].max()]
        
        try:
            ncount_fig = px.violin(
                adata.obs,
                y="nCount_RNA",
                box=True,
                points=False,
                title="nCount_RNA",
            )
            ncount_fig.update_layout(template="plotly_white", showlegend=False, xaxis_title="", yaxis_title="Counts")

            ncount_fig.add_hline(y=min_counts[0], line_dash="dash", line_color="red", annotation_text="Min Counts", annotation_position="top left")
            ncount_fig.add_hline(y=min_counts[1], line_dash="dash", line_color="red", annotation_text="Max Counts", annotation_position="top right")
            return ncount_fig
        except Exception:
            logger.exception("Failed to update nCount violin figure.")
            return px.violin(title="nCount_RNA")
    
    @app.callback(
        Output("nfeature-violin", "figure"),
        Input("min-features-slider", "value")
    )
    def update_nfeature_violin(min_features):
        adata = get_working_dataset()
        state = get_state_store()
        if adata is None:
            logger.warning("nFeature violin update requested without dataset.")
            return px.violin(title="nFeature_RNA")

        if not min_features or len(min_features) != 2:
            stored = state.get("nFeature_RNA", {})
            stored_min = stored.get("min")
            stored_max = stored.get("max")
            if stored_min is not None and stored_max is not None:
                min_features = [stored_min, stored_max]
            else:
                min_features = [adata.obs["nFeature_RNA"].min(), adata.obs["nFeature_RNA"].max()]
        
        try:
            nfeature_fig = px.violin(
                adata.obs,
                y="nFeature_RNA",
                box=True,
                points=False,
                title="nFeature_RNA",
            )
            nfeature_fig.update_layout(template="plotly_white", showlegend=False, xaxis_title="", yaxis_title="Features")

            nfeature_fig.add_hline(y=min_features[0], line_dash="dash", line_color="red", annotation_text="Min Features", annotation_position="top left")
            nfeature_fig.add_hline(y=min_features[1], line_dash="dash", line_color="red", annotation_text="Max Features", annotation_position="top right")
            return nfeature_fig
        except Exception:
            logger.exception("Failed to update nFeature violin figure.")
            return px.violin(title="nFeature_RNA")
    
    @app.callback(
        Output("percent-mt-violin", "figure"),
        Input("max-percent-mt-slider", "value")
    )
    def update_percent_mt_violin(max_percent_mt):
        adata = get_working_dataset()
        state = get_state_store()

        if adata is None:
            logger.warning("percent.mt violin update requested without dataset.")
            return px.violin(title="Percent Mitochondrial Genes")

        if max_percent_mt is None:
            max_percent_mt = state.get("percent_mt")
            if max_percent_mt is None:
                max_percent_mt = adata.obs["percent.mt"].max()
        
        try:
            percent_mt_fig = px.violin(
                adata.obs,
                y="percent.mt",
                box=True,
                points=False,
                title="Percent Mitochondrial Genes",
            )
            percent_mt_fig.update_layout(template="plotly_white", showlegend=False, xaxis_title="", yaxis_title="Percent")

            percent_mt_fig.add_hline(y=max_percent_mt, line_dash="dash", line_color="red", annotation_text="Max Percent MT", annotation_position="top right")
            return percent_mt_fig
        except Exception:
            logger.exception("Failed to update percent.mt violin figure.")
            return px.violin(title="Percent Mitochondrial Genes")
    
    @app.callback(
        Output("qc-plot-container", "children"),
        Input("apply-qc-filters-button", "n_clicks"),
        State("min-counts-slider", "value"),
        State("min-features-slider", "value"),
        State("max-percent-mt-slider", "value")
    )
    def apply_qc_filters(n_clicks, min_counts, min_features, max_percent_mt):
        if n_clicks == 0:
            return no_update
        
        adata = get_working_dataset()
        if adata is None:
            logger.warning("Apply QC filters clicked without dataset.")
            return "No dataset available to apply QC filters."

        try:
            if not min_counts or len(min_counts) != 2:
                raise ValueError("Invalid nCount range provided.")
            if not min_features or len(min_features) != 2:
                raise ValueError("Invalid nFeature range provided.")
            if max_percent_mt is None:
                raise ValueError("Invalid percent.mt threshold provided.")

            adata = adata[adata.obs["nCount_RNA"] >= min_counts[0]]
            adata = adata[adata.obs["nCount_RNA"] <= min_counts[1]]
            adata = adata[adata.obs["nFeature_RNA"] >= min_features[0]]
            adata = adata[adata.obs["nFeature_RNA"] <= min_features[1]]
            adata = adata[adata.obs["percent.mt"] <= max_percent_mt]

            set_working_dataset(adata)
            update_state_store(
                nCount_RNA={"min": float(min_counts[0]), "max": float(min_counts[1])},
                nFeature_RNA={"min": float(min_features[0]), "max": float(min_features[1])},
                percent_mt=float(max_percent_mt),
            )
            logger.info(
                "Applied QC filters: nCount=%s, nFeature=%s, max_percent_mt=%s. Remaining cells: %s",
                min_counts,
                min_features,
                max_percent_mt,
                adata.n_obs,
            )
            result = create_qc_plots(adata)
            cache_results(**{"qc-temp-container": result, "qc-plot-container": result})
            return result
        except Exception:
            logger.exception("Failed to apply QC filters.")
            return html.Div("An error occurred while applying QC filters.")
