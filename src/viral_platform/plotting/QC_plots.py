import plotly.express as px
from dash import dcc, html
from viral_platform.state.dataset_store import sync_state_with_dataset

def create_qc_plots(adata):
    if adata is None:
        return "No dataset available for QC plotting."
    
    required_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    missing_cols = [col for col in required_cols if col not in adata.obs]
    if missing_cols:
        return f"Missing QC metrics after calculation: {', '.join(missing_cols)}"

    adata.obs["nCount_RNA"] = adata.obs["total_counts"]
    adata.obs["nFeature_RNA"] = adata.obs["n_genes_by_counts"]
    adata.obs["percent.mt"] = adata.obs["pct_counts_mt"]
    sync_state_with_dataset(adata)

    ncount_fig = px.violin(
        adata.obs,
        y="nCount_RNA",
        box=True,
        points=False,
        title="nCount_RNA",
    )
    ncount_fig.update_layout(template="plotly_white", showlegend=False, xaxis_title="", yaxis_title="Counts")

    nfeature_fig = px.violin(
        adata.obs,
        y="nFeature_RNA",
        box=True,
        points=False,
        title="nFeature_RNA",
    )
    nfeature_fig.update_layout(template="plotly_white", showlegend=False, xaxis_title="", yaxis_title="Features")

    percent_mt_fig = px.violin(
        adata.obs,
        y="percent.mt",
        box=True,
        points=False,
        title="percent.mt",
    )
    percent_mt_fig.update_layout(template="plotly_white", showlegend=False, xaxis_title="", yaxis_title="Percent")

    return html.Div([
        dcc.Graph(figure=ncount_fig, id="ncount-violin"),
        dcc.RangeSlider(
            id="min-counts-slider",
            min=adata.obs["nCount_RNA"].min(),
            max=adata.obs["nCount_RNA"].max(),
            step=100,
            allowCross=False,
            value=[adata.obs["nCount_RNA"].min(), adata.obs["nCount_RNA"].max()],
        ),
        dcc.Graph(figure=nfeature_fig, id="nfeature-violin"),
        dcc.RangeSlider(
            id="min-features-slider",
            min=adata.obs["nFeature_RNA"].min(),
            max=adata.obs["nFeature_RNA"].max(),
            step=10,
            allowCross=False,
            value=[adata.obs["nFeature_RNA"].min(), adata.obs["nFeature_RNA"].max()],
        ),
        dcc.Graph(figure=percent_mt_fig, id="percent-mt-violin"),
        dcc.Slider(
            id="max-percent-mt-slider",
            min=adata.obs["percent.mt"].min(),
            max=adata.obs["percent.mt"].max(),
            value=adata.obs["percent.mt"].max(),
        ),
        html.Button("Apply QC Filters", id="apply-qc-filters-button", n_clicks=0, style={"marginTop": "20px"}),
        html.Hr(),
        html.H4("Doublet Detection (Scrublet)"),
        dcc.Checklist(
            id="doublet-filter-checklist",
            options=[{"label": " Remove predicted doublets", "value": "filter"}],
            value=[],
        ),
        html.Button("Run Scrublet", id="run-scrublet-button", n_clicks=0, style={"marginTop": "10px"}),
        dcc.Loading(html.Div(id="scrublet-result")),
    ],
    id="qc-plot-container")

