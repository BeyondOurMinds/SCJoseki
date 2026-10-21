from dash import html, dcc
from viral_platform.state.dataset_store import get_cached_result
import dash_bootstrap_components as dbc

def create_pathway_enrichment_panel():
    return html.Div(
        style={
            "backgroundColor": "#eff7ff",
            "padding": "20px",
            "borderRadius": "5px",
            "border": "1px solid #000000",
            "margin": "0 0 20px 0",
        },
        children=[
             # title and run button
            dbc.Row([
                dbc.Col(
                    html.Div(
                        [
                            html.H4(
                                ["Pathway Enrichment Analysis"],
                                style={"display": "inline-flex", "alignItems": "center", "marginBottom": "2px"},
                            ),
                            html.P(
                                "Perform pathway enrichment analysis",
                                style={"color": "#6c757d", "marginBottom": "16px", "fontSize": "14px"},
                            ),
                        ]
                    ),
                ),
                dbc.Col(
                    dbc.Button("Run Pathway Enrichment", id="run-pathway-enrichment-button", n_clicks=0, color="primary", className="mb-3"),
                    width="auto"
                ),
                dbc.Col(
                    dbc.Button("Advanced", id="pathway-enrichment-advanced-options-button", n_clicks=0, color="secondary", className="mb-3"),
                    width="auto"
                ),
            ]),

            dbc.Collapse(
                id="pathway-enrichment-advanced-options-collapse",
                is_open=False,
                children=[
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    html.P("Adjusted p-value / FDR cutoff"),
                                    dcc.Input(
                                        id="pathway-enrichment-pvalue-cutoff-input",
                                        type="number",
                                        min=0,
                                        max=1,
                                        step=any,
                                        value=0.05,
                                        style={"width": "100%"},
                                    ),
                                ],
                                md=4,
                            ),
                            dbc.Col(
                                [
                                    html.P("Absolute log2 fold-change cutoff (ORA)"),
                                    dcc.Input(
                                        id="pathway-enrichment-logfc-cutoff-input",
                                        type="number",
                                        min=0,
                                        max=10,
                                        step=0.1,
                                        value=1.0,
                                        style={"width": "100%"},
                                    ),
                                ],
                                md=4,
                            ),
                            dbc.Col(
                                dbc.Button(
                                    "Reset Advanced Defaults",
                                    id="pathway-enrichment-advanced-reset-button",
                                    n_clicks=0,
                                    color="secondary",
                                    outline=True,
                                    className="mt-4",
                                ),
                                md=4,
                            ),
                        ],
                        style={"marginBottom": "12px"},
                    ),
                ],
            ),

            # parameter selection
            dbc.Row([
                dbc.Col(
                    html.Div(
                        [
                            html.P("Select Cell Type"),
                            dcc.Dropdown(
                                id="pathway-enrichment-celltype-dropdown",
                                options=[],
                                placeholder="Select Cell Type",
                                searchable=True,
                            ),
                        ]
                    ),
                ),
                dbc.Col(
                    html.Div(
                        [
                            html.P("Select Analysis"),
                            dcc.Dropdown(
                                id="pathway-enrichment-analysis-dropdown",
                                options=[({"label": "Over Representation Analysis (ORA)", "value": "ORA"}),
                                          ({"label": "Gene Set Enrichment Analysis (GSEA)", "value": "GSEA"})],
                                placeholder="Select Analysis",
                                searchable=True,
                            ),
                        ]
                    ),
                ),
                dbc.Col(
                    html.Div(
                        [
                            html.P("Select Gene Set"),
                            dcc.Dropdown(
                                id="pathway-enrichment-gene-set-dropdown",
                                options=[
                                    ({"label": "GO", "value": "GO_Biological_Process_2026"}),
                                ],
                                placeholder="Select Gene Set",
                                searchable=True,
                            ),
                        ]
                    ),
                ),
            ]),

            # loading signal
            dcc.Loading(
                type="default",
                children=html.Div(id="pathway-enrichment-loading-signal", style={"display": "none"}),
            ),

            # results container
            dbc.Tabs(
                id="pathway-enrichment-tabs",
                active_tab="table-tab",
                children=[
                    dbc.Tab(
                        label="Results Table",
                        tab_id="table-tab",
                        children=html.Div(
                            id="pathway-enrichment-results-container", 
                            children=get_cached_result("pathway-enrichment-results-container", ["No pathway enrichment results yet. Run the analysis to see results here."])
                        ),
                    ),
                    dbc.Tab(
                        label="Dot Plot",
                        tab_id="dot-plot-tab",
                        children=html.Div(
                            id="pathway-enrichment-dot-plot-container", 
                            children=get_cached_result("pathway-enrichment-dot-plot-container", ["No dot plot yet. Run the analysis to see results here."])
                        ),
                    ),
                    dbc.Tab(
                        label="Bar Plot",
                        tab_id="bar-plot-tab",
                        children=html.Div(
                            id="pathway-enrichment-bar-plot-container", 
                            children=get_cached_result("pathway-enrichment-bar-plot-container", ["No bar plot yet. Run the analysis to see results here."])
                        ),
                    ),
                ],
            ),
            

        ], id="pathway-enrichment-panel"
    )