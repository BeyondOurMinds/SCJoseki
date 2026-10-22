from dash import html, dcc
import dash_bootstrap_components as dbc
from viral_platform.state.dataset_store import get_cached_result

def create_viral_burden_panel():
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
                                ["2. Viral Burden Analysis"],
                                style={"display": "inline-flex", "alignItems": "center", "marginBottom": "2px"},
                            ),
                            html.P(
                                "Calculate the viral burden from the detected viral genes",
                                style={"color": "#6c757d", "marginBottom": "16px", "fontSize": "14px"},
                            ),
                        ]
                    ),
                ),
                dbc.Col(
                    dbc.Button("Run Viral Burden Analysis", id="run-viral-burden-analysis-button", n_clicks=0, color="primary", className="mb-3"),
                    width="auto"
                ),
                dbc.Col(
                    dbc.Button("Advanced Options", id="viral-burden-advanced-options-button", n_clicks=0, color="secondary", className="mb-3"),
                    width="auto"
                ),
            ]),

            dbc.Collapse(
                id="viral-burden-advanced-options-collapse",
                is_open=False,
                children=[
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Dropdown(
                                    id="viral-burden-celltype-column-dropdown",
                                    options=[
                                        {"label": "Upload a dataset to select cell type column", "value": ""},
                                    ],
                                    placeholder="Select Cell Type Column (optional)",
                                    searchable=True,
                                    clearable=True,
                                ),
                                md=4,
                            ),
                            dbc.Col(
                                dcc.Dropdown(
                                    id="viral-burden-condition-column-dropdown",
                                    options=[
                                        {"label": "Upload a dataset to select condition column", "value": ""},
                                    ],
                                    placeholder="Select Condition Column (optional)",
                                    searchable=True,
                                    clearable=True,
                                ),
                                md=4,
                            ),
                            dbc.Col(
                                dcc.Dropdown(
                                    id="viral-burden-sample-column-dropdown",
                                    options=[
                                        {"label": "Upload a dataset to select sample column", "value": ""},
                                    ],
                                    placeholder="Select Sample Column (optional)",
                                    searchable=True,
                                    clearable=True,
                                ),
                                md=4,
                            ),
                        ],
                        style={"marginBottom": "12px"},
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    html.P("Infected cell threshold (viral counts)"),
                                    dcc.Input(
                                        id="viral-burden-infected-threshold-input",
                                        type="number",
                                        min=0,
                                        step=1,
                                        value=0,
                                        style={"width": "100%"},
                                    ),
                                ],
                                md=3,
                            ),
                            dbc.Col(
                                [
                                    html.P("Association min cells per host gene"),
                                    dcc.Input(
                                        id="viral-burden-association-min-cells-input",
                                        type="number",
                                        min=1,
                                        max=1000,
                                        step=1,
                                        value=10,
                                        style={"width": "100%"},
                                    ),
                                ],
                                md=3,
                            ),
                            dbc.Col(
                                [
                                    html.P("Association correlation cutoff"),
                                    dcc.Input(
                                        id="viral-burden-association-corr-threshold-input",
                                        type="number",
                                        min=0,
                                        max=1,
                                        step=0.01,
                                        value=0.3,
                                        style={"width": "100%"},
                                    ),
                                ],
                                md=3,
                            ),
                            dbc.Col(
                                [
                                    html.P("Association FDR cutoff"),
                                    dcc.Input(
                                        id="viral-burden-association-fdr-threshold-input",
                                        type="number",
                                        min=0,
                                        max=1,
                                        step="any",
                                        value=0.05,
                                        style={"width": "100%"},
                                    ),
                                ],
                                md=3,
                            ),
                        ],
                        style={"marginBottom": "12px"},
                    ),
                    dbc.Button(
                        "Reset Advanced Defaults",
                        id="viral-burden-advanced-reset-button",
                        n_clicks=0,
                        color="secondary",
                        outline=True,
                        className="mb-3",
                    ),
                ],
            ),

            # loading signal
            dcc.Loading(
                type="default",
                children=html.Div(id="viral-burden-loading-signal", style={"display": "none"}),
            ),

            dbc.Tabs(
                id="viral-burden-output-tabs",
                active_tab="viral-burden-results-tab",
                children=[
                    dbc.Tab(
                        label="Results Table",
                        tab_id="viral-burden-results-tab",
                        children=html.Div(
                            id="viral-burden-results-container",
                            children=get_cached_result("viral-burden-results-container", [
                                html.P(
                                    "No viral burden results yet. Run the analysis to see results here.",
                                    style={"color": "#6c757d", "fontSize": "14px"},
                                )
                            ]),
                            style={"padding": "15px"},
                        ),
                    ),
                    dbc.Tab(
                        label="Infection UMAP",
                        tab_id="viral-burden-infection-umap-tab",
                        children=html.Div(
                            id="viral-burden-infection-umap-container",
                            children=get_cached_result("viral-burden-infection-umap-container", "Infection UMAP output will appear here."),
                            style={"padding": "15px"},
                        ),
                    ),
                    dbc.Tab(
                        label="Viral Burden UMAP",
                        tab_id="viral-burden-umap-tab",
                        children=html.Div(
                            id="viral-burden-umap-container",
                            children=get_cached_result("viral-burden-umap-container", "Viral burden UMAP output will appear here."),
                            style={"padding": "15px"},
                        ),
                    ),
                    dbc.Tab(
                        label="Violin Plots",
                        tab_id="viral-burden-violin-tab",
                        children=html.Div(
                            id="viral-burden-violin-container",
                            children=get_cached_result("viral-burden-violin-container", "Violin plot output will appear here."),
                            style={"padding": "15px"},
                        ),
                    ),
                    dbc.Tab(
                        label="Cell Type Infection Fraction",
                        tab_id="viral-burden-celltype-fraction-tab",
                        children=html.Div(
                            id="viral-burden-celltype-fraction-container",
                            children=get_cached_result("viral-burden-celltype-fraction-container", "Cell type infection fraction output will appear here."),
                            style={"padding": "15px"},
                        ),
                    ),
                ],
            ),
            # Viral Burden Associations
            html.Div(
                id="viral-burden-associations-container",
                children=[
                    dbc.Row([
                        dbc.Col(
                            children=[
                                html.Div(
                                    id="viral-burden-associations-title-container",
                                    children=get_cached_result("viral-burden-associations-results-container", [
                                        html.H4(
                                            ["2.1 Viral Burden Association"],
                                            style={"display": "inline-flex", "alignItems": "center", "marginBottom": "2px"},
                                        ),
                                        html.P(
                                            "Find viral burden associations or correlations",
                                            style={"color": "#6c757d", "marginBottom": "16px", "marginTop": "8px", "fontSize": "14px"},
                                        ),
                                    ]),
                                ),
                            ],
                        ),
                        dbc.Col(
                            dbc.Button("Run Viral Burden Association", id="run-viral-burden-association-button", n_clicks=0, color="primary", className="mb-3"),
                            width="auto"
                        ),
                        
                        
                    ]),

                    # loading signal
                    dcc.Loading(
                        type="default",
                        children=html.Div(id="viral-burden-associations-loading-signal", style={"display": "none"}),
                    ),

                    # results container
                    dbc.Row([
                        dbc.Col(
                            children=[
                                html.Div(
                                    id="viral-burden-associations-results-container",
                                    children=get_cached_result("viral-burden-associations-significant-results-container", [
                                        html.P(
                                            "No viral burden associations results yet. Run the analysis to see results here.",
                                        ),
                                    ]),
                                ),
                                html.Div(
                                    id="viral-burden-associations-significant-results-container",
                                    children=[
                                        html.P(
                                            "",
                                        ),
                                    ],
                                ),
                            ],
                        )
                    ]),
                ],
                style={"marginTop": "20px", "padding": "10px", "border": "1px solid #dee2e6", "borderRadius": "5px"},
            )
        ]
    )
