from dash import html, dcc
import dash_bootstrap_components as dbc
from viral_platform.state.dataset_store import get_cached_result

def create_host_virus_interaction_panel():
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
                                ["4. Host-Virus Interaction Analysis"],
                                style={"display": "inline-flex", "alignItems": "center", "marginBottom": "2px"},
                            ),
                            html.P(
                                "Identify host genes associated with viral burden",
                                style={"color": "#6c757d", "marginBottom": "16px", "fontSize": "14px"},
                            ),
                        ]
                    ),
                ),
                dbc.Col(
                    dbc.Button("Run Host-Virus Interaction Analysis", id="run-host-virus-interaction-analysis-button", n_clicks=0, color="primary", className="mb-3"),
                    width="auto"
                ),
                dbc.Col(
                    dbc.Button("Advanced Options", id="host-virus-advanced-options-button", n_clicks=0, color="secondary", className="mb-3"),
                    width="auto"
                )
            ]),
            dbc.Collapse(
                id="host-virus-advanced-options-collapse",
                is_open=False,
                children=[
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    html.P("Minimum expressing cells per host gene"),
                                    dcc.Input(
                                        id="host-virus-min-cells-input",
                                        type="number",
                                        min=1,
                                        max=1000,
                                        step=1,
                                        value=10,
                                        style={"width": "100%"},
                                    ),
                                ],
                                md=4,
                            ),
                            dbc.Col(
                                [
                                    html.P("Adjusted p-value significance cutoff"),
                                    dcc.Input(
                                        id="host-virus-adj-p-cutoff-input",
                                        type="number",
                                        min=0,
                                        max=1,
                                        step="any",
                                        value=0.05,
                                        style={"width": "100%"},
                                    ),
                                ],
                                md=4,
                            ),
                            dbc.Col(
                                [
                                    html.P("Absolute correlation cutoff"),
                                    dcc.Input(
                                        id="host-virus-corr-cutoff-input",
                                        type="number",
                                        min=0,
                                        max=1,
                                        step=0.01,
                                        value=0.15,
                                        style={"width": "100%"},
                                    ),
                                ],
                                md=4,
                            ),
                        ],
                        style={"marginBottom": "12px"},
                    ),
                    dbc.Button(
                        "Reset Advanced Defaults",
                        id="host-virus-advanced-reset-button",
                        n_clicks=0,
                        color="secondary",
                        outline=True,
                        className="mb-3",
                    ),
                ],
            ),
            dbc.Row([
                dcc.Dropdown(
                    id="host-virus-interaction-dropdown",
                    options=[],
                )
            ]),
            # loading signal
            dcc.Loading(
                type="default",
                children=html.Div(id="host-virus-interaction-loading-signal", style={"display": "none"}),
            ),
            # summary box
            html.Div(
                id="host-virus-interaction-summary-container",
                children=get_cached_result("host-virus-interaction-summary-container", [
                    html.P(
                        "No host-virus interaction results yet. Run the analysis to see results here.",
                        style={"color": "#6c757d", "fontSize": "14px", "marginTop": "10px"},
                    )
                ]),
            ),
            # results container
            html.Div(
                id="host-virus-interaction-results-container",
                children=get_cached_result("host-virus-interaction-results-container", [
                    html.P(
                        "",
                        style={"color": "#6c757d", "fontSize": "14px"},
                    )
                ]),
            ),
        ]
    )
