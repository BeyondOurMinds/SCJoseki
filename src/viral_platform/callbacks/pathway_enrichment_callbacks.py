from dash import Output, Input, State, dash_table, no_update, dcc
import plotly.express as px
import numpy as np
import textwrap
from viral_platform.state.dataset_store import cache_results, get_state_store
from viral_platform.analysis.PathwayEnrichment import run_pathway_enrichment
from viral_platform.utils.runtime_tracker import get_tracker
import logging

logger = logging.getLogger(__name__)

_pathway_tracker = get_tracker("Pathway Enrichment")


def create_pe_table(enrichment_results):
    table = dash_table.DataTable(
            id="pathway-enrichment-results-table",
            columns=[{"name": col, "id": col} for col in enrichment_results.columns],
            data=enrichment_results.to_dict("records"),
            style_table={"overflowX": "auto"},
            style_cell={"textAlign": "left"},
            page_size=10,
        )
    return table

def _build_ora_dotplot(results_df):
    plot_df = results_df.copy()
    plot_df["Gene Ratio"] = (
        plot_df["Overlap"]
        .str.split("/")
        .apply(lambda x: int(x[0]) / int(x[1]))
    )
    plot_df["Gene Count"] = (
        plot_df["Overlap"]
        .str.split("/")
        .str[0]
        .astype(int)
    )
    plot_df["Wrapped Term"] = plot_df["Term"].apply(
        lambda x: "<br>".join(textwrap.wrap(str(x), width=40))
    )
    plot_df["Wrapped Genes"] = plot_df["Genes"].apply(
        lambda x: "<br>".join(textwrap.wrap(str(x), width=45))
    )
    plot_df["-log10(adj p)"] = -np.log10(plot_df["Adjusted P-value"])
    plot_df = plot_df.sort_values(
        "Adjusted P-value"
    )
    plot_df = plot_df.head(20)
    fig = px.scatter(
        plot_df,
        x="Gene Ratio",
        y="Wrapped Term",
        size="Gene Count",
        color="-log10(adj p)",
        hover_name="Term",
        hover_data={
            "Adjusted P-value":":.2e",
            "Gene Count":True,
            "Wrapped Genes":True,
            "Genes":False,
            "Combined Score":":.2f",
        },
    )
    fig.update_yaxes(
        autorange="reversed",
        categoryorder="array",
        categoryarray=plot_df["Wrapped Term"].tolist(),
    )
    max_lines = plot_df["Wrapped Term"].str.count("<br>").max() + 1
    height = max(600, len(plot_df) * (25 + 12 * max_lines))
    fig.update_layout(
        template="plotly_white",
        title="Pathway Enrichment Dot Plot (ORA)",
        xaxis_title="Gene Ratio",
        yaxis_title="Pathway Term",
        height=height,
        margin=dict(l=300, r=40, t=60, b=60),
    )
    return fig

def _build_gsea_dotplot(results_df, fdr_cutoff=0.05):
    plot_df = results_df.copy()
    plot_df["Gene Count"] = (
        plot_df["Lead_genes"]
        .str.split(";")
        .str.len()
    )
    plot_df["Wrapped Term"] = plot_df["Term"].apply(
        lambda x: "<br>".join(textwrap.wrap(str(x), width=40))
    )
    plot_df["Wrapped Lead Genes"] = plot_df["Lead_genes"].apply(
        lambda x: "<br>".join(textwrap.wrap(str(x), width=45))
    )
    plot_df["-log10(FDR)"] = -np.log10(plot_df["FDR q-val"])
    plot_df = plot_df[
        plot_df["FDR q-val"] <= fdr_cutoff
    ].copy()

    plot_df["absNES"] = plot_df["NES"].abs()

    plot_df = plot_df.sort_values(
        ["absNES", "FDR q-val"],
        ascending=[False, True]
    )
    plot_df = plot_df.head(20)
    fig = px.scatter(
        plot_df,
        x="NES",
        y="Wrapped Term",
        size="Gene Count",
        color="NES",
        hover_name="Term",
        hover_data={
            "FDR q-val":":.2e",
            "NES":":.3f",
            "ES":":.3f",
            "Gene Count":True,
            "Wrapped Lead Genes":True,
            "Lead_genes":False,
        },
    )
    fig.update_yaxes(
        autorange="reversed",
        categoryorder="array",
        categoryarray=plot_df["Wrapped Term"].tolist(),
    )
    max_lines = plot_df["Wrapped Term"].str.count("<br>").max() + 1
    height = max(600, len(plot_df) * (25 + 12 * max_lines))
    fig.update_layout(
        template="plotly_white",
        title="Pathway Enrichment Dot Plot (GSEA)",
        xaxis_title="NES",
        yaxis_title="Pathway Term",
        height=height,
        margin=dict(l=350, r=40, t=60, b=60),
    )
    return fig

def _build_ora_barplot(results_df):
    plot_df = results_df.copy()

    plot_df["Wrapped Term"] = plot_df["Term"].apply(
        lambda x: "<br>".join(textwrap.wrap(str(x), width=40))
    )

    plot_df["Wrapped Genes"] = plot_df["Genes"].apply(
        lambda x: "<br>".join(textwrap.wrap(str(x), width=45))
    )

    plot_df["-log10(adj p)"] = -np.log10(plot_df["Adjusted P-value"])

    plot_df = plot_df.sort_values("Adjusted P-value").head(20)

    fig = px.bar(
        plot_df,
        x="-log10(adj p)",
        y="Wrapped Term",
        orientation="h",
        color="Combined Score",
        hover_name="Term",
        hover_data={
            "Adjusted P-value": ":.2e",
            "Combined Score": ":.2f",
            "Wrapped Genes": True,
            "Genes": False,
        },
        color_continuous_scale="Viridis",
    )
    fig.update_yaxes(
        autorange="reversed",
        categoryorder="array",
        categoryarray=plot_df["Wrapped Term"].tolist(),
    )
    
    height = max(600, 45 * len(plot_df))

    fig.update_layout(
        template="plotly_white",
        title="Pathway Enrichment Bar Plot (ORA)",
        xaxis_title="-log10(Adjusted P-value)",
        yaxis_title="Pathway Term",
        height=height,
        margin=dict(l=300, r=40, t=60, b=60),
    )
    return fig

def _build_gsea_barplot(results_df, fdr_cutoff=0.05):
    plot_df = results_df.copy()

    plot_df["Wrapped Term"] = plot_df["Term"].apply(
        lambda x: "<br>".join(textwrap.wrap(str(x), width=40))
    )

    plot_df["Wrapped Lead Genes"] = plot_df["Lead_genes"].apply(
        lambda x: "<br>".join(textwrap.wrap(str(x), width=45))
    )

    plot_df["-log10(FDR)"] = -np.log10(plot_df["FDR q-val"])

    plot_df["absNES"] = plot_df["NES"].abs()

    plot_df = plot_df[plot_df["FDR q-val"] <= fdr_cutoff].copy()
    plot_df = plot_df.sort_values(["absNES", "FDR q-val"], ascending=[False, True]).head(20)

    fig = px.bar(
        plot_df,
        x="NES",
        y="Wrapped Term",
        orientation="h",
        color="NES",
        hover_name="Term",
        hover_data={
            "FDR q-val":":.2e",
            "NES":":.3f",
            "ES":":.3f",
            "Wrapped Lead Genes": True,
            "Lead_genes": False,
        },
        color_continuous_scale="RdBu_r",
    )
    fig.update_yaxes(
        autorange="reversed",
        categoryorder="array",
        categoryarray=plot_df["Wrapped Term"].tolist(),
    )
    
    height = max(600, 45 * len(plot_df))

    fig.update_layout(
        template="plotly_white",
        title="Pathway Enrichment Bar Plot (GSEA)",
        xaxis_title="Normalised Enrichment Score (NES)",
        yaxis_title="Pathway Term",
        height=height,
        margin=dict(l=300, r=40, t=60, b=60),
    )
    return fig


def register_pathway_enrichment_callbacks(app):
    @app.callback(
        Output("pathway-enrichment-advanced-options-collapse", "is_open"),
        Input("pathway-enrichment-advanced-options-button", "n_clicks"),
        State("pathway-enrichment-advanced-options-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_pathway_advanced_options(n_clicks, is_open):
        if not n_clicks:
            return is_open
        return not is_open

    @app.callback(
        Output("pathway-enrichment-pvalue-cutoff-input", "value"),
        Output("pathway-enrichment-logfc-cutoff-input", "value"),
        Input("pathway-enrichment-advanced-reset-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_pathway_advanced_defaults(n_clicks):
        if not n_clicks:
            return no_update, no_update
        return 0.05, 1.0

    @app.callback(
        Output("pathway-enrichment-celltype-dropdown", "options"),
        Output("pathway-enrichment-celltype-dropdown", "value"),
        Input("active-dataset-version", "data"),
    )
    def populate_pathway_enrichment_celltype_dropdown(dataset_version):
        """
        Populate the cell type dropdown for pathway enrichment based on available DE results.
        """
        # Access the global state to get the DE results by cell type
        state = get_state_store()
        de_results_by_celltype = state.get("DE_results", {}).get("results_by_celltype", {})
        
        # Extract cell types from the keys of the DE results dictionary
        cell_types = list(de_results_by_celltype.keys())
        
        # Create options for the dropdown
        options = [{"label": ct, "value": ct} for ct in cell_types]
        
        # Set a default value if there are available cell types
        default_value = cell_types[0] if cell_types else None
        
        return options, default_value
    
    @app.callback(
        Output("pathway-enrichment-loading-signal", "children"),
        Output("pathway-enrichment-results-container", "children"),
        Output("pathway-enrichment-dot-plot-container", "children"),
        Output("pathway-enrichment-bar-plot-container", "children"),
        Input("run-pathway-enrichment-button", "n_clicks"),
        State("pathway-enrichment-celltype-dropdown", "value"),
        State("pathway-enrichment-analysis-dropdown", "value"),
        State("pathway-enrichment-gene-set-dropdown", "value"),
        State("pathway-enrichment-pvalue-cutoff-input", "value"),
        State("pathway-enrichment-logfc-cutoff-input", "value"),
    )
    @_pathway_tracker.track
    def run_pathway_enrichment_callback(
        n_clicks,
        celltype,
        method,
        gene_set,
        pvalue_cutoff,
        logfc_cutoff,
    ):
        """
        Run pathway enrichment analysis based on the selected cell type, method, and gene set.
        """
        if n_clicks is None or n_clicks == 0:
            return no_update, no_update, no_update, no_update
        
        # Access the global state to get the DE results by cell type
        state = get_state_store()
        de_results_by_celltype = state.get("DE_results", {}).get("results_by_celltype", {})
        
        # Check if DE results for the selected cell type are available
        if celltype not in de_results_by_celltype:
            logger.warning(f"No DE results available for cell type: {celltype}")
            return no_update, f"No DE results available for cell type: {celltype}", f"No DE results available for cell type: {celltype}", f"No DE results available for cell type: {celltype}"
        
        # Get the DE results for the selected cell type
        de_results = de_results_by_celltype[celltype]
        
        # Run pathway enrichment analysis
        with _pathway_tracker.phase("analysis"):
            enrichment_payload = run_pathway_enrichment(
                de_results,
                method,
                gene_set,
                pvalue_cutoff=pvalue_cutoff if pvalue_cutoff is not None else 0.05,
                logfc_cutoff=logfc_cutoff if logfc_cutoff is not None else 1.0,
            )
        if enrichment_payload is None:
            logger.warning("Pathway enrichment analysis returned no results.")
            return no_update, "No pathway enrichment results found.", "No pathway enrichment results", "No pathway enrichment results"

        enrichment_results, enr = enrichment_payload
        
        if enrichment_results is None or enrichment_results.empty:
            logger.warning("Pathway enrichment analysis returned no results.")
            return no_update, "No pathway enrichment results found.", "No pathway enrichment results", "No pathway enrichment results"

        
        # Cache the results for future use
        with _pathway_tracker.phase("processing"):
            table = create_pe_table(enrichment_results)
            if method == "ORA":
                dotplot_fig = _build_ora_dotplot(enrichment_results)
                dotplot_graph = dcc.Graph(figure=dotplot_fig)
                barplot_fig = _build_ora_barplot(enrichment_results)
                barplot_graph = dcc.Graph(figure=barplot_fig)
            elif method == "GSEA":
                dotplot_fig = _build_gsea_dotplot(
                    enrichment_results,
                    fdr_cutoff=pvalue_cutoff if pvalue_cutoff is not None else 0.05,
                )
                dotplot_graph = dcc.Graph(figure=dotplot_fig)
                barplot_fig = _build_gsea_barplot(
                    enrichment_results,
                    fdr_cutoff=pvalue_cutoff if pvalue_cutoff is not None else 0.05,
                )
                barplot_graph = dcc.Graph(figure=barplot_fig)


        with _pathway_tracker.phase("visualization"):
            cache_results(**{
                "pathway-enrichment-results-container": table,
                "pathway-enrichment-dot-plot-container": dotplot_graph,
                "pathway-enrichment-bar-plot-container": barplot_graph
            })

        return "done", table, dotplot_graph, barplot_graph