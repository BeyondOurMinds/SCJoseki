import gseapy as gp
import pandas as pd
import logging

logger = logging.getLogger(__name__)

def run_pathway_enrichment(de_results, method, gene_set="GO_Biological_Process_2023", pvalue_cutoff=0.05, logfc_cutoff=1.0):
    """
    Perform pathway enrichment analysis using GSEApy.

    Parameters
    ----------
    de_results : pd.DataFrame
        DataFrame containing differential expression results with 'gene' and 'log2FoldChange' columns.
    method : str
        The enrichment method to use ('GSEA', 'ORA', etc.).
    gene_set : str
        The gene set to use for enrichment analysis.
    pvalue_cutoff : float
        The p-value cutoff for filtering significant results.
    logfc_cutoff : float
        The log fold change cutoff for filtering significant results.
    """

    # Perform enrichment analysis using GSEApy
    try:
        if method == 'GSEA':
           
            # Prepare gene list for GSEA (ranked by logFC)
            gene_list = de_results[['gene', 'stat']].dropna().drop_duplicates(subset='gene').sort_values(by='stat', ascending=False)

            if gene_list is None or len(gene_list) == 0:
                logger.warning("No genes available for enrichment analysis after filtering.")
                return None
            
            enr = gp.prerank(rnk=gene_list, gene_sets=gene_set, outdir=None)
            results_df = enr.res2d
            print(results_df.columns.tolist())
            print(results_df.head())
            # gseapy may return stats as object dtype
            results_df['FDR q-val'] = pd.to_numeric(results_df['FDR q-val'], errors='coerce')
            logger.info("GSEA FDR cutoff=%s; min FDR=%s", pvalue_cutoff, results_df['FDR q-val'].min())
            significant_results = results_df[results_df['FDR q-val'] <= pvalue_cutoff]
        elif method == 'ORA':

            # Filter DE results based on logFC and p-value cutoffs
            filtered_results = de_results[
                (de_results['log2FoldChange'].abs() >= logfc_cutoff) &
                (de_results['padj'] <= pvalue_cutoff)
            ]

            if filtered_results.empty:
                logger.warning("No significant genes found after filtering with logFC and p-value cutoffs.")
                return None
            
            # Prepare gene list for enrichment analysis
            gene_list = filtered_results['gene'].dropna().drop_duplicates().tolist()

            if gene_list is None or len(gene_list) == 0:
                logger.warning("No genes available for enrichment analysis after filtering.")
                return None

            enr = gp.enrichr(gene_list=gene_list, gene_sets=gene_set, organism='human', outdir=None)
            results_df = enr.results
            print(results_df.columns.tolist())
            print(results_df.head())
            # Filter results based on p-value cutoff
            significant_results = results_df[results_df['Adjusted P-value'] <= pvalue_cutoff]
        else:
            logger.error("Unsupported enrichment method: %s", method)
            return None
    except Exception as e:
        logger.error("Error during enrichment analysis: %s", str(e))
        return None

    

    if significant_results.empty:
        logger.warning("No significant pathways found after enrichment analysis.")
        return None

    return significant_results, enr