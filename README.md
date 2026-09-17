# <img width="100" height="100" alt="SCJosekiLogo" src="https://github.com/user-attachments/assets/8b281b9a-40d7-464c-a0d4-66f6eb123f48" /> SCJoseki: A Viral Single-Cell RNA-seq Analysis Platform

## About
SCJoseki is a user-friendly Plotly Dash web application for the analysis of single-cell RNA-seq datasets with additional support for viral datasets.

The application contains many of the standard modular workflows found in scRNA-seq pipelines, while addressing the current gaps in viral analysis functionality.

<img width="1121" height="510" alt="FinalNewPlotBabyBoy" src="https://github.com/user-attachments/assets/5345cb3b-79e5-4611-929a-f60566f59440" />


An overview of the application and its modules can be found on the first overview page whenever you load the application. 
Currently, the application only supports H5AD files or zipped 10x Genomics Cell Ranger files.


<img width="1896" height="908" alt="Screenshot 2026-08-26 211800" src="https://github.com/user-attachments/assets/545abd10-0954-464d-a166-d4fd4d31a942" />


Almost every module offers some form of data visualisation, such as violin, bar, volcano, UMAP, and network plots

<img width="1892" height="913" alt="Screenshot 2026-08-26 214142" src="https://github.com/user-attachments/assets/8f5effd3-8f46-496b-b42c-397af7ed759d" />


## Testing
The dataset used for the development and testing of SCJoseki can be found within the assets folder, called TestEBVDataset.h5ad, 

### Example Results

Using the testing dataset, the following results were produced for the differential expression/pseudobulk, pathway enrichment, and viral analysis modules:

#### DE & Pseudobulk Analysis
<img width="1657" height="2738" alt="VolanoReplace" src="https://github.com/user-attachments/assets/ad91b67b-aea0-4aca-98d1-76361118e168" />

Results from DE analysis across cell types **A)** Differential expression analysis of memory B cells, showing the corresponding volcano plot and heatmap of the top differentially expressed genes **B)** Differential expression analysis of Age-associated B cells, showing the corresponding volcano plot and heatmap of the top differentially expressed genes **C)** Heatmap showing the expression of top differentially expressed genes across the identified cell populations to compare expression patterns between populations. The volcano plots display log2FC against -log10(p-adj), with significantly upregulated (red) and downregulated (blue) genes distinguished from non-significant (gray) genes. Heatmaps show the top 20 differentially expressed genes, where the cell specific ones are coloured by Z-score and the cross-population is coloured by average expression. 

#### Pathway Enrichment
<img width="1343" height="1890" alt="Figure4" src="https://github.com/user-attachments/assets/2e3b7fe9-a456-48cc-8406-33e17cd83d9d" />

Pathway enrichment analysis using over-representation analysis (ORA) and gene set enrichment analysis (GSEA). **A)** ORA bar plot and dot plot. Bar plot shows the most significantly enriched pathways, with bar length representing gene ratio and colour representing -log10padj value. Dot plot shows enriched pathways, with the dot position representing the gene ratio and colour representing -log10padj value. **B)** GSEA bar and dot plots. Bar plots pathway enrichment ranked by NES, with bar length representing normalised enrichment score (NES) and colour representing -log10padj value. Dot plot shows the most enriched pathways, with dot position representing NES and colour representing -log10padj value.

#### Viral Burden
<img width="1400" height="1175" alt="Figure5" src="https://github.com/user-attachments/assets/304f9697-ff85-44da-9d11-1356bc5ebeb2" />

UMAP visualisation of cell type clustering, viral infection, and viral burden. **A)** UMAP of the analysed single-cell dataset coloured according to cell type annotations predicted using CellTypist majority voting. **B)** UMAP showing infection status, with infected cells (green) and non-infected cells (gray) based on detection of viral reads. **C)** UMAP coloured by continuous viral burden. Each point represents an individual cell, with cells positioned according to their transcriptional similarity in the UMAP embedding

<img width="1400" height="1000" alt="Figure6" src="https://github.com/user-attachments/assets/469d2985-67a6-4bff-b831-67799a907fbc" />

Distribution of EBV infection and viral burden across cell populations, conditions, and samples. **A)** Proportion of cells classified as EBV infected within each major cell population. **B)** Viral burden percentage between experimental conditions. **C)** Viral burden percentage across major cell populations. D) Viral burden across individual samples, identified by their GEO accession numbers




## Installation
SCJoseki can be installed in 3 different ways:

1. The EXE available within the GitHub releases.
2. Through Docker
```bash
docker pull beyondourminds/scjoseki:latest
```
3. Via PIP
   a. Clone the repo:
   ```bash
   git clone https://github.com/BeyondOurMinds/SCJoseki 
   ```
   b. Switch to the current directory:
   ```bash
   cd SCJoseki
   ```
   c. Create a virtual environment:
   ```bash
   python.exe -m venv .venv
   ```
   d. Activate the virtual environment:
   
   **Windows (PowerShell):**
    ```powershell
    .venv\Scripts\Activate.ps1
    ```
    
    **Windows (Command Prompt):**
    ```cmd
    .venv\Scripts\activate.bat
    ```
    
    **Linux/Mac:**
    ```bash
    source .venv/bin/activate
    ```

    e. Install the package:
    ```bash
    pip install .
    ```
    

## How To Run
1. **EXE File**

Double-click the file, and it runs!

2. **Docker**

To run with Docker, follow the functional example below. -v mounts a location for exported files. The code below will create a folder, if one does not already exist, in your current working directory called "exports", where the outputs will be saved. To choose/create a different export folder, update the path before the semicolon (:). 
```bash
docker run --rm -p 8050:8050 -e CHROME_BIN=/usr/bin/chromium -v "${PWD}/exports:/exports" scjoseki
```
Once run, navigate to https://localhost:8050 to see the running app.

3. Pip

Within the package's root directory, run:
```bash
scjoseki
```

## Prerequisites:
Python < 3.14

