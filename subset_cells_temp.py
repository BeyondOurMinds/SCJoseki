"""Temporary script: randomly subset a dataset to N cells and save as h5ad. Delete when done."""
from pathlib import Path

import numpy as np

INPUT_PATH = r"C:\Users\jtspy\OneDrive\Desktop\Bioinformatics\ViralDatasets\EpsteinBarr\GSE165080_HC_LCL_annotated_compatible.h5ad"  # placeholder: .zip (10x) or .h5ad
N_CELLS = 10000  # placeholder: number of cells to keep
OUTPUT_PATH = None  # None -> "<input stem>_subset<N>.h5ad" next to the input
SEED = 0

# --- Loader: keep ONE of the two blocks below, comment out the other ---

# Option A: 10x loader from scjoseki (zip of 10x files)
# from viral_platform.io.loaders import load_file_from_path
# adata = load_file_from_path(INPUT_PATH)

# Option B: standard h5ad loading
import scanpy as sc
adata = sc.read_h5ad(INPUT_PATH)

# --- Subset ---
if N_CELLS >= adata.n_obs:
    print(f"Requested {N_CELLS} cells but dataset has {adata.n_obs}; keeping all cells.")
    subset = adata.copy()
else:
    rng = np.random.default_rng(SEED)
    idx = np.sort(rng.choice(adata.n_obs, size=N_CELLS, replace=False))
    subset = adata[idx].copy()

import anndata
anndata.settings.allow_write_nullable_strings = True  # pandas string dtypes fail to write otherwise

out = Path(OUTPUT_PATH) if OUTPUT_PATH else Path(INPUT_PATH).with_name(f"{Path(INPUT_PATH).stem}_subset{subset.n_obs}.h5ad")
subset.write_h5ad(out)
print(f"Saved {subset.n_obs} cells x {subset.n_vars} genes to {out}")
