import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

NOTEBOOK = Path("/p1/zulab_users/jtian/my_jupyter/IRI-seq/IRIS/demo-FY1_3-pipeline.ipynb")


def run_cell(env, cell_number):
    data = json.loads(NOTEBOOK.read_text())
    source = "".join(data["cells"][cell_number - 1]["source"])
    print(f"\n===== running notebook cell {cell_number} =====", flush=True)
    exec(compile(source, f"{NOTEBOOK}#cell{cell_number}", "exec"), env)


env = {"__name__": "__main__"}

# Configuration and imports.
run_cell(env, 3)

# Reuse existing full-pipeline outputs; do not rerun FASTQ/STAR/bead reconstruction.
env["RUN_RECONSTRUCTION"] = False
env["adata_full_path"] = env["ADATA_FOLDER"] / "adata_full.h5ad"
env["spatial_rmdup_path"] = env["DEDUPLICATE_SPATIAL"] / f"{env['BEAD_SAMPLE']}.spatial.csv.gz"

# New metadata selection by row-index prefix.
run_cell(env, 15)

# Recompute expression UMAP table with the new annotation mapping.
run_cell(env, 17)

# Define spatial helper functions without recomputing the spatial UMAP.
run_cell(env, 23)

# Re-draw expression/spatial/two-panel figures with the new annotations.
run_cell(env, 25)

# Re-run validation so the common-barcode counts use the same metadata selection.
run_cell(env, 27)

print("\nRerun finished.", flush=True)
