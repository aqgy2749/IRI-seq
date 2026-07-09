import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.spatial import procrustes
from sklearn.neighbors import NearestNeighbors

SAMPLE_ID = "FY1_3"
INPUT_DIR = Path("/p2/zulab/jtian/data/IRISeq/demo_FY1_3/data")
OUTPUT_DIR = Path("/p2/zulab/jtian/data/IRISeq/demo_FY1_3/output/demo-FY1_3-pipeline")

metadata_path = INPUT_DIR / "GSE270383_meta_data.csv.gz"
processed_count_mtx = INPUT_DIR / "GSE270383_count.mtx.gz"
processed_barcodes_tsv = INPUT_DIR / "GSE270383_barcodes.tsv.gz"
processed_genes_tsv = INPUT_DIR / "GSE270383_genes.tsv.gz"
adata_path = OUTPUT_DIR / "cDNA/Adata/adata_full.h5ad"
spatial_umap_path = OUTPUT_DIR / "connections/20231221_FY1_3_connection_spatial_reconstruction_umap.csv"
validation_dir = OUTPUT_DIR / "validation"
summary_json = validation_dir / "validation_summary.json"
summary_md = validation_dir / "validation_summary.md"


def strip_sample_prefix(value):
    s = str(value)
    if "." in s:
        s = s.split(".", 1)[1]
    if "-" in s:
        head, tail = s.rsplit("-", 1)
        if tail.isdigit():
            s = head
    return s


def compact_barcode(value):
    return strip_sample_prefix(value).replace("-", "")


def safe_corr(x, y, method="pearson"):
    x = pd.Series(x, dtype="float64")
    y = pd.Series(y, dtype="float64")
    mask = x.notna() & y.notna()
    if mask.sum() < 3:
        return np.nan
    return x[mask].corr(y[mask], method=method)


validation_dir.mkdir(parents=True, exist_ok=True)
metrics = {}
if summary_json.exists():
    metrics.update(json.loads(summary_json.read_text()))

metadata_all = pd.read_csv(metadata_path, index_col=0)
metadata_index_mask = metadata_all.index.astype(str).str.contains(SAMPLE_ID, regex=False)
metadata_fy = metadata_all[metadata_index_mask].copy()
metadata_fy["receiver_barcode"] = metadata_fy.index.map(compact_barcode)
metadata_fy = metadata_fy.drop_duplicates("receiver_barcode", keep="first")
fy_compact_set = set(metadata_fy["receiver_barcode"].astype(str))

reproduced = sc.read_h5ad(adata_path)
reproduced.obs["receiver_barcode"] = reproduced.obs_names.map(compact_barcode)
reproduced_totals = pd.Series(
    np.asarray(reproduced.X.sum(axis=1)).ravel(),
    index=reproduced.obs["receiver_barcode"].astype(str),
)
reproduced_genes = pd.Index(reproduced.var_names.astype(str))

processed_barcodes = pd.read_csv(processed_barcodes_tsv, sep="\t", header=None)[0].astype(str)
processed_genes_df = pd.read_csv(processed_genes_tsv, sep="\t", header=None)
processed_gene_names = processed_genes_df.iloc[:, -1].astype(str)
processed_compact = processed_barcodes.map(compact_barcode)
fy_col_mask = processed_compact.isin(fy_compact_set).to_numpy()

metrics["metadata_match_source"] = f"index contains {SAMPLE_ID}"
metrics["metadata_index_contains_FY1_3_rows"] = int(metadata_index_mask.sum())
metrics["metadata_FY1_3_after_barcode_dedup"] = int(metadata_fy.shape[0])
metrics["processed_total_barcodes"] = int(processed_barcodes.shape[0])
metrics["processed_FY1_3_barcodes_from_metadata"] = int(fy_col_mask.sum())
metrics["reproduced_cDNA_barcodes"] = int(reproduced.n_obs)
metrics["reproduced_cDNA_genes"] = int(reproduced.n_vars)
metrics["processed_reproduced_gene_name_overlap"] = int(len(set(reproduced_genes) & set(processed_gene_names)))

if fy_col_mask.sum() > 0:
    processed_mtx = mmread(processed_count_mtx).tocsr()
    processed_fy = processed_mtx[:, fy_col_mask]
    processed_totals = pd.Series(
        np.asarray(processed_fy.sum(axis=0)).ravel(),
        index=processed_compact[fy_col_mask].to_numpy().astype(str),
    ).groupby(level=0).sum()
    common = reproduced_totals.index.intersection(processed_totals.index)
    metrics["cDNA_common_barcodes_for_total_UMI"] = int(len(common))
    metrics["cDNA_total_UMI_pearson"] = safe_corr(reproduced_totals.loc[common], processed_totals.loc[common], "pearson")
    metrics["cDNA_total_UMI_spearman"] = safe_corr(reproduced_totals.loc[common], processed_totals.loc[common], "spearman")
    pd.DataFrame(
        {
            "receiver_barcode": common,
            "reproduced_total_UMI": reproduced_totals.loc[common].to_numpy(),
            "processed_total_UMI": processed_totals.loc[common].to_numpy(),
        }
    ).to_csv(validation_dir / "cDNA_total_UMI_reproduced_vs_processed.csv", index=False)

spatial_df = pd.read_csv(spatial_umap_path)
spatial_df["receiver_barcode"] = spatial_df["receiver_barcode"].astype(str)
ref_coords = metadata_fy[["receiver_barcode", "UMAP1_spatial", "UMAP2_spatial", "Annotation"]].copy()
compare = spatial_df.merge(ref_coords, on="receiver_barcode", how="inner")
compare = compare.dropna(subset=["spatial_UMAP1", "spatial_UMAP2", "UMAP1_spatial", "UMAP2_spatial"])
metrics["spatial_common_barcodes_for_coordinate_validation"] = int(compare.shape[0])
if compare.shape[0] >= 3:
    ref = compare[["UMAP1_spatial", "UMAP2_spatial"]].astype(float).to_numpy()
    rep = compare[["spatial_UMAP1", "spatial_UMAP2"]].astype(float).to_numpy()
    _, rep_aligned, disparity = procrustes(ref, rep)
    ref_scaled = (ref - ref.mean(axis=0)) / ref.std(axis=0)
    metrics["spatial_procrustes_disparity"] = float(disparity)
    metrics["spatial_procrustes_axis1_pearson"] = safe_corr(ref_scaled[:, 0], rep_aligned[:, 0], "pearson")
    metrics["spatial_procrustes_axis2_pearson"] = safe_corr(ref_scaled[:, 1], rep_aligned[:, 1], "pearson")
    k = min(10, compare.shape[0] - 1)
    if k >= 1:
        ref_nn = NearestNeighbors(n_neighbors=k + 1).fit(ref).kneighbors(return_distance=False)[:, 1:]
        rep_nn = NearestNeighbors(n_neighbors=k + 1).fit(rep).kneighbors(return_distance=False)[:, 1:]
        overlap = [len(set(a).intersection(b)) / k for a, b in zip(ref_nn, rep_nn)]
        metrics["spatial_knn10_mean_overlap"] = float(np.mean(overlap))
    compare.to_csv(validation_dir / "spatial_coordinates_reproduced_vs_metadata.csv", index=False)

summary_json.write_text(json.dumps(metrics, indent=2, sort_keys=True, allow_nan=True))
summary_md.write_text(
    "# FY1_3 validation summary\n\n"
    + "\n".join(f"- **{key}**: {value}" for key, value in sorted(metrics.items()))
    + "\n"
)

print(json.dumps(metrics, indent=2, sort_keys=True, allow_nan=True))
