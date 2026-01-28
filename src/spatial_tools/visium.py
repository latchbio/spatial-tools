"""Convert 10x Visium HD SpaceRanger output to AnnData (h5ad) format.

This module expects a SpaceRanger output directory with the following structure:

    spaceranger_output/
    └── binned_outputs/
        ├── square_002um/
        │   ├── filtered_feature_bc_matrix/
        │   │   ├── barcodes.tsv.gz
        │   │   ├── features.tsv.gz
        │   │   └── matrix.mtx.gz
        │   ├── raw_feature_bc_matrix/
        │   │   ├── barcodes.tsv.gz
        │   │   ├── features.tsv.gz
        │   │   └── matrix.mtx.gz
        │   ├── spatial/
        │   │   ├── tissue_positions.parquet  (or .csv or tissue_positions_list.csv)
        │   │   ├── tissue_hires_image.png
        │   │   ├── tissue_lowres_image.png
        │   │   └── scalefactors_json.json
        │   └── analysis/  (optional)
        │       ├── pca/
        │       │   └── gene_expression_10_components/
        │       │       ├── projection.csv
        │       │       ├── components.csv
        │       │       └── variance.csv
        │       ├── umap/
        │       │   └── gene_expression_2_components/
        │       │       └── projection.csv
        │       ├── clustering/
        │       │   └── gene_expression_graphclust/
        │       │       └── clusters.csv
        │       └── diffexp/
        │           └── gene_expression_graphclust/
        │               └── differential_expression.csv
        ├── square_008um/  (same structure)
        └── square_016um/  (same structure)

Example usage:
    python -m spatial_tools visium /path/to/spaceranger_output
    python -m spatial_tools visium /path/to/spaceranger_output --bin-name square_002um
    python -m spatial_tools visium /path/to/spaceranger_output --matrix-type raw -o output.h5ad

Reference:
    https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview

Test data on Latch:
    latch://38438.account/visium/colon
"""

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from PIL import Image

from .args import VisiumArguments
from .ctx import Ctx

def read_tissue_positions(spatial_dir: Path) -> pd.DataFrame:
    """Load tissue positions, preferring Parquet, then CSV, then old CSV list."""
    pq = spatial_dir / "tissue_positions.parquet"
    if pq.exists():
        tp = pd.read_parquet(pq)
        return tp.set_index("barcode")

    csv = spatial_dir / "tissue_positions.csv"
    if csv.exists():
        tp = pd.read_csv(csv)
        return tp.set_index("barcode")

    csv_list = spatial_dir / "tissue_positions_list.csv"
    if csv_list.exists():
        tp = pd.read_csv(
            csv_list,
            header=None,
            names=[
                "barcode",
                "in_tissue",
                "array_row",
                "array_col",
                "pxl_row_in_fullres",
                "pxl_col_in_fullres",
            ],
        )
        return tp.set_index("barcode")

    raise FileNotFoundError(f"No tissue_positions file found in {spatial_dir}")


def load_images_and_scalefactors(spatial_dir: Path):
    """Load hires/lowres images and scalefactors into dicts."""
    images = {}
    hires_path = spatial_dir / "tissue_hires_image.png"
    lowres_path = spatial_dir / "tissue_lowres_image.png"

    if hires_path.exists():
        images["hires"] = np.array(Image.open(hires_path))
    if lowres_path.exists():
        images["lowres"] = np.array(Image.open(lowres_path))

    scalefactors = {}
    sf_path = spatial_dir / "scalefactors_json.json"
    if sf_path.exists():
        with open(sf_path) as f:
            scalefactors = json.load(f)

    return images, scalefactors


def attach_pca(adata: ad.AnnData, pca_dir: Path, prefix: str = "X_pca"):
    """
    Attach SpaceRanger PCA projections and metadata.
    Any subdir under pca_dir that has projection.csv is used.
    """
    if not pca_dir.exists():
        return

    adata.uns.setdefault("spaceranger_pca", {})

    for subdir in sorted(pca_dir.iterdir()):
        if not subdir.is_dir():
            continue

        proj_path = subdir / "projection.csv"
        if not proj_path.exists():
            continue

        proj = pd.read_csv(proj_path, index_col=0)
        proj = proj.reindex(adata.obs_names)
        key = f"{prefix}_{subdir.name}"
        adata.obsm[key] = proj.to_numpy()

        meta = {}
        comps_path = subdir / "components.csv"
        var_path = subdir / "variance.csv"
        feats_path = subdir / "features_selected.csv"
        disp_path = subdir / "dispersion.csv"

        if comps_path.exists():
            meta["components"] = pd.read_csv(comps_path, index_col=0)
        if var_path.exists():
            meta["variance"] = pd.read_csv(var_path, index_col=0)
        if feats_path.exists():
            meta["features_selected"] = pd.read_csv(feats_path, index_col=0)
        if disp_path.exists():
            meta["dispersion"] = pd.read_csv(disp_path, index_col=0)

        if meta:
            adata.uns["spaceranger_pca"][subdir.name] = meta


def attach_umap(adata: ad.AnnData, umap_dir: Path, prefix: str = "X_umap"):
    """
    Attach SpaceRanger UMAP projections.
    Any subdir under umap_dir that has projection.csv is used.
    """
    if not umap_dir.exists():
        return

    for subdir in sorted(umap_dir.iterdir()):
        if not subdir.is_dir():
            continue

        proj_path = subdir / "projection.csv"
        if not proj_path.exists():
            continue

        proj = pd.read_csv(proj_path, index_col=0)
        proj = proj.reindex(adata.obs_names)
        key = f"{prefix}_{subdir.name}"
        adata.obsm[key] = proj.to_numpy()


def attach_clustering(adata: ad.AnnData, clustering_dir: Path):
    """
    Attach SpaceRanger clustering assignments as columns in adata.obs.
    Any subdir under clustering_dir that has clusters.csv is used.
    """
    if not clustering_dir.exists():
        return

    for method_dir in sorted(clustering_dir.iterdir()):
        if not method_dir.is_dir():
            continue

        clusters_path = method_dir / "clusters.csv"
        if not clusters_path.exists():
            continue

        df = pd.read_csv(clusters_path, index_col=0)
        if "cluster" in df.columns:
            cluster_series = df["cluster"]
        else:
            cluster_series = df.iloc[:, 0]

        cluster_series = cluster_series.reindex(adata.obs_names)
        col_name = f"spaceranger_{method_dir.name}"
        adata.obs[col_name] = cluster_series.astype("category")


def attach_diffexp(adata: ad.AnnData, diffexp_dir: Path):
    """
    Attach diffexp tables into adata.uns['spaceranger_diffexp'].
    Any subdir under diffexp_dir that has differential_expression.csv is used.
    """
    if not diffexp_dir.exists():
        return

    diffexp_store = {}
    for method_dir in sorted(diffexp_dir.iterdir()):
        if not method_dir.is_dir():
            continue

        de_path = method_dir / "differential_expression.csv"
        if not de_path.exists():
            continue

        df = pd.read_csv(de_path)
        diffexp_store[method_dir.name] = df

    if diffexp_store:
        adata.uns.setdefault("spaceranger_diffexp", {})
        adata.uns["spaceranger_diffexp"] = diffexp_store


def build_visium_hd_adata_for_bin(
    bin_dir: Path,
    matrix_type: str = "filtered",
    load_analysis: bool = True,
) -> ad.AnnData:
    """
    Build an AnnData for a given bin directory, for example:
    visium_hd_data/binned_outputs/square_008um

    matrix_type: 'filtered' or 'raw' → choose filtered_feature_bc_matrix/ or raw_feature_bc_matrix/
    load_analysis: if False, skip PCA/UMAP/clustering/diffexp from analysis/
    """
    bin_dir = bin_dir.resolve()

    if matrix_type == "filtered":
        mtx_dir = bin_dir / "filtered_feature_bc_matrix"
    elif matrix_type == "raw":
        mtx_dir = bin_dir / "raw_feature_bc_matrix"
    else:
        raise ValueError("matrix_type must be 'filtered' or 'raw'")

    spatial_dir = bin_dir / "spatial"
    analysis_dir = bin_dir / "analysis"

    if not mtx_dir.exists():
        raise FileNotFoundError(f"Missing {mtx_dir.name} in {bin_dir}")
    if not spatial_dir.exists():
        raise FileNotFoundError(f"Missing spatial directory in {bin_dir}")

    # counts
    adata = sc.read_10x_mtx(
        mtx_dir,
        var_names="gene_symbols",
        make_unique=True,
    )

    # spatial metadata
    tp = read_tissue_positions(spatial_dir)
    adata.obs = adata.obs.join(tp, how="left")

    if {"pxl_row_in_fullres", "pxl_col_in_fullres"}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[
            ["pxl_row_in_fullres", "pxl_col_in_fullres"]
        ].to_numpy(dtype=float)
    else:
        raise ValueError(
            f"Pixel coordinate columns missing in tissue positions for {bin_dir}"
        )

    images, scalefactors = load_images_and_scalefactors(spatial_dir)

    library_id = bin_dir.name
    adata.uns.setdefault("spatial", {})
    adata.uns["spatial"][library_id] = {
        "images": images,
        "scalefactors": scalefactors,
    }

    # analysis outputs
    if load_analysis and analysis_dir.exists():
        attach_pca(adata, analysis_dir / "pca")
        attach_umap(adata, analysis_dir / "umap")
        attach_clustering(adata, analysis_dir / "clustering")
        attach_diffexp(adata, analysis_dir / "diffexp")

    adata.var_names_make_unique()
    return adata


def populate_ctx(ctx: Ctx, args: VisiumArguments) -> None:
    """Convert 10x Visium HD SpaceRanger bin (square_XXXum) to h5ad.
    
    See module docstring for expected SpaceRanger directory structure.
    """
    spaceranger_dir = args.input.resolve()
    bin_dir = spaceranger_dir / "binned_outputs" / args.bin_name

    ctx.log(f"Bin name: {args.bin_name}")
    ctx.log(f"Matrix type: {args.matrix_type}")
    if args.no_analysis:
        ctx.log("Will skip analysis outputs due to --no-analysis")
    ctx.log()

    if not bin_dir.exists():
        raise FileNotFoundError(f"Bin directory not found: {bin_dir}")

    ctx.log(f"Loading data from bin: {bin_dir}")
    adata = build_visium_hd_adata_for_bin(
        bin_dir=bin_dir,
        matrix_type=args.matrix_type,
        load_analysis=not args.no_analysis,
    )

    ctx.log(f"Loaded AnnData: {adata.n_obs:,} observations, {adata.n_vars:,} variables")

    output_resolved = args.output.resolve()
    auto_name = f"{spaceranger_dir.name}_{args.bin_name}_{args.matrix_type}.h5ad"

    if output_resolved.suffix == ".h5ad":
        out_path = output_resolved
    else:
        output_resolved.mkdir(parents=True, exist_ok=True)
        out_path = output_resolved / auto_name

    ctx.log(f"Saving h5ad to: {out_path}")
    adata.write_h5ad(out_path)
    ctx.log(f"Saved h5ad: {out_path}")
