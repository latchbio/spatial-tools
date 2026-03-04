from tifffile import TiffFile, imread as tiff_imread
import json
import numpy as np
from typing import cast
import duckdb
import scanpy
import pandas as pd
import tarfile
import re

from .args import XeniumArguments
from .ctx import Ctx
from .rect import Fov
from .vec2 import Vec2
from .write_tile import Tile, _ome_tiff_to_2d

# https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs


def populate_ctx(ctx: Ctx, args: XeniumArguments) -> None:
    meta = json.loads((args.input / "experiment.xenium").read_text())

    ctx.mm_per_px = meta["pixel_size"] / 1000
    ctx.log(f"{ctx.mm_per_px} mm/px ({1 / ctx.mm_per_px} px/mm)")
    ctx.log()

    if not args.no_rescale:
        ctx.log("Loading images, calculating quantiles")
    else:
        ctx.log("Listing images")

    # Newer Xenium runs can omit `morphology_mip_filepath` and only provide
    # `morphology_filepath` + `morphology_focus_filepath`. Be flexible and fall
    # back to `morphology_filepath` as the MIP-equivalent when needed.
    images = meta.get("images", {})
    image_specs: list[tuple[str, str]] = []

    if "morphology_mip_filepath" in images:
        image_specs.append(("mip", images["morphology_mip_filepath"]))
    elif "morphology_filepath" in images:
        image_specs.append(("mip", images["morphology_filepath"]))

    if "morphology_focus_filepath" in images:
        image_specs.append(("focus", images["morphology_focus_filepath"]))

    if not image_specs:
        raise KeyError(
            "No usable morphology image paths found in experiment.xenium['images']"
        )

    fov: Fov | None = None
    for cat, subpath in image_specs:
        img_p = args.input / subpath
        ctx.log(f"{cat}: {subpath}")
        # Use same 2D reduction as write_tile so quantiles and tiling match.
        data = np.asarray(tiff_imread(img_p))
        ctx.log(f"  Raw shape: {data.shape} dtype: {data.dtype}")
        data_2d = _ome_tiff_to_2d(data)
        ctx.log(
            f"  2D shape: {data_2d.shape} min: {data_2d.min()} max: {data_2d.max()}"
        )
        height_2d, width_2d = data_2d.shape[0], data_2d.shape[1]

        if not args.no_rescale:
            cur_np = np.quantile(data_2d, [0.05, 0.95])
            ctx.quantiles[cat] = Vec2(int(cur_np[0]), int(cur_np[1]))
            ctx.log(f"  Data range: {ctx.quantiles[cat].x}-{ctx.quantiles[cat].y}")

        if fov is not None:
            assert Vec2(width_2d, height_2d) == fov.size_px
            fov.paths[cat] = img_p
        else:
            fov = Fov(
                ctx=ctx,
                paths={cat: img_p},
                id="all",
                pos_mm=Vec2(0, 0),
                size_px=Vec2(width_2d, height_2d),
            )

    ctx.add_fov(fov)


def generate_extras(ctx: Ctx, args: XeniumArguments) -> None:
    slide = ctx.fovs[0]
    tile_px = Tile(ctx=ctx, z=0, pos_idx=Vec2(0, 0)).resolution()

    mm_per_tile_px = (slide.size_mm.x - slide.pos_mm.x) / tile_px.x

    ctx.log("Writing cell boundaries")
    if not args.dryrun:
        (args.output / "cell_boundaries.duckdb").unlink(missing_ok=True)

    con = duckdb.connect(
        args.output / "cell_boundaries.duckdb" if not args.dryrun else None
    )
    try:
        con.sql(
            f"""
            install spatial;
            load spatial;
            drop table if exists cell_boundaries;
            create table
                cell_boundaries
            as
            select
                cell_id
                    as EntityID,
                ST_MakeLine(
                    list(
                        ST_Point(
                            (vertex_x / 1000 - $x_mm) / $mm_per_tile_px,
                            -(vertex_y / 1000 - $y_mm) / $mm_per_tile_px
                        )
                    )
                )
                    as Geometry
            from
                read_parquet($path)
                data
            group by
                cell_id
            """,
            params={
                "path": str(args.input / "cell_boundaries.parquet"),
                "x_mm": slide.pos_mm.x,
                "y_mm": slide.pos_mm.y,
                "mm_per_tile_px": mm_per_tile_px,
            },
        )
    finally:
        con.close()

    ctx.log("Writing transcripts")
    if not args.dryrun:
        (args.output / "transcripts.duckdb").unlink(missing_ok=True)

    con = duckdb.connect(
        args.output / "transcripts.duckdb" if not args.dryrun else None
    )
    try:
        con.sql(
            f"""
            drop table if exists final_transcripts;
            create table
                final_transcripts
            as
            select
                cell_id::text
                    as cell_id,
                feature_name::text
                    as target,
                (x_location / 1000 - $x_mm) / $mm_per_tile_px
                    as global_x,
                -(y_location / 1000 - $y_mm) / $mm_per_tile_px
                    as global_y,
                -1
                    as cell_comp
            from
                read_parquet($path)
                data
            """,
            params={
                "path": str(args.input / "transcripts.parquet"),
                "x_mm": slide.pos_mm.x,
                "y_mm": slide.pos_mm.y,
                "mm_per_tile_px": mm_per_tile_px,
            },
        )
    finally:
        con.close()

    ctx.log("Loading the cell feature matrix")
    adata = scanpy.read_10x_h5(args.input / "cell_feature_matrix.h5")

    ctx.log("  Adding spatial coordinates and cell information")
    data = duckdb.sql(
        """
        select
            cell_id::text
                as cell_id,
            (x_centroid / 1000 - $x_mm) / $mm_per_tile_px
                as x,
            -(y_centroid / 1000 - $y_mm) / $mm_per_tile_px
                as y,
            transcript_counts,
            cell_area,
            nucleus_area
        from
            read_parquet($path)
            data
        """,
        params={
            "path": str(args.input / "cells.parquet"),
            "x_mm": slide.pos_mm.x,
            "y_mm": slide.pos_mm.y,
            "mm_per_tile_px": mm_per_tile_px,
        },
    ).fetchnumpy()
    # Check order
    assert list(adata.obs_names) == list(data["cell_id"])
    adata.obsm["Spatial"] = np.column_stack((data["x"], data["y"]))
    adata.obs["Transcript Counts"] = data["transcript_counts"]
    adata.obs["Cell Area"] = data["cell_area"]
    adata.obs["Nucelus Area"] = data["nucleus_area"]

    analysis_tar = args.input / "analysis.tar.gz"
    analysis_dir = args.input / "analysis"

    if analysis_tar.exists():
        ctx.log("  Adding embeddings from analysis.tar.gz")

        with tarfile.open(analysis_tar) as analysis:
            ctx.log("    - PCA")
            data = pd.read_csv(
                analysis.extractfile(
                    "analysis/pca/gene_expression_10_components/projection.csv"
                ),
                index_col="Barcode",
            )
            data = data.reindex(adata.obs_names)
            adata.obsm["PCA"] = np.column_stack((data["PC-1"], data["PC-2"]))

            ctx.log("    - TSNE")
            data = pd.read_csv(
                analysis.extractfile(
                    "analysis/tsne/gene_expression_2_components/projection.csv"
                ),
                index_col="Barcode",
            )
            data = data.reindex(adata.obs_names)
            adata.obsm["TSNE"] = np.column_stack((data["TSNE-1"], data["TSNE-2"]))

            ctx.log("    - UMAP")
            data = pd.read_csv(
                analysis.extractfile(
                    "analysis/umap/gene_expression_2_components/projection.csv"
                ),
                index_col="Barcode",
            )
            data = data.reindex(adata.obs_names)
            adata.obsm["UMAP"] = np.column_stack((data["UMAP-1"], data["UMAP-2"]))

            clustering_re = re.compile(
                r"^analysis/clustering/gene_expression_([^/]+)/clusters.csv$"
            )
            kmeans_name_re = re.compile(r"^kmeans_(\d+)_clusters$")

            ctx.log("  Adding clusterings")
            for x in analysis.getnames():
                m = clustering_re.match(x)
                if m is None:
                    continue

                name = m.group(1)
                if name == "graphclust":
                    name = "Graph-based"

                m = kmeans_name_re.match(name)
                if m is not None:
                    name = f"K={m.group(1)} K-means Gene Expressions K-means"

                print(f"    - {name}")
                data = pd.read_csv(analysis.extractfile(x), index_col="Barcode")
                data = data.reindex(adata.obs_names)
                adata.obs[name] = pd.Categorical(data["Cluster"])

    elif analysis_dir.is_dir():
        ctx.log("  Adding embeddings from analysis/ directory")

        # PCA
        pca_csv = analysis_dir / "pca/gene_expression_10_components/projection.csv"
        if pca_csv.exists():
            ctx.log("    - PCA")
            data = pd.read_csv(pca_csv, index_col="Barcode")
            data = data.reindex(adata.obs_names)
            adata.obsm["PCA"] = np.column_stack((data["PC-1"], data["PC-2"]))

        # TSNE
        tsne_csv = analysis_dir / "tsne/gene_expression_2_components/projection.csv"
        if tsne_csv.exists():
            ctx.log("    - TSNE")
            data = pd.read_csv(tsne_csv, index_col="Barcode")
            data = data.reindex(adata.obs_names)
            adata.obsm["TSNE"] = np.column_stack((data["TSNE-1"], data["TSNE-2"]))

        # UMAP
        umap_csv = analysis_dir / "umap/gene_expression_2_components/projection.csv"
        if umap_csv.exists():
            ctx.log("    - UMAP")
            data = pd.read_csv(umap_csv, index_col="Barcode")
            data = data.reindex(adata.obs_names)
            adata.obsm["UMAP"] = np.column_stack((data["UMAP-1"], data["UMAP-2"]))

        # Clustering
        clustering_re = re.compile(
            r"^clustering/gene_expression_([^/]+)/clusters\.csv$"
        )
        kmeans_name_re = re.compile(r"^kmeans_(\d+)_clusters$")

        ctx.log("  Adding clusterings")
        for f in analysis_dir.rglob("clustering/gene_expression_*/clusters.csv"):
            rel = f.relative_to(analysis_dir).as_posix()
            m = clustering_re.match(rel)
            if m is None:
                continue

            name = m.group(1)
            if name == "graphclust":
                name = "Graph-based"

            m = kmeans_name_re.match(name)
            if m is not None:
                name = f"K={m.group(1)} K-means Gene Expressions K-means"

            print(f"    - {name}")
            data = pd.read_csv(f, index_col="Barcode")
            data = data.reindex(adata.obs_names)
            adata.obs[name] = pd.Categorical(data["Cluster"])

    else:
        ctx.log(
            "  analysis.tar.gz or analysis/ not found; skipping embeddings and clustering"
        )

    if not args.dryrun:
        ctx.log("  Saving")
        adata.write(args.output / "cell_by_gene.h5ad")

    ctx.log()
