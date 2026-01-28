import json
import time
from pathlib import Path

import duckdb
import numpy as np
import squidpy as sq  # type: ignore
from scipy.sparse import issparse

from spatial_tools.ctx import ctx


def get_transcripts(dataset_dir: Path) -> duckdb.DuckDBPyConnection:
    images_dir = dataset_dir / "images"
    manifest_path = images_dir / "manifest.json"
    transcripts_csv = dataset_dir / "detected_transcripts.csv"

    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found at {manifest_path}")
    if not transcripts_csv.exists():
        raise FileNotFoundError(f"Transcript CSV not found at {transcripts_csv}")

    with manifest_path.open("r", encoding="utf-8") as f:
        manifest = json.load(f)

    microns_per_pixel = float(manifest["microns_per_pixel"])
    bbox_microns = manifest.get("bbox_microns", [0.0, 0.0, 0.0, 0.0])
    min_x_microns, min_y_microns = float(bbox_microns[0]), float(bbox_microns[1])

    print("Copying detected_transcripts.csv into DuckDB")
    start_time = time.time()

    con = duckdb.connect(database=":memory:")
    con.execute(
        """
        create or replace table transcripts as
        -- Load the raw CSV exactly as-is first; we will perform casting in a later
        -- query so that we **only** use an existing "cell_id" column (if present)
        -- rather than incorrectly deriving one from the barcode.
        select *
        from read_csv_auto(? , sample_size=-1)
        """,
        [str(transcripts_csv)],
    )

    row = con.execute("select count(*) as cnt from transcripts").fetchone()
    total_transcripts = row[0] if row is not None else 0

    print(
        f"Loaded {total_transcripts:,} transcripts in {time.time() - start_time:.1f}s"
    )

    px_per_micron = 1.0 / microns_per_pixel

    slide_width_px = int(manifest["mosaic_width_pixels"])
    slide_height_px = int(manifest["mosaic_height_pixels"])

    scale0 = min(
        ctx.viewport_size_px.x / slide_width_px,
        ctx.viewport_size_px.y / slide_height_px,
    )

    # Determine whether the CSV already contains a 'cell_id' column produced by
    # VPT's `partition_transcripts` step. We prefer that because it encodes the
    # *EntityID* assigned during segmentation and is globally unique.

    columns_info = con.execute("PRAGMA table_info('transcripts')").fetchall()
    column_names = {row[1] for row in columns_info}

    if "cell_id" not in column_names:
        raise RuntimeError(
            "detected_transcripts.csv does not contain a 'cell_id' column. "
            "Run the Vizgen PartitionTranscripts step (vpt.partition_transcripts) "
            "before building the h5ad, or fall back to the pre-computed cell_by_gene.csv pipeline."
        )

    con.execute(
        """
            create or replace table final_transcripts as
            select
                cast(coalesce(fov, 0) as integer)                   as fov,
                cast(cell_id as bigint)                             as cell_id,
                cast(gene as text)                                  as target,
                -- X/Y in *slide* pixel space (origin = top-left)
                ( (cast(global_x as double) - ?) * ? ) * ?          as global_x,
                -1 * ( (cast(global_y as double) - ?) * ? ) * ?     as global_y
            from transcripts
        """,
        [min_x_microns, px_per_micron, scale0, min_y_microns, px_per_micron, scale0],
    )

    return con


def create_h5ad_from_files(dataset_dir: Path, output_path: Path) -> None:
    counts_file = "cell_by_gene.csv"
    meta_file = "cell_metadata.csv"
    transformation_file = "micron_to_mosaic_pixel_transform.csv"

    adata = sq.read.vizgen(
        path=dataset_dir,
        counts_file=counts_file,
        meta_file=meta_file,
        transformation_file=transformation_file,
    )

    for key, val in list(adata.obsm.items()):
        if not isinstance(val, np.ndarray):
            adata.obsm[key] = np.asarray(val, dtype=np.float32)

    if issparse(adata.X):
        adata.obs["total_transcripts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    else:
        adata.obs["total_transcripts"] = adata.X.sum(axis=1)

    if "spatial" in adata.obsm:
        spatial_arr = adata.obsm.pop("spatial")
        adata.obsm = {"spatial": spatial_arr, **adata.obsm}

    if "spatial" in adata.obsm:
        adata.obsm["spatial"][:, 1] *= -1

    if "spatial" in adata.obsm:
        images_dir = dataset_dir / "images"
        manifest_path = images_dir / "manifest.json"
        with manifest_path.open("r", encoding="utf-8") as f:
            manifest = json.load(f)

        microns_per_pixel = float(manifest["microns_per_pixel"])
        bbox_microns = manifest.get("bbox_microns", [0.0, 0.0, 0.0, 0.0])
        min_x_microns, min_y_microns = float(bbox_microns[0]), float(bbox_microns[1])

        px_per_micron = 1.0 / microns_per_pixel
        slide_width_px = int(manifest["mosaic_width_pixels"])
        slide_height_px = int(manifest["mosaic_height_pixels"])

        scale0 = min(
            ctx.viewport_size_px.x / slide_width_px,
            ctx.viewport_size_px.y / slide_height_px,
        )

        spatial_coords = adata.obsm["spatial"].copy()
        spatial_coords[:, 0] = (
            (spatial_coords[:, 0] - min_x_microns) * px_per_micron * scale0
        )
        spatial_coords[:, 1] = (
            (spatial_coords[:, 1] + min_y_microns) * px_per_micron * scale0
        )

        adata.obsm["spatial"] = spatial_coords

    print(f"Writing {output_path}")
    adata.write_h5ad(output_path, compression="gzip")
    print(
        f"Saved h5ad with {adata.n_obs:,} cells x {adata.n_vars:,} genes -- "
        f"{output_path.stat().st_size / 1_048_576:.1f} MiB"
    )
