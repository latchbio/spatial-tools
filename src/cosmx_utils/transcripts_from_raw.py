import time
from pathlib import Path
from typing import TYPE_CHECKING

import anndata as ad
import duckdb
import numpy as np
import pandas as pd
import tqdm
from scipy.sparse import coo_matrix

if TYPE_CHECKING:
    from spatial_tools.ctx import Ctx


def get_transcripts(slide_dir: Path, ctx: "Ctx") -> duckdb.DuckDBPyConnection:
    start_time = time.time()
    print(f"Have {len(ctx.fovs)} FOVs")

    transcript_files = list(
        slide_dir.glob("**/*__complete_code_cell_target_call_coord.csv")
    )
    print(f"Found {len(transcript_files)} transcript CSV files")
    if len(transcript_files) == 0:
        raise ValueError("No transcript CSV files found")

    print("Copying into DuckDB...")
    con = duckdb.connect(database=":memory:")
    _ = con.execute("""
        create table transcripts (
            fov integer,
            cell_id integer,
            x double,
            y double,
            target text,
            cell_comp text
        )
    """)

    for csv_file in tqdm.tqdm(transcript_files, desc="Copying into DuckDB"):
        _ = con.execute(f"""
            insert into transcripts
            select
                cast(fov as integer)       as fov,
                cast(CellId as integer)    as cell_id,
                cast(x as double)          as x,
                cast(y as double)          as y,
                target                     as target,
                CellComp                   as cell_comp
            from read_csv('{csv_file}', header=true, sample_size=-1)
        """)

    print(f"Copying into DuckDB took {time.time() - start_time:.2f} seconds")

    start_time = time.time()
    print("Processing FOV data...")
    fov_data = {
        "fov": [fov.id for fov in ctx.fovs],
        "x_mm": [fov.pos_mm.x for fov in ctx.fovs],
        "y_mm": [fov.pos_mm.y for fov in ctx.fovs],
    }

    _ = con.execute("""
        create table fov_data (
            fov integer,
            x_mm double,
            y_mm double
        )
    """)

    for fov, x_mm, y_mm in zip(
        fov_data["fov"], fov_data["x_mm"], fov_data["y_mm"], strict=True
    ):
        _ = con.execute("insert into fov_data values (?, ?, ?)", [fov, x_mm, y_mm])

    # Global coordinates in same space as stitched PMTiles: origin = ctx.slide.pos_mm (top-left), scale0
    px_per_mm = 1 / ctx.mm_per_px
    scale0 = min(
        ctx.viewport_size_px.x / ctx.slide.size_px.x,
        ctx.viewport_size_px.y / ctx.slide.size_px.y,
    )
    slide_origin_x = ctx.slide.pos_mm.x
    slide_origin_y = ctx.slide.pos_mm.y
    # X is fine as-is; Y is negated to flip orientation
    _ = con.execute(
        """
        create table final_transcripts as
        select
            t.fov,
            t.cell_id,
            t.target,
            t.cell_comp,
            ((f.x_mm - ?) * ? + t.x) * ? as global_x,
            -((f.y_mm - ?) * ? + t.y) * ? as global_y
        from
            transcripts t
        join
            fov_data f
        on
            t.fov = f.fov
    """,
        [slide_origin_x, px_per_mm, scale0, slide_origin_y, px_per_mm, scale0],
    )
    print(f"Processing FOV data took {time.time() - start_time:.2f} seconds")

    return con


def create_h5ad(con: duckdb.DuckDBPyConnection, output_path: Path) -> None:
    start_time = time.time()

    print("Creating counts matrix...")
    counts_df = con.execute("""
        select
            concat(fov, '_', cell_id) as unique_cell_id,
            target as gene,
            count(*) as count
        from final_transcripts
        group by fov, cell_id, target
    """).fetchdf()
    n_cells = counts_df["unique_cell_id"].nunique()
    n_genes = counts_df["gene"].nunique()
    print(f"Found {n_cells} cells and {n_genes} genes")
    print("Building sparse counts matrix...")

    # Build a sparse cell x gene matrix directly instead of a huge dense pivot
    cells_c = counts_df["unique_cell_id"].astype("category")
    genes_c = counts_df["gene"].astype("category")

    row_idx = cells_c.cat.codes.to_numpy()
    col_idx = genes_c.cat.codes.to_numpy()
    data = counts_df["count"].to_numpy(dtype=np.int32)

    counts_matrix = coo_matrix(
        (data, (row_idx, col_idx)),
        shape=(cells_c.cat.categories.size, genes_c.cat.categories.size),
        dtype=np.int32,
    ).tocsr()

    cells = cells_c.cat.categories.astype(str).tolist()
    genes = genes_c.cat.categories.astype(str).tolist()
    print(f"Creating counts matrix took {time.time() - start_time:.2f} seconds")

    start_time = time.time()
    print("Getting cell metadata...")
    obs_df = con.execute("""
        select
            concat(fov, '_', cell_id) as unique_cell_id,
            fov,
            cell_id,
            count(*) as total_transcripts
        from final_transcripts
        group by fov, cell_id
        order by unique_cell_id
    """).fetchdf()
    obs_df = obs_df.set_index("unique_cell_id").loc[cells].reset_index()

    spatial_df = con.execute("""
        select
            concat(fov, '_', cell_id) as unique_cell_id,
            avg(global_x) as x,
            avg(global_y) as y
        from final_transcripts
        group by fov, cell_id
        order by unique_cell_id
    """).fetchdf()
    spatial_df = spatial_df.set_index("unique_cell_id").loc[cells].reset_index()
    spatial_coords = spatial_df[["x", "y"]].to_numpy()
    print(f"Getting metadata took {time.time() - start_time:.2f} seconds")

    start_time = time.time()
    print("Creating AnnData object...")
    obs = pd.DataFrame(
        {
            "fov": obs_df["fov"].to_numpy(),
            "cell_id": obs_df["cell_id"].to_numpy(),
            "total_transcripts": obs_df["total_transcripts"].to_numpy(),
        },
        index=obs_df["unique_cell_id"].astype(str),
    )
    obs.index.name = "unique_cell_id"
    var = pd.DataFrame(index=genes)
    var.index.name = "gene"
    adata = ad.AnnData(X=counts_matrix, obs=obs, var=var)
    adata.obsm["spatial"] = spatial_coords
    adata.strings_to_categoricals()
    adata.uns["name"] = "CosMx transcripts"

    print(f"Saving to {output_path}...")
    adata.write_h5ad(output_path, compression="gzip")
    print(f"Created h5ad file with {adata.n_obs} cells and {adata.n_vars} genes")
    print(f"Total processing took {time.time() - start_time:.2f} seconds")
