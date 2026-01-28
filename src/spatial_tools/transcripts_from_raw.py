import time
from pathlib import Path

import anndata as ad
import duckdb
import numpy as np
import pandas as pd
import tqdm

from .ctx import ctx


def get_transcripts(slide_dir: Path) -> duckdb.DuckDBPyConnection:
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

    # global pixel coordinates
    px_per_mm = 1 / ctx.mm_per_px
    fov_height = ctx.fovs[0].size_px.y
    fov_width = ctx.fovs[0].size_px.x

    dash = (fov_height / fov_width) != 1

    # compute min/max values
    tl_mm_result = con.execute(
        """
        select
            case when ? then -max(x_mm) else min(y_mm) end as x,
            case when ? then min(y_mm) else -max(x_mm) end as y
        from fov_data
    """,
        [dash, dash],
    ).fetchone()

    if tl_mm_result is None:
        raise ValueError("Failed to compute min/max values from FOV data")

    tl_mm = (tl_mm_result[0], tl_mm_result[1])

    # scale into the same viewport as the stitched FOVs in pmtiles
    scale0 = min(
        ctx.viewport_size_px.x / ctx.slide.size_px.x,
        ctx.viewport_size_px.y / ctx.slide.size_px.y,
    )

    _ = con.execute(
        """
        create table final_transcripts as
        select
            t.fov,
            t.cell_id,
            t.target,
            t.cell_comp,
            (case
                when ? then t.x + f.y_mm * ? - ? * ?
                else t.x - f.x_mm * ? - ? * ?
            end) * ? as global_x,
            (case
                when ? then -(t.y - f.x_mm * ? - ? * ?)
                else -(t.y + f.y_mm * ? - ? * ?)
            end) * ? as global_y
        from
            transcripts t
        join
            fov_data f
        on
            t.fov = f.fov
    """,
        [
            dash,
            px_per_mm,
            tl_mm[1],
            px_per_mm,
            px_per_mm,
            tl_mm[1],
            px_per_mm,
            scale0,
            dash,
            px_per_mm,
            tl_mm[0],
            px_per_mm,
            px_per_mm,
            tl_mm[0],
            px_per_mm,
            scale0,
        ],
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
    print(
        f"Found {counts_df['unique_cell_id'].nunique()} cells and {counts_df['gene'].nunique()} genes"
    )
    print("Pivoting to create counts matrix...")
    counts_pivot = counts_df.pivot_table(
        index="unique_cell_id",
        columns="gene",
        values="count",
        fill_value=0,
        aggfunc="sum",
    )
    cells = counts_pivot.index.tolist()
    genes = counts_pivot.columns.tolist()
    counts_matrix = counts_pivot.to_numpy(dtype=np.int32)
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
