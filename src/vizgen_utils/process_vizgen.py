import argparse
import json
import time
from math import ceil, log2
from pathlib import Path

import duckdb
import geopandas as gpd  # type: ignore
import numpy as np
from PIL import Image
from shapely.affinity import scale as shp_scale
from shapely.affinity import translate as shp_translate

from spatial_tools.ctx import ctx
from spatial_tools.pmtiles import write_pmtiles
from spatial_tools.vec2 import Vec2
from spatial_tools.write_tile import iterate_tiles

from .transcripts_from_csv import create_h5ad_from_files as create_h5ad_vizgen
from .transcripts_from_csv import get_transcripts

Image.MAX_IMAGE_PIXELS = None  # TIFFs can be huge, and PIL detects a decompression bomb on many ~400mb examples otherwise

# CLAHE (Contrast-Limited Adaptive Histogram Equalization) greatly improves
# contrast in regions that are otherwise flat. It matches Vizgen's own
# post-processing pipeline.

import cv2  # type: ignore


def process_vizgen_stain(
    *,
    stain: str,
    mosaic_path: Path,
    slide_size_px: Vec2[int],
    output_dir: Path,
    max_zoom: int,
    z_plane: int,
) -> None:
    print(f"\n>>> Generating PMTiles for stain '{stain}' from {mosaic_path.name}")

    stitched_dir = output_dir / f"stitched_{stain.replace(' ', '_').lower()}_z{z_plane}"
    stitched_dir.mkdir(parents=True, exist_ok=True)

    with Image.open(mosaic_path) as mosaic_img:
        # Determine intensity range of the whole mosaic so every tile is scaled consistently
        extrema = mosaic_img.getextrema()
        # Pillow returns (min,max) for single-band images. Some type stubs mark it
        # as a tuple-of-tuple possibility; handle that defensively.
        if isinstance(extrema[0], tuple):
            extrema = extrema[0]  # type: ignore[assignment]
        glob_min = int(extrema[0])  # type: ignore[arg-type]
        glob_max = int(extrema[1])  # type: ignore[arg-type]
        # Avoid zero division
        if glob_max == glob_min:
            glob_max = glob_min + 1

        for tile in iterate_tiles(max_zoom):
            tile_size_px = tile.size_px

            left = int(tile.pos_idx.x * tile_size_px.x)
            upper = int(tile.pos_idx.y * tile_size_px.y)
            right = int(min(left + tile_size_px.x, slide_size_px.x))
            lower = int(min(upper + tile_size_px.y, slide_size_px.y))

            if right <= left or lower <= upper:
                continue

            region = mosaic_img.crop((left, upper, right, lower))

            # Rescale 16-bit → 8-bit using floating-point math to avoid integer
            # truncation (which previously produced all-black tiles).
            arr16 = np.array(region, dtype=np.uint16)
            rng = glob_max - glob_min
            if rng <= 0:
                arr8 = np.zeros_like(arr16, dtype=np.uint8)
            else:
                arr_scaled = (arr16.astype(np.float32) - glob_min) / rng

                arr8 = (arr_scaled * 255).clip(0, 255).astype(np.uint8)

                # Apply CLAHE (same defaults as Vizgen tools)
                clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
                arr8 = clahe.apply(arr8)

            gray = Image.fromarray(arr8, mode="L")
            region = gray.convert("RGB")

            region = region.resize(
                tile.resolution().to_tuple(), resample=Image.Resampling.LANCZOS
            )

            out_path = stitched_dir / f"{tile.z}/{tile.pos_idx.x}-{tile.pos_idx.y}.webp"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            region.save(out_path, format="webp")

    pmtiles_path = output_dir / f"{stain.replace(' ', '_').lower()}_z{z_plane}.pmtiles"
    write_pmtiles(stitched_dir, pmtiles_path, mirror_x=False)


def process_vizgen_data(
    *,
    dataset_dir: Path,
    output_dir: Path,
    max_zoom: int | None = None,
    generate_pmtiles: bool = False,
    generate_for_stain: str | None = None,
    generate_transcripts: bool = False,
    generate_h5ad: bool = False,
    generate_boundaries: bool = False,
) -> None:
    if generate_pmtiles:
        images_dir = dataset_dir / "images"
        manifest_path = images_dir / "manifest.json"

        if not manifest_path.exists():
            raise FileNotFoundError(f"Manifest not found at {manifest_path}")

        with manifest_path.open("r", encoding="utf-8") as f:
            manifest = json.load(f)

        microns_per_pixel = float(manifest["microns_per_pixel"])
        slide_width_px = int(manifest["mosaic_width_pixels"])
        slide_height_px = int(manifest["mosaic_height_pixels"])

        ctx.mm_per_px = microns_per_pixel / 1_000  # convert to mm
        ctx.slide.pos_mm = Vec2(0.0, 0.0)  # origin
        ctx.slide.size_px = Vec2(slide_width_px, slide_height_px)

        print(
            f"Slide size: {slide_width_px}x{slide_height_px}px • {ctx.mm_per_px * 1e3:.3f}μm/px"
        )

        if max_zoom is None:
            viewports_per_slide = ctx.slide.size_px.pw_div(ctx.viewport_size_px)
            max_zoom = ceil(
                max(0, log2(viewports_per_slide.x), log2(viewports_per_slide.y))
            )

        print(f"Using max zoom level: {max_zoom}")

        pyramid_entries: list[dict[str, object]] = manifest.get("mosaic_files", [])
        if len(pyramid_entries) == 0:
            raise ValueError("No 'mosaic_files' entries found")

        stains_to_entries = {}
        for entry in pyramid_entries:
            stain = str(entry["stain"])
            stains_to_entries.setdefault(stain, []).append(entry)

        # todo(aidan): support multiple z-planes?
        for stain, entries in stains_to_entries.items():
            for entry in entries:
                entry_z = int(entry.get("z", -1))
                if entry_z == -1:
                    print(f"Skipping stain '{stain}': z-plane {entry_z} not specified")
                    continue

                if (
                    generate_for_stain is not None
                    and stain.lower() != generate_for_stain.lower()
                ):
                    print(f"Skipping stain '{stain}': not {generate_for_stain}")
                    continue

                mosaic_path = images_dir / str(entry["file_name"])
                if not mosaic_path.exists():
                    print(
                        f"Skipping stain '{stain}': expected mosaic image {mosaic_path} missing."
                    )
                    continue

                start = time.monotonic()
                print(f"Processing '{stain} at z={entry_z}'")
                process_vizgen_stain(
                    stain=stain,
                    mosaic_path=mosaic_path,
                    slide_size_px=Vec2(slide_width_px, slide_height_px),
                    output_dir=output_dir,
                    max_zoom=max_zoom,
                    z_plane=entry_z,
                )
                print(
                    f"Completed '{stain}' at z={entry_z} in {(time.monotonic() - start):.2f}s"
                )

    if generate_transcripts:
        print("\nProcessing transcripts")
        transcript_start = time.monotonic()

        con = get_transcripts(dataset_dir)

        db_path = output_dir / "transcripts.duckdb"
        print(f"Writing DuckDB database to {db_path}")
        persistent_con = duckdb.connect(str(db_path))

        con.execute(f"attach '{db_path}' as persistent_db")
        con.execute(
            "create or replace table persistent_db.final_transcripts as select * from final_transcripts"
        )
        con.execute("detach persistent_db")
        persistent_con.close()

        print(f"DuckDB created in {time.monotonic() - transcript_start:.1f}s")

    if generate_h5ad:
        print("\nGenerating h5ad using squidpy")
        h5ad_path = output_dir / "cell_by_gene.h5ad"
        create_h5ad_vizgen(dataset_dir, h5ad_path)

    if generate_boundaries:
        print("\nGenerating scaled cell boundaries GeoParquet")
        boundaries_path = dataset_dir / "cell_boundaries.parquet"
        images_dir = dataset_dir / "images"
        manifest_path = images_dir / "manifest.json"

        if not boundaries_path.exists():
            raise FileNotFoundError(
                f"cell_boundaries.parquet not found at {boundaries_path}"
            )
        if not manifest_path.exists():
            raise FileNotFoundError(f"manifest.json not found at {manifest_path}")
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

        gdf = gpd.read_parquet(boundaries_path)
        geom_col = "Geometry" if "Geometry" in gdf.columns else "geometry"

        gdf[geom_col] = gdf[geom_col].apply(
            lambda g: shp_translate(g, xoff=-min_x_microns, yoff=-min_y_microns)
        )
        gdf[geom_col] = gdf[geom_col].apply(
            lambda g: shp_scale(
                g,
                xfact=px_per_micron * scale0,
                yfact=-px_per_micron * scale0,
                origin=(0, 0),
            )
        )

        out_path = output_dir / "cell_boundaries_slide_px.parquet"
        gdf.to_parquet(out_path, engine="pyarrow")
        print(f"Saved boundaries to {out_path}")

        duckdb_path = output_dir / "cell_boundaries.duckdb"
        con_duck = duckdb.connect(str(duckdb_path))

        con_duck.execute("install spatial; load spatial;")
        con_duck.execute(
            "create table cell_boundaries as select * from parquet_scan(?)",
            (str(out_path),),
        )
        con_duck.close()
        print(f"Saved boundaries DuckDB to {duckdb_path}")

    print("\nProcessing complete!")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Process Vizgen dataset to generate PMTiles and transcript data."
    )
    parser.add_argument(
        "dataset_dir",
        type=Path,
        help="Path to the root Vizgen dataset directory (contains images/)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Destination directory (default: ./output)",
    )
    parser.add_argument(
        "--max-zoom",
        type=int,
        default=None,
        help="Maximum zoom level (auto-calculated if omitted)",
    )
    parser.add_argument(
        "--generate-pmtiles",
        action="store_true",
        default=False,
        help="Generate PMTiles for all stains",
    )
    parser.add_argument(
        "--generate-for-stain",
        type=str,
        default=None,
        help="Generate PMTiles for a specific stain",
    )
    parser.add_argument(
        "--generate-transcripts",
        action="store_true",
        default=False,
        help="Generate transcripts DuckDB database from detected_transcripts.csv",
    )
    parser.add_argument(
        "--generate-h5ad",
        action="store_true",
        default=False,
        help="Generate h5ad using cell_by_gene.csv and cell_metadata.csv",
    )
    parser.add_argument(
        "--generate-boundaries",
        action="store_true",
        default=False,
        help="Generate scaled cell boundaries GeoParquet",
    )

    args = parser.parse_args()

    process_vizgen_data(
        dataset_dir=args.dataset_dir,
        output_dir=args.output_dir,
        max_zoom=args.max_zoom,
        generate_pmtiles=args.generate_pmtiles,
        generate_for_stain=args.generate_for_stain,
        generate_transcripts=args.generate_transcripts,
        generate_h5ad=args.generate_h5ad,
        generate_boundaries=args.generate_boundaries,
    )


if __name__ == "__main__":
    main()
