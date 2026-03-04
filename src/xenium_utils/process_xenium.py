import argparse
import sys
import time
from io import TextIOWrapper
from pathlib import Path
from typing import cast

from spatial_tools.args import XeniumArguments
from spatial_tools.ctx import Ctx
from spatial_tools.pmtiles import write_pmtiles
from spatial_tools.vec2 import Vec2
from spatial_tools.write_tile import Tile, iterate_tiles, write_tile
from spatial_tools.xenium import generate_extras, populate_ctx


def process_xenium_data(
    *,
    xenium_dir: Path,
    output_dir: Path,
    max_zoom: int | None = None,
    generate_mip_pmtiles: bool = False,
    generate_focus_pmtiles: bool = False,
    generate_transcripts: bool = False,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    any_generate = (
        generate_mip_pmtiles or generate_focus_pmtiles or generate_transcripts
    )
    if not any_generate:
        print(
            "No output generated: pass at least one of "
            "--generate-mip-pmtiles, --generate-focus-pmtiles, "
            "--generate-transcripts"
        )
        return

    args = XeniumArguments()
    args.input = xenium_dir
    args.output = output_dir
    args.no_rescale = False
    args.dryrun = False

    ctx = Ctx(args=args, log_f=cast(TextIOWrapper, sys.stdout))

    populate_ctx(ctx, args)
    ctx.complete()

    if max_zoom is not None:
        ctx._max_z = max_zoom

    print(f"\nSlide @ {ctx.slide} ({ctx.slide.size_px.size_str('px')})")
    print(
        f"FOVs {ctx.fovs[0].size_mm.size_str('mm')} ({ctx.fovs[0].size_px.size_str('px')})"
    )
    for fov in ctx.fovs:
        print(f"{fov.id.rjust(len(str(len(ctx.fovs))))} @ {fov.pos_str()}")
    print()

    print(f"Max zoom level: {ctx.max_z}")
    for z in range(ctx.max_z + 1):
        tile = Tile(ctx=ctx, z=z, pos_idx=Vec2(0, 0))
        if z == 0:
            print(f"Tile resolution: {tile.resolution().size_str('px')}")
        print(
            f"{z}: {tile.size_mm.size_str('mm')} | {tile.scale:.2f}x zoom,"
            f" FOVs {tile.fov_size_spx(ctx.fovs[0]).size_str('px')}"
        )
    print()

    if generate_mip_pmtiles:
        print("Generating MIP PMTiles...")
        start = time.monotonic()
        stitched_dir = output_dir / "stitched_mip"
        stitched_dir.mkdir(exist_ok=True)

        for tile in iterate_tiles(ctx):
            write_tile(ctx, tile, stitched_dir, category="mip", color_mode="I;16")

        pmtiles_path = output_dir / "mip.pmtiles"
        write_pmtiles(ctx, stitched_dir, pmtiles_path)
        print(f"MIP PMTiles generated in {time.monotonic() - start:.2f}s")

    if generate_focus_pmtiles:
        print("\nGenerating Focus PMTiles...")
        start = time.monotonic()
        stitched_dir = output_dir / "stitched_focus"
        stitched_dir.mkdir(exist_ok=True)

        for tile in iterate_tiles(ctx):
            write_tile(ctx, tile, stitched_dir, category="focus", color_mode="I;16")

        pmtiles_path = output_dir / "focus.pmtiles"
        write_pmtiles(ctx, stitched_dir, pmtiles_path)
        print(f"Focus PMTiles generated in {time.monotonic() - start:.2f}s")

    if generate_transcripts:
        print("\nGenerating transcripts, boundaries, and h5ad...")
        start = time.monotonic()
        generate_extras(ctx, args)
        print(
            f"Transcripts/boundaries/h5ad generated in {time.monotonic() - start:.2f}s"
        )

    print("\nProcessing complete!")
    if generate_mip_pmtiles:
        print(f"  MIP PMTiles: {output_dir / 'mip.pmtiles'}")
    if generate_focus_pmtiles:
        print(f"  Focus PMTiles: {output_dir / 'focus.pmtiles'}")
    if generate_transcripts:
        print(f"  Cell boundaries: {output_dir / 'cell_boundaries.duckdb'}")
        print(f"  Transcripts: {output_dir / 'transcripts.duckdb'}")
        print(f"  H5ad: {output_dir / 'cell_by_gene.h5ad'}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Process 10x Xenium outs directory to generate PMTiles and transcript data"
    )
    parser.add_argument(
        "xenium_dir",
        type=Path,
        help="Path to Xenium outs directory (contains experiment.xenium)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Output directory (default: ./output)",
    )
    parser.add_argument(
        "--max-zoom",
        type=int,
        default=None,
        help="Maximum zoom level for PMTiles (default: auto-calculated from slide size)",
    )
    parser.add_argument(
        "--generate-mip-pmtiles",
        action="store_true",
        default=False,
        help="Generate MIP PMTiles",
    )
    parser.add_argument(
        "--generate-focus-pmtiles",
        action="store_true",
        default=False,
        help="Generate focus PMTiles",
    )
    parser.add_argument(
        "--generate-transcripts",
        action="store_true",
        default=False,
        help="Generate transcripts DuckDB, cell boundaries DuckDB, and cell_by_gene h5ad",
    )

    args = parser.parse_args()
    process_xenium_data(
        xenium_dir=args.xenium_dir,
        output_dir=args.output_dir,
        max_zoom=args.max_zoom,
        generate_mip_pmtiles=args.generate_mip_pmtiles,
        generate_focus_pmtiles=args.generate_focus_pmtiles,
        generate_transcripts=args.generate_transcripts,
    )


if __name__ == "__main__":
    main()
