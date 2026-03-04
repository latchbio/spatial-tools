import argparse
import json
import sys
import time
from io import TextIOWrapper
from pathlib import Path
from typing import cast

from PIL import Image
from PIL.TiffTags import TAGS_V2

from spatial_tools.args import BaseArguments
from spatial_tools.ctx import Ctx
from spatial_tools.pmtiles import write_pmtiles
from spatial_tools.rect import Fov
from spatial_tools.vec2 import Vec2
from spatial_tools.write_tile import Tile, iterate_tiles, write_tile

from .transcripts_from_raw import create_h5ad, get_transcripts


def process_cosmx_data(
    *,
    cell_stats_dir: Path,
    transcript_dir: Path,
    output_dir: Path,
    max_zoom: int | None = None,
    generate_morphology_pmtiles: bool = False,
    generate_composite_pmtiles: bool = False,
    generate_overlay_pmtiles: bool = False,
    generate_transcripts: bool = False,
    morphology_dir: Path | None = None,
    composite_dir: Path | None = None,
    overlay_dir: Path | None = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    any_generate = (
        generate_morphology_pmtiles
        or generate_composite_pmtiles
        or generate_overlay_pmtiles
        or generate_transcripts
    )
    if not any_generate:
        print(
            "No output generated: pass at least one of "
            "--generate-morphology-pmtiles, --generate-composite-pmtiles, "
            "--generate-overlay-pmtiles, --generate-transcripts"
        )
        return

    morph_2d_dir = (
        morphology_dir
        if morphology_dir is not None
        else cell_stats_dir / "Morphology2D"
    )
    composite_dir_resolved = (
        composite_dir if composite_dir is not None else cell_stats_dir / "CellComposite"
    )
    overlay_dir_resolved = (
        overlay_dir if overlay_dir is not None else cell_stats_dir / "CellOverlay"
    )

    image_dirs = {
        "morphology": morph_2d_dir,
        "composite": composite_dir_resolved,
        "overlay": overlay_dir_resolved,
    }

    base_args = BaseArguments()
    base_args.output = output_dir
    ctx = Ctx(args=base_args, log_f=cast(TextIOWrapper, sys.stdout))
    ctx.mirrored_x = True

    min_pos_mm = Vec2(float("inf"), float("inf"))
    max_extent_mm = Vec2(float("-inf"), float("-inf"))

    def find_image(dir_p: Path, fov_id: str | int) -> Path:
        fov_id_int = int(fov_id)
        matching_files = []
        for f in dir_p.glob("*_F*.*"):
            stem = f.stem
            if "_F" not in stem:
                continue
            suffix = stem.split("_F")[-1]
            try:
                if int(suffix) == fov_id_int:
                    matching_files.append(f)
            except ValueError:
                continue
        if len(matching_files) != 1:
            raise FileNotFoundError(
                f"0 or > 1 file found for FOV {fov_id} in {dir_p} (found {len(matching_files)} with _F{{number}} in name)"
            )
        return matching_files[0]

    for p in morph_2d_dir.glob("*_F*.TIF"):
        with Image.open(p) as img:
            tags = {TAGS_V2[k].name: v for k, v in img.tag_v2.items()}
            img_desc = json.loads(tags["ImageDescription"])

            mm_per_pix = img_desc["PixelSize_um"] / img_desc["Magnification"] / 1000
            if ctx._mm_per_px is not None and ctx.mm_per_px != mm_per_pix:
                raise ValueError(
                    f"Inconsistent mm_per_px: {mm_per_pix} vs {ctx.mm_per_px}"
                )
            ctx.mm_per_px = mm_per_pix

            cur = Fov(
                ctx=ctx,
                paths={"morphology": p},
                id=str(img_desc["Fov"]),
                pos_mm=Vec2(img_desc["X_mm"], img_desc["Y_mm"]),
                size_px=Vec2.from_tuple(img.size),
            )
            for k, dir_p in image_dirs.items():
                cur.paths[k] = find_image(dir_p, cur.id)

            min_pos_mm = min_pos_mm.pw_min(cur.pos_mm)
            max_extent_mm = max_extent_mm.pw_max(cur.pos_mm + cur.size_mm)

            ctx.fovs.append(cur)
            ctx.fovs_by_id[cur.id] = cur

    if not ctx.fovs:
        raise ValueError("No FOV images found")

    ctx.slide.pos_mm = min_pos_mm
    ctx.slide.size_px = (max_extent_mm - min_pos_mm) / ctx.mm_per_px

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
            f"{z}: {tile.size_mm.size_str('mm')} | {tile.scale:.2f}x zoom, FOVs {tile.fov_size_spx(ctx.fovs[0]).size_str('px')}"
        )
    print()

    if generate_morphology_pmtiles:
        print("Generating Morphology PMTiles...")
        start = time.monotonic()
        stitched_dir = output_dir / "stitched_morphology"
        stitched_dir.mkdir(exist_ok=True)

        for tile in iterate_tiles(ctx):
            write_tile(
                ctx, tile, stitched_dir, category="morphology", color_mode="I;16"
            )

        pmtiles_path = output_dir / "morphology.pmtiles"
        write_pmtiles(ctx, stitched_dir, pmtiles_path)
        print(f"Morphology PMTiles generated in {time.monotonic() - start:.2f}s")

    if generate_composite_pmtiles:
        print("\nGenerating Composite PMTiles...")
        start = time.monotonic()
        stitched_dir = output_dir / "stitched_composite"
        stitched_dir.mkdir(exist_ok=True)

        for tile in iterate_tiles(ctx):
            write_tile(ctx, tile, stitched_dir, category="composite")

        pmtiles_path = output_dir / "composite.pmtiles"
        write_pmtiles(ctx, stitched_dir, pmtiles_path)
        print(f"Composite PMTiles generated in {time.monotonic() - start:.2f}s")

    if generate_overlay_pmtiles:
        print("\nGenerating Overlay PMTiles...")
        start = time.monotonic()
        stitched_dir = output_dir / "stitched_overlay"
        stitched_dir.mkdir(exist_ok=True)

        for tile in iterate_tiles(ctx):
            write_tile(ctx, tile, stitched_dir, category="overlay")

        pmtiles_path = output_dir / "overlay.pmtiles"
        write_pmtiles(ctx, stitched_dir, pmtiles_path)
        print(f"Overlay PMTiles generated in {time.monotonic() - start:.2f}s")

    if generate_transcripts:
        print("\nProcessing transcripts...")
        start = time.monotonic()
        transcripts_con = get_transcripts(transcript_dir, ctx)

        db_path = output_dir / "transcripts.duckdb"
        transcripts_con.execute(f"attach '{db_path}' AS persistent_db")
        transcripts_con.execute(
            "create or replace table persistent_db.final_transcripts as select * from final_transcripts"
        )
        transcripts_con.execute("detach persistent_db")

        print(f"Transcripts generated in {time.monotonic() - start:.2f}s")

        print("\nGenerating h5ad file...")
        h5ad_start = time.monotonic()
        h5ad_path = output_dir / "transcripts.h5ad"
        create_h5ad(transcripts_con, h5ad_path)
        print(f"H5ad file generated in {time.monotonic() - h5ad_start:.2f}s")

    print("\nProcessing complete!")
    if generate_morphology_pmtiles:
        print(f"Morphology PMTiles output: {output_dir / 'morphology.pmtiles'}")
    if generate_composite_pmtiles:
        print(f"Composite PMTiles output: {output_dir / 'composite.pmtiles'}")
    if generate_overlay_pmtiles:
        print(f"Overlay PMTiles output: {output_dir / 'overlay.pmtiles'}")
    if generate_transcripts:
        print(f"Transcripts output: {output_dir / 'transcripts.duckdb'}")
        print(f"H5ad output: {output_dir / 'transcripts.h5ad'}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Process CosMx output directory to generate PMTiles and transcript data"
    )
    parser.add_argument(
        "cell_stats_dir", type=Path, help="Path to CosMx CellStats directory"
    )
    parser.add_argument(
        "transcript_dir", type=Path, help="Path to CosMx transcript directory"
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
        "--generate-morphology-pmtiles",
        action="store_true",
        default=False,
        help="Generate morphology PMTiles",
    )
    parser.add_argument(
        "--generate-composite-pmtiles",
        action="store_true",
        default=False,
        help="Generate composite PMTiles",
    )
    parser.add_argument(
        "--generate-overlay-pmtiles",
        action="store_true",
        default=False,
        help="Generate overlay PMTiles",
    )
    parser.add_argument(
        "--generate-transcripts",
        action="store_true",
        default=False,
        help="Generate transcripts",
    )
    parser.add_argument(
        "--morphology-dir",
        type=Path,
        default=None,
        help="Custom path to Morphology2D directory (default: cell_stats_dir/Morphology2D)",
    )
    parser.add_argument(
        "--composite-dir",
        type=Path,
        default=None,
        help="Custom path to CellComposite directory (default: cell_stats_dir/CellComposite)",
    )
    parser.add_argument(
        "--overlay-dir",
        type=Path,
        default=None,
        help="Custom path to CellOverlay directory (default: cell_stats_dir/CellOverlay)",
    )

    args = parser.parse_args()

    if args.max_zoom is not None and args.max_zoom > 5:
        print("Fix header encoding in pmtiles.py to use zoom > 5")

    process_cosmx_data(
        cell_stats_dir=args.cell_stats_dir,
        transcript_dir=args.transcript_dir,
        output_dir=args.output_dir,
        max_zoom=args.max_zoom,
        generate_morphology_pmtiles=args.generate_morphology_pmtiles,
        generate_composite_pmtiles=args.generate_composite_pmtiles,
        generate_overlay_pmtiles=args.generate_overlay_pmtiles,
        generate_transcripts=args.generate_transcripts,
        morphology_dir=args.morphology_dir,
        composite_dir=args.composite_dir,
        overlay_dir=args.overlay_dir,
    )


if __name__ == "__main__":
    main()
