import argparse
import json
import time
from math import ceil, log2
from pathlib import Path

import duckdb
from PIL import Image
from PIL.TiffTags import TAGS_V2

from .ctx import ctx
from .pmtiles import write_pmtiles
from .rect import Fov
from .transcripts_from_raw import create_h5ad, get_transcripts
from .vec2 import Vec2
from .write_tile import Tile, iterate_tiles, write_tile


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
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    min_pos_mm = Vec2(float("inf"), float("inf"))
    max_extent_mm = Vec2(float("-inf"), float("-inf"))

    morph_2d_dir = cell_stats_dir / "Morphology2D"
    composite_dir = cell_stats_dir / "CellComposite"
    overlay_dir = cell_stats_dir / "CellOverlay"

    image_dirs = {
        "morphology": morph_2d_dir,
        "composite": composite_dir,
        "overlay": overlay_dir,
    }

    def find_image(dir_p: Path, id: str) -> Path:
        fov_pattern = f"*_F{id.zfill(3)}.*"
        matching_files = list(dir_p.glob(fov_pattern))
        if len(matching_files) != 1:
            raise FileNotFoundError(
                f"0 or > 1 file found for FOV {fov.id} in {dir_p} matching pattern {fov_pattern}"
            )

        return matching_files[0]

    for p in morph_2d_dir.glob("*_F*.TIF"):
        with Image.open(p) as img:
            tags = {TAGS_V2[k].name: v for k, v in img.tag_v2.items()}
            img_desc = json.loads(tags["ImageDescription"])

            # Verify mm_per_px is consistent
            mm_per_pix = img_desc["PixelSize_um"] / img_desc["Magnification"] / 1000
            assert ctx.mm_per_px in {-1, mm_per_pix}, (
                f"Inconsistent mm_per_px: {mm_per_pix} vs {ctx.mm_per_px}"
            )
            ctx.mm_per_px = mm_per_pix

            cur = Fov(
                paths={"morphology": p},
                id=img_desc["Fov"],
                pos_mm=Vec2(img_desc["X_mm"], img_desc["Y_mm"]),
                size_px=Vec2.from_tuple(img.size),
            )
            for k, dir_p in image_dirs.items():
                cur.paths[k] = find_image(dir_p, cur.id)

            min_pos_mm = min_pos_mm.pw_min(cur.pos_mm)
            max_extent_mm = max_extent_mm.pw_max(cur.pos_mm + cur.size_mm)

            ctx.fovs.append(cur)

    if not ctx.fovs:
        raise ValueError("No FOV images found")

    ctx.slide.pos_mm = min_pos_mm
    ctx.slide.size_px = (max_extent_mm - min_pos_mm) / ctx.mm_per_px

    ctx.complete()

    print(f"\nSlide @ {ctx.slide} ({ctx.slide.size_px.size_str('px')})")
    print(
        f"FOVs {ctx.fovs[0].size_mm.size_str('mm')} ({ctx.fovs[0].size_px.size_str('px')})"
    )
    for fov in ctx.fovs:
        print(f"{fov.id.rjust(len(str(len(ctx.fovs))))} @ {fov.pos_str()}")
    print()

    if max_zoom is None:
        viewports_per_slide = ctx.slide.size_px.pw_div(ctx.viewport_size_px)
        max_zoom = ceil(
            max(0, log2(viewports_per_slide.x), log2(viewports_per_slide.y))
        )

    print(f"Max zoom level: {max_zoom}")
    for z in range(max_zoom + 1):
        tile = Tile(z=z, pos_idx=Vec2(0, 0))
        if z == 0:
            print(f"Tile resolution: {tile.resolution().size_str('px')}")
        print(
            f"{z}: {tile.size_mm.size_str('mm')} | {tile.scale:.2f}x zoom, FOVs {tile.fov_size_spx(ctx.fovs[0]).size_str('px')}"
        )
    print()

    # Generate PMTiles for morphology images
    if generate_morphology_pmtiles:
        print("Generating Morphology PMTiles...")
        start = time.monotonic()
        stitched_dir = output_dir / "stitched_morphology"
        stitched_dir.mkdir(exist_ok=True)

        for tile in iterate_tiles(max_zoom):
            write_tile(tile, stitched_dir, category="morphology", color_mode="I;16")

        pmtiles_path = output_dir / "morphology.pmtiles"
        write_pmtiles(stitched_dir, pmtiles_path)
        print(f"Morphology PMTiles generated in {time.monotonic() - start:.2f}s")

    # Generate PMTiles for composite images
    if generate_composite_pmtiles:
        print("\nGenerating Composite PMTiles...")
        start = time.monotonic()
        stitched_dir = output_dir / "stitched_composite"
        stitched_dir.mkdir(exist_ok=True)

        for tile in iterate_tiles(max_zoom):
            write_tile(tile, stitched_dir, category="composite")

        pmtiles_path = output_dir / "composite.pmtiles"
        write_pmtiles(stitched_dir, pmtiles_path)
        print(f"Composite PMTiles generated in {time.monotonic() - start:.2f}s")

    # Generate PMTiles for overlay images
    if generate_overlay_pmtiles:
        print("\nGenerating Overlay PMTiles...")
        start = time.monotonic()
        stitched_dir = output_dir / "stitched_overlay"
        stitched_dir.mkdir(exist_ok=True)

        for tile in iterate_tiles(max_zoom):
            write_tile(tile, stitched_dir, category="overlay")

        pmtiles_path = output_dir / "overlay.pmtiles"
        write_pmtiles(stitched_dir, pmtiles_path)
        print(f"Overlay PMTiles generated in {time.monotonic() - start:.2f}s")

    if generate_transcripts:
        print("\nProcessing transcripts...")
        start = time.monotonic()
        transcripts_con = get_transcripts(transcript_dir)

        db_path = output_dir / "transcripts.duckdb"
        persistent_con = duckdb.connect(str(db_path))

        transcripts_con.execute(f"attach '{db_path}' AS persistent_db")
        transcripts_con.execute(
            "create table persistent_db.final_transcripts as select * from final_transcripts"
        )
        transcripts_con.execute("detach persistent_db")

        persistent_con.close()
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
        default=5,
        type=int,
        help="Maximum zoom level for PMTiles (if not specified, will be calculated based on slide size)",
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

    args = parser.parse_args()

    if args.max_zoom > 5:
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
    )


if __name__ == "__main__":
    main()
