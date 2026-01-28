import sys
import time
from io import TextIOWrapper
from pathlib import Path

from PIL import Image

from .args import get_args
from .ctx import Ctx
from .pmtiles import write_pmtiles
from .utils import format_duration
from .write_tile import iterate_tiles, write_tile

Image.MAX_IMAGE_PIXELS = None

assert isinstance(sys.stdout, TextIOWrapper)
sys.stdout.reconfigure(line_buffering=True)
assert isinstance(sys.stderr, TextIOWrapper)
sys.stderr.reconfigure(line_buffering=True)


def main() -> None:
    args = get_args()

    args.output.mkdir(parents=True, exist_ok=True)

    log_p = args.output / "log.txt"
    with log_p.open("w") as log_f:
        ctx = Ctx(args=args, log_f=log_f)

        ctx.log(f"Input: {args.input}")
        ctx.log(f"Output: {args.output}")

        if args.command != "visium":
            if args.dryrun:
                ctx.log("⚠️ Will not write any output due to --dryrun")

            if args.command == "element" and args.single_well is not None:
                ctx.log(f"⚠️ Will only process well {args.single_well} due to --single-well")

            if args.single_category is not None:
                ctx.log(
                    f"⚠️ Will only process category {args.single_category} due to --single-category"
                )

        ctx.log()
        match args.command:
            case "visium":
                from .visium import populate_ctx  # noqa: PLC0415

                populate_ctx(ctx, args)
                return

            case "element":
                from .element import populate_ctx  # noqa: PLC0415

                populate_ctx(ctx, args)

            case "svs":
                from .svs import populate_ctx  # noqa: PLC0415

                populate_ctx(ctx, args)

            case "xenium":
                from .xenium import populate_ctx  # noqa: PLC0415

                populate_ctx(ctx, args)

            case x:
                raise RuntimeError(f"implementation error: unhandled command: {x!r}")

        ctx.complete()

        match args.command:
            case "xenium":
                from .xenium import generate_extras  # noqa: PLC0415

                generate_extras(ctx, args)

            case x:
                ...

        start = time.monotonic()

        for category in ctx.categories:
            if args.single_category is not None and category != args.single_category:
                continue

            cat_p = ctx.viewports_p / category
            if not (ctx.viewports_p / category).is_dir():
                time_per_z: dict[int, float] = {}
                start_cat = time.monotonic()

                last_z = 0
                start_z = time.monotonic()

                def finish_z_level(z: int) -> None:
                    nonlocal start_z

                    duration = time.monotonic() - start_z
                    time_per_z[z] = duration
                    ctx.log(f">>> z={z} time: {format_duration(duration)}")
                    ctx.log()

                    start_z = time.monotonic()

                quants = ctx.quantiles.get(category)
                if quants is not None:
                    ctx.log(f"{category} {{ Data range: {quants.x}-{quants.y}")
                else:
                    ctx.log(f"{category} {{")

                for tile in iterate_tiles(ctx):
                    if tile.z != last_z:
                        finish_z_level(tile.z)
                        last_z = tile.z

                    if not args.dryrun:
                        write_tile(
                            ctx,
                            tile,
                            cat_p,
                            category=category,
                            # todo(maximsmol): auto-detect this based on source FOV images
                            # todo(maximsmol): should be different for each category for e.g. CosMX
                            color_mode="I;16" if args.command in {"element", "xenium"} else "RGB",
                            rescale=quants,
                        )
                finish_z_level(last_z)
                ctx.log(
                    f"}} {category} total: {format_duration(time.monotonic() - start_cat)}"
                )

                for z, d in time_per_z.items():
                    ctx.log(f"z={z}: {d:.1f}s")
                ctx.log()
            else:
                ctx.log(f"{category} already stitched")
                ctx.log()

            start_cat = time.monotonic()
            ctx.log("Writing PMTiles")
            if not args.dryrun:
                write_pmtiles(
                    ctx,
                    cat_p,
                    (args.output / "pmtiles" / category).with_suffix(".pmtiles"),
                    max_zoom=ctx.max_z,
                    verbose=args.print_pmtiles_progress,
                )
            ctx.log(f"  Done in {format_duration(time.monotonic() - start_cat)}")
            ctx.log()

        ctx.log(f"Total duration: {format_duration(time.monotonic() - start)}")
        ctx.log()


if __name__ == "__main__":
    main()
