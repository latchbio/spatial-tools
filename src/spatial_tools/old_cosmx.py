import json
import time
from argparse import ArgumentParser
from math import ceil, log2
from pathlib import Path

from PIL import Image
from PIL.TiffTags import TAGS_V2

from .ctx import ctx
from .rect import Fov
from .vec2 import Vec2
from .write_tile import Tile, iterate_tiles, write_tile

ctx.mirrored_x = True

ctx.out_p = Path("stitched")
ctx.log_p = ctx.out_p / "info.txt"
ctx.log_p_f = ctx.log_p.open("w")

argp = ArgumentParser()
_ = argp.add_argument("--max-z", type=int, required=False)
args = argp.parse_args()

exts = {".jpg", ".jpeg", ".png", ".tiff", ".tif"}

root_p = Path("/root/NormalLiverFiles")
dir_p = root_p / "CellStatsDir/Morphology2D/"

ctx.log(f"Input: {dir_p}")

n_fovs = 0

for p in dir_p.iterdir():
    if not any(p.name.lower().endswith(x) for x in exts):
        continue

    with Image.open(p) as img:
        tags = {TAGS_V2[k].name: v for k, v in img.tag_v2.items()}

    img_desc = json.loads(tags["ImageDescription"])
    ctx.mm_per_px = img_desc["PixelSize_um"] / img_desc["Magnification"] / 1000

    n_fovs = img_desc["NFov"]
    ctx.log(
        f"{img_desc['NFov']} FOVs, {ctx.mm_per_px} mm/px ({1 / ctx.mm_per_px} px/mm)"
    )

    break

for p in dir_p.iterdir():
    if not any(p.name.lower().endswith(x) for x in exts):
        continue

    # todo(maximsmol): support non-tiff
    with Image.open(p) as img:
        tags = {TAGS_V2[k].name: v for k, v in img.tag_v2.items()}

        img_desc = json.loads(tags["ImageDescription"])
        cur = Fov(
            paths={"morphology": p},
            id=img_desc["Fov"],
            pos_mm=Vec2(img_desc["X_mm"], img_desc["Y_mm"]),
            size_px=Vec2.from_tuple(img.size),
        )

    my_mm_per_px = img_desc["PixelSize_um"] / img_desc["Magnification"] / 1000
    assert my_mm_per_px == ctx.mm_per_px

    ctx.add_fov(cur)

ctx.complete()

ctx.log()
ctx.log(f"Slide @ {ctx.slide} ({ceil(ctx.slide.size_px).size_str('px')})")
ctx.log(
    f"FOVs {ctx.fovs[0].size_mm.size_str('mm')} ({ctx.fovs[0].size_px.size_str('px')})"
)

for x in ctx.fovs:
    ctx.log(
        f"{x.id.rjust(len(str(n_fovs)))} @ {x.pos_str()}: {x.paths['morphology'].name}"
    )
ctx.log()

viewports_per_slide = ctx.slide.size_px.pw_div(ctx.viewport_size_px)
max_z = ceil(max(0, log2(viewports_per_slide.x), log2(viewports_per_slide.y)))

ctx.log(f"Max zoom level: {max_z}")
for z in range(max_z + 1):
    tile = Tile(z=z, pos_idx=Vec2(0, 0))

    if z == 0:
        print(f"Tile resolution: {tile.resolution().size_str('px')}")

    ctx.log(
        f"{z}: {tile.size_mm.size_str('mm')} | {tile.scale:.2f}x zoom, FOVs {tile.fov_size_spx(ctx.fovs[0]).size_str('px')}"
    )
ctx.log()

# >>>

if args.max_z is not None:
    ctx.log(f">>> Stopping at z={args.max_z} as requested <<<")
    max_z = int(args.max_z)
    ctx.log()

start_first = time.monotonic()

time_per_z: dict[int, float] = {}

last_z = 0
start = time.monotonic()


def finish_z_level() -> None:
    global start, last_z

    duration = time.monotonic() - start
    time_per_z[tile.z] = duration
    ctx.log(f">>> z={tile.z} time: {duration:.1f}s")
    ctx.log()

    start = time.monotonic()
    last_z = tile.z


for tile in iterate_tiles(max_z):
    if tile.z != last_z:
        finish_z_level()

    write_tile(tile, ctx.out_p, category="morphology")
finish_z_level()

ctx.log()
ctx.log(f"Total runtime: {time.monotonic() - start_first:.1f}s")
for z, d in time_per_z.items():
    ctx.log(f"z={z}: {d:.1f}s")
