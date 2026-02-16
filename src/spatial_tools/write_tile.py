import time
from collections.abc import Generator
from contextlib import contextmanager
from copy import copy
from dataclasses import dataclass
from math import ceil, floor
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import PIL
from PIL import Image
import pyvips

from .rect import Fov, Rect
from .vec2 import Vec2

if TYPE_CHECKING:
    from .ctx import Ctx


@dataclass(init=False)
class Tile(Rect):
    # note(maximsmol): redundant with the parent class but pyright is being dumb
    size_px: Vec2[float]
    pos_mm: Vec2[float]

    def __init__(self, *, ctx: "Ctx", z: int, pos_idx: Vec2[int]) -> None:
        super().__init__(
            ctx=ctx, pos_mm=copy(ctx.slide.pos_mm), size_px=copy(ctx.slide.size_px)
        )

        self.z: int = z
        self.pos_idx: Vec2[int] = pos_idx

        self.size_px /= int(2**z)
        self.pos_mm += pos_idx.pw_mul(self.size_mm)

        pw_scale = ctx.viewport_size_px.pw_div(self.size_px)
        self.scale: float = min(pw_scale.x, pw_scale.y)

    def resolution(self) -> Vec2[int]:
        return ceil(self.size_px * self.scale)

    def fov_pos_spx(self, fov: Fov) -> Vec2[int]:
        res = floor((fov.pos_mm - self.pos_mm) * self.scale / self.ctx.mm_per_px)
        mirrored = self.resolution() - self.fov_size_spx(fov) - res

        if self.ctx.mirrored_x:
            res.x = mirrored.x
        if self.ctx.mirrored_y:
            res.y = mirrored.y

        return res

    def fov_size_spx(self, fov: Fov) -> Vec2[int]:
        return ceil(fov.size_px * self.scale)


def iterate_tiles(ctx: "Ctx") -> Generator[Tile]:
    for z in range(ctx.max_z + 1):
        # todo(maximsmol): pyright wot
        l = 2**z
        assert isinstance(l, int)

        for x in range(l):
            for y in range(l):
                yield Tile(ctx=ctx, z=z, pos_idx=Vec2(x, y))


# todo(maximsmol): webp-compress every FOV first?
# todo(maximsmol): doing this inside-out is probably better
# i.e. loop over each fov and add it to each zoom tile that requires it
# todo(maximsmol): parallelize?
# todo(maximsmol): iterate over Z levels first and cache the open FOV images?
# todo(maximsmol): stitch from high Z to low, reusing previous levels as thumbnails


def _get_cached_fov_img(tile: Tile, fov: Fov, *, category: str) -> Image.Image:
    img = pyvips.Image.new_from_file(fov.paths[category])

    target_size = tile.fov_size_spx(fov)
    img = img.resize(
        target_size.x / img.width,
        vscale=target_size.y / img.height,
        kernel=pyvips.Kernel.NEAREST,
    )

    return img


def write_tile(
    ctx: "Ctx",
    tile: Tile,
    out_dir: Path,
    *,
    category: str,
    color_mode: Literal["I;16", "RGB"] = "RGB",
    rescale: Vec2[int] | None = None,
) -> None:
    start = time.monotonic()
    ctx.log(f"z={tile.z} x={tile.pos_idx.x} y={tile.pos_idx.y} @ {tile}")

    res_size = tile.resolution()
    channel = pyvips.Image.black(res_size.x, res_size.y)
    if color_mode == "I;16":
        res = channel.cast(pyvips.BandFormat.USHORT)
    else:
        res = channel.bandjoin([channel, channel])
    # todo(maximsmol): support color_mode

    total_fovs = 0
    for fov in ctx.fovs:
        if not tile.overlaps(fov):
            continue
        total_fovs += 1

    fov_pct = total_fovs / len(ctx.fovs)
    ctx.log(f"  Using {total_fovs}/{len(ctx.fovs)} FOVs ({fov_pct * 100:.2f}%)")

    i = 0
    for fov in ctx.fovs:
        if not tile.overlaps(fov):
            continue

        box_pos_spx = tile.fov_pos_spx(fov)

        progress_pct = i / total_fovs
        ctx.log(
            f"  {i}/{total_fovs} ({progress_pct * 100:.2f}%): {box_pos_spx.x}px, {box_pos_spx.y}px <- FOV {fov.id} @ {fov.pos_str()}"
        )
        i += 1

        fov_img = _get_cached_fov_img(tile, fov, category=category)
        if fov_img.hasalpha():
            fov_img = fov_img.flatten()

        # todo(maximsmol): use arrayjoin?
        res = res.insert(fov_img, box_pos_spx.x, box_pos_spx.y)
        del fov_img

    # todo(maximsmol): implement
    # if color_mode == "I;16" and rescale is not None:
    #     lo, hi = rescale.to_tuple()
    #     res = res.convert("I").point(
    #         [(x - lo) / (hi - lo) * 255 for x in range(256 * 256)], "L"
    #     )

    res_p = out_dir / f"{tile.z}/{tile.pos_idx.x}-{tile.pos_idx.y}.webp"
    res_p.parent.mkdir(parents=True, exist_ok=True)

    res.write_to_file(res_p)
    del res

    ctx.log(f"  Time {time.monotonic() - start:.1f}s")
