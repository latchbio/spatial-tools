import time
from collections.abc import Generator
from contextlib import contextmanager
from copy import copy
from dataclasses import dataclass
from math import ceil, floor
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np
import PIL
from PIL import Image

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


# todo(maximsmol): try pyvips


# todo(maximsmol): webp-compress every FOV first?
# todo(maximsmol): doing this inside-out is probably better
# i.e. loop over each fov and add it to each zoom tile that requires it
# todo(maximsmol): parallelize?
# todo(maximsmol): iterate over Z levels first and cache the open FOV images?
# todo(maximsmol): stitch from high Z to low, reusing previous levels as thumbnails

_cached_fov_imgs: dict[str, Image.Image] = {}
_cached_zoom: int | None = None
_cached_category: str | None = None


def _ome_tiff_to_2d(data: np.ndarray) -> np.ndarray:
    """Reduce OME-TIFF / multi-dim array to 2D (Y, X) for display."""
    data = np.asarray(data)
    if data.ndim == 2:
        return data
    if data.ndim == 1:
        return data.reshape(1, -1)
    # 3+ D: use the two largest dimensions as spatial (Y, X), max over the rest.
    shape = data.shape
    axes_by_size = sorted(range(data.ndim), key=lambda i: shape[i], reverse=True)
    keep_axes = tuple(axes_by_size[:2])
    max_axes = tuple(i for i in range(data.ndim) if i not in keep_axes)
    out = data.max(axis=max_axes)
    # Ensure order (Y, X) = (height, width); keep_axes may be (2,1) for (Z,Y,X).
    if keep_axes[0] > keep_axes[1]:
        out = np.swapaxes(out, 0, 1)
    return out


@contextmanager
def load_image(path: Path) -> Generator[Image.Image]:
    try:
        with Image.open(path) as img:
            yield img
            return
    except PIL.UnidentifiedImageError:
        import tifffile  # noqa: PLC0415

        # Fall back to tifffile for OME-TIFF and other complex TIFFs.
        data = tifffile.imread(path)
        data = _ome_tiff_to_2d(data)
        if np.issubdtype(data.dtype, np.floating):
            data = (np.clip(data, 0, 1) * 65535).astype(np.uint16)
        elif data.dtype != np.uint16 and np.issubdtype(data.dtype, np.integer):
            data = np.clip(data, 0, 65535).astype(np.uint16)
        # Use mode "I" (32-bit signed int) so paste() works reliably.
        img = Image.fromarray(data.astype(np.int32), mode="I")
        yield img
        return


def _get_cached_fov_img(tile: Tile, fov: Fov, *, category: str) -> Image.Image:
    global _cached_zoom, _cached_category

    if _cached_zoom != tile.z or _cached_category != category:
        _cached_fov_imgs.clear()
        _cached_zoom = tile.z
        _cached_category = category

    img = _cached_fov_imgs.get(fov.id)
    if img is not None:
        return img

    image_path = fov.paths[category]

    with load_image(image_path) as fov_img:
        fov_img = fov_img.resize(
            tile.fov_size_spx(fov).to_tuple(), resample=Image.Resampling.NEAREST
        )

        _cached_fov_imgs[fov.id] = fov_img.copy()

    return _cached_fov_imgs[fov.id]


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

    # Use "I" (32-bit int) internally for 16-bit sources so paste() is reliable.
    canvas_mode = "I" if color_mode == "I;16" else color_mode
    with Image.new(canvas_mode, tile.resolution().to_tuple()) as res:
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
            res.paste(fov_img, box=box_pos_spx.to_tuple())

        if color_mode == "I;16":
            lo, hi = rescale.to_tuple() if rescale is not None else (0, 65535)
            if hi <= lo:
                lo, hi = 0, 65535
            lut = [
                int(max(0, min(255, (x - lo) / (hi - lo) * 255)))
                for x in range(256 * 256)
            ]
            res = res.point(lut, "L")

        # Always save as RGB so viewers that expect 3-channel tiles work correctly.
        if res.mode != "RGB":
            res = res.convert("RGB")

        res_p = out_dir / f"{tile.z}/{tile.pos_idx.x}-{tile.pos_idx.y}.webp"
        res_p.parent.mkdir(parents=True, exist_ok=True)
        res.save(res_p, format="webp")

        ctx.log(f"  Time {time.monotonic() - start:.1f}s")
