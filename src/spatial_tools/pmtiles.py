import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

from .ctx import Ctx
from .utils import to_si


# https://github.com/protomaps/PMTiles/blob/26c857ff404f76dc0363f849b3b63113a00a0805/python/pmtiles/pmtiles/tile.py#L19C1-L43C15
def rotate(n: int, x: int, y: int, rx: int, ry: int) -> tuple[int, int]:
    if ry == 0:
        if rx != 0:
            x = n - 1 - x
            y = n - 1 - y
        x, y = y, x
    return x, y


def zxy_to_tileid(z: int, x: int, y: int) -> int:
    if z > 31:
        raise OverflowError("tile zoom exceeds 64-bit limit")
    if x > (1 << z) - 1 or y > (1 << z) - 1:
        raise ValueError("tile x/y outside zoom level bounds")

    acc = ((1 << (z * 2)) - 1) // 3
    a = z - 1
    while a >= 0:
        s = 1 << a
        rx = s & x
        ry = s & y
        acc += ((3 * rx) ^ ry) << a
        (x, y) = rotate(s, x, y, rx, ry)
        a -= 1
    return acc


@dataclass(kw_only=True)
class Tile:
    x: int
    y: int
    z: int
    offset: int
    len: int
    src: Path

    @property
    def tile_id(self) -> int:
        return zxy_to_tileid(self.z, self.x, self.y)

    def encode(self, prev: "Tile") -> bytes:
        return b"".join([
            encode_varint(self.tile_id - prev.tile_id),
            encode_varint(1),
            encode_varint(self.len),
            encode_varint(0)
            if self.offset == prev.offset + prev.len
            else encode_varint(self.offset + 1),
        ])


def encode_varint(x: int) -> bytes:
    res = b""

    while True:
        msb = 0b1000_0000 if x > 0b0111_1111 else 0
        res += (msb | (x & 0b0111_1111)).to_bytes(1, "little")

        x >>= 7

        if msb == 0:
            break

    return res


def encode_pos(pos: tuple[int, int]) -> bytes:
    return pos[0].to_bytes(4, "little", signed=True) + pos[1].to_bytes(
        4, "little", signed=True
    )


# https://github.com/protomaps/PMTiles/blob/e232df55745642b39f0cc1edfc85df2633bce29c/spec/v3/spec.md
def write_pmtiles(
    ctx: Ctx,
    input_dir: Path,
    output_path: Path,
    *,
    max_zoom: int | None = None,
    mirror_x: bool = False,
    verbose: bool = False,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("wb") as f:
        header_len = 127

        tile_entries: list[Tile] = []

        coords_re = re.compile(r"^(\d+)-(\d+)\.webp$", re.IGNORECASE)
        zoom_ps: list[Path] = []
        for zoom in input_dir.iterdir():
            try:
                z = int(zoom.name)
            except ValueError:
                continue

            if max_zoom is not None and z > max_zoom:
                continue

            zoom_ps.append(zoom)
        zoom_ps.sort(key=lambda x: x.name)

        for zoom in zoom_ps:
            z = int(zoom.name)

            old_l = len(tile_entries)
            total_size = 0
            for img in zoom.iterdir():
                m = coords_re.match(img.name)
                assert m is not None
                x_orig = int(m.group(1))
                y = int(m.group(2))

                # coxmx tile grid is produced with x increasing from right to left
                # vizgen tile grid is produced with x increasing from left to right
                if mirror_x:
                    max_x = (1 << z) - 1  # 2**z - 1
                    x_final = max_x - x_orig
                else:
                    x_final = x_orig

                size = img.stat().st_size

                tile_entries.append(
                    Tile(x=x_final, y=y, z=z, offset=0, len=size, src=img)
                )
                total_size += size

            ctx.log(f"z={z}: {len(tile_entries) - old_l} ({to_si(total_size)}B)")

        tile_entries.sort(key=lambda x: x.tile_id)
        ctx.log(
            f">>> Total: {len(tile_entries)} ({to_si(sum(t.len for t in tile_entries))}B)"
        )

        prev = None
        for t in tile_entries:
            p = prev.offset + prev.len if prev is not None else 0
            t.offset = p
            prev = t

        # assert len(tile_entries) >= 1
        root_dir = encode_varint(len(tile_entries))
        prev = None
        for t in tile_entries:
            p = prev.tile_id if prev is not None else 0
            root_dir += encode_varint(t.tile_id - p)
            prev = t

        for _ in tile_entries:
            root_dir += encode_varint(1)

        for t in tile_entries:
            root_dir += encode_varint(t.len)

        prev = None
        for t in tile_entries:
            p = prev.offset + prev.len if prev is not None else -1
            root_dir += (
                encode_varint(0) if t.offset == p else encode_varint(t.offset + 1)
            )
            prev = t

        scale0 = min(
            ctx.viewport_size_px.x / ctx.slide.size_px.x,
            ctx.viewport_size_px.y / ctx.slide.size_px.y,
        )

        slide_w = math.ceil(ctx.slide.size_px.x * scale0)
        slide_h = math.ceil(ctx.slide.size_px.y * scale0)
        meta = json.dumps({"slide_px": {"width": slide_w, "height": slide_h}}).encode(
            "utf-8"
        )
        meta_off = header_len
        meta_len = len(meta)

        root_dir_off = meta_off + meta_len
        root_dir_len = len(root_dir)

        assert root_dir_off + root_dir_len <= 16384

        leaf_off = 0
        leaf_len = 0

        tiles_off = root_dir_off + root_dir_len
        tiles_len = sum(x.len for x in tile_entries)

        num_tiles_before_rle = len(tile_entries)
        num_tile_entries = len(tile_entries)
        num_images = len(tile_entries)

        mul = 10000000
        min_pos = (-180 * mul, -85 * mul)
        max_pos = (180 * mul, 85 * mul)

        _ = f.write(b"PMTiles")
        _ = f.write(b"\x03")  # version

        _ = f.write(root_dir_off.to_bytes(8, "little"))
        _ = f.write(root_dir_len.to_bytes(8, "little"))

        _ = f.write(meta_off.to_bytes(8, "little"))
        _ = f.write(meta_len.to_bytes(8, "little"))

        _ = f.write(leaf_off.to_bytes(8, "little"))
        _ = f.write(leaf_len.to_bytes(8, "little"))

        _ = f.write(tiles_off.to_bytes(8, "little"))
        _ = f.write(tiles_len.to_bytes(8, "little"))

        _ = f.write(num_tiles_before_rle.to_bytes(8, "little"))
        _ = f.write(num_tile_entries.to_bytes(8, "little"))
        _ = f.write(num_images.to_bytes(8, "little"))

        _ = f.write(b"\x00")  # clustered
        _ = f.write(b"\x01")  # meta not compressed
        _ = f.write(b"\x01")  # tiles not compressed
        _ = f.write(b"\x04")  # tiles are webp

        _ = f.write(min(t.z for t in tile_entries).to_bytes(1, "little"))  # min zoom
        _ = f.write(max(t.z for t in tile_entries).to_bytes(1, "little"))  # max zoom

        _ = f.write(encode_pos(min_pos))
        _ = f.write(encode_pos(max_pos))

        _ = f.write((0).to_bytes(1, "little"))  # center zoom

        center_pos = (
            min_pos[0] + (max_pos[0] - min_pos[0]) // 2,
            min_pos[1] + (max_pos[1] - min_pos[1]) // 2,
        )
        _ = f.write(encode_pos(center_pos))  # center pos

        _ = f.write(meta)
        _ = f.write(root_dir)

        for t in tile_entries:
            if verbose:
                ctx.log(f"z={t.z} x={t.x} y={t.y}")
            with t.src.open("rb") as ft:
                while True:
                    chunk = ft.read(4096)
                    if len(chunk) == 0:
                        break

                    _ = f.write(chunk)
