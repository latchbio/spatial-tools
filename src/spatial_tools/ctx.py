import sys
from dataclasses import dataclass
from io import TextIOWrapper
from math import ceil, log2
from pathlib import Path

from .args import BaseArguments
from .rect import Fov, Rect
from .vec2 import Vec2
from .write_tile import Tile


@dataclass(init=False)
class Ctx:
    args: BaseArguments
    log_f: TextIOWrapper

    _mm_per_px: float | None = None
    _categories: list[str] | None = None
    _max_z: int | None = None

    slide: Rect
    fovs: list[Fov]
    fovs_by_id: dict[str, Fov]
    quantiles: dict[str, Vec2[int]]

    mirrored_x: bool = False
    mirrored_y: bool = False

    viewport_size_px: Vec2[int]

    def __init__(self, *, args: BaseArguments, log_f: TextIOWrapper) -> None:
        self.args = args
        self.log_f = log_f
        self.slide = Rect(
            ctx=self,
            pos_mm=Vec2(float("inf"), float("inf")),
            size_px=Vec2(float("-inf"), float("-inf")),
        )

        self.fovs = []
        self.fovs_by_id = {}
        self.quantiles = {}

        self.viewport_size_px = Vec2(1920, 1080)

    @property
    def viewports_p(self) -> Path:
        return self.args.output / "viewports"

    @property
    def categories(self) -> list[str]:
        if self._categories is None:
            raise RuntimeError(
                "implementation error: ctx.categories accessed before ctx.complete()"
            )
        return self._categories

    @property
    def max_z(self) -> int:
        if self._max_z is None:
            raise RuntimeError(
                "implementation error: ctx.max_z accessed before ctx.complete()"
            )
        return self._max_z

    @property
    def mm_per_px(self) -> float:
        if self._mm_per_px is None:
            raise RuntimeError("implementation error: ctx.mm_per_px never set")
        return self._mm_per_px

    @mm_per_px.setter
    def mm_per_px(self, x: float) -> None:
        self._mm_per_px = x

    def add_fov(self, x: Fov) -> None:
        self.fovs.append(x)
        self.fovs_by_id[x.id] = x
        self.slide.expand_to_contain(x)

    def complete(self) -> None:
        self.fovs.sort(key=lambda x: x.pos_mm.y)
        self.fovs.sort(key=lambda x: x.pos_mm.x)

        self._categories = sorted(self.fovs[0].paths.keys())
        self.log()
        self.log("Categories:")
        for cat in self.categories:
            self.log(f"- {cat}")

        for fov in self.fovs:
            assert set(fov.paths.keys()) == set(self.categories)

        n_fovs = len(self.fovs)

        self.log()
        self.log(f"{n_fovs} FOVs")
        self.log(f"Slide @ {self.slide} ({ceil(self.slide.size_px).size_str('px')})")
        self.log(
            f"FOVs {self.fovs[0].size_mm.size_str('mm')} ({self.fovs[0].size_px.size_str('px')})"
        )

        for x in self.fovs:
            self.log(f"{x.id.rjust(len(str(n_fovs)))} @ {x.pos_str()}")
            if self.args.print_fov_paths:
                for k, p in x.paths.items():
                    self.log(f"  {k}: {p.relative_to(self.args.input)}")
        self.log()

        viewports_per_slide = self.slide.size_px.pw_div(self.viewport_size_px)
        self._max_z = ceil(
            max(0, log2(viewports_per_slide.x), log2(viewports_per_slide.y))
        )

        self.log(f"Max zoom level: {self.max_z}")
        for z in range(self.max_z + 1):
            tile = Tile(ctx=self, z=z, pos_idx=Vec2(0, 0))

            if z == 0:
                self.log(f"Tile resolution: {tile.resolution().size_str('px')}")

            self.log(
                f"{z}: {tile.size_mm.size_str('mm')} | {tile.scale:.2f}x zoom, FOVs {tile.fov_size_spx(self.fovs[0]).size_str('px')}"
            )
        self.log()

        if self.args.max_z is not None:
            if self.max_z <= self.args.max_z:
                self.log(">>> Calculated z level below command-line limit <<<")
            else:
                self.log(f">>> Stopping at z={self.args.max_z} as requested <<<")
                self._max_z = self.args.max_z
            self.log()

    def log(self, *args: object, sep: str = " ", end: str = "\n") -> None:
        line = sep.join(str(x) for x in args) + end
        _ = self.log_f.write(line)
        self.log_f.flush()

        _ = sys.stdout.write(line)
