from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Self

from typing_extensions import override

from .vec2 import Vec2

if TYPE_CHECKING:
    from .ctx import Ctx


@dataclass(kw_only=True)
class Rect:
    ctx: "Ctx"
    pos_mm: Vec2[float]
    size_px: Vec2[float]

    @property
    def size_mm(self) -> Vec2[float]:
        return self.size_px * self.ctx.mm_per_px

    def overlaps(self, that: "Rect") -> bool:
        self_extent = self.pos_mm + self.size_mm
        that_extent = that.pos_mm + that.size_mm
        return self.pos_mm < that_extent and that.pos_mm < self_extent

    def pos_str(self) -> str:
        return f"{self.pos_mm.x:.1f}mm, {self.pos_mm.y:.1f}mm"

    @override
    def __str__(self) -> str:
        return f"{self.pos_str()} | {self.size_mm.size_str('mm')}"

    def expand_to_contain(self, x: Self) -> None:
        self.pos_mm = self.pos_mm.pw_min(x.pos_mm)
        self.size_px = self.size_px.pw_max(
            (x.pos_mm + x.size_mm - self.pos_mm) / self.ctx.mm_per_px
        )


@dataclass(kw_only=True)
class Fov(Rect):
    id: str
    paths: dict[str, Path]
