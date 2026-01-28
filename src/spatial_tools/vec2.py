from dataclasses import dataclass
from math import ceil, floor
from typing import Any, Generic, TypeVar, overload

T = TypeVar("T", bound=int | float)
X = TypeVar("X", bound=int | float)


@dataclass
class Vec2(Generic[T]):
    x: T
    y: T

    def size_str(self, unit: str) -> str:
        if isinstance(self.x, int):
            return f"{self.x}{unit} x {self.y}{unit}"

        return f"{self.x:.1f}{unit} x {self.y:.1f}{unit}"

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, Vec2):
            return NotImplemented

        return self.x < other.x and self.y < other.y

    @classmethod
    def from_tuple(cls, that: tuple[X, X]) -> "Vec2[X]":
        # pyright is being an idiot
        return cls(that[0], that[1])

    def to_tuple(self) -> tuple[T, T]:
        return (self.x, self.y)

    @overload
    def pw_min(self, other: "Vec2[int]") -> "Vec2[T]": ...
    @overload
    def pw_min(self, other: "Vec2[float]") -> "Vec2[float]": ...
    def pw_min(self, other: "Vec2[Any]") -> "Vec2[Any]":
        return Vec2(min(self.x, other.x), min(self.y, other.y))

    @overload
    def pw_max(self, other: "Vec2[int]") -> "Vec2[T]": ...
    @overload
    def pw_max(self, other: "Vec2[float]") -> "Vec2[float]": ...
    def pw_max(self, other: "Vec2[Any]") -> "Vec2[Any]":
        return Vec2(max(self.x, other.x), max(self.y, other.y))

    @overload
    def pw_mul(self, other: "Vec2[int]") -> "Vec2[T]": ...
    @overload
    def pw_mul(self, other: "Vec2[float]") -> "Vec2[float]": ...
    def pw_mul(self, other: "Vec2[Any]") -> "Vec2[Any]":
        return Vec2(self.x * other.x, self.y * other.y)

    # Vec2[int | float] does not work because it collapses to
    # Vec2[float] and Vec2 is invariant
    @overload
    def pw_div(self, other: "Vec2[int]") -> "Vec2[float]": ...
    @overload
    def pw_div(self, other: "Vec2[float]") -> "Vec2[float]": ...
    def pw_div(self, other: "Vec2[Any]") -> "Vec2[float]":
        return Vec2(self.x / other.x, self.y / other.y)

    def __floor__(self) -> "Vec2[int]":
        return Vec2(floor(self.x), floor(self.y))

    def __ceil__(self) -> "Vec2[int]":
        return Vec2(ceil(self.x), ceil(self.y))

    def __mul__(self, other: object) -> "Vec2[T]":
        if not isinstance(other, int | float):
            return NotImplemented

        return Vec2(self.x * other, self.y * other)

    def __truediv__(self, other: object) -> "Vec2[float]":
        if not isinstance(other, int | float):
            return NotImplemented

        return Vec2(self.x / other, self.y / other)

    def __floordiv__(self, other: object) -> "Vec2[T]":
        if not isinstance(other, int | float):
            return NotImplemented

        return Vec2(self.x // other, self.y // other)

    def __add__(self, other: object) -> "Vec2[T]":
        if not isinstance(other, Vec2):
            return NotImplemented

        return Vec2(self.x + other.x, self.y + other.y)

    def __sub__(self, other: object) -> "Vec2[T]":
        if not isinstance(other, Vec2):
            return NotImplemented

        return Vec2(self.x - other.x, self.y - other.y)

    def __neg__(self) -> "Vec2[T]":
        return Vec2(-self.x, -self.y)
