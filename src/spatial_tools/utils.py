from datetime import timedelta


def to_si(x: float, *, base: int = 1000) -> str:
    for prefix in ["", "K", "M", "G", "T", "P"]:
        if x < base:
            return f"{x:.2f}{prefix}"

        x /= float(base)

    return f"{x}Y"


def format_duration(x: timedelta | float) -> str:
    parts = {"h": 60 * 60, "m": 60}

    if isinstance(x, timedelta):
        x = x.total_seconds()

    leftover = x

    res: list[str] = []
    for name, size in parts.items():
        cur, leftover = divmod(leftover, size)
        if cur > 0:
            res.append(f"{cur}{name}")

    res.append(f"{leftover:.1f}s")

    return " ".join(res)
