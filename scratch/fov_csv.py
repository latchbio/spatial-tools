from pathlib import Path
from typing import TypedDict

root_p = Path("/root/NormalLiverFiles")
from csv import DictReader


class FovRow(TypedDict):
    Slide: int
    X_mm: float
    Y_mm: float
    Z_mm: float
    ZOffset_mm: float
    ROI: int
    FOV: int
    Order: str | None


fov_rows: dict[int, FovRow] = {}
with (root_p / "RunSummary/latest.fovs.csv").open("r", encoding="utf-8") as f:
    for row in DictReader(
        f,
        # fieldnames from the Napari plugin's _stitch.py
        fieldnames=[
            "Slide",
            "X_mm",
            "Y_mm",
            "Z_mm",
            "ZOffset_mm",
            "ROI",
            "FOV",
            "Order",
        ],
    ):
        data: dict[str, object] = {}
        for k, t in FovRow.__annotations__.items():
            cur = row[k]
            if cur is not None:
                cur = cur.strip()

            if t is int or t is float:
                cur = t(cur)

            data[k] = cur

        fov_rows[data["FOV"]] = data
