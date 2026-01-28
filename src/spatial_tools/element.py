import json
import re
from typing import Literal, TypedDict, cast

import numpy as np
from PIL import Image

from .args import ElementArguments
from .ctx import Ctx
from .rect import Fov, Rect
from .vec2 import Vec2


class ImageInfo(TypedDict):
    ImageWidth: int
    ImageHeight: int
    PixelSizeUm: float


class ElementTile(TypedDict):
    Name: str
    XMillimeters: float
    YMillimeters: float


class ElementWell(TypedDict):
    WellLocation: str
    Tiles: list[ElementTile]


class RunParameters(TypedDict):
    FileVersion: Literal["6.1.0"]
    RunName: str
    RunID: str
    RunType: str
    RunDescription: str
    Date: str
    OperatorName: str
    InstrumentName: str
    ImageInfo: ImageInfo
    Wells: list[ElementWell]


tile_name_re = re.compile(
    r"""
^
(
    (?P<prefix>
        CP\d{2}
    )
    _
)?
(?P<name>[^_]+)
(
    _
    (?P<suffix>
        [A-Za-z\-]+
    )
)?
\.
tif
$
""",
    re.VERBOSE | re.IGNORECASE,
)


def populate_ctx(ctx: Ctx, args: ElementArguments) -> None:
    run = cast(
        RunParameters, json.loads((args.input / "RunParameters.json").read_text())
    )

    assert run["FileVersion"] == "6.1.0"
    ctx.log(f"{run['RunName']} #{run['RunID']} ({run['RunType']})")
    ctx.log(run["RunDescription"])
    ctx.log(run["Date"])
    ctx.log(f"by {run['OperatorName']} on {run['InstrumentName']}")

    fov_w = run["ImageInfo"]["ImageWidth"]
    fov_h = run["ImageInfo"]["ImageHeight"]

    ctx.mm_per_px = run["ImageInfo"]["PixelSizeUm"] / 1000
    ctx.log(f"{ctx.mm_per_px} mm/px ({1 / ctx.mm_per_px} px/mm)")

    ctx.log()
    ctx.log(f"{len(run['Wells'])} wells")

    for well in run["Wells"]:
        if args.single_well is not None and well["WellLocation"] != args.single_well:
            continue

        tiles_by_name = {x["Name"]: x for x in well["Tiles"]}
        ctx.log(f"{well['WellLocation']}: {len(tiles_by_name)} tiles")

        well_rect = Rect(
            ctx=ctx,
            pos_mm=Vec2(float("inf"), float("inf")),
            size_px=Vec2(float("-inf"), float("-inf")),
        )

        well_p = args.input / "Projection" / f"Well{well['WellLocation']}"
        for p in well_p.iterdir():
            m = tile_name_re.match(p.name)
            assert m is not None

            name = m.group("name")
            category = m.group("suffix")
            if args.single_category is not None and category != args.single_category:
                continue

            tile = tiles_by_name[name]

            if not (ctx.viewports_p / category).is_dir():
                # if the category is already stitched, we want to skip the expensive image operations
                # since their results will not be used anyway
                #
                # we still want to collect the FOVs to construct a valid ctx.slide in case the pmtiles step runs
                with Image.open(p) as img:
                    assert img.width == fov_w
                    assert img.height == fov_h

                    if not args.no_rescale:
                        old = ctx.quantiles.setdefault(category, Vec2(65_535, 0))

                        ctx.log(f"- {p.relative_to(well_p)}")
                        data = np.array(img)
                        cur_np = np.quantile(data, [0.05, 0.95])
                        cur = Vec2(cast(int, cur_np[0]), cast(int, cur_np[1]))

                        ctx.quantiles[category] = Vec2(
                            min(old.x, cur.x), max(old.x, cur.y)
                        )

            cur = ctx.fovs_by_id.get(name)
            if cur is None:
                cur = Fov(
                    ctx=ctx,
                    paths={},
                    id=name,
                    # todo(maximsmol): figure out how this interacts with the mirroring
                    # in Ctx and in pmtiles
                    # pos_mm=Vec2(tile["XMillimeters"], tile["YMillimeters"])
                    pos_mm=Vec2(tile["YMillimeters"], -tile["XMillimeters"]),
                    size_px=Vec2(fov_w, fov_h),
                )
                ctx.add_fov(cur)

            cur.paths[category] = p
            well_rect.expand_to_contain(cur)

        ctx.log(f"@ {well_rect}")

        if args.print_well_fovs:
            tile_names = sorted(tiles_by_name.keys())
            for t in tile_names:
                fov = ctx.fovs_by_id[t]
                ctx.log(f"- {t} {fov}")
