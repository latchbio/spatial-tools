from tifffile import TiffFile

from .args import SvsArguments
from .ctx import Ctx
from .rect import Fov
from .vec2 import Vec2


def populate_ctx(ctx: Ctx, args: SvsArguments) -> None:
    ctx.mm_per_px = 1
    ctx.log(f"{ctx.mm_per_px} mm/px ({1 / ctx.mm_per_px} px/mm)")
    ctx.log()

    with TiffFile(args.input) as img:
        sizes = img.pages[0].sizes
        ctx.add_fov(
            Fov(
                ctx=ctx,
                paths={"main": args.input},
                id="all",
                pos_mm=Vec2(0, 0),
                size_px=Vec2(sizes["width"], sizes["height"]),
            )
        )
