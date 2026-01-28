from argparse import ArgumentParser
from pathlib import Path
from typing import Literal, cast


class BaseArguments:
    input: Path = Path()
    output: Path = Path()
    print_fov_paths: bool = False
    max_z: int | None = None
    single_category: str | None = None
    no_rescale: bool = False
    dryrun: bool = False
    print_pmtiles_progress: bool = False


class CosmxArguments(BaseArguments):
    command: Literal["cosmx"] = "cosmx"


class ElementArguments(BaseArguments):
    command: Literal["element"] = "element"
    print_well_fovs: bool = False
    single_well: str | None = None


class SvsArguments(BaseArguments):
    command: Literal["svs"] = "svs"


class XeniumArguments(BaseArguments):
    command: Literal["xenium"] = "xenium"


class VisiumArguments(BaseArguments):
    command: Literal["visium"] = "visium"
    bin_name: str = "square_008um"
    matrix_type: str = "filtered"
    no_analysis: bool = False


def get_args() -> CosmxArguments | ElementArguments | SvsArguments | XeniumArguments | VisiumArguments:
    argp = ArgumentParser()

    subp = argp.add_subparsers(dest="command")
    cosmx_p = subp.add_parser("cosmx")
    element_p = subp.add_parser("element")
    svs_p = subp.add_parser("svs")
    xenium_p = subp.add_parser("xenium")
    visium_p = subp.add_parser("visium")

    # Base arguments for commands that use the Ctx pattern
    for parser in [cosmx_p, element_p, svs_p, xenium_p, visium_p]:
        _ = parser.add_argument("input", type=Path, help="Input data")
        _ = parser.add_argument(
            "--output", default=Path("stitched"), type=Path, help="Output location"
        )
        _ = parser.add_argument(
            "--print-fov-paths",
            action="store_true",
            default=False,
            help="List the paths of the images for each FOV",
        )
        _ = parser.add_argument(
            "--max-z", type=int, required=False, help="Only stitch up to a given Z-level"
        )
        _ = parser.add_argument(
            "--single-category",
            type=str,
            required=False,
            help="Only output a specific category",
        )
        _ = parser.add_argument(
            "--no-rescale",
            action="store_true",
            default=False,
            help="Clip HDR images rather than rescale the value range (faster but some images can be unreadable)",
        )
        _ = parser.add_argument(
            "--dryrun", action="store_true", default=False, help="Do not write output files"
        )
        _ = parser.add_argument(
            "--print-pmtiles-progress",
            action="store_true",
            default=False,
            help="Print progress logs while writing pmtiles",
        )

    _ = element_p.add_argument(
        "--single-well", type=str, required=False, help="Only output a specific well"
    )
    _ = element_p.add_argument(
        "--print-well-fovs",
        action="store_true",
        default=False,
        help="List FOVs for each well",
    )

    _ = visium_p.add_argument(
        "--bin-name",
        default="square_008um",
        help="Bin directory name under binned_outputs/ (default: square_008um). Examples: square_002um, square_016um.",
    )
    _ = visium_p.add_argument(
        "--matrix-type",
        choices=["filtered", "raw"],
        default="filtered",
        help="Which matrix to load (default: filtered). 'filtered' uses filtered_feature_bc_matrix/, 'raw' uses raw_feature_bc_matrix/.",
    )
    _ = visium_p.add_argument(
        "--no-analysis",
        action="store_true",
        help="Skip loading analysis outputs (PCA, UMAP, clustering, differential expression). By default, these are loaded if present.",
    )

    return cast(
        CosmxArguments | ElementArguments | SvsArguments | XeniumArguments | VisiumArguments,
        cast(object, argp.parse_args()),
    )
