# Latch Spatial Widget Data Format

## Introduction

For spatial data visualization, the widget requires a directory in LData which has a predefined layout. This directory can contain:

1. Background images
2. Cell boundary database
3. Transcripts database

This information is overlayed on top of the cell coordinates in the H5/AnnData file. Typically there are multiple sets of coordinates—one each for the various dimensionality reduction techniques and other analyses, and one set of spatial coordinates. The widget does not have a way to recognize whether the coordinate set is spatial, so the user must select the correct coordinate set for spatial overlays to align.

## Converting Machine Outputs

Each spatial machine outputs data in different formats depending on the technology used. Latch has support for converting CosMX, Visium, Xenium, and VisGen data into the common format that the widget understands. These scripts are available upon request.

The rest of this document outlines the conversion process for users implementing support for new spatial technologies. An appendix is provided which includes high-level descriptions of how each officially supported technology was implemented.

## Background Images

### Overview

A Latch spatial directory may contain multiple background images. Typically they represent different stains of the slide, different z-levels/slices through the tissue, or different focus settings on the microscope.

The name of the file determines the name of the layer in the widget, so they should be descripive e.g. `Maxium Intensity Projection.pmtiles` instead of `2273917474_mip_out_final.pmtiles`

### PMTiles

To support quick zooming and panning via streaming from LData (backed by AWS S3), the files are stored in [a pyramidal format called PMTiles,](https://github.com/protomaps/PMTiles/tree/main/spec/v3) which was originally designed for streaming geographic data.

PMTiles is not well-supported by editing programs and does not have official libraries for most languages. It is recommended to use [the official Python library](https://pypi.org/project/pmtiles/) for both reading and writing. Files can be previewed or checked using [the official online viewer.](https://pmtiles.io/) If using another language, note that the format is quite simple to read and write without a library. By far the most difficult part is generating the Z-curve indices, however [the reference implementation has easily protable code.](https://github.com/protomaps/PMTiles/blob/26c857ff404f76dc0363f849b3b63113a00a0805/python/pmtiles/pmtiles/tile.py#L19C1-L43C15)

_Note that the Latch spatial widget does not support PMTiles leaf directories._ Please let us know if this seriously limits the resolution of your images.

---

Generating the file is done step-by-step, starting at the top, lowest, 0th zoom level (least zoomed-in). And proceeding down, towards higher zoom levels until the entire resolution of the source image is represented. Each level consists of a grid of individual images, known as tiles, which can be placed side-by-side to create the full picture of the slide. At each level, the tiles cover a smaller physical area, decreasing by a factor of two in both length and width but they keep the same digital size in pixels. This means that each level of zoom increases the resolution (pixels per mm of physical space) by a factor of four.

For example, consider a 10mm × 10mm microscope slide. Say the microscope produces a photograph, where each cell (assume 50μm = 0.05mm size) is covered by at least 200px. This means the entire image is 10mm × 100 pixels / 0.05mm = 40,000 pixels on each side.

This is too large to display without loss of detail on laptop screens, which we can safely but arbitrarily assume are at least Full HD (1,920px × 1,080px), and take that as our viewport dimensions. To fit the slide on the screen, we have to shrink by 40,000px / 1,080px (the smaller of the two dimensions) = 37.1x on each side. We will create new, more detailed zoom levels until this scaling factor is 1x or lower, meaning we are stretching the image rather than shrinking it, so no detail is lost.

The following is the structure of the resulting PMTiles image:

| Z-Level | Zoom Factor      | Tile Size                                 | Total Size                         |
| ------: | ---------------- | ----------------------------------------- | ---------------------------------- |
|       0 | 1x (1 tile)      | 10mm × 10mm (1,920px × 1,080px)           | 10mm × 10mm (1920px × 1080px)      |
|       1 | 2x (4 tiles)     | 5mm × 5mm (1,920px × 1,080px)             | 10mm × 10mm (3840px × 2160px)      |
|       2 | 4x (16 tiles)    | 2.5mm × 2.5mm (1,920px × 1,080px)         | 10mm × 10mm (7,680px × 4,320px)    |
|       3 | 8x (64 tiles)    | 1.25mm × 1.25mm (1,920px × 1,080px)       | 10mm × 10mm (15,360px × 8,640px)   |
|       4 | 16x (256 tiles)  | 0.625mm × 0.625mm (1,920px × 1,080px)     | 10mm × 10mm (30,720px × 17,280px)  |
|       5 | 32x (1024 tiles) | 0.3125mm × 0.3125mm (1,920px × 1,080px)   | 10mm × 10mm (61,440px × 34,560px)  |
|       6 | 64x (4096 tiles) | 0.15625mm × 0.15625mm (1,920px × 1,080px) | 10mm × 10mm (122,880px × 69,120px) |

Note that the max zoom level can be computed as `max( ceil(log_2(slide_width / viewport_width)), ceil(log_2(slide_height / viewport_height)) )`

Latch images always use the 1,920px × 1,080px viewport size.

### Geographical Coordinates

Since the PMTiles format was originally designed for map data, it requires fields that are nonsensical for flat images like slide photographs. These have been chosen such that existing viewers (like [pmtiles.io](https://pmtiles.io)) accurately display Latch PMTiles:

```python
mul = 10000000

min_pos = (-180 * mul, -85 * mul)
max_pos = (180 * mul, 85 * mul)
center_pos = (
    min_pos[0] + (max_pos[0] - min_pos[0]) // 2,
    min_pos[1] + (max_pos[1] - min_pos[1]) // 2,
)
```

### Resolution and Offset

When generating spatial background images, the most important step is to determine the appropriate resolution and translation of the image relative to the cell coordinates from the H5/AnnData file.

For the resolution, typically a dedicated metadata file will specify a `mm_per_px` field or similar (see appendix for examples). In some cases the original data comes in TIFF files which contain the resolution field in TIFF tags.

For the slide position, the same metadata will sometimes contain a pair of coordinates e.g. `X_mm` and `Y_mm`. Many formats, however, implicitly pre-align the top-left corner of the slide image with the physical origin i.e. coordinates (0, 0), and thus have no position field anywhere. Explicit position information is most common when, like CosMX, the source format uses a tiled representation stored in directories of loose image files rather than a container like OME TIFF.

---

The Latch widget will always place the top-left corner of the PMTiles image at the origin, so if the source format differs from this, you must offset all the spatial data correspondingly.

The resolution data is used throughout the process of converting spatial data.

### Rescaling/Contrast Correction

Many technologies produce images in a `Int16` format rather than RGB. The Latch widget does not currently support this directly and requires ahead-of-time conversion to 8-bit color. In most cases, linearly shrinking the entire `Int16` range will produce images that are too dark and lose too much contrast to be usable. It is recommended to at least clamp the data to its 10th-90th percentile range before doing the mapping, though more advanced techniques are possible like non-linear mapping functions or [CLAHE.](https://en.wikipedia.org/wiki/Adaptive_histogram_equalization) Multiple images in one dataset should typically be normalized to the same range (take the percentile statistics across all images rather than each one individually) to avoid misleading the viewer about their relative brightness.

Note that if you intend to visually compare different datasets you should also normalize them all to the same data range.

Our intent is to eventually support configurable contrast correction and color mapping in the widget itself. If your data would seriously benefit from this feature, let us know.

### Conversion Process Outline

1. Determine what format the source images are stored in. Common options are directories full of TIFF files, OME TIFF, Zarr arrays.

2. Find the resolution, dimensions, and position of each source image. Verify whether the resolution is the same for all images. Record the dimensions and position of each image in both physical units (mm) and in digital units (px), using the resolution to translate.

3. Calculate the overall dimensions and position of the slide. Typically this is computed by taking the minimum and maximum position across all source images. Note the offset between the top-left corner of the slide and the origin. Compensate for this offset in all output data.

4. Store 1,920px × 1,080px as the viewport size. Alternatively, choose an arbitrary viewport size, somewhere between 0.5x to 1.5x the size of the surface on which the image will be displayed.

5. Determine the maximum required zoom level: `max_z = max( ceil(log_2(slide_width / viewport_width)), ceil(log_2(slide_height / viewport_height)) )`

6. Clamp maximum zoom to 5, since levels higher than that require PMTiles leaf directory support and will not work with the widget.

7. Iterating through each zoom level, collect all relevant source image files by calculating whether any part of the source image overlaps the tile. If you only have one source image (e.g. using OME TIFF), this step is trivial.

   i. Resize each source image to the current zoom level
   ii. Composite each output tile by cropping and positioning the scaled source images.
   iii. Encode each image as WebP. Other formats may work depending on browser support. The widget renders the result using `<img>` tags.

8. Write a PMTiles file containing all the composited tiles. The widget requires that the pixel size of the top zoom level is specified in the PMTiles metadata section in the `slide_px` key as `{width: number; height: number}`. This can be calculated as follows:

   ```python
   # note: all sizes in pixels
   scale0 = min(
       viewport_w / slide_w,
       viewport_h / slide_h,
   )

   slide_vw = math.ceil(slide_w * scale0)
   slide_vh = math.ceil(slide_h * scale0)
   meta = {"slide_px": {"width": slide_vw, "height": slide_vh}}
   ```

### Pseudocode

```python
from pathlib import Path

from pmtiles.tile import zxy_to_tileid, tileid_to_zxy, TileType, Compression
from pmtiles.writer import Writer

viewport_w = 1920
viewport_h = 1080

with Path("Maximum Intensity Projection.pmtiles".open("wb") as f:
    w = Writer(f)

    metadata = load_metadata()
    sources = load_source_data(mm_per_px=metadata.mm_per_px)

    slide_x = 0
    slide_y = 0
    slide_x1 = 0
    slide_y1 = 0
    for img in sources:
        slide_x = min(img.x)
        slide_y = min(img.y)

        slide_x1 = max(img.x)
        slide_y1 = max(img.y)

    slide_w = slide_x1 - slide_x
    slide_h = slide_y1 - slide_y

    # to fit the slide in the viewport, how much do we need to scale by?
    slide_vscale = min(
        viewport_w / slide_w,
        viewport_h / slide_h,
    )

    # slide-shaped rectangle (same aspect-ratio) that fits in the viewport
    slide_vw = ceil(slide_w * slide_vscale)
    slide_vh = ceil(slide_h * slide_vscale)

    max_z = ceil(
      max(0, log2(slide_w / viewport_w), log2(slide_h / viewport_h))
    )
    max_z = min(max_z, 5)

    for z in range(max_z):
        n_tiles = 2 ** z

        tile_w = slide_w // n_tiles
        tile_h = slide_h // n_tiles

        # remember that all tiles are the same size, determined by the viewport
        #
        # to fit the tile in the viewport, how much do we need to scale by?
        # this is always <= 1.0 and increases at each tile level
        tile_vscale = min(
          viewport_w / tile_w,
          viewport_h / tile_h
        )

        # tile-shaped rectangle (same aspect-ratio) that fits in the viewport
        tile_vw = ceil(tile_w * tile_vscale)
        tile_vh = ceil(tile_h * tile_vscale)

        tiles = []
        for x in range(n_tiles):
            tile_x = slide_x + tile_w * x

            for y in range(n_tiles):
                tile_y = slide_y + tile_h * y

                img = composite(
                  [
                    {
                      "image": src,
                      # position in the tile viewport
                      # this can intentionally be negative or overflow the viewport
                      # because sometimes only a little piece of the source image
                      # is in the tile
                      "x": (src.x - tile_x) * tile_vscale,
                      "y": (src.y - tile_y) * tile_vscale,
                      #
                      "w": src.w * tile_vscale,
                      "h": src.h * tile_vscale
                    } for src in sources
                    if src.overlaps(
                      x=tile_x, y=tile_y,
                      w=tile_w, y=tile_w
                    )
                  ]
                )

                tiles.append(
                  Tile(
                    z=z,
                    x=x,
                    y=y,
                    img=img
                  )
                )

        # tiles must be written in increasing order of tileid
        tiles.sort(key=lambda x: zxy_to_tileid(tile.z, tile.x, tile.y))
        for t in tiles:
            w.write_tile(t.img)

    mul = 10000000

    min_pos = (-180 * mul, -85 * mul)
    max_pos = (180 * mul, 85 * mul)
    center_pos = (
        min_pos[0] + (max_pos[0] - min_pos[0]) // 2,
        min_pos[1] + (max_pos[1] - min_pos[1]) // 2,
    )

    writer.finalize(
        {
            "tile_type": TileType.WEBP,
            "tile_compression": Compression.NONE,
            "min_zoom": 0,
            "max_zoom": max_z,
            "min_lon_e7": int(min_pos[0]),
            "min_lat_e7": int(min_pos[1]),
            "max_lon_e7": int(max_pos[0]),
            "max_lat_e7": int(max_pos[1]),
            "center_zoom": 0,
            "center_lon_e7": center_pos[0],
            "center_lat_e7": center_pos[1],
        },
        {"slide_px": {"width": slide_vw, "height": slide_vh}},
    )

```

## Cell Boundaries

The widget looks for an optional `cell_boundaries.duckdb` file, which comes from [the DuckDB database.](https://duckdb.org/)

The schema is as follows:

```sql
create table cell_boundaries as (
  EntityID text not null,
  Geometry geometry not null
);
```

The `EntityID` is the cell ID from the AnnData/H5 file. The `Geometry` column contains [spatial data,](https://duckdb.org/docs/stable/core_extensions/spatial/functions) usually created with [`ST_MakeLine`.](https://duckdb.org/docs/stable/core_extensions/spatial/functions#st_makeline) Make sure you load the `spatial` extension using

```sql
install spatial; load spatial;
```

The coordinates should be in pixels and the axes should align with the slide's axes in scale and offset.

## Transcripts

The widget looks for an optional `transcripts.duckdb` file, which comes from [the DuckDB database.](https://duckdb.org/)

The schema is as follows:

```sql
create table final_transcripts as (
  cell_id text not null,
  feature_name text not null,
  global_x float not null,
  global_y float not null
);
```

The `cell_id` is the cell ID from the AnnData/H5 file. `feature_name` is the name of the transcript. `global_x` and `global_y` are the coordinates of the transcript.

The coordinates should be in pixels and the axes should align with the slide's axes in scale and offset.
