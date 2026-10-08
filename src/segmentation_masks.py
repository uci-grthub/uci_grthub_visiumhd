#!/usr/bin/env python3
"""
Rasterise Space Ranger's cell and nucleus segmentation polygons into label masks.

Space Ranger (>= 4.0, with an H&E image) writes the segmentation only as GeoJSON
polygons in ``outs/segmented_outputs/``; the raster it segmented on is not kept. This
draws each polygon back onto the full-resolution H&E grid: pixel value = ``cell_id``
(uint32), 0 = background. A cell and its nucleus share the same ``cell_id``, which is
the integer in the matrix barcode ``cellid_000000001-1``.

The grid is the ``--image`` passed to ``spaceranger count`` (the microscope TIFF), in
whose full-resolution pixels the polygons are given. Only its shape is read. A pixel
is labelled when its centre lies inside the polygon (``skimage.draw.polygon``), so
cells that are adjacent in the GeoJSON stay adjacent in the mask; where two polygons
claim a pixel, the later ``cell_id`` wins and the overlap is counted in the log.

Masks are tiled, zlib-compressed BigTIFFs with the pixel size from
``microns_per_pixel`` in ``scalefactors_json.json`` stored as the resolution.
"""

import argparse
import json
import os

import numpy as np
import tifffile
from skimage.draw import polygon as draw_polygon


def image_shape(path):
    """(height, width) of the full-resolution page, without reading pixels."""
    with tifffile.TiffFile(path) as tif:
        shape = tif.series[0].shape
        axes = tif.series[0].axes
    return shape[axes.index("Y")], shape[axes.index("X")]


def rasterise(geojson_path, shape):
    """Label mask from one segmentation GeoJSON, plus per-mask counts for the log."""
    with open(geojson_path) as handle:
        features = json.load(handle)["features"]

    mask = np.zeros(shape, dtype=np.uint32)
    empty, overlapped = 0, 0
    for feature in features:
        cell_id = int(feature["properties"]["cell_id"])
        geometry = feature["geometry"]
        if geometry["type"] == "Polygon":
            rings = [geometry["coordinates"][0]]
        elif geometry["type"] == "MultiPolygon":
            rings = [part[0] for part in geometry["coordinates"]]
        else:
            raise ValueError(f"cell_id {cell_id}: unsupported geometry "
                             f"{geometry['type']}")
        drawn = 0
        for ring in rings:
            xy = np.asarray(ring, dtype=float)
            # GeoJSON is (x, y) = (column, row); skimage wants rows first. Polygon
            # holes (later rings) are not used by Space Ranger and are ignored.
            rows, cols = draw_polygon(xy[:, 1], xy[:, 0], shape=shape)
            overlapped += int(np.count_nonzero(mask[rows, cols]))
            mask[rows, cols] = cell_id
            drawn += len(rows)
        if drawn == 0:
            empty += 1
    return mask, len(features), empty, overlapped


def write_mask(mask, path, microns_per_pixel):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    pixels_per_cm = 1e4 / microns_per_pixel
    tifffile.imwrite(path, mask, bigtiff=True, tile=(512, 512), compression="zlib",
                     resolution=(pixels_per_cm, pixels_per_cm), resolutionunit="CENTIMETER",
                     metadata=None)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--image", required=True,
                        help="full-resolution H&E TIFF passed to spaceranger --image")
    parser.add_argument("--scalefactors", required=True,
                        help="segmented_outputs/spatial/scalefactors_json.json")
    parser.add_argument("--cell-geojson", required=True)
    parser.add_argument("--nucleus-geojson", required=True)
    parser.add_argument("--cell-mask", required=True, help="cell label TIFF to write")
    parser.add_argument("--nucleus-mask", required=True, help="nucleus label TIFF to write")
    args = parser.parse_args()

    shape = image_shape(args.image)
    with open(args.scalefactors) as handle:
        microns_per_pixel = float(json.load(handle)["microns_per_pixel"])
    print(f"Grid {shape[0]} x {shape[1]} px from {args.image}, "
          f"{microns_per_pixel:.4f} um per px")

    for name, source, target in (("cell", args.cell_geojson, args.cell_mask),
                                 ("nucleus", args.nucleus_geojson, args.nucleus_mask)):
        mask, n, empty, overlapped = rasterise(source, shape)
        labels = np.unique(mask)
        print(f"{name}: {n} polygons, {labels.size - (labels[0] == 0)} labels in the "
              f"mask, {empty} polygons too small to cover a pixel centre, "
              f"{overlapped} pixels claimed twice, "
              f"{np.count_nonzero(mask) / mask.size:.3%} of the image labelled")
        write_mask(mask, target, microns_per_pixel)
        print(f"Wrote {target}")


if __name__ == "__main__":
    main()
