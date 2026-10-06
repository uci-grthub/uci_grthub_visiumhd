#!/usr/bin/env python3
"""
Extract H&E image features for the bins of one capture area.

Space Ranger registers the high-resolution microscope image and reports bin positions as
``pxl_row_in_fullres`` / ``pxl_col_in_fullres`` in that image's own coordinate frame, so
``obsm['spatial']`` indexes straight into the TIFF.

That correspondence is verified before any feature is computed, from the scale factors
rather than by eye: ``tissue_hires_image.png`` is the full-resolution frame scaled by
``tissue_hires_scalef``, so dividing its dimensions by that factor recovers the size the
microscope image must have. If it disagrees with the TIFF on disk, the run was registered
against a different image and the features would be meaningless.

Separately, some bins genuinely fall outside the image: the capture area is 6.5 mm across
(~13,250 full-resolution pixels) while the microscope images here are smaller than that,
and Space Ranger's tissue detection over-called on several capture areas, flagging bins
well beyond the imaged region -- A1_TMA4 has bins at x = -2056. Those bins are dropped
with a count, not treated as an error. Running ``src/squidpy_make_adata.py`` with a
correct core annotation removes most of them beforehand.

These microscope images are pyramidal TIFFs -- ``TMA 4.tif`` holds a 8337x9826 full
resolution page and a 2779x3275 downsample. Handed the path, ``sq.im.ImageContainer``
stacks the two levels as a Z-dimension and then rejects the library id, so the
full-resolution page is read explicitly instead.

Features are summary statistics and Haralick texture per bin, over a window ``--spot-scale``
times the bin diameter.
"""

import argparse
import os
import re

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import tifffile


def load_image_container(path, level=0):
    """Read one pyramid level of a TIFF into an ImageContainer.

    Passing the path straight to ImageContainer stacks every page as a Z-dimension,
    because these are pyramidal TIFFs; level 0 is the full-resolution page.
    """
    with tifffile.TiffFile(path) as handle:
        n_pages = len(handle.pages)
        shapes = [tuple(page.shape) for page in handle.pages]
    print(f"{os.path.basename(path)}: {n_pages} page(s) {shapes}; using level {level}")
    image = tifffile.imread(path, key=level)
    return sq.im.ImageContainer(image)


def check_alignment(adata, container, run_id, tolerance=0.02):
    """Verify the image is the one Space Ranger registered, using the scale factors.

    ``tissue_hires_image.png`` is the full-resolution frame scaled by
    ``tissue_hires_scalef``, so its dimensions divided by that factor give the size the
    microscope image must have. This is exact, unlike comparing against where the bins
    happen to land.
    """
    height, width = container.shape
    library = adata.uns.get("spatial", {}).get(run_id, {})
    factors = library.get("scalefactors", {})
    hires = library.get("images", {}).get("hires")
    scalef = factors.get("tissue_hires_scalef")

    if hires is None or not scalef:
        print("Warning: no hires image or tissue_hires_scalef in uns['spatial']; "
              "cannot verify that this is the registered image")
        return

    hires_height, hires_width = np.asarray(hires).shape[:2]
    expected_width = hires_width / scalef
    expected_height = hires_height / scalef
    print(f"image {width} x {height} px; scale factors imply "
          f"{expected_width:.0f} x {expected_height:.0f}")

    if (abs(width - expected_width) > expected_width * tolerance
            or abs(height - expected_height) > expected_height * tolerance):
        raise SystemExit(
            f"{run_id}: the image is {width} x {height} but the Space Ranger scale "
            f"factors imply {expected_width:.0f} x {expected_height:.0f}. This is not "
            "the microscope image this capture area was registered against, and the "
            "features would come from the wrong frame."
        )


def drop_offimage_bins(adata, container, spot_scale, spot_diameter):
    """Remove bins whose crop window leaves the image.

    The 6.5 mm capture area is wider than these microscope images, and tissue detection
    over-called on several capture areas, so a handful of bins sit outside the picture
    entirely. Cropping there yields garbage rather than an error, so they go.
    """
    height, width = container.shape
    margin = 0.5 * spot_scale * spot_diameter
    xs, ys = adata.obsm["spatial"][:, 0], adata.obsm["spatial"][:, 1]
    inside = ((xs - margin >= 0) & (xs + margin < width)
              & (ys - margin >= 0) & (ys + margin < height))
    dropped = int((~inside).sum())
    if dropped:
        print(f"Dropping {dropped} of {adata.n_obs} bins whose {margin:.0f} px crop "
              "window falls outside the image")
    if not inside.any():
        raise SystemExit(
            "every bin's crop window falls outside the image; check that this capture "
            "area was registered against this microscope image"
        )
    return adata[inside].copy()


def resolve_image(run_id, image, image_dir, pattern, regex):
    """Explicit --image wins; otherwise derive the TIFF name from the run id."""
    if image:
        return image
    if not (image_dir and pattern and regex):
        raise SystemExit("pass --image, or all of --image-dir, --image-pattern and "
                         "--image-regex")
    match = re.search(regex, run_id)
    if not match:
        raise SystemExit(f"--image-regex {regex!r} does not match run id {run_id!r}")
    path = os.path.join(image_dir, pattern.format(token=match.group(1), sample=run_id))
    if not os.path.exists(path):
        raise SystemExit(f"no microscope image at {path} for {run_id}")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="clustered h5ad")
    parser.add_argument("--run-id", required=True, help="capture area to featurise")
    parser.add_argument("--output-csv", required=True, help="feature table to write")
    parser.add_argument("--image", help="microscope TIFF for this capture area")
    parser.add_argument("--image-dir", help="directory of microscope TIFFs")
    parser.add_argument("--image-pattern", default="TMA {token}.tif",
                        help="filename pattern; {token} comes from --image-regex")
    parser.add_argument("--image-regex", default=r"TMA\s*(\d+)$",
                        help="captures {token} from the run id")
    parser.add_argument("--image-level", type=int, default=0,
                        help="pyramid level to read; 0 is full resolution")
    parser.add_argument("--library-key", default="run_id")
    parser.add_argument("--features", nargs="+", default=["summary", "texture"])
    parser.add_argument("--scale", type=float, default=0.5,
                        help="downscale factor applied to the crop before featurising")
    parser.add_argument("--spot-scale", type=float, default=2.0,
                        help="crop size as a multiple of the bin diameter")
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    if args.library_key not in adata.obs:
        raise SystemExit(f"{args.input} has no {args.library_key} column")
    subset = adata[adata.obs[args.library_key] == args.run_id].copy()
    if subset.n_obs == 0:
        raise SystemExit(
            f"{args.run_id} is not present in {args.input}; available: "
            f"{sorted(adata.obs[args.library_key].unique())}"
        )
    del adata
    print(f"{args.run_id}: {subset.n_obs} bins")

    path = resolve_image(args.run_id, args.image, args.image_dir,
                         args.image_pattern, args.image_regex)
    print(f"Reading {path}")
    container = load_image_container(path, args.image_level)
    check_alignment(subset, container, args.run_id)

    factors = subset.uns.get("spatial", {}).get(args.run_id, {}).get("scalefactors", {})
    spot_diameter = factors.get("spot_diameter_fullres", 32.0)
    subset = drop_offimage_bins(subset, container, args.spot_scale, spot_diameter)

    sq.im.calculate_image_features(
        subset, container, features=args.features, key_added="image_features",
        spot_scale=args.spot_scale, scale=args.scale, n_jobs=args.n_jobs,
        show_progress_bar=False,
        # the subset still carries every capture area in uns['spatial']
        library_id=args.run_id,
    )

    # calculate_image_features stores a DataFrame already indexed by barcode
    frame = pd.DataFrame(subset.obsm["image_features"]).copy()
    frame.index = subset.obs_names
    frame.index.name = "barcode"
    frame.insert(0, "run_id", args.run_id)
    os.makedirs(os.path.dirname(args.output_csv) or ".", exist_ok=True)
    frame.to_csv(args.output_csv)
    print(f"Wrote {frame.shape[0]} bins x {frame.shape[1] - 1} features to "
          f"{args.output_csv}")


if __name__ == "__main__":
    main()
