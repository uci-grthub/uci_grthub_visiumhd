#!/usr/bin/env python3
"""
Combine per-capture-area Space Ranger metrics_summary.csv files into one table.

spaceranger aggr does not accept Visium HD runs (Space Ranger 3.1.2 and 4.0
both refuse them at PARSE_CSV), so the per-run metrics are concatenated here
instead. Every capture area writes the same 38 columns, prefixed with the run
id taken from the output directory name.
"""

import argparse
import os

import pandas as pd


def run_id_from_path(path):
    """output/spaceranger/<run id>/outs/metrics_summary.csv -> <run id>"""
    return os.path.basename(os.path.dirname(os.path.dirname(path)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", nargs="+", required=True,
                        help="per-run metrics_summary.csv files")
    parser.add_argument("--output", required=True, help="combined CSV to write")
    args = parser.parse_args()

    frames = []
    for path in sorted(args.metrics):
        run_id = run_id_from_path(path)
        try:
            frame = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            print(f"Warning: empty metrics file, skipping: {path}")
            continue
        frame.insert(0, "run_id", run_id)
        frames.append(frame)

    if not frames:
        raise SystemExit("No metrics files could be read")

    combined = pd.concat(frames, ignore_index=True, sort=False)

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    combined.to_csv(args.output, index=False)
    print(f"Wrote {combined.shape[0]} capture areas x {combined.shape[1]} "
          f"columns to {args.output}")


if __name__ == "__main__":
    main()
