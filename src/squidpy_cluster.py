#!/usr/bin/env python3
"""
Normalise, integrate across capture areas and cluster the combined Visium HD object.

Batch correction is Harmony on two batch variables, ``sample_id`` and ``patient_id``:
one batch per biological sample (22 of them against the experiment's 24 TMA cores --
the two ``A1_TMA6`` cores of each condition are technical replicates of the same
patient 12 W0 sample and share a ``sample_id``) and, on top, one batch per patient.
``sample_id`` is built in ``src/squidpy_make_adata.py`` as
``<patient_id>_<condition>_W<timepoint_weeks>``. Because every patient occupies exactly
one capture area, the ``patient_id`` correction is also a capture-area correction:
it removes the between-patient/between-section component of variation on top of the
between-sample component.

The cost, which is not small: every ``sample_id`` batch is a single condition x
timepoint, so ``condition`` and ``timepoint_weeks`` are perfectly nested inside the
batch keys, and adding ``patient_id`` flattens the remaining between-patient structure.
Harmony aligns cluster composition across batches, so it attenuates genuine
normal-vs-pigmented and W0-vs-W16 composition differences along with the technical
ones. Three things make that visible rather than assumed away:

* Harmony only rewrites the PCA embedding. ``rank_genes_groups`` here and the pyDESeq2
  pseudobulk path in ``src/squidpy_de.py`` read ``.raw`` / ``layers['counts']``, so
  expression-level differential expression is untouched. Only clustering, UMAP and any
  differential-abundance reading are affected.
* ``X_pca`` (uncorrected) is kept next to ``X_pca_harmony``.
* ``--theta`` softens the correction, and ``--crosstab-dir`` writes the cluster-by-batch
  and cluster-by-condition tables. A cluster_by_condition table that is ~50/50 in every
  cluster means the contrast has been flattened; lower ``--theta`` or compare against
  ``--no-harmony``.

``--no-harmony`` clusters on uncorrected PCA, for that comparison or if harmonypy is
unavailable.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

#: obs columns worth a UMAP panel, when present.
UMAP_KEYS = ("sample_id", "run_id", "condition", "timepoint_weeks", "patient",
             "core_id", "tma_position", "broad_cell_type", "coarse_cell_type",
             "total_counts")

#: annotations the cluster-composition crosstabs are written for, when present
CROSSTAB_KEYS = ("condition", "timepoint_weeks", "run_id", "core_id", "patient")


def preprocess(adata, n_top_genes, batch_key):
    """Normalise, log, select HVGs and run PCA. Raw counts stay in layers['counts'].

    ``batch_key`` is the HVG batch, which is None by default and is *not* the Harmony
    batch: per-batch selection lets a near-empty section contribute noise genes to the
    2000, and several sections here run at a median of a few counts per bin.
    """
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata

    # seurat flavor on log-normalised data; batch-aware only when asked for
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes,
        batch_key=batch_key if batch_key and batch_key in adata.obs else None,
    )
    print(f"HVG selection batch_key: {batch_key or 'none (whole object)'}")
    print(f"{int(adata.var['highly_variable'].sum())} highly variable genes")

    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=50, svd_solver="arpack")
    return adata


def integrate(adata, batch_keys, theta=2.0, max_iter_harmony=20):
    """Harmony on one or more batch keys; returns the obsm key to build neighbors from.

    ``theta`` is the diversity-clustering penalty: higher corrects harder. Every
    ``sample_id`` batch here is a single condition x timepoint, so the default 2.0 is
    already enough to pull conditions together -- lower it if cluster_by_condition.csv
    comes back uniform. A single ``theta`` is applied to every batch key.
    """
    import harmonypy

    missing = [k for k in batch_keys if k not in adata.obs]
    if missing:
        print(f"Batch key(s) {missing} not in obs; skipping Harmony")
        return "X_pca"
    counts = {k: adata.obs[k].nunique() for k in batch_keys}
    if all(n < 2 for n in counts.values()):
        print(f"No batch key with 2+ levels among {batch_keys}; skipping Harmony")
        return "X_pca"
    print(f"Harmony on {batch_keys} "
          f"({', '.join(f'{k}: {n} batches' for k, n in counts.items())}), "
          f"theta={theta}")
    # harmonypy directly rather than sce.pp.harmony_integrate: scanpy 1.11 transposes
    # Z_corr unconditionally, which is right for harmonypy 0.x (PCs x cells) and wrong
    # for 2.x (cells x PCs)
    pcs = adata.obsm["X_pca"]
    result = harmonypy.run_harmony(pcs, adata.obs, batch_keys, theta=theta,
                                   max_iter_harmony=max_iter_harmony)
    corrected = np.asarray(result.Z_corr)
    if corrected.shape != pcs.shape:
        corrected = corrected.T
    adata.obsm["X_pca_harmony"] = corrected
    return "X_pca_harmony"


def report_graph_components(adata, crosstab_dir):
    """Connected components of the kNN graph, which bound the cluster count from below.

    Leiden returns at least one cluster per component, so a fragmented graph inflates the
    cluster count no matter the resolution. Near-empty bins are what fragment it: a run
    with 32k bins under 10 counts produced 110 components and 273 clusters, 213 of them
    holding under 100 bins and each drawn from a single sample.
    """
    from scipy.sparse.csgraph import connected_components

    n_components, labels = connected_components(
        adata.obsp["connectivities"], directed=False
    )
    sizes = np.bincount(labels)
    small = sizes < 100
    print(f"kNN graph: {n_components} connected component(s); largest holds "
          f"{sizes.max()} of {adata.n_obs} bins ({sizes.max() / adata.n_obs:.3f})")
    if small.sum():
        print(f"  {int(small.sum())} component(s) under 100 bins, holding "
              f"{int(sizes[small].sum())} bins -- these force extra clusters")
    if crosstab_dir:
        os.makedirs(crosstab_dir, exist_ok=True)
        path = os.path.join(crosstab_dir, "graph_components.csv")
        pd.DataFrame({"component": np.arange(n_components), "n_bins": sizes}).to_csv(
            path, index=False)
        print(f"Wrote {path}")
    return n_components


def write_crosstabs(adata, cluster_key, batch_keys, crosstab_dir):
    """Cluster composition by batch and by each biological annotation.

    Two readings, in opposite directions. By each batch key: a cluster that is the
    property of a single sample means Harmony under-corrected. By ``condition``: a table
    that is uniform across every cluster means it over-corrected and the contrast is
    gone from the embedding. Twenty-two batches do not fit in a log line, so these go to
    CSV; only the batch tables are printed, and only their extremes.
    """
    keys = list(batch_keys) + [k for k in CROSSTAB_KEYS
                               if k in adata.obs and k not in batch_keys]
    for key in keys:
        crosstab = pd.crosstab(adata.obs[cluster_key], adata.obs[key],
                               normalize="index")
        if crosstab_dir:
            os.makedirs(crosstab_dir, exist_ok=True)
            path = os.path.join(crosstab_dir, f"cluster_by_{key}.csv")
            crosstab.round(5).to_csv(path)
            print(f"Wrote {path}")
    for batch_key in batch_keys:
        batch_table = pd.crosstab(adata.obs[cluster_key], adata.obs[batch_key],
                                  normalize="index")
        worst = batch_table.max(axis=1).sort_values(ascending=False)
        print(f"\nLargest single-{batch_key} share per cluster (1.0 = unintegrated):")
        print(worst.round(3).head(10).to_string())


def plot_umaps(adata, cluster_key, figure_dir):
    """One UMAP per annotation, so batch effects are visible next to the clusters."""
    os.makedirs(figure_dir, exist_ok=True)
    keys = [cluster_key] + [k for k in UMAP_KEYS if k in adata.obs]
    for key in keys:
        figure, axis = plt.subplots(figsize=(6, 5))
        sc.pl.umap(adata, color=key, ax=axis, show=False, frameon=False, size=2)
        path = os.path.join(figure_dir, f"umap_{key}.pdf")
        figure.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(figure)
        print(f"Wrote {path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="combined h5ad")
    parser.add_argument("--output", required=True, help="clustered h5ad to write")
    parser.add_argument("--figure-dir", required=True, help="directory for UMAP panels")
    parser.add_argument("--markers-csv", help="rank_genes_groups table to write")
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--resolution", type=float, default=1.0)
    parser.add_argument("--batch-key", nargs="+", default=["sample_id"],
                        help="obs column(s) Harmony integrates on, space-separated")
    parser.add_argument("--hvg-batch-key", default=None,
                        help="batch key for highly_variable_genes; omit to select over "
                             "the whole object, which keeps near-empty sections from "
                             "contributing noise genes")
    parser.add_argument("--theta", type=float, default=2.0,
                        help="Harmony diversity penalty; lower corrects more gently")
    parser.add_argument("--max-iter-harmony", type=int, default=20)
    parser.add_argument("--crosstab-dir",
                        help="directory for the cluster-composition crosstabs")
    parser.add_argument("--cluster-key", default="leiden")
    parser.add_argument("--no-harmony", dest="use_harmony", action="store_false",
                        help="cluster on uncorrected PCA")
    parser.add_argument("--subsample", type=int,
                        help="cluster on at most this many bins (for the 8 um object)")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    print(f"Loaded {adata.n_obs} bins x {adata.n_vars} genes")

    if args.subsample and adata.n_obs > args.subsample:
        sc.pp.subsample(adata, n_obs=args.subsample, random_state=args.seed)
        print(f"Subsampled to {adata.n_obs} bins")

    adata = preprocess(adata, args.n_top_genes, args.hvg_batch_key)

    use_rep = "X_pca"
    if args.use_harmony:
        use_rep = integrate(adata, args.batch_key, args.theta, args.max_iter_harmony)

    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs,
                    use_rep=use_rep, random_state=args.seed)
    report_graph_components(adata, args.crosstab_dir)
    sc.tl.umap(adata, random_state=args.seed)
    sc.tl.leiden(adata, resolution=args.resolution, key_added=args.cluster_key,
                 flavor="igraph", n_iterations=2, directed=False,
                 random_state=args.seed)
    counts = adata.obs[args.cluster_key].value_counts().sort_index()
    print(f"{len(counts)} clusters: {counts.to_dict()}")

    sc.tl.rank_genes_groups(adata, args.cluster_key, method="wilcoxon")
    if args.markers_csv:
        markers = sc.get.rank_genes_groups_df(adata, group=None)
        os.makedirs(os.path.dirname(args.markers_csv) or ".", exist_ok=True)
        markers.to_csv(args.markers_csv, index=False)
        print(f"Wrote {args.markers_csv}")

    plot_umaps(adata, args.cluster_key, args.figure_dir)

    write_crosstabs(adata, args.cluster_key, args.batch_key, args.crosstab_dir)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    adata.write_h5ad(args.output, compression="gzip")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
