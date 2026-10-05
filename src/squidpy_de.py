#!/usr/bin/env python3
"""
Pseudobulk differential expression for one comparison with pyDESeq2.

Bins (or cells) are not independent replicates, so testing them against each other
inflates p-values by orders of magnitude. The unit of replication here is the biopsy
core: raw counts are summed per core, once over the whole core and
once per cluster, and pyDESeq2 fits a negative-binomial GLM across cores. The two TMA6
technical replicates share a ``core_index`` and are summed into one sample.

The comparison is a JSON object from ``squidpy.de.comparisons`` in config.yaml:

    {"factor": "sample_type", "test": "melasma", "reference": "normal",
     "subset": {"diagnosis": "melasma", "timepoint_weeks": 0}, "paired": true}

``subset`` restricts the cores first, then only cores at ``test`` or ``reference`` of
``factor`` are kept. A paired comparison fits ``~ patient_id + factor``: the patient term
absorbs patient and capture area (the same thing in this design), and patients seen at
only one level of the factor are dropped because they carry no within-patient contrast.
An unpaired comparison fits ``~ factor``. log2FoldChange is test over reference.

Counts come from combined.h5ad (raw counts, every gene), because clustered.h5ad keeps
only the highly variable genes; only the cluster labels are read from clustered.h5ad,
matched by obs name.

The clinical columns (patient_id, diagnosis, sample_type) are joined from
sample_metadata.csv by core_id here rather than read from obs, so adding them to the
metadata does not require rebuilding the AnnData objects.
"""

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

METADATA_COLUMNS = ["core_id", "core_index", "patient_id", "diagnosis", "sample_type",
                    "timepoint_weeks", "de_group"]
RESULT_COLUMNS = ["gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]


def attach_metadata(adata, sample_metadata):
    """Replace the per-core obs columns with the ones in sample_metadata, by core_id."""
    metadata = pd.read_csv(sample_metadata)
    missing = [c for c in METADATA_COLUMNS if c not in metadata.columns]
    if missing:
        raise SystemExit(f"{sample_metadata} lacks {missing}; regenerate it with "
                         "src/make_sample_metadata.py --clinical-xlsx")
    metadata = metadata[METADATA_COLUMNS].drop_duplicates("core_id").set_index("core_id")
    core_id = adata.obs["core_id"].astype(str)
    unknown = sorted(set(core_id) - set(metadata.index))
    if unknown:
        raise SystemExit(f"core_id(s) {unknown} are not in {sample_metadata}")
    for column in metadata.columns:
        # strings throughout: the formula treats every factor, timepoint included, as
        # categorical, and the contrast levels arrive from JSON as strings
        adata.obs[column] = core_id.map(metadata[column].astype(str)).to_numpy()
    return adata


def attach_clusters(adata, clustered_path, cluster_key):
    """Copy cluster labels from the clustered object; unmatched obs get none."""
    import anndata as ad

    labels = ad.read_h5ad(clustered_path, backed="r").obs[cluster_key]
    adata.obs[cluster_key] = labels.astype(str).reindex(adata.obs_names).to_numpy()
    missing = adata.obs[cluster_key].isna()
    print(f"{int((~missing).sum())} of {adata.n_obs} obs carry a {cluster_key} label")
    return adata


def pseudobulk(adata, keys, min_obs):
    """Sum counts over keys; returns (counts genes-in-columns frame, per-sample obs)."""
    layer = "counts" if "counts" in adata.layers else None
    if layer is None:
        print("no layers['counts']; summing X, which must hold raw counts")
    summed = sc.get.aggregate(adata, by=keys, func="sum", layer=layer)
    counts = pd.DataFrame(
        np.asarray(summed.layers["sum"]), index=summed.obs_names, columns=summed.var_names
    ).round().astype(int)

    n_obs = adata.obs.groupby(keys, observed=True).size().rename("n_obs").reset_index()
    per_sample = summed.obs[keys].astype(str).reset_index(names="pseudobulk")
    per_sample = per_sample.merge(n_obs.astype({k: str for k in keys}), on=keys, how="left")
    sample_columns = [c for c in METADATA_COLUMNS if c not in keys and c != "core_id"]
    first = adata.obs.groupby(keys, observed=True)[sample_columns].first().reset_index()
    per_sample = per_sample.merge(first.astype({k: str for k in keys}), on=keys, how="left")
    per_sample = per_sample.set_index("pseudobulk").loc[counts.index]
    per_sample["total_counts"] = counts.sum(axis=1)
    per_sample["kept"] = per_sample["n_obs"] >= min_obs
    return counts, per_sample


def select_samples(per_sample, comparison):
    """Apply subset, factor levels and pairing; returns the per-sample frame with 'used'."""
    factor = comparison["factor"]
    test, reference = str(comparison["test"]), str(comparison["reference"])
    used = per_sample["kept"].copy()
    for column, value in comparison.get("subset", {}).items():
        used &= per_sample[column].astype(str) == str(value)
    used &= per_sample[factor].isin([test, reference])
    if comparison.get("paired", False):
        levels = per_sample[used].groupby("patient_id")[factor].nunique()
        unpaired = levels[levels < 2].index
        if len(unpaired):
            print(f"  paired: dropping patient(s) {list(unpaired)} seen at one level only")
        used &= ~per_sample["patient_id"].isin(unpaired)
    per_sample = per_sample.copy()
    per_sample["used"] = used
    return per_sample


def run_deseq(counts, per_sample, comparison, min_samples, min_gene_counts, n_cpus):
    """Fit pyDESeq2 on the used samples; None when a group is too small."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    factor = comparison["factor"]
    test, reference = str(comparison["test"]), str(comparison["reference"])
    samples = per_sample[per_sample["used"]]
    sizes = samples[factor].value_counts()
    if sizes.get(test, 0) < min_samples or sizes.get(reference, 0) < min_samples:
        print(f"  skipped: {test} n={sizes.get(test, 0)}, {reference} n={sizes.get(reference, 0)}"
              f" (need {min_samples} each)")
        return None

    matrix = counts.loc[samples.index]
    matrix = matrix.loc[:, matrix.sum(axis=0) >= min_gene_counts]
    design = f"~ patient_id + {factor}" if comparison.get("paired", False) else f"~ {factor}"
    print(f"  {design}: {test} n={sizes[test]} vs {reference} n={sizes[reference]}, "
          f"{matrix.shape[1]} genes")

    metadata = samples[["patient_id", factor]].astype("category")
    dds = DeseqDataSet(counts=matrix, metadata=metadata, design=design,
                       inference=DefaultInference(n_cpus=n_cpus), quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=[factor, test, reference],
                       inference=DefaultInference(n_cpus=n_cpus), quiet=True)
    stats.summary()
    result = stats.results_df.rename_axis("gene").reset_index()
    result["n_test"] = sizes[test]
    result["n_reference"] = sizes[reference]
    return result.sort_values("pvalue")


def plot_volcano(result, title, path, alpha, n_labels=15):
    figure, axis = plt.subplots(figsize=(6, 5))
    if result is not None and len(result):
        frame = result.dropna(subset=["padj"])
        significant = frame["padj"] < alpha
        y = -np.log10(frame["pvalue"].clip(lower=1e-300))
        axis.scatter(frame["log2FoldChange"], y, s=4, c=np.where(significant, "C3", "0.6"),
                     linewidths=0)
        for _, row in frame[significant].head(n_labels).iterrows():
            axis.annotate(row["gene"], (row["log2FoldChange"], -np.log10(max(row["pvalue"], 1e-300))),
                          fontsize=6)
        title += f" ({int(significant.sum())} genes padj < {alpha})"
    else:
        axis.text(0.5, 0.5, "not tested", ha="center", transform=axis.transAxes)
    axis.axvline(0, color="0.8", linewidth=0.5)
    axis.set_xlabel("log2 fold change")
    axis.set_ylabel("-log10 p")
    axis.set_title(title, fontsize=9)
    figure.savefig(path, bbox_inches="tight")
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", required=True, help="combined h5ad (raw counts)")
    parser.add_argument("--clustered", required=True, help="clustered h5ad (labels)")
    parser.add_argument("--sample-metadata", required=True, help="metadata/sample_metadata.csv")
    parser.add_argument("--name", required=True, help="comparison name")
    parser.add_argument("--comparison", required=True, help="comparison as a JSON object")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample-key", default="core_index")
    parser.add_argument("--cluster-key", default="leiden")
    parser.add_argument("--min-obs-per-pseudobulk", type=int, default=50)
    parser.add_argument("--min-samples-per-group", type=int, default=2)
    parser.add_argument("--min-gene-counts", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()

    comparison = json.loads(args.comparison)
    os.makedirs(args.outdir, exist_ok=True)
    print(f"{args.name}: {comparison}")

    adata = sc.read_h5ad(args.counts)
    adata = attach_metadata(adata, args.sample_metadata)
    adata = attach_clusters(adata, args.clustered, args.cluster_key)

    print("whole core")
    counts, per_sample = pseudobulk(adata, [args.sample_key], args.min_obs_per_pseudobulk)
    per_sample = select_samples(per_sample, comparison)
    whole = run_deseq(counts, per_sample, comparison, args.min_samples_per_group,
                      args.min_gene_counts, args.n_jobs)
    if whole is None:
        raise SystemExit(f"{args.name}: too few cores for the whole-core test")
    whole.to_csv(os.path.join(args.outdir, "whole_core.csv"), index=False)
    plot_volcano(whole, f"{args.name}, whole core", os.path.join(args.outdir, "volcano.pdf"),
                 args.alpha)
    samples = [per_sample.assign(**{args.cluster_key: "whole_core"})]

    clustered_results = []
    keys = [args.sample_key, args.cluster_key]
    labelled = adata[adata.obs[args.cluster_key].notna()]
    counts, per_cluster = pseudobulk(labelled, keys, args.min_obs_per_pseudobulk)
    for cluster, block in per_cluster.groupby(args.cluster_key, sort=False):
        print(f"{args.cluster_key} {cluster}")
        # pairing is per cluster: a patient may have the cluster at one level only
        block = select_samples(block, comparison)
        result = run_deseq(counts, block, comparison, args.min_samples_per_group,
                           args.min_gene_counts, args.n_jobs)
        if result is not None:
            clustered_results.append(result.assign(**{args.cluster_key: cluster}))
        samples.append(block)

    if clustered_results:
        clustered = pd.concat(clustered_results, ignore_index=True)
    else:
        clustered = pd.DataFrame(columns=RESULT_COLUMNS + ["n_test", "n_reference",
                                                           args.cluster_key])
    clustered.to_csv(os.path.join(args.outdir, "per_cluster.csv"), index=False)
    pd.concat(samples).to_csv(os.path.join(args.outdir, "samples.csv"))

    print(f"whole core: {int((whole['padj'] < args.alpha).sum())} genes padj < {args.alpha}")
    if len(clustered):
        print((clustered[clustered["padj"] < args.alpha]
               .groupby(args.cluster_key).size().rename("n_padj_lt_alpha").to_string()))
    print(f"Wrote {args.outdir}")


if __name__ == "__main__":
    main()
