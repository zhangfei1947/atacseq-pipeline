"""Exploratory distal-peak to promoter-peak correlations across bulk samples."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


counts = pd.read_csv(snakemake.input.counts, sep="\t")
manifest = pd.read_csv(snakemake.input.manifest, sep="\t")
annotation = pd.read_csv(snakemake.input.annotation, sep="\t")
enabled = bool(snakemake.params.enabled)
minimum = int(snakemake.params.minimum_samples)
status = {"enabled": enabled, "sample_count": len(manifest), "minimum_samples": minimum}
Path(snakemake.output.links).parent.mkdir(parents=True, exist_ok=True)

if not enabled or len(manifest) < minimum:
    pd.DataFrame(columns=["distal_peak_id", "promoter_peak_id", "gene_id", "distance", "spearman_rho", "pvalue", "padj"]).to_csv(
        snakemake.output.links, sep="\t", index=False
    )
    status["status"] = "disabled" if not enabled else "ineligible"
    status["reason"] = "module_disabled" if not enabled else "fewer_than_20_independent_samples"
    Path(snakemake.output.status).write_text(json.dumps(status, indent=2) + "\n")
else:
    sample_ids = manifest.sample_id.tolist()
    matrix = counts.set_index("peak_id")[sample_ids].astype(float)
    lib = matrix.sum(axis=0)
    logcpm = np.log2(matrix.divide(lib, axis=1) * 1e6 + 1)
    annotation["peak_id"] = annotation["peak_id"].astype(str)
    promoters = annotation[annotation["annotation"].astype(str).str.startswith("Promoter")].copy()
    records = []
    max_distance = int(snakemake.params.max_distance)
    coords = counts.set_index("peak_id")[["chrom", "start_1based", "end"]]
    for _, promoter in promoters.iterrows():
        promoter_id = promoter["peak_id"]
        if promoter_id not in logcpm.index:
            continue
        chrom = coords.loc[promoter_id, "chrom"]
        center = (int(coords.loc[promoter_id, "start_1based"]) - 1 + int(coords.loc[promoter_id, "end"])) // 2
        local = coords[(coords.chrom == chrom)].copy()
        local["center"] = ((local.start_1based - 1 + local.end) // 2).astype(int)
        local = local[(local.center - center).abs() <= max_distance]
        for distal_id, row in local.iterrows():
            if distal_id == promoter_id:
                continue
            rho, pvalue = spearmanr(logcpm.loc[distal_id], logcpm.loc[promoter_id])
            records.append([distal_id, promoter_id, promoter.get("geneId", ""), int(row.center - center), rho, pvalue])
    links = pd.DataFrame(records, columns=["distal_peak_id", "promoter_peak_id", "gene_id", "distance", "spearman_rho", "pvalue"])
    if len(links):
        order = np.argsort(links.pvalue.values)
        ranked = links.pvalue.values[order] * len(links) / np.arange(1, len(links) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        links["padj"] = np.nan
        links.loc[order, "padj"] = np.minimum(ranked, 1.0)
        links = links[(links.padj <= 0.05) & (links.spearman_rho.abs() >= 0.6)]
    else:
        links["padj"] = []
    links.to_csv(snakemake.output.links, sep="\t", index=False)
    status.update({"status": "exploratory_complete", "reported_links": len(links), "warning": "ATAC-only correlation is not a validated enhancer-gene interaction."})
    Path(snakemake.output.status).write_text(json.dumps(status, indent=2) + "\n")

