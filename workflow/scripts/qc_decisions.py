"""Create non-destructive sample QC and module-eligibility decisions."""

import gzip
import json
import re
from pathlib import Path

import pandas as pd


def line_count_gzip(path):
    with gzip.open(path, "rt") as handle:
        return sum(1 for _ in handle)


def picard_duplication(path):
    lines = Path(path).read_text().splitlines()
    for idx, line in enumerate(lines):
        if line.startswith("LIBRARY\t") and idx + 1 < len(lines):
            header = line.split("\t")
            values = lines[idx + 1].split("\t")
            row = dict(zip(header, values))
            try:
                return float(row["PERCENT_DUPLICATION"])
            except (KeyError, ValueError):
                return None
    return None


def flagstat_rates(path):
    text = Path(path).read_text()
    total_match = re.search(r"^(\d+) \+ \d+ in total", text, re.M)
    mapped_match = re.search(r"^(\d+) \+ \d+ mapped", text, re.M)
    proper_match = re.search(r"^(\d+) \+ \d+ properly paired", text, re.M)
    total = int(total_match.group(1)) if total_match else 0
    return (
        int(mapped_match.group(1)) / total if total and mapped_match else None,
        int(proper_match.group(1)) / total if total and proper_match else None,
    )


def recursive_numbers(value, output=None):
    output = {} if output is None else output
    if isinstance(value, dict):
        for key, child in value.items():
            normalized = re.sub(r"[^a-z0-9]", "", str(key).lower())
            if isinstance(child, (int, float)):
                output[normalized] = float(child)
            recursive_numbers(child, output)
    elif isinstance(value, list):
        for child in value:
            recursive_numbers(child, output)
    return output


def first_metric(metrics, candidates):
    for candidate in candidates:
        key = re.sub(r"[^a-z0-9]", "", candidate.lower())
        if key in metrics:
            return metrics[key]
    return None


manifest = pd.read_csv(snakemake.input.manifest, sep="\t", dtype=str, keep_default_na=False)
sample_ids = list(snakemake.params.sample_ids)
thresholds = dict(snakemake.params.thresholds)

cor = pd.read_csv(snakemake.input.correlations, sep="\t", comment="#", index_col=0)
cor.index = [Path(str(x)).stem.replace(".analysis", "") for x in cor.index]
cor.columns = [Path(str(x)).stem.replace(".analysis", "") for x in cor.columns]

rows = []
for idx, sample in enumerate(sample_ids):
    reasons = []
    usable = line_count_gzip(snakemake.input.fragments[idx])
    frip_tab = pd.read_csv(snakemake.input.frip[idx], sep="\t")
    frip = float(frip_tab.frip.iloc[0])
    duplicate = picard_duplication(snakemake.input.picard[idx])
    mapped, proper = flagstat_rates(snakemake.input.flagstat[idx])
    ataqv = recursive_numbers(json.loads(Path(snakemake.input.ataqv[idx]).read_text()))
    tss = first_metric(ataqv, ["tss_enrichment", "tss_enrichment_score"])
    mito_reads = first_metric(ataqv, ["mitochondrial_reads", "mitochondrial_read_count"])
    total_reads = first_metric(ataqv, ["total_reads", "total_read_count"])
    mito = mito_reads / total_reads if mito_reads is not None and total_reads else None

    meta = manifest.loc[manifest.sample_id == sample].iloc[0]
    peers = manifest.loc[(manifest.replicate_group == meta.replicate_group) & (manifest.sample_id != sample), "sample_id"].tolist()
    peer_cor = [float(cor.loc[sample, peer]) for peer in peers if sample in cor.index and peer in cor.columns]
    replicate_cor = min(peer_cor) if peer_cor else None

    if usable < int(thresholds["usable_fragments_warn"]): reasons.append("low_usable_fragments")
    elif usable < int(thresholds["usable_fragments_good"]): reasons.append("usable_fragments_warning")
    if frip < float(thresholds["frip_warn"]): reasons.append("low_frip")
    elif frip < float(thresholds["frip_good"]): reasons.append("frip_warning")
    if mapped is not None and mapped < 0.60: reasons.append("low_alignment_rate")
    if proper is not None and proper < 0.60: reasons.append("low_proper_pair_rate")
    if duplicate is not None and duplicate > 0.50: reasons.append("high_duplicate_fraction")
    if mito is not None and mito > 0.40: reasons.append("high_mitochondrial_fraction")
    if replicate_cor is not None and replicate_cor < float(thresholds["replicate_correlation_warn"]): reasons.append("low_replicate_correlation")

    severe = {"low_usable_fragments", "low_frip", "low_alignment_rate", "low_proper_pair_rate", "high_duplicate_fraction", "high_mitochondrial_fraction", "low_replicate_correlation"}
    qc_status = "REVIEW" if severe.intersection(reasons) else ("WARN" if reasons else "PASS")
    technical = "VALID" if usable > 0 else "INVALID"
    core_ready = technical == "VALID" and qc_status != "REVIEW"
    rows.append({
        "analysis_id": snakemake.wildcards.analysis,
        "sample_id": sample,
        "technical_status": technical,
        "qc_status": qc_status,
        "analysis_include": meta.analysis_include,
        "reasons": ";".join(reasons),
        "usable_nuclear_fragments": usable,
        "alignment_rate": mapped,
        "proper_pair_rate": proper,
        "duplicate_fraction": duplicate,
        "mitochondrial_fraction": mito,
        "frip": frip,
        "tss_enrichment": tss,
        "replicate_correlation_min": replicate_cor,
        "peak_calling_eligible": technical == "VALID",
        "differential_eligible": core_ready,
        "motif_eligible": technical == "VALID",
        "footprinting_eligible": core_ready and usable >= int(thresholds["usable_fragments_warn"]) and frip >= float(thresholds["frip_warn"]),
        "recommended_action": "review_before_primary_analysis" if qc_status == "REVIEW" else "include",
    })

Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(snakemake.output[0], sep="\t", index=False, na_rep="NA")
