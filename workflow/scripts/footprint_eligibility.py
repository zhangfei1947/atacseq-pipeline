import json
from pathlib import Path

import pandas as pd


qc = pd.read_csv(snakemake.input.qc, sep="\t", dtype=str, keep_default_na=False)
manifest = pd.read_csv(snakemake.input.manifest, sep="\t", dtype=str, keep_default_na=False)
factor = snakemake.params.factor
numerator = str(snakemake.params.numerator)
denominator = str(snakemake.params.denominator)
minimum_sites = int(snakemake.params.minimum_sites)

groups = {
    numerator: manifest.loc[manifest[factor] == numerator, "sample_id"].tolist(),
    denominator: manifest.loc[manifest[factor] == denominator, "sample_id"].tolist(),
}
replicate_ok = all(len(samples) >= 2 for samples in groups.values())
selected = groups[numerator] + groups[denominator]
qc_selected = qc[qc.sample_id.isin(selected)]
qc_ok = len(qc_selected) == len(selected) and all(qc_selected.footprinting_eligible.str.lower() == "true")

try:
    fimo = pd.read_csv(snakemake.input.fimo, sep="\t", comment="#")
    motif_column = "motif_id" if "motif_id" in fimo.columns else fimo.columns[0]
    site_counts = fimo[motif_column].value_counts()
    eligible_motifs = site_counts[site_counts >= minimum_sites]
except pd.errors.EmptyDataError:
    eligible_motifs = pd.Series(dtype=int)

Path(snakemake.output.status).parent.mkdir(parents=True, exist_ok=True)
eligible_motifs.rename("accessible_site_count").to_csv(snakemake.output.motifs, sep="\t", header=True)
eligible = replicate_ok and qc_ok and len(eligible_motifs) > 0
status = {
    "eligible": eligible,
    "replicates_per_level": {key: len(value) for key, value in groups.items()},
    "replicate_requirement_met": replicate_ok,
    "all_samples_footprinting_qc_eligible": qc_ok,
    "minimum_accessible_sites_per_motif": minimum_sites,
    "eligible_motif_count": len(eligible_motifs),
    "interpretation": "computational footprint inference; not direct TF binding evidence",
}
Path(snakemake.output.status).write_text(json.dumps(status, indent=2) + "\n")
