"""Build majority/strict loci from sample peaks gated by pairwise IDR support."""

import json
import math
from pathlib import Path


def read_bed(path):
    rows = []
    if not Path(path).exists():
        return rows
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            rows.append((fields[0], int(fields[1]), int(fields[2]), fields))
    return rows


def overlaps(interval, rows):
    chrom, start, end = interval
    return any(c == chrom and s < end and e > start for c, s, e, _ in rows)


samples = list(snakemake.params.samples)
sample_peaks = [read_bed(path) for path in snakemake.input.peaks]
idr_peaks = [read_bed(path) for path in snakemake.input.idr]

all_intervals = sorted(
    [(c, s, e) for peaks in sample_peaks for c, s, e, _ in peaks],
    key=lambda row: (row[0], row[1], row[2]),
)
merged = []
for chrom, start, end in all_intervals:
    if merged and merged[-1][0] == chrom and start < merged[-1][2]:
        merged[-1] = (chrom, merged[-1][1], max(end, merged[-1][2]))
    else:
        merged.append((chrom, start, end))

n = len(samples)
majority_n = math.ceil(n / 2)
majority = []
strict = []
for locus in merged:
    support = sum(overlaps(locus, peaks) for peaks in sample_peaks)
    idr_supported = n == 1 or any(overlaps(locus, peaks) for peaks in idr_peaks)
    if idr_supported and support >= majority_n:
        majority.append((*locus, support))
    if n >= 2 and idr_supported and support == n:
        strict.append((*locus, support))

for path, rows, label in (
    (snakemake.output.majority, majority, "majority"),
    (snakemake.output.strict, strict, "strict"),
):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        for idx, (chrom, start, end, support) in enumerate(rows, 1):
            handle.write(f"{chrom}\t{start}\t{end}\t{label}_{idx}\t{support}\t.\n")

status = {
    "replicate_count": n,
    "samples": samples,
    "pairwise_idr_comparisons": len(idr_peaks),
    "idr_threshold": float(snakemake.params.threshold),
    "majority_required_replicates": majority_n,
    "majority_peak_count": len(majority),
    "strict_peak_count": len(strict),
    "interpretation": "no_biological_idr" if n == 1 else ("direct_pair" if n == 2 else "all_pairwise"),
}
Path(snakemake.output.status).write_text(json.dumps(status, indent=2) + "\n")

