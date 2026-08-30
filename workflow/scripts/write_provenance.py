import hashlib
import json
import os
import shutil
import subprocess
import importlib.metadata
from datetime import datetime, timezone
from pathlib import Path

import yaml


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


Path(snakemake.output.methods).parent.mkdir(parents=True, exist_ok=True)
images = yaml.safe_load(Path(snakemake.input.images).read_text())
try:
    git_commit = subprocess.run(["git", "rev-parse", "HEAD"], check=True, text=True, capture_output=True).stdout.strip()
except subprocess.SubprocessError:
    git_commit = "uncommitted"

methods = f"""# Methods

Paired-end bulk ATAC-seq technical units were adapter-trimmed with Cutadapt
{images['cutadapt']['version']} (minimum retained length {snakemake.params.parameters['trimming']['minimum_length']} nt)
and aligned end-to-end to {snakemake.params.project['genome_build']} with Bowtie2
{images['bowtie2']['version']} using `--very-sensitive --no-mixed --no-discordant`
and maximum insert size {snakemake.params.parameters['alignment']['max_insert']} bp.
Technical units were merged per biological sample. Picard
{images['picard']['version']} marked duplicates; optical duplicate detection was
recorded as unavailable because SRA-converted read names may not preserve tile
coordinates. Analysis BAMs retained proper, primary, non-supplementary,
nonduplicate nuclear pairs at MAPQ >= {snakemake.params.parameters['alignment']['mapq']}
and excluded configured blacklist regions.

Peaks were called with MACS3 {images['macs3']['version']} using paired-fragment
mode: `macs3 callpeak -f BAMPE -q {snakemake.params.parameters['peaks']['qvalue']}
--call-summits --keep-dup all`. No shift/extsize override or control library was
used. Relaxed MACS3 `p={snakemake.params.parameters['peaks']['idr_peak_pvalue']}`
lists were generated specifically for IDR. Pairwise IDR {images['idr']['version']}
used `--rank p.value` and global IDR threshold
{snakemake.params.parameters['peaks']['idr_threshold']}. Two replicates
used a direct pair; three or four used every pair; pseudoreplicates were not
created. Majority loci required support from at least half (rounded up) of
biological replicates and an IDR-supported pair; strict loci required every
replicate.

Reproducible summits were represented as nonoverlapping
{snakemake.params.parameters['peaks']['da_window']}-bp windows. Paired fragments
were counted once with featureCounts {images['subread']['version']}. DESeq2
{images['deseq2']['version']} used raw counts, the configured full-rank design,
FDR <= {snakemake.params.parameters['differential']['fdr']} and absolute log2 fold
change >= {snakemake.params.parameters['differential']['abs_log2fc']}. DESeq2
compositional normalization does not demonstrate an absolute genome-wide
accessibility shift.

Known and de novo motif enrichment used AME and STREME from MEME Suite
{images['meme']['version']} on summit-centered
{snakemake.params.parameters['peaks']['motif_window']}-bp sequences with
GC-decile-matched tested-peak backgrounds. FIMO scanned the tested universe.
Eligible contrasts optionally used TOBIAS {images['tobias']['version']} for Tn5
bias correction and footprint inference. Footprints are computational evidence,
not direct TF-binding measurements. ChIPseeker {images['chipseeker']['version']}
provided positional peak-to-gene annotation; optional REDfly overlap and
eligibility-gated bulk ATAC promoter correlations were kept as separate evidence.
"""
Path(snakemake.output.methods).write_text(methods)

Path(snakemake.output.config).write_text(yaml.safe_dump(dict(snakemake.params.full_config), sort_keys=False))
shutil.copyfile(snakemake.input.samples, snakemake.output.samples)
shutil.copyfile(snakemake.input.analyses, snakemake.output.analyses)

locked = Path("containers/locked.tsv")
provenance = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "project": dict(snakemake.params.project),
    "git_commit": git_commit,
    "snakemake_version": importlib.metadata.version("snakemake"),
    "config_sha256": sha256(snakemake.output.config),
    "samples_sha256": sha256(snakemake.input.samples),
    "analyses_sha256": sha256(snakemake.input.analyses),
    "reference": dict(snakemake.params.reference),
    "parameters": dict(snakemake.params.parameters),
    "modules": dict(snakemake.params.modules),
    "container_images": images,
    "sif_lock": locked.read_text().splitlines() if locked.exists() else "not_prefetched_in_this_checkout",
    "note": "Exact shell commands remain in rule logs and Snakemake job logs.",
}
Path(snakemake.output.provenance).write_text(json.dumps(provenance, indent=2) + "\n")
