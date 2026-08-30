# Drosophila bulk ATAC-seq pipeline

A provenance-first Snakemake workflow for paired-end bulk ATAC-seq. Version 1
targets *Drosophila melanogaster* FlyBase r6.68, runs on a Linux workstation or
TAMU HPRC Grace, and accepts local `fastq.gz` files (including SRA data converted
beforehand).

## Scope

- technical-unit trimming/QC and Bowtie2 alignment;
- biological-sample merging, duplicate marking, nuclear/MAPQ/blacklist filtering;
- fragment and Tn5-cutsite files, CPM coverage and cutsite tracks;
- sample and pooled MACS3 `BAMPE` peaks;
- practical 2–4 replicate pairwise IDR, strict and majority peak sets;
- fixed-width peak universe and fragment-level featureCounts;
- ATAC-only differential accessibility with DESeq2;
- known and de novo motif enrichment plus FIMO scanning;
- eligibility-gated TOBIAS footprinting;
- positional peak-to-gene annotation and optional REDfly evidence;
- MultiQC, machine-readable QC decisions, methods and provenance.

RNA-seq, ChIP-seq and CUT&Tag integration are intentionally outside v1. Stable
sample metadata, peak BEDs, raw count matrices, normalized matrices and gene
annotation tables are retained as future integration interfaces.

## Quick start

```bash
micromamba create -f environment.yaml
micromamba activate atacseq-control
cp config/config.example.yaml config/config.yaml
cp config/samples.example.tsv config/samples.tsv
cp config/analyses.example.yaml config/analyses.yaml
snakemake --configfile config/config.yaml --profile profiles/local --dry-run
snakemake --configfile config/config.yaml --profile profiles/local
```

Edit all three copied files and add the r6.68 reference bundle described in
[`docs/reference-bundle.md`](docs/reference-bundle.md) before a real run.

On Grace, submit the small controller job from a login node:

```bash
sbatch bin/submit-grace-controller.sh /scratch/user/project/config/config.yaml
```

The controller runs on a compute node. It submits rule jobs through
`cluster-generic`; it does not require the Snakemake Slurm executor plugin.

See [`docs/design.md`](docs/design.md), [`docs/grace.md`](docs/grace.md), and
[`docs/outputs.md`](docs/outputs.md) for the frozen v1 decisions.
