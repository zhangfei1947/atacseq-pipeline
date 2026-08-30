# Changelog

## Unreleased

- Make peak-to-gene annotation compatible with Bioconductor 3.22 by deriving
  nearest-gene/TSS annotations directly from the configured GTF. This avoids
  the now-defunct `GenomicFeatures::makeTxDbFromGFF()` call and keeps gene
  identifiers tied to the selected reference release.
- Write AME and STREME results in complete text formats because the pinned
  MEME 5.5.9 BioContainer does not include the HTML web-template data files.
- Load SciPy only when optional bulk co-accessibility is eligible to run, so
  disabled or sample-size-ineligible analyses still emit their status files in
  lean workflow environments.
- Allow an explicit `autosomes_file` for ataqv QC, with configuration validation and a backward-compatible fallback to `nuclear_contigs`.
- Retry transient Singularity registry failures and identify each image while populating the digest-pinned cache.
- Coordinate-sort Tn5 cut sites before `bedtools genomecov`, which requires each chromosome to occur in one sequential block.

## 0.1.0 — 2026-08-29

- Initial paired-end bulk ATAC-seq workflow for FlyBase r6.68.
- Local and TAMU Grace `cluster-generic` execution profiles.
- ATAC-only DA, practical IDR, motif, gated footprinting and peak-to-gene modules.
