# Frozen v1 design

## Data model

- `unit_id`: one paired FASTQ technical unit; aligned independently.
- `sample_id`: one biological sample; technical units merge here.
- `replicate_group`: biological replicates expected to share a condition.
- `analysis_id`: a reusable sample set with zero or more contrasts.

One BioProject/paper should normally be one project configuration. A sample may
belong to multiple analysis sets. Technical units never count as biological
replicates.

## Alignment and analysis BAM

Cutadapt performs adapter removal with a retained-read minimum of 20 nt and no
aggressive quality trimming. Bowtie2 uses end-to-end `--very-sensitive`, paired
concordant alignment, `--no-mixed`, `--no-discordant`, and maximum insert 2,000.
Dovetail pairs are disabled by default and configurable.

Technical-unit BAMs merge by biological sample. Picard marks duplicates before
the analysis BAM is filtered to proper pairs, MAPQ >=30, primary,
non-supplementary, non-QC-fail, nonduplicate nuclear alignments. Blacklisted
fragments are removed. Mitochondrial reads remain in pre-filter metrics but are
excluded from analysis. SRA headers without optical-coordinate fields produce
an unavailable optical-duplicate metric rather than a fabricated value.

## Peaks and reproducibility

MACS3 is called as `callpeak -f BAMPE -q 0.01 --call-summits --keep-dup all`.
There is no ATAC shift/extsize override because BAMPE uses inferred fragments;
duplicate filtering occurs upstream. Exact invocation is captured in logs and
provenance.

IDR uses relaxed MACS3 `p=0.1` input lists, p-value ranking, and a 0.05 global
threshold. Primary reported peaks remain `q=0.01`. Two replicates use one direct
pair. Three or four replicates use every pair. Pseudoreplicates are not created
by default. A locus must overlap an IDR-supported pair and peaks from at least
ceil(n/2) replicates for the majority set, or every replicate for the strict set.
Pooled peaks are retained as an oracle/discovery output, not as evidence of
replicate reproducibility. Groups outside 2–4 replicates receive a documented
skip status instead of an invented procedure.

## Differential accessibility

The analysis universe merges group-reproducible summits into nonoverlapping
500-bp windows. Original MACS3 peaks and 200-bp motif windows are retained.
Paired fragments are counted once using featureCounts paired-fragment mode.
DESeq2 uses raw counts; the default significance rule is FDR <=0.05 and
|log2FC| >=1. Results preserve tested, significant, opening and closing tables.
VST data support PCA/correlation; shrunken effect estimates support ranking and
plots, while unshrunk estimates and p-values remain the test record.

No spike-in/global-shift correction is implied. A report warning explains that
DESeq2's compositional normalization cannot establish a genome-wide absolute
accessibility shift.

## Module gates

QC never silently deletes a sample. It emits `technical_status`, `qc_status`,
reasons, a recommended action, and module-specific eligibility. The frozen
analysis manifest records the actual inclusion decision.

- DA: at least two independent replicates per condition and a full-rank design;
  three or more are preferred.
- motif: independent of RNA-seq; <50 foreground peaks skip, 50–199 is
  exploratory, >=200 is standard.
- footprinting: requires the analysis to opt in, >=2 replicates/condition,
  acceptable cutsite/QC evidence, and >=100 accessible sites per motif. Results
  are computational footprint inference, not direct binding proof.
- positional peak-to-gene annotation always runs. Optional REDfly overlap runs
  when provided. Bulk ATAC-only co-accessibility is off by default and requires
  >=20 independent samples before an exploratory promoter-correlation workflow
  is permitted.

## Stable future-integration files

Sample metadata snapshots, fragment BED, consensus/universe BED, raw counts,
normalized matrices, DA tables, motif results, gene annotation and provenance
are stable interfaces. They use explicit FlyBase build and gene identifiers so
future RNA-seq, ChIP-seq or CUT&Tag integration does not require reprocessing
FASTQs.
