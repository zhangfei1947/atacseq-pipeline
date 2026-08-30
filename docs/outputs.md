# Output contract

All paths are rooted at `output_dir`.

| Directory | Stable contents |
|---|---|
| `qc/` | FastQC, alignment/duplicate/fragment metrics, replicate correlation, QC decisions and MultiQC. |
| `bam/` | coordinate-sorted unit BAMs, marked sample BAMs and filtered analysis BAMs. |
| `fragments/` | one-record-per-pair fragments and Tn5 cut sites. |
| `tracks/` | CPM fragment bigWigs, cutsite bigWigs and optional TOBIAS corrected tracks. |
| `peaks/sample/` | original sample MACS3 peaks and summits. |
| `peaks/groups/` | pooled, pairwise IDR, majority and strict reproducible peaks. |
| `analyses/<id>/` | frozen manifest, peak universe, raw count matrix and analysis-set QC. |
| `differential/<contrast>/` | complete DESeq2 table, opening/closing BED and plots. |
| `motifs/<contrast>/` | AME/STREME/FIMO results and eligibility record. |
| `footprints/<contrast>/` | eligibility, bias-corrected tracks and TOBIAS results when eligible. |
| `annotation/` | peak-to-gene positional annotation and optional REDfly overlaps. |
| `report/` | MultiQC report, methods, configuration snapshot and provenance JSON/TSV. |

Intermediate technical-unit trimmed FASTQs are reproducible working products,
not stable integration interfaces. BAMs and stable tables are never removed by
the workflow.

