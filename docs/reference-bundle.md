# Reference bundle contract

Reference files are external to Git and configured under `reference`.

Required for a complete run:

| Key | Contract |
|---|---|
| `fasta` | FlyBase r6.68 genomic FASTA; contig names must match every file below. |
| `bowtie2_index` | Bowtie2 index prefix, not an individual `.bt2` filename. |
| `chrom_sizes` | Two-column contig and length file derived from the FASTA. |
| `nuclear_contigs` | One allowed nuclear contig per line. |
| `blacklist_bed` | Sorted BED3+ regions excluded after duplicate filtering. May be an explicitly empty file. |
| `gtf` | FlyBase r6.68 annotation with stable FlyBase gene IDs. |
| `tss_bed` | BED6 TSS positions derived from the selected gene/transcript policy. |
| `motif_database` | MEME-format motif collection; source and version must be documented. |
| `effective_genome_size` | MACS3 numeric effective genome size, recorded in provenance. |

Optional:

- `redfly_bed`: versioned REDfly regulatory elements in the same coordinate build.
- `autosomes`: comma-separated contigs for ataqv complexity summaries.

`bin/validate-config` checks existence, FASTA/index completeness, contig
compatibility, sort order and basic formats before analysis. Reference creation
itself is deliberately outside v1 so that curated files can be supplied later.
