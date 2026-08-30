from pathlib import Path

import pandas as pd


counts = pd.read_csv(snakemake.input.counts, sep="\t", comment="#")
fixed = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
sample_columns = [column for column in counts.columns if column not in fixed]
sample_ids = list(snakemake.params.sample_ids)
if len(sample_columns) != len(sample_ids):
    raise ValueError(f"featureCounts returned {len(sample_columns)} BAM columns for {len(sample_ids)} samples")
counts = counts[["Geneid", "Chr", "Start", "End", "Length"] + sample_columns]
counts.columns = ["peak_id", "chrom", "start_1based", "end", "length"] + sample_ids
Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)
counts.to_csv(snakemake.output[0], sep="\t", index=False)
