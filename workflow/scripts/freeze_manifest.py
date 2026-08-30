from pathlib import Path

import pandas as pd


units = pd.read_csv(snakemake.input.samples, sep="\t", dtype=str, keep_default_na=False)
samples = units.drop_duplicates("sample_id").set_index("sample_id")
selected = samples.loc[list(snakemake.params.sample_ids)].reset_index()
selected.insert(0, "analysis_id", snakemake.wildcards.analysis)
selected.insert(1, "genome_build", snakemake.params.genome)
selected.insert(2, "analysis_include", "yes")
Path(snakemake.output[0]).parent.mkdir(parents=True, exist_ok=True)
selected.to_csv(snakemake.output[0], sep="\t", index=False)

