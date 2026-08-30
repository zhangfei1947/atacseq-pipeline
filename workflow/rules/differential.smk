def groups_for_analysis(analysis_id):
    return sorted(SAMPLES_DF.loc[samples_for_analysis(analysis_id), "replicate_group"].unique())


def analysis_majority_peaks(wildcards):
    return [f"{OUT}/peaks/groups/{group}.majority.bed" for group in groups_for_analysis(wildcards.analysis)]


def analysis_group_summits(wildcards):
    return [f"{OUT}/peaks/groups/{group}.pooled.summits.bed" for group in groups_for_analysis(wildcards.analysis)]


def analysis_bams(wildcards):
    return [f"{OUT}/bam/sample/{sample}.analysis.bam" for sample in samples_for_analysis(wildcards.analysis)]


def contrast_spec(wildcards):
    return ANALYSES[wildcards.analysis]["contrasts"][wildcards.contrast]


rule analysis_manifest:
    input:
        samples=config["samples"],
        analyses=config["analyses"],
    output:
        f"{OUT}/analyses/{{analysis}}.manifest.tsv"
    params:
        sample_ids=lambda wc: samples_for_analysis(wc.analysis),
        genome=config["project"]["genome_build"],
    script:
        "../scripts/freeze_manifest.py"


rule peak_universe:
    input:
        peaks=analysis_majority_peaks,
        summits=analysis_group_summits,
        sizes=REF["chrom_sizes"],
    output:
        bed=f"{OUT}/analyses/{{analysis}}.peak_universe.bed",
        motif=f"{OUT}/analyses/{{analysis}}.motif_windows.bed",
        saf=f"{OUT}/analyses/{{analysis}}.peak_universe.saf",
    params:
        da_window=PARAMS["peaks"]["da_window"],
        motif_window=PARAMS["peaks"]["motif_window"],
    resources:
        mem_mb=4000,
        runtime=60
    script:
        "../scripts/build_peak_universe.py"


rule featurecounts_fragments:
    input:
        bams=analysis_bams,
        saf=f"{OUT}/analyses/{{analysis}}.peak_universe.saf",
    output:
        raw=temp(f"{OUT}/analyses/{{analysis}}.featureCounts.txt")
    log:
        f"{OUT}/logs/featurecounts/{{analysis}}.log"
    threads: 8
    resources:
        mem_mb=12000,
        runtime=360
    container:
        CONTAINERS["subread"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.raw}) $(dirname {log})
        featureCounts -T {threads} -F SAF -p --countReadPairs -B -C \
            -a {input.saf} -o {output.raw} {input.bams} > {log} 2>&1
        """


rule counts_matrix:
    input:
        counts=f"{OUT}/analyses/{{analysis}}.featureCounts.txt",
        manifest=f"{OUT}/analyses/{{analysis}}.manifest.tsv",
    output:
        f"{OUT}/analyses/{{analysis}}.counts.tsv"
    params:
        sample_ids=lambda wc: samples_for_analysis(wc.analysis)
    script:
        "../scripts/featurecounts_to_matrix.py"


rule differential_accessibility:
    input:
        counts=f"{OUT}/analyses/{{analysis}}.counts.tsv",
        manifest=f"{OUT}/analyses/{{analysis}}.manifest.tsv",
    output:
        results=f"{OUT}/differential/{{analysis}}/{{contrast}}/results.tsv",
        opening=f"{OUT}/differential/{{analysis}}/{{contrast}}/opening.bed",
        closing=f"{OUT}/differential/{{analysis}}/{{contrast}}/closing.bed",
        normalized=f"{OUT}/differential/{{analysis}}/{{contrast}}/normalized_counts.tsv",
        pca=f"{OUT}/differential/{{analysis}}/{{contrast}}/pca.pdf",
        summary=f"{OUT}/differential/{{analysis}}/{{contrast}}/summary.json",
    params:
        design=lambda wc: ANALYSES[wc.analysis]["design"],
        factor=lambda wc: contrast_spec(wc)["factor"],
        numerator=lambda wc: contrast_spec(wc)["numerator"],
        denominator=lambda wc: contrast_spec(wc)["denominator"],
        fdr=PARAMS["differential"]["fdr"],
        lfc=PARAMS["differential"]["abs_log2fc"],
        min_count=PARAMS["differential"]["min_count"],
    log:
        f"{OUT}/logs/deseq2/{{analysis}}.{{contrast}}.log"
    threads: 4
    resources:
        mem_mb=16000,
        runtime=360
    container:
        CONTAINERS["deseq2"]["uri"]
    script:
        "../scripts/deseq2_atac.R"

