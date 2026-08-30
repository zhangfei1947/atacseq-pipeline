def contrast_level_samples(analysis_id, contrast_id, role):
    spec = ANALYSES[analysis_id]["contrasts"][contrast_id]
    level = str(spec[role])
    factor = spec["factor"]
    selected = samples_for_analysis(analysis_id)
    return [sample for sample in selected if str(SAMPLES_DF.loc[sample, factor]) == level]


def contrast_group_bam(wildcards, role):
    samples = contrast_level_samples(wildcards.analysis, wildcards.contrast, role)
    groups = sorted(SAMPLES_DF.loc[samples, "replicate_group"].unique())
    if len(groups) != 1:
        raise ValueError(
            f"Footprinting requires one replicate_group per contrast level; "
            f"{wildcards.analysis}/{wildcards.contrast}/{role} has {groups}"
        )
    return f"{OUT}/bam/groups/{groups[0]}.pooled.bam"


rule footprint_eligibility:
    input:
        qc=f"{OUT}/analyses/{{analysis}}.qc_decisions.tsv",
        manifest=f"{OUT}/analyses/{{analysis}}.manifest.tsv",
        fimo=f"{OUT}/motifs/{{analysis}}/fimo/fimo.tsv",
    output:
        status=f"{OUT}/footprints/{{analysis}}/{{contrast}}/eligibility.json",
        motifs=f"{OUT}/footprints/{{analysis}}/{{contrast}}/eligible_motifs.tsv",
    params:
        factor=lambda wc: contrast_spec(wc)["factor"],
        numerator=lambda wc: contrast_spec(wc)["numerator"],
        denominator=lambda wc: contrast_spec(wc)["denominator"],
        minimum_sites=100,
    script:
        "../scripts/footprint_eligibility.py"


rule tobias_footprinting:
    input:
        numerator=lambda wc: contrast_group_bam(wc, "numerator"),
        denominator=lambda wc: contrast_group_bam(wc, "denominator"),
        eligibility=f"{OUT}/footprints/{{analysis}}/{{contrast}}/eligibility.json",
        motif_eligibility=f"{OUT}/footprints/{{analysis}}/{{contrast}}/eligible_motifs.tsv",
        peaks=f"{OUT}/analyses/{{analysis}}.peak_universe.bed",
        genome=REF["fasta"],
        blacklist=REF["blacklist_bed"],
        motifs=REF["motif_database"],
    output:
        f"{OUT}/footprints/{{analysis}}/{{contrast}}/done"
    params:
        outdir=lambda wc: f"{OUT}/footprints/{wc.analysis}/{wc.contrast}",
        numerator_name=lambda wc: str(contrast_spec(wc)["numerator"]),
        denominator_name=lambda wc: str(contrast_spec(wc)["denominator"]),
    log:
        f"{OUT}/logs/tobias/{{analysis}}.{{contrast}}.log"
    threads: 8
    resources:
        mem_mb=32000,
        runtime=1440
    container:
        CONTAINERS["tobias"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        if ! grep -q '"eligible": true' {input.eligibility}; then
            printf 'skipped: eligibility criteria not met\n' > {output}
            exit 0
        fi
        TOBIAS ATACorrect --bam {input.numerator} --genome {input.genome} --peaks {input.peaks} \
            --blacklist {input.blacklist} --outdir {params.outdir}/corrected \
            --prefix numerator --cores {threads} >> {log} 2>&1
        TOBIAS ATACorrect --bam {input.denominator} --genome {input.genome} --peaks {input.peaks} \
            --blacklist {input.blacklist} --outdir {params.outdir}/corrected \
            --prefix denominator --cores {threads} >> {log} 2>&1
        TOBIAS ScoreBigwig --signal {params.outdir}/corrected/numerator_corrected.bw \
            --regions {input.peaks} --output {params.outdir}/numerator_footprints.bw \
            --cores {threads} >> {log} 2>&1
        TOBIAS ScoreBigwig --signal {params.outdir}/corrected/denominator_corrected.bw \
            --regions {input.peaks} --output {params.outdir}/denominator_footprints.bw \
            --cores {threads} >> {log} 2>&1
        TOBIAS BINDetect --motifs {input.motifs} \
            --signals {params.outdir}/numerator_footprints.bw {params.outdir}/denominator_footprints.bw \
            --cond-names {params.numerator_name} {params.denominator_name} \
            --genome {input.genome} --peaks {input.peaks} --outdir {params.outdir}/bindetect \
            --prefix bindetect --cores {threads} >> {log} 2>&1
        printf 'complete: computational footprint inference\n' > {output}
        """

