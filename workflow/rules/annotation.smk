def optional_redfly(wildcards):
    path = REF.get("redfly_bed", "")
    return [path] if MODULES.get("redfly", False) and path else []


rule annotate_peaks:
    input:
        peaks=f"{OUT}/analyses/{{analysis}}.peak_universe.bed",
        gtf=REF["gtf"],
        redfly=optional_redfly,
    output:
        genes=f"{OUT}/annotation/{{analysis}}.peak_to_gene.tsv",
        redfly=f"{OUT}/annotation/{{analysis}}.redfly.tsv",
    params:
        has_redfly=lambda wc: bool(optional_redfly(wc))
    log:
        f"{OUT}/logs/annotation/{{analysis}}.log"
    resources:
        mem_mb=16000,
        runtime=360
    container:
        CONTAINERS["chipseeker"]["uri"]
    script:
        "../scripts/annotate_peaks.R"


rule bulk_coaccessibility:
    input:
        counts=f"{OUT}/analyses/{{analysis}}.counts.tsv",
        annotation=f"{OUT}/annotation/{{analysis}}.peak_to_gene.tsv",
        manifest=f"{OUT}/analyses/{{analysis}}.manifest.tsv",
    output:
        links=f"{OUT}/annotation/{{analysis}}.bulk_coaccessibility.tsv",
        status=f"{OUT}/annotation/{{analysis}}.bulk_coaccessibility_status.json",
    params:
        enabled=MODULES.get("bulk_coaccessibility", False),
        minimum_samples=20,
        max_distance=100000,
    resources:
        mem_mb=16000,
        runtime=360
    script:
        "../scripts/bulk_coaccessibility.py"

