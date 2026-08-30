def direction_bed(wildcards):
    return f"{OUT}/differential/{wildcards.analysis}/{wildcards.contrast}/{wildcards.direction}.bed"


rule universe_fasta:
    input:
        bed=f"{OUT}/analyses/{{analysis}}.motif_windows.bed",
        fasta=REF["fasta"],
    output:
        f"{OUT}/motifs/{{analysis}}.universe.fa"
    resources:
        mem_mb=2000,
        runtime=60
    container:
        CONTAINERS["bedtools"]["uri"]
    shell:
        r"""
        mkdir -p $(dirname {output})
        bedtools getfasta -name -fi {input.fasta} -bed {input.bed} -fo {output}
        """


rule motif_foreground_fasta:
    input:
        sig=direction_bed,
        windows=f"{OUT}/analyses/{{analysis}}.motif_windows.bed",
        fasta=REF["fasta"],
    output:
        bed=f"{OUT}/motifs/{{analysis}}/{{contrast}}/{{direction}}.foreground.bed",
        fasta=f"{OUT}/motifs/{{analysis}}/{{contrast}}/{{direction}}.foreground.fa",
    resources:
        mem_mb=2000,
        runtime=60
    container:
        CONTAINERS["bedtools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bed})
        bedtools intersect -wa -a {input.windows} -b {input.sig} > {output.bed}
        if [[ -s {output.bed} ]]; then
            bedtools getfasta -name -fi {input.fasta} -bed {output.bed} -fo {output.fasta}
        else
            : > {output.fasta}
        fi
        """


rule motif_matched_background:
    input:
        foreground=f"{OUT}/motifs/{{analysis}}/{{contrast}}/{{direction}}.foreground.fa",
        universe=f"{OUT}/motifs/{{analysis}}.universe.fa",
    output:
        fasta=f"{OUT}/motifs/{{analysis}}/{{contrast}}/{{direction}}.background.fa",
        status=f"{OUT}/motifs/{{analysis}}/{{contrast}}/{{direction}}.eligibility.json",
    params:
        exploratory=PARAMS["motif"]["minimum_exploratory_peaks"],
        standard=PARAMS["motif"]["minimum_standard_peaks"],
    script:
        "../scripts/match_motif_background.py"


rule fimo_universe:
    input:
        fasta=f"{OUT}/motifs/{{analysis}}.universe.fa",
        motifs=REF["motif_database"],
    output:
        done=f"{OUT}/motifs/{{analysis}}/fimo.done",
        tsv=f"{OUT}/motifs/{{analysis}}/fimo/fimo.tsv",
    params:
        outdir=lambda wc: f"{OUT}/motifs/{wc.analysis}/fimo",
        threshold=PARAMS["motif"]["fimo_threshold"],
    log:
        f"{OUT}/logs/meme/{{analysis}}.fimo.log"
    resources:
        mem_mb=12000,
        runtime=720
    container:
        CONTAINERS["meme"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.done}) $(dirname {log})
        fimo --oc {params.outdir} --thresh {params.threshold} {input.motifs} {input.fasta} > {log} 2>&1
        printf 'complete\n' > {output.done}
        """


rule motif_enrichment:
    input:
        foreground=f"{OUT}/motifs/{{analysis}}/{{contrast}}/{{direction}}.foreground.fa",
        background=f"{OUT}/motifs/{{analysis}}/{{contrast}}/{{direction}}.background.fa",
        eligibility=f"{OUT}/motifs/{{analysis}}/{{contrast}}/{{direction}}.eligibility.json",
        fimo=f"{OUT}/motifs/{{analysis}}/fimo.done",
        motifs=REF["motif_database"],
    output:
        f"{OUT}/motifs/{{analysis}}/{{contrast}}/{{direction}}.done"
    params:
        outbase=lambda wc: f"{OUT}/motifs/{wc.analysis}/{wc.contrast}/{wc.direction}",
        minimum=PARAMS["motif"]["minimum_exploratory_peaks"],
    log:
        f"{OUT}/logs/meme/{{analysis}}.{{contrast}}.{{direction}}.log"
    resources:
        mem_mb=12000,
        runtime=720
    container:
        CONTAINERS["meme"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outbase} $(dirname {log})
        n=$(grep -c '^>' {input.foreground} || true)
        if [[ "$n" -lt {params.minimum} ]]; then
            printf 'skipped: only %s foreground peaks\n' "$n" > {output}
            exit 0
        fi
        ame --oc {params.outbase}/ame --control {input.background} {input.foreground} {input.motifs} >> {log} 2>&1
        streme --oc {params.outbase}/streme --p {input.foreground} --n {input.background} >> {log} 2>&1
        printf 'complete: %s foreground peaks\n' "$n" > {output}
        """
