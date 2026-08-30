PAIR_LOOKUP = {
    (group, token): (a, b)
    for group, pairs in GROUP_PAIRS.items()
    for a, b, token in pairs
}


def idr_sample_peak(wildcards, which):
    sample = PAIR_LOOKUP[(wildcards.group, wildcards.pair)][which]
    return f"{OUT}/peaks/idr_input/sample/{sample}.narrowPeak"


def idr_pair_files(wildcards):
    return [f"{OUT}/peaks/groups/idr/{wildcards.group}/{token}.idr.narrowPeak" for _, _, token in GROUP_PAIRS[wildcards.group]]


def group_sample_peaks(wildcards):
    return [f"{OUT}/peaks/sample/{sample}.peaks.narrowPeak" for sample in samples_for_group(wildcards.group)]


rule pairwise_idr:
    input:
        p1=lambda wc: idr_sample_peak(wc, 0),
        p2=lambda wc: idr_sample_peak(wc, 1),
        oracle=f"{OUT}/peaks/idr_input/groups/{{group}}.narrowPeak",
    output:
        peaks=f"{OUT}/peaks/groups/idr/{{group}}/{{pair}}.idr.narrowPeak",
        plot=f"{OUT}/peaks/groups/idr/{{group}}/{{pair}}.idr.narrowPeak.png",
    params:
        threshold=PARAMS["peaks"]["idr_threshold"]
    log:
        f"{OUT}/logs/idr/{{group}}/{{pair}}.log"
    resources:
        mem_mb=8000,
        runtime=240
    container:
        CONTAINERS["idr"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.peaks}) $(dirname {log})
        idr --samples {input.p1} {input.p2} --peak-list {input.oracle} \
            --input-file-type narrowPeak --output-file-type narrowPeak --rank p.value \
            --idr-threshold {params.threshold} --soft-idr-threshold {params.threshold} \
            --plot --output-file {output.peaks} --log-output-file {log}
        """


rule reproducible_group_peaks:
    input:
        peaks=group_sample_peaks,
        idr=idr_pair_files,
    output:
        majority=f"{OUT}/peaks/groups/{{group}}.majority.bed",
        strict=f"{OUT}/peaks/groups/{{group}}.strict.bed",
        status=f"{OUT}/peaks/groups/{{group}}.idr_status.json",
    params:
        samples=lambda wc: samples_for_group(wc.group),
        threshold=PARAMS["peaks"]["idr_threshold"],
    resources:
        mem_mb=4000,
        runtime=60
    script:
        "../scripts/build_reproducible_peaks.py"
