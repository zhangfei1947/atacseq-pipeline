def group_bams(wildcards):
    return [f"{OUT}/bam/sample/{sample}.analysis.bam" for sample in samples_for_group(wildcards.group)]


rule macs3_sample:
    input:
        f"{OUT}/bam/sample/{{sample}}.analysis.bam"
    output:
        peak=f"{OUT}/peaks/sample/raw/{{sample}}_peaks.narrowPeak",
        summit=f"{OUT}/peaks/sample/raw/{{sample}}_summits.bed",
        xls=f"{OUT}/peaks/sample/raw/{{sample}}_peaks.xls",
    params:
        outdir=lambda wc: f"{OUT}/peaks/sample/raw",
        genome=REF["effective_genome_size"],
        q=PARAMS["peaks"]["qvalue"],
    log:
        f"{OUT}/logs/macs3/sample/{{sample}}.log"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=240
    container:
        CONTAINERS["macs3"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        macs3 callpeak -t {input} -f BAMPE -g {params.genome} -n {wildcards.sample} \
            --outdir {params.outdir} -q {params.q} --call-summits --keep-dup all \
            > {log} 2>&1
        """


rule filter_sample_peaks:
    input:
        peak=f"{OUT}/peaks/sample/raw/{{sample}}_peaks.narrowPeak",
        summit=f"{OUT}/peaks/sample/raw/{{sample}}_summits.bed",
        blacklist=REF["blacklist_bed"],
    output:
        peak=f"{OUT}/peaks/sample/{{sample}}.peaks.narrowPeak",
        summit=f"{OUT}/peaks/sample/{{sample}}.summits.bed",
    resources:
        mem_mb=2000,
        runtime=60
    container:
        CONTAINERS["bedtools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.peak})
        bedtools intersect -v -a {input.peak} -b {input.blacklist} > {output.peak}
        bedtools intersect -v -a {input.summit} -b {input.blacklist} > {output.summit}
        """


rule macs3_sample_idr_input:
    input:
        f"{OUT}/bam/sample/{{sample}}.analysis.bam"
    output:
        f"{OUT}/peaks/idr_input/raw/sample/{{sample}}_peaks.narrowPeak"
    params:
        outdir=lambda wc: f"{OUT}/peaks/idr_input/raw/sample",
        genome=REF["effective_genome_size"],
        p=PARAMS["peaks"]["idr_peak_pvalue"],
    log:
        f"{OUT}/logs/macs3/idr_input/sample/{{sample}}.log"
    threads: 2
    resources:
        mem_mb=8000,
        runtime=240
    container:
        CONTAINERS["macs3"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        macs3 callpeak -t {input} -f BAMPE -g {params.genome} -n {wildcards.sample} \
            --outdir {params.outdir} -p {params.p} --call-summits --keep-dup all \
            > {log} 2>&1
        """


rule filter_sample_idr_input:
    input:
        peak=f"{OUT}/peaks/idr_input/raw/sample/{{sample}}_peaks.narrowPeak",
        blacklist=REF["blacklist_bed"],
    output:
        f"{OUT}/peaks/idr_input/sample/{{sample}}.narrowPeak"
    resources:
        mem_mb=2000,
        runtime=60
    container:
        CONTAINERS["bedtools"]["uri"]
    shell:
        "mkdir -p $(dirname {output}) && bedtools intersect -v -a {input.peak} -b {input.blacklist} > {output}"


rule merge_group_bam:
    input:
        group_bams
    output:
        temp(f"{OUT}/bam/groups/{{group}}.pooled.bam")
    threads: 4
    resources:
        mem_mb=8000,
        runtime=180
    container:
        CONTAINERS["samtools"]["uri"]
    shell:
        r"""
        mkdir -p $(dirname {output})
        samtools merge -@ {threads} -f -o {output} {input}
        """


rule macs3_group_pooled:
    input:
        f"{OUT}/bam/groups/{{group}}.pooled.bam"
    output:
        peak=f"{OUT}/peaks/groups/raw/{{group}}_peaks.narrowPeak",
        summit=f"{OUT}/peaks/groups/raw/{{group}}_summits.bed",
        xls=f"{OUT}/peaks/groups/raw/{{group}}_peaks.xls",
    params:
        outdir=lambda wc: f"{OUT}/peaks/groups/raw",
        genome=REF["effective_genome_size"],
        q=PARAMS["peaks"]["qvalue"],
    log:
        f"{OUT}/logs/macs3/groups/{{group}}.log"
    threads: 2
    resources:
        mem_mb=12000,
        runtime=360
    container:
        CONTAINERS["macs3"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        macs3 callpeak -t {input} -f BAMPE -g {params.genome} -n {wildcards.group} \
            --outdir {params.outdir} -q {params.q} --call-summits --keep-dup all \
            > {log} 2>&1
        """


rule filter_group_peaks:
    input:
        peak=f"{OUT}/peaks/groups/raw/{{group}}_peaks.narrowPeak",
        summit=f"{OUT}/peaks/groups/raw/{{group}}_summits.bed",
        blacklist=REF["blacklist_bed"],
    output:
        peak=f"{OUT}/peaks/groups/{{group}}.pooled.narrowPeak",
        summit=f"{OUT}/peaks/groups/{{group}}.pooled.summits.bed",
    resources:
        mem_mb=2000,
        runtime=60
    container:
        CONTAINERS["bedtools"]["uri"]
    shell:
        r"""
        bedtools intersect -v -a {input.peak} -b {input.blacklist} > {output.peak}
        bedtools intersect -v -a {input.summit} -b {input.blacklist} > {output.summit}
        """


rule macs3_group_idr_oracle:
    input:
        f"{OUT}/bam/groups/{{group}}.pooled.bam"
    output:
        f"{OUT}/peaks/idr_input/raw/groups/{{group}}_peaks.narrowPeak"
    params:
        outdir=lambda wc: f"{OUT}/peaks/idr_input/raw/groups",
        genome=REF["effective_genome_size"],
        p=PARAMS["peaks"]["idr_peak_pvalue"],
    log:
        f"{OUT}/logs/macs3/idr_input/groups/{{group}}.log"
    threads: 2
    resources:
        mem_mb=12000,
        runtime=360
    container:
        CONTAINERS["macs3"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        macs3 callpeak -t {input} -f BAMPE -g {params.genome} -n {wildcards.group} \
            --outdir {params.outdir} -p {params.p} --call-summits --keep-dup all \
            > {log} 2>&1
        """


rule filter_group_idr_oracle:
    input:
        peak=f"{OUT}/peaks/idr_input/raw/groups/{{group}}_peaks.narrowPeak",
        blacklist=REF["blacklist_bed"],
    output:
        f"{OUT}/peaks/idr_input/groups/{{group}}.narrowPeak"
    resources:
        mem_mb=2000,
        runtime=60
    container:
        CONTAINERS["bedtools"]["uri"]
    shell:
        "mkdir -p $(dirname {output}) && bedtools intersect -v -a {input.peak} -b {input.blacklist} > {output}"
