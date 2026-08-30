def qc_bams(wildcards):
    return [f"{OUT}/bam/sample/{sample}.analysis.bam" for sample in samples_for_analysis(wildcards.analysis)]


def qc_ataqv_files(wildcards):
    return [f"{OUT}/qc/ataqv/{sample}.json" for sample in samples_for_analysis(wildcards.analysis)]


def qc_frip_files(wildcards):
    return [f"{OUT}/qc/frip/{sample}.tsv" for sample in samples_for_analysis(wildcards.analysis)]


def qc_fragment_files(wildcards):
    return [f"{OUT}/fragments/{sample}.fragments.bed.gz" for sample in samples_for_analysis(wildcards.analysis)]


def qc_picard_files(wildcards):
    return [f"{OUT}/qc/duplicates/{sample}.picard_metrics.txt" for sample in samples_for_analysis(wildcards.analysis)]


def qc_flagstat_files(wildcards):
    return [f"{OUT}/qc/alignment/{sample}.flagstat.txt" for sample in samples_for_analysis(wildcards.analysis)]


rule sample_flagstat:
    input:
        f"{OUT}/bam/sample/{{sample}}.marked.bam"
    output:
        f"{OUT}/qc/alignment/{{sample}}.flagstat.txt"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=60
    container:
        CONTAINERS["samtools"]["uri"]
    shell:
        "mkdir -p $(dirname {output}) && samtools flagstat -@ {threads} {input} > {output}"


rule sample_frip:
    input:
        fragments=f"{OUT}/fragments/{{sample}}.fragments.bed.gz",
        peaks=f"{OUT}/peaks/sample/{{sample}}.peaks.narrowPeak",
    output:
        f"{OUT}/qc/frip/{{sample}}.tsv"
    resources:
        mem_mb=3000,
        runtime=60
    container:
        CONTAINERS["bedtools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output})
        total=$(gzip -dc {input.fragments} | wc -l)
        in_peaks=$(bedtools intersect -u -a {input.fragments} -b {input.peaks} | wc -l)
        awk -v t="$total" -v p="$in_peaks" 'BEGIN{{OFS="\t"; print "total_fragments","fragments_in_peaks","frip"; print t,p,(t?p/t:0)}}' > {output}
        """


rule ataqv:
    input:
        bam=f"{OUT}/bam/sample/{{sample}}.marked.bam",
        bai=f"{OUT}/bam/sample/{{sample}}.marked.bai",
        peaks=f"{OUT}/peaks/sample/{{sample}}.peaks.narrowPeak",
        tss=REF["tss_bed"],
        blacklist=REF["blacklist_bed"],
        autosomes=REF["nuclear_contigs"],
    output:
        json=f"{OUT}/qc/ataqv/{{sample}}.json",
        text=f"{OUT}/qc/ataqv/{{sample}}.txt",
    params:
        mito=REF.get("mitochondrial_contigs", ["mitochondrion_genome"])[0]
    threads: 4
    resources:
        mem_mb=12000,
        runtime=360
    container:
        CONTAINERS["ataqv"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.json})
        ataqv --threads {threads} --peak-file {input.peaks} --tss-file {input.tss} \
            --excluded-region-file {input.blacklist} --autosomal-reference-file {input.autosomes} \
            --mitochondrial-reference-name {params.mito} --name {wildcards.sample} \
            --ignore-read-groups --metrics-file {output.json} fly {input.bam} > {output.text}
        """


rule multibam_summary:
    input:
        bams=qc_bams
    params:
        labels=lambda wc: " ".join(samples_for_analysis(wc.analysis))
    output:
        npz=f"{OUT}/qc/replicates/{{analysis}}.bins.npz",
        raw=f"{OUT}/qc/replicates/{{analysis}}.raw_counts.tsv",
    log:
        f"{OUT}/logs/deeptools/{{analysis}}.multibam.log"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=360
    container:
        CONTAINERS["deeptools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.npz}) $(dirname {log})
        multiBamSummary bins --bamfiles {input.bams} --labels {params.labels} \
            --binSize 1000 --outFileName {output.npz} --outRawCounts {output.raw} \
            --numberOfProcessors {threads} > {log} 2>&1
        """


rule replicate_plots:
    input:
        f"{OUT}/qc/replicates/{{analysis}}.bins.npz"
    output:
        correlation=f"{OUT}/qc/replicates/{{analysis}}.spearman.pdf",
        matrix=f"{OUT}/qc/replicates/{{analysis}}.spearman.tsv",
        pca=f"{OUT}/qc/replicates/{{analysis}}.pca.pdf",
    resources:
        mem_mb=6000,
        runtime=120
    container:
        CONTAINERS["deeptools"]["uri"]
    shell:
        r"""
        plotCorrelation -in {input} --corMethod spearman --skipZeros --whatToPlot heatmap \
            --plotFile {output.correlation} --outFileCorMatrix {output.matrix}
        plotPCA -in {input} -o {output.pca}
        """


rule qc_decisions:
    input:
        manifest=f"{OUT}/analyses/{{analysis}}.manifest.tsv",
        ataqv=qc_ataqv_files,
        frip=qc_frip_files,
        fragments=qc_fragment_files,
        picard=qc_picard_files,
        flagstat=qc_flagstat_files,
        correlations=f"{OUT}/qc/replicates/{{analysis}}.spearman.tsv",
    output:
        f"{OUT}/analyses/{{analysis}}.qc_decisions.tsv"
    params:
        thresholds=PARAMS["qc"],
        sample_ids=lambda wc: samples_for_analysis(wc.analysis),
    script:
        "../scripts/qc_decisions.py"
