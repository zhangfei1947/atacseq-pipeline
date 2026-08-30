def unit_fastq(wildcards, mate):
    return UNITS_DF.loc[wildcards.unit, f"fastq_r{mate}"]


def unit_adapter_mode(wildcards):
    value = UNITS_DF.loc[wildcards.unit, "adapter_mode"].strip()
    return value or PARAMS["trimming"]["adapter_mode"]


def sample_unit_bams(wildcards):
    return [f"{OUT}/bam/unit/{unit}.sorted.bam" for unit in units_for_sample(wildcards.sample)]


rule fastqc_raw:
    input:
        r1=lambda wc: unit_fastq(wc, 1),
        r2=lambda wc: unit_fastq(wc, 2),
    output:
        r1=f"{OUT}/qc/fastqc/raw/{{unit}}_R1_fastqc.html",
        r2=f"{OUT}/qc/fastqc/raw/{{unit}}_R2_fastqc.html",
    log:
        f"{OUT}/logs/fastqc/raw/{{unit}}.log"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=30
    container:
        CONTAINERS["fastqc"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.r1}) $(dirname {log})
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        fastqc --threads {threads} --outdir "$tmp" {input.r1} >> {log} 2>&1
        mv "$tmp"/*_fastqc.html {output.r1}
        rm -f "$tmp"/*_fastqc.zip
        fastqc --threads {threads} --outdir "$tmp" {input.r2} >> {log} 2>&1
        mv "$tmp"/*_fastqc.html {output.r2}
        """


rule trim_adapters:
    input:
        r1=lambda wc: unit_fastq(wc, 1),
        r2=lambda wc: unit_fastq(wc, 2),
    output:
        r1=f"{OUT}/fastq/trimmed/{{unit}}_R1.fastq.gz",
        r2=f"{OUT}/fastq/trimmed/{{unit}}_R2.fastq.gz",
    params:
        mode=unit_adapter_mode,
        a1=PARAMS["trimming"]["adapter_r1"],
        a2=PARAMS["trimming"]["adapter_r2"],
        minlen=PARAMS["trimming"]["minimum_length"],
    log:
        f"{OUT}/logs/cutadapt/{{unit}}.log"
    threads: 4
    resources:
        mem_mb=4000,
        runtime=120
    container:
        CONTAINERS["cutadapt"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.r1}) $(dirname {log})
        if [[ "{params.mode}" == "none" ]]; then
            cp {input.r1} {output.r1}
            cp {input.r2} {output.r2}
            printf 'adapter_mode\tnone\n' > {log}
        elif [[ "{params.mode}" == "auto" || "{params.mode}" == "explicit" ]]; then
            cutadapt --cores {threads} --minimum-length {params.minlen} --pair-filter=any \
                -a {params.a1} -A {params.a2} \
                -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log} 2>&1
        else
            echo "Unknown adapter_mode: {params.mode}" >&2
            exit 2
        fi
        """


rule fastqc_trimmed:
    input:
        r1=f"{OUT}/fastq/trimmed/{{unit}}_R1.fastq.gz",
        r2=f"{OUT}/fastq/trimmed/{{unit}}_R2.fastq.gz",
    output:
        r1=f"{OUT}/qc/fastqc/trimmed/{{unit}}_R1_fastqc.html",
        r2=f"{OUT}/qc/fastqc/trimmed/{{unit}}_R2_fastqc.html",
    log:
        f"{OUT}/logs/fastqc/trimmed/{{unit}}.log"
    threads: 2
    resources:
        mem_mb=2000,
        runtime=30
    container:
        CONTAINERS["fastqc"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.r1}) $(dirname {log})
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        fastqc --threads {threads} --outdir "$tmp" {input.r1} >> {log} 2>&1
        mv "$tmp"/*_fastqc.html {output.r1}
        rm -f "$tmp"/*_fastqc.zip
        fastqc --threads {threads} --outdir "$tmp" {input.r2} >> {log} 2>&1
        mv "$tmp"/*_fastqc.html {output.r2}
        """


rule bowtie2_align:
    input:
        r1=f"{OUT}/fastq/trimmed/{{unit}}_R1.fastq.gz",
        r2=f"{OUT}/fastq/trimmed/{{unit}}_R2.fastq.gz",
        index=lambda wc: REF["bowtie2_index"] + ".1.bt2",
    output:
        temp(f"{OUT}/bam/unit/{{unit}}.sam")
    params:
        index=REF["bowtie2_index"],
        max_insert=PARAMS["alignment"]["max_insert"],
        dovetail="--dovetail" if PARAMS["alignment"].get("allow_dovetail", False) else "",
        sample=lambda wc: UNITS_DF.loc[wc.unit, "sample_id"],
    log:
        f"{OUT}/logs/bowtie2/{{unit}}.log"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=360
    container:
        CONTAINERS["bowtie2"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output}) $(dirname {log})
        bowtie2 --very-sensitive --end-to-end --no-mixed --no-discordant \
            --maxins {params.max_insert} {params.dovetail} \
            --rg-id {wildcards.unit} --rg SM:{params.sample} --rg PL:ILLUMINA \
            -x {params.index} -1 {input.r1} -2 {input.r2} -p {threads} \
            -S {output} 2> {log}
        """


rule sort_unit_bam:
    input:
        f"{OUT}/bam/unit/{{unit}}.sam"
    output:
        bam=f"{OUT}/bam/unit/{{unit}}.sorted.bam",
        bai=f"{OUT}/bam/unit/{{unit}}.sorted.bam.bai",
    log:
        f"{OUT}/logs/samtools/{{unit}}.sort.log"
    threads: 4
    resources:
        mem_mb=6000,
        runtime=180
    container:
        CONTAINERS["samtools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {log})
        samtools sort -@ {threads} -o {output.bam} {input} 2> {log}
        samtools index -@ {threads} {output.bam} {output.bai}
        """


rule merge_sample_bam:
    input:
        sample_unit_bams
    output:
        temp(f"{OUT}/bam/sample/{{sample}}.merged.bam")
    log:
        f"{OUT}/logs/samtools/{{sample}}.merge.log"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=180
    container:
        CONTAINERS["samtools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output}) $(dirname {log})
        samtools merge -@ {threads} -f -o {output} {input} 2> {log}
        """


rule mark_duplicates:
    input:
        f"{OUT}/bam/sample/{{sample}}.merged.bam"
    output:
        bam=f"{OUT}/bam/sample/{{sample}}.marked.bam",
        bai=f"{OUT}/bam/sample/{{sample}}.marked.bai",
        metrics=f"{OUT}/qc/duplicates/{{sample}}.picard_metrics.txt",
    log:
        f"{OUT}/logs/picard/{{sample}}.log"
    resources:
        mem_mb=12000,
        runtime=360
    container:
        CONTAINERS["picard"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam}) $(dirname {output.metrics}) $(dirname {log})
        picard -Xmx10g MarkDuplicates I={input} O={output.bam} M={output.metrics} \
            REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
            READ_NAME_REGEX=null CREATE_INDEX=true 2> {log}
        """


rule filter_analysis_bam:
    input:
        bam=f"{OUT}/bam/sample/{{sample}}.marked.bam",
        contigs=REF["nuclear_contigs"],
    output:
        temp(f"{OUT}/bam/sample/{{sample}}.nuclear.bam")
    params:
        mapq=PARAMS["alignment"]["mapq"]
    log:
        f"{OUT}/logs/samtools/{{sample}}.filter.log"
    threads: 4
    resources:
        mem_mb=6000,
        runtime=180
    container:
        CONTAINERS["samtools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        contigs=$(paste -sd' ' {input.contigs})
        samtools view -@ {threads} -b -f 2 -F 3852 -q {params.mapq} \
            -o {output} {input.bam} $contigs 2> {log}
        """


rule remove_blacklist:
    input:
        bam=f"{OUT}/bam/sample/{{sample}}.nuclear.bam",
        blacklist=REF["blacklist_bed"],
    output:
        temp(f"{OUT}/bam/sample/{{sample}}.analysis.unsorted.bam")
    log:
        f"{OUT}/logs/bedtools/{{sample}}.blacklist.log"
    resources:
        mem_mb=4000,
        runtime=120
    container:
        CONTAINERS["bedtools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output}) $(dirname {log})
        bedtools intersect -v -abam {input.bam} -b {input.blacklist} > {output} 2> {log}
        """


rule finalize_analysis_bam:
    input:
        f"{OUT}/bam/sample/{{sample}}.analysis.unsorted.bam"
    output:
        bam=f"{OUT}/bam/sample/{{sample}}.analysis.bam",
        bai=f"{OUT}/bam/sample/{{sample}}.analysis.bam.bai",
    log:
        f"{OUT}/logs/samtools/{{sample}}.finalize.log"
    threads: 4
    resources:
        mem_mb=6000,
        runtime=180
    container:
        CONTAINERS["samtools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        samtools sort -@ {threads} -o {output.bam} {input} 2> {log}
        samtools index -@ {threads} {output.bam} {output.bai}
        samtools quickcheck -v {output.bam}
        """


rule namesort_analysis_bam:
    input:
        f"{OUT}/bam/sample/{{sample}}.analysis.bam"
    output:
        temp(f"{OUT}/bam/sample/{{sample}}.analysis.namesort.bam")
    threads: 4
    resources:
        mem_mb=6000,
        runtime=180
    container:
        CONTAINERS["samtools"]["uri"]
    shell:
        "samtools sort -n -@ {threads} -o {output} {input}"


rule fragments:
    input:
        f"{OUT}/bam/sample/{{sample}}.analysis.namesort.bam"
    output:
        f"{OUT}/fragments/{{sample}}.fragments.bed.gz"
    log:
        f"{OUT}/logs/bedtools/{{sample}}.fragments.log"
    resources:
        mem_mb=4000,
        runtime=180
    container:
        CONTAINERS["bedtools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output}) $(dirname {log})
        bedtools bamtobed -bedpe -i {input} 2> {log} | \
            awk 'BEGIN{{OFS="\t"}} $1==$4 && $2>=0 {{s=($2<$5?$2:$5); e=($3>$6?$3:$6); if(e>s) print $1,s,e,$7,1}}' | \
            gzip -c > {output}
        """


rule cutsites:
    input:
        f"{OUT}/fragments/{{sample}}.fragments.bed.gz"
    output:
        f"{OUT}/fragments/{{sample}}.cutsites.bed.gz"
    resources:
        mem_mb=1000,
        runtime=60
    shell:
        r"""
        set -euo pipefail
        gzip -dc {input} | \
            awk 'BEGIN{{OFS="\t"}} {{l=$2+4; r=$3-5; if(l>=0) print $1,l,l+1,$4,"+"; if(r>=0) print $1,r,r+1,$4,"-"}}' | \
            LC_ALL=C sort -k1,1 -k2,2n | gzip -c > {output}
        """


rule fragment_bigwig:
    input:
        bam=f"{OUT}/bam/sample/{{sample}}.analysis.bam",
        bai=f"{OUT}/bam/sample/{{sample}}.analysis.bam.bai",
        blacklist=REF["blacklist_bed"],
    output:
        f"{OUT}/tracks/{{sample}}.fragments.cpm.bw"
    log:
        f"{OUT}/logs/deeptools/{{sample}}.bamcoverage.log"
    threads: 4
    resources:
        mem_mb=8000,
        runtime=180
    container:
        CONTAINERS["deeptools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output}) $(dirname {log})
        bamCoverage -b {input.bam} -o {output} --binSize 10 --normalizeUsing CPM \
            --extendReads --blackListFileName {input.blacklist} -p {threads} 2> {log}
        """


rule cutsite_bedgraph:
    input:
        cuts=f"{OUT}/fragments/{{sample}}.cutsites.bed.gz",
        sizes=REF["chrom_sizes"],
    output:
        temp(f"{OUT}/tracks/{{sample}}.cutsites.cpm.bedGraph")
    resources:
        mem_mb=4000,
        runtime=120
    container:
        CONTAINERS["bedtools"]["uri"]
    shell:
        r"""
        set -euo pipefail
        n=$(gzip -dc {input.cuts} | wc -l)
        test "$n" -gt 0
        scale=$(awk -v n="$n" 'BEGIN{{print 1000000/n}}')
        bedtools genomecov -bg -scale "$scale" -i {input.cuts} -g {input.sizes} | \
            sort -k1,1 -k2,2n > {output}
        """


rule cutsite_bigwig:
    input:
        bg=f"{OUT}/tracks/{{sample}}.cutsites.cpm.bedGraph",
        sizes=REF["chrom_sizes"],
    output:
        f"{OUT}/tracks/{{sample}}.cutsites.cpm.bw"
    resources:
        mem_mb=1000,
        runtime=30
    container:
        CONTAINERS["bedgraphtobigwig"]["uri"]
    shell:
        "bedGraphToBigWig {input.bg} {input.sizes} {output}"
