def report_qc_inputs(wildcards):
    paths = []
    paths += [f"{OUT}/qc/fastqc/raw/{unit}_R1_fastqc.html" for unit in UNITS]
    paths += [f"{OUT}/qc/fastqc/raw/{unit}_R2_fastqc.html" for unit in UNITS]
    paths += [f"{OUT}/qc/fastqc/trimmed/{unit}_R1_fastqc.html" for unit in UNITS]
    paths += [f"{OUT}/qc/fastqc/trimmed/{unit}_R2_fastqc.html" for unit in UNITS]
    paths += [f"{OUT}/qc/alignment/{sample}.flagstat.txt" for sample in SAMPLES]
    paths += [f"{OUT}/qc/duplicates/{sample}.picard_metrics.txt" for sample in SAMPLES]
    paths += [f"{OUT}/qc/ataqv/{sample}.json" for sample in SAMPLES]
    paths += [f"{OUT}/analyses/{analysis}.qc_decisions.tsv" for analysis in ANALYSIS_IDS]
    return paths


def provenance_dependencies(wildcards):
    paths = [f"{OUT}/analyses/{analysis}.counts.tsv" for analysis in ANALYSIS_IDS]
    paths += [f"{OUT}/peaks/groups/{group}.idr_status.json" for group in GROUPS]
    if MODULES.get("differential", True):
        paths += [f"{OUT}/differential/{a}/{c}/summary.json" for a, c in CONTRASTS]
    return paths


rule multiqc:
    input:
        report_qc_inputs
    output:
        html=f"{OUT}/report/multiqc.html",
        data=directory(f"{OUT}/report/multiqc_data"),
    log:
        f"{OUT}/logs/multiqc.log"
    resources:
        mem_mb=8000,
        runtime=180
    container:
        CONTAINERS["multiqc"]["uri"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.html}) $(dirname {log})
        multiqc {OUT} --force --filename $(basename {output.html}) \
            --outdir $(dirname {output.html}) --data-dir --no-ansi > {log} 2>&1
        """


rule methods_and_provenance:
    input:
        dependencies=provenance_dependencies,
        samples=config["samples"],
        analyses=config["analyses"],
        images=config["containers"],
    output:
        methods=f"{OUT}/report/methods.md",
        provenance=f"{OUT}/report/provenance.json",
        config=f"{OUT}/report/config.snapshot.yaml",
        samples=f"{OUT}/report/samples.snapshot.tsv",
        analyses=f"{OUT}/report/analyses.snapshot.yaml",
    params:
        project=config["project"],
        reference=REF,
        parameters=PARAMS,
        modules=MODULES,
        full_config=config,
    script:
        "../scripts/write_provenance.py"
