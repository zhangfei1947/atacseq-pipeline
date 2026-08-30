#!/usr/bin/env python3
"""Submit one Snakemake jobscript to Slurm and print only its job ID."""

import os
import re
import subprocess
import sys
from pathlib import Path

from snakemake.utils import read_job_properties


def minutes_to_slurm(value):
    minutes = max(1, int(float(value)))
    days, remainder = divmod(minutes, 24 * 60)
    hours, mins = divmod(remainder, 60)
    prefix = f"{days}-" if days else ""
    return f"{prefix}{hours:02d}:{mins:02d}:00"


jobscript = Path(sys.argv[-1]).resolve()
props = read_job_properties(str(jobscript))
resources = props.get("resources", {})
threads = max(1, int(props.get("threads", 1)))
mem_mb = max(256, int(resources.get("mem_mb", 2000)))
runtime = resources.get("runtime", 60)
rule = props.get("rule", "job")
jobid = props.get("jobid", "na")
name = re.sub(r"[^A-Za-z0-9_.-]", "_", f"atac.{rule}.{jobid}")[:128]

Path("logs/slurm").mkdir(parents=True, exist_ok=True)
command = [
    "sbatch",
    "--parsable",
    f"--job-name={name}",
    f"--cpus-per-task={threads}",
    f"--mem={mem_mb}M",
    f"--time={minutes_to_slurm(runtime)}",
    f"--output=logs/slurm/{name}.%j.out",
    f"--error=logs/slurm/{name}.%j.err",
]
account = str(resources.get("slurm_account", "")).strip()
partition = str(resources.get("slurm_partition", "")).strip()
if account:
    command.append(f"--account={account}")
if partition:
    command.append(f"--partition={partition}")
command.append(str(jobscript))

result = subprocess.run(command, check=True, text=True, capture_output=True)
print(result.stdout.strip().split(";")[0])
