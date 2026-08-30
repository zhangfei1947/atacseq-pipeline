#!/usr/bin/env python3
"""Translate Slurm states to cluster-generic's running/success/failed API."""

import subprocess
import sys


jobid = sys.argv[-1]
success = {"COMPLETED"}
failed = {
    "BOOT_FAIL", "CANCELLED", "DEADLINE", "FAILED", "NODE_FAIL",
    "OUT_OF_MEMORY", "PREEMPTED", "REVOKED", "TIMEOUT",
}
running = {
    "CONFIGURING", "COMPLETING", "PENDING", "REQUEUED", "RESIZING",
    "RUNNING", "SUSPENDED", "STAGE_OUT",
}

query = subprocess.run(
    ["sacct", "-X", "-j", jobid, "--format=State", "--noheader", "--parsable2"],
    text=True, capture_output=True, check=False,
)
states = [line.strip().split("|")[0].split("+")[0].split()[0] for line in query.stdout.splitlines() if line.strip()]
if any(state in failed for state in states):
    print("failed")
elif states and all(state in success for state in states):
    print("success")
elif any(state in running for state in states):
    print("running")
else:
    queued = subprocess.run(
        ["squeue", "-h", "-j", jobid, "-o", "%T"],
        text=True, capture_output=True, check=False,
    )
    print("running" if queued.stdout.strip() else "failed")
