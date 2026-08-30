"""Deterministically match background sequences by GC-decile."""

import hashlib
import json
from collections import defaultdict
from pathlib import Path


def read_fasta(path):
    records = []
    name = None
    seq = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq).upper()))
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        records.append((name, "".join(seq).upper()))
    return records


def gc_bin(seq):
    valid = [base for base in seq if base in "ACGT"]
    if not valid:
        return 0
    return min(9, int(10 * sum(base in "GC" for base in valid) / len(valid)))


fg = read_fasta(snakemake.input.foreground)
universe = read_fasta(snakemake.input.universe)
fg_ids = {name.split("::")[0] for name, _ in fg}
available = defaultdict(list)
for name, seq in universe:
    if name.split("::")[0] not in fg_ids:
        available[gc_bin(seq)].append((name, seq))
for key in available:
    available[key].sort(key=lambda row: hashlib.sha256(row[0].encode()).hexdigest())

chosen = []
need = defaultdict(int)
for _, seq in fg:
    need[gc_bin(seq)] += 1
for key, count in need.items():
    chosen.extend(available[key][:count])

# If a sparse bin is exhausted, fill deterministically from all unused peaks.
target = len(fg)
used = {name for name, _ in chosen}
if len(chosen) < target:
    remainder = sorted(
        [row for rows in available.values() for row in rows if row[0] not in used],
        key=lambda row: hashlib.sha256(row[0].encode()).hexdigest(),
    )
    chosen.extend(remainder[: target - len(chosen)])

Path(snakemake.output.fasta).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.output.fasta, "w") as handle:
    for name, seq in chosen:
        handle.write(f">{name}\n{seq}\n")

n = len(fg)
exploratory = int(snakemake.params.exploratory)
standard = int(snakemake.params.standard)
status = "standard" if n >= standard else ("exploratory" if n >= exploratory else "skip")
Path(snakemake.output.status).write_text(json.dumps({
    "foreground_sequences": n,
    "background_sequences": len(chosen),
    "status": status,
    "background_method": "GC-decile matched within tested peak universe",
}, indent=2) + "\n")

