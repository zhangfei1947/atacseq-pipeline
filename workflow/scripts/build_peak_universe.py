"""Select strongest nonoverlapping fixed-width summit-centered windows."""

from pathlib import Path


def read_bed(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                f = line.rstrip().split("\t")
                rows.append((f[0], int(f[1]), int(f[2]), f))
    return rows


sizes = {}
with open(snakemake.input.sizes) as handle:
    for line in handle:
        chrom, size = line.split()[:2]
        sizes[chrom] = int(size)

allowed = [read_bed(path) for path in snakemake.input.peaks]
allowed = [row for group in allowed for row in group]

def in_allowed(chrom, pos):
    return any(c == chrom and start <= pos < end for c, start, end, _ in allowed)


candidates = []
for path in snakemake.input.summits:
    for chrom, start, end, fields in read_bed(path):
        pos = (start + end) // 2
        score = float(fields[4]) if len(fields) > 4 else 0.0
        if chrom in sizes and in_allowed(chrom, pos):
            candidates.append((score, chrom, pos))

# A safe fallback for unusually sparse/failed summit output.
if not candidates:
    for chrom, start, end, fields in allowed:
        score = float(fields[4]) if len(fields) > 4 else 0.0
        candidates.append((score, chrom, (start + end) // 2))

da_width = int(snakemake.params.da_window)
motif_width = int(snakemake.params.motif_window)
accepted = []
for score, chrom, pos in sorted(candidates, reverse=True):
    start = max(0, pos - da_width // 2)
    end = min(sizes[chrom], start + da_width)
    start = max(0, end - da_width)
    if end - start != da_width:
        continue
    if any(c == chrom and s < end and e > start for c, s, e, _, _ in accepted):
        continue
    accepted.append((chrom, start, end, score, pos))

accepted.sort(key=lambda row: (row[0], row[1]))
Path(snakemake.output.bed).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.output.bed, "w") as bed, open(snakemake.output.motif, "w") as motif, open(snakemake.output.saf, "w") as saf:
    saf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    for idx, (chrom, start, end, score, summit) in enumerate(accepted, 1):
        peak_id = f"peak_{idx:07d}"
        bed.write(f"{chrom}\t{start}\t{end}\t{peak_id}\t{score:g}\t.\n")
        mstart = max(0, summit - motif_width // 2)
        mend = min(sizes[chrom], mstart + motif_width)
        mstart = max(0, mend - motif_width)
        motif.write(f"{chrom}\t{mstart}\t{mend}\t{peak_id}\t{score:g}\t.\n")
        saf.write(f"{peak_id}\t{chrom}\t{start + 1}\t{end}\t.\n")

