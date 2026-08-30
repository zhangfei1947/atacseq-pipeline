#!/bin/bash
#SBATCH --job-name=atac-controller
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/slurm/controller.%j.out
#SBATCH --error=logs/slurm/controller.%j.err
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: sbatch $0 /absolute/path/to/config.yaml" >&2
    exit 2
fi

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
config_file="$(realpath "$1")"
mkdir -p "$repo_dir/logs/slurm"
cd "$repo_dir"

module load Singularity
eval "$(micromamba shell hook --shell bash)"
micromamba activate atacseq-control

container_cache="${SCRATCH:?SCRATCH is not defined}/atacseq-container-cache"
export TMPDIR="$SCRATCH/atacseq-tmp"
export SINGULARITY_CACHEDIR="$SCRATCH/.singularity"
mkdir -p "$TMPDIR" "$SINGULARITY_CACHEDIR" "$container_cache"
snakemake --profile profiles/grace --configfile "$config_file" \
    --apptainer-prefix "$container_cache"
