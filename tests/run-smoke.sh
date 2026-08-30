#!/usr/bin/env bash
set -euo pipefail
bin/validate-config tests/config/config.yaml --check-files
snakemake --configfile tests/config/config.yaml --cores 1 --dry-run
snakemake --configfile tests/config/config.yaml --lint

