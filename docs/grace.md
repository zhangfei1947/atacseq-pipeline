# TAMU HPRC Grace

## Control environment

Create the small micromamba environment once. Analysis programs are not
installed into it.

```bash
module load WebProxy
micromamba create -f environment.yaml
```

## Containers

Analysis rules use versioned OCI/Biocontainer URIs declared in
`containers/images.yaml`. On Grace, load Singularity and WebProxy on a compute
node, then run:

```bash
bin/pull-containers "$SCRATCH/atacseq-container-cache"
```

The generated `containers/locked.tsv` records URI, local SIF and SHA256. Keep
the cache in scratch. Runtime provenance records the same values. Containers are
never committed to Git.

## Controller and child jobs

`bin/submit-grace-controller.sh` starts Snakemake in a small compute-node job.
The Grace profile uses:

- `executor: cluster-generic`;
- `profiles/grace/slurm-submit.py` to map rule threads, memory and runtime;
- `profiles/grace/slurm-status.py` to query `sacct` and fall back to `squeue`;
- `scancel` for cancellation.

Set account/partition in `profiles/grace/config.yaml` or pass them as Snakemake
resource overrides. No free-space or quota gate is implemented; project work is
assumed to reside in the user's 10-TB scratch allocation.

The controller script accepts an absolute config path:

```bash
sbatch --account=ACCOUNT bin/submit-grace-controller.sh \
  /scratch/USER/PROJECT/config/config.yaml
```
