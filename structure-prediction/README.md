# Structure Prediction

Retrieve or predict protein structures for candidate protein sequences. The
default pipeline searches PDB and AlphaFold DB templates first, downloads
direct structures, then optionally runs MODELLER homology modeling and
user-specified Phase B prediction.

## Required tools (user-configurable)

You can use your own installed tools and databases. This skill does not assume
fixed server paths or a fixed user:

- DIAMOND (`diamond`)
- A local PDB seqres FASTA and/or `.dmnd` database
- A local AFDB full-sequence FASTA and/or `.dmnd` database
- MODELLER in a conda environment (license required); set the environment name
  with `STRUCTURE_PRED_MODELLER_ENV`, for example
  `export STRUCTURE_PRED_MODELLER_ENV=modeller`
- A Phase B predictor only if the user requests Phase B, for example NVIDIA NIM
  OpenFold2 or ESMFold

If you use conda and want a minimal Python helper env:

```bash
conda env create -f environment.yml
conda activate structure-prediction
```

## Environment variables

```bash
export STRUCTURE_PRED_ROOT=/path/to/workdir
export STRUCTURE_PRED_AFDB_FASTA=/path/to/afdb_sequences.fasta
export STRUCTURE_PRED_AFDB_MIN_PLDDT=60.0
export STRUCTURE_PRED_MODELLER_ENV=meta
```

`STRUCTURE_PRED_MODELLER_ENV` defaults to `meta` for backward compatibility
with existing examples, but should be changed to your own MODELLER conda
environment name.

## Quick run

```bash
bash /Users/lily/.codex/skills/structure-prediction/scripts/run_structure_pipeline.sh   --pdb-db /path/to/pdb.dmnd   --afdb-db /path/to/afdb.dmnd   --threads 44   --model-source pdb
```

See `SKILL.md` for the full stage-by-stage workflow, quality gates, provenance,
and resumability requirements.
