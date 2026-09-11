#!/bin/bash
#SBATCH --job-name=x15_si_thz
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#
# One run directory per job: submit with the directory as the only argument,
#     for d in runs/*/; do sbatch run_lomonosov.sbatch "$d"; done
# Each directory is self-contained, so the jobs cannot clobber one another.
set -euo pipefail
RUNDIR="${1:?usage: sbatch run_lomonosov.sbatch <run-directory>}"
SALMON="${SALMON:-salmon}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd "$RUNDIR"
"$SALMON" < Si_prim_epm_gs.inp  > gs.log 2>&1
"$SALMON" < Si_prim_sbe_rt.inp  > rt.log 2>&1
