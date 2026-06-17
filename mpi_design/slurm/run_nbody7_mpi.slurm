#!/bin/bash
# ===========================================================================
# NBODY7 internal-MPI (Path B) -- Slurm production launcher
# ===========================================================================
# Submit from the gpudyn-slurm gateway:
#     sbatch /users/sambaran/updated-nbody7/mpi_design/slurm/run_nbody7_mpi.slurm
#
# What this runs: the standalone MPI binary that parallelises the O(N^2) force
# evaluations + start-up over MPI ranks with the copy algorithm (every rank
# holds the full identical state -> bit-identical to the serial run for any np;
# see mpi_design/README.md and 01_data_distribution.md).
#
# Defaults below reproduce the VERIFIED multi-GPU config on gpudyn1
# (results_multigpu_throughput_gpudyn1_2026-06-17.log): mpi-gpu, np=2, strict
# one-rank-per-GPU on the 2x RTX 2080. Pure-MPI (OMP_NUM_THREADS=1) and the
# rank-0 I/O guard (production default) are enforced below.
#
# Partition -> binary (architecture-gated CUDA build):
#     turing  (gpudyn1/2, 2x RTX 2080, sm_75) -> GPU2/run_turing/...
#     ampere  (gpudyn3,   1x A40,      sm_86)  -> GPU2/run_ampere/...
# selected automatically from $SLURM_JOB_PARTITION.
# ---------------------------------------------------------------------------
#SBATCH --job-name=nbody7-mpi
#SBATCH --partition=turing
#SBATCH --nodelist=gpudyn1
#SBATCH --nodes=1
#SBATCH --ntasks=2
# GPU ranks need spare CPU cores for the CUDA driver/runtime HOST threads even
# though OMP_NUM_THREADS=1: under the Slurm cpuset, --cpus-per-task=1 crams each
# rank's integrator thread + CUDA helper threads onto one core and kernel
# launch/sync latency balloons (job 259: 437s vs 273s bare-metal, same np2).
# 8 cores/rank = the whole 16-core turing node (you already hold both GPUs).
# For mpi-cpu set this back to 1 (no GPU host threads to feed).
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
# Slurm job log (script echoes + timing summary). Absolute path on the gpudyn1
# local disk -- KEEP IN SYNC WITH $RUNDIR below. Physics output goes to
# $RUNDIR/run.out and $RUNDIR/err.out (the validated convention), not here.
#SBATCH --output=/vol/gpudyn1/data1/sambaran/Test_nb7_mpi/slurm-%x-%j.out
# ===========================================================================

set -u

# ============================ USER CONFIG ==================================
REPO=/users/sambaran/updated-nbody7                       # NFS home (shared)
RUNDIR=/vol/gpudyn1/data1/sambaran/Test_nb7_mpi           # gpudyn1 local disk
VARIANT=mpi-gpu              # mpi-gpu (default) | mpi-cpu
INPUT=input                 # unit-5 parameter file inside $RUNDIR
# ---- optional knobs (uncomment / export to enable) ----
# export NBODY_PHASE_TIMERS=1     # per-phase wall timers on rank 0 stderr
# export NBODY_GPU_RANK_BIND=0    # opt out of rank->GPU pinning (mpi-gpu only)
# export NBODY_RANK0_IO=0         # disable rank-0 I/O guard (per-rank-dir mode!)
# ===========================================================================
#
# Switching to the CPU build:  set VARIANT=mpi-cpu, drop "--gres" above, and
# raise --ntasks to the core count (turing 16, ampere 32). Everything else is
# identical -- the launch recipe is variant-agnostic.
#
# For mpi-gpu DO NOT oversubscribe the GPUs: keep --ntasks <= #GPUs on the node
# (np > nGPU serializes kernel launches and is slower -- see the throughput log).

# --- partition -> architecture-matched binary directory --------------------
case "${SLURM_JOB_PARTITION:-turing}" in
    turing) BINDIR=run_turing ;;   # sm_75
    ampere) BINDIR=run_ampere ;;   # sm_86
    *) echo "ERROR: no binary mapped for partition '${SLURM_JOB_PARTITION}'"; exit 1 ;;
esac
BIN="$REPO/GPU2/$BINDIR/nbody7b.$VARIANT"
NP=${SLURM_NTASKS:-1}

# --- pre-flight sanity checks ----------------------------------------------
[ -x "$BIN" ]        || { echo "ERROR: binary not found / not executable: $BIN"; exit 1; }
command -v mpirun >/dev/null || { echo "ERROR: mpirun not on PATH"; exit 1; }
cd "$RUNDIR"         || { echo "ERROR: cannot cd to RUNDIR=$RUNDIR"; exit 1; }
[ -f "$INPUT" ]      || { echo "ERROR: missing unit-5 input '$INPUT' in $RUNDIR"; exit 1; }
[ -f fort.10 ]       || { echo "ERROR: missing lowercase 'fort.10' in $RUNDIR (gfortran unit-10)"; exit 1; }

# --- pure-MPI: REQUIRED for the bit-identical guarantee --------------------
export OMP_NUM_THREADS=1

echo "=========================================================================="
echo " NBODY7 Path B MPI run"
echo "   job          : ${SLURM_JOB_NAME} (${SLURM_JOB_ID:-interactive})"
echo "   node(s)      : ${SLURM_JOB_NODELIST:-$(hostname)}   partition=${SLURM_JOB_PARTITION:-?}"
echo "   variant      : $VARIANT      binary=$BIN"
echo "   ranks (np)   : $NP   gres=${SLURM_JOB_GRES:-none}   CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "   rundir       : $RUNDIR     input(unit5)=$INPUT"
echo "   started      : $(date '+%F %T')"
echo "=========================================================================="

# --- launch ----------------------------------------------------------------
# `bash -c 'exec BIN < INPUT'` runs the stdin redirect INSIDE each rank, so every
# rank gets unit 5 = the input file (mpirun forwards its own stdin to rank 0
# only; NBODY7 also does lazy mid-run unit-5 reads). --bind-to none + OMP=1 is
# the validated recipe. Rank->GPU pinning is automatic via OMPI_COMM_WORLD_LOCAL_RANK.
SECONDS=0
mpirun --bind-to none \
       -x OMP_NUM_THREADS=1 \
       -np "$NP" \
       bash -c 'exec "'"$BIN"'" < "'"$INPUT"'"' \
       > run.out 2> err.out
rc=$?

echo "=========================================================================="
echo "   finished     : $(date '+%F %T')   wall=${SECONDS}s   mpirun_rc=$rc"
if grep -q 'END RUN' run.out 2>/dev/null; then
    echo "   status       : reached END RUN"
else
    echo "   status       : DID NOT reach END RUN -- inspect $RUNDIR/{run,err}.out"
fi
# rank->device placement banner (PROFILE build prints it to stderr at GPU open)
grep -E 'device [0-9]+:|Rank-0 I/O guard' err.out run.out 2>/dev/null | sed 's/^/   /'
echo "=========================================================================="
exit $rc
