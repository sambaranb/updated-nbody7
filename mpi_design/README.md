# Running NBODY7 with internal MPI (Path B)

The `mpi-cpu` build variant parallelises the O(N²) force evaluations (regular +
irregular) and the O(N²) start-up (FPOLY0 + FPOLY2) over MPI ranks with a
copy-algorithm decomposition: every rank holds the full, identical particle
state, so the results are **bit-identical to the serial run for any number of
ranks** (pure-MPI mode, see caveats below). KS/chain/ARCHAIN/BSE and all
diagnostics stay replicated. Design notes: `01_data_distribution.md`;
validation harnesses and result logs: `poc_validation/`.

## Build

```bash
conda activate Amuse-env            # Mac; on Linux hosts with system-wide
                                    # Open MPI (e.g. gpudyn3) no env needed
cd GPU2
make mpi-cpu CXX="$CXX"             # -> GPU2/run_versions/nbody7b.mpi-cpu
```

CUDA variant (Path B G0; needs nvcc + an NVIDIA cuda-samples checkout for
`helper_cuda.h`):

```bash
make clean                          # REQUIRED when switching variants:
                                    # targets alias gpunb/gpupot/cnbint/
                                    # gpuirr .o names via cp -p
make mpi-gpu SDK_PATH=/usr/local/cuda-samples CUDA_PATH=/usr/local/cuda
                                    # -> GPU2/run_versions/nbody7b.mpi-gpu
```

Each rank binds to one GPU (`local_rank % nGPU`; opt-out
`NBODY_GPU_RANK_BIND=0`). On both CUDA targets (`gpu` and `mpi-gpu`) the
CPU-side irregular-force pair defaults to the portable auto-vectorized
sources (3x faster than the hand-SSE pair on AVX-512 Xeons); `GPUIRR=sse`
restores the hand-written SSE pair for pre-AVX Intel CPUs (bit-identity
holds within a pair, not across them).

(Mac toolchain notes — conda clang for `-fopenmp`, system linker via
`-B/usr/bin`, `-lc++` — are in `01_data_distribution.md` §9.)

## Run: one shared directory (production default)

Put the input files in a single run directory, exactly as for a serial run
(`Fort.10`/`fort.10` initial conditions, `input_bse` if BSE is on, the
parameter input file):

```bash
cd test                             # the run directory
export OMP_NUM_THREADS=1            # pure-MPI: required for bit-identical results
mpirun -np 8 bash -c 'exec ../path/to/nbody7b.mpi-cpu < input' > run.out 2> err.out
```

**Why `bash -c 'exec ... < input'` instead of `mpirun ... < input`:** every rank
replicates the integrator state, so every rank must read the *full* input —
including the main parameter file on unit 5 (stdin), and NBODY7 also performs
lazy mid-run unit-5 reads (the ksint/chain CLIGHT lines). `mpirun` connects its
own stdin to **rank 0 only** by default, so a plain `< input` redirect makes
ranks 1..N-1 fail on their first `READ (5,...)`. The `bash -c 'exec ... <
input'` form performs the redirect *inside each rank's process*, giving every
rank unit 5 = the regular file `test/input`. (`mpirun --stdin all` is the
launcher-side alternative where supported.)

What happens in the shared directory:

- **Reads — every rank.** All ranks concurrently read the same `Fort.10`,
  `input_bse`, restart `fort.1`, and the input file. Concurrent reads of
  ordinary files are harmless; nothing to configure.
- **Writes — rank 0 only.** The broad rank-0 I/O guard (on by default for
  np > 1; banner `Rank-0 I/O guard active` in `run.out`) connects stdout and
  every output unit of ranks > 0 to `/dev/null`. Rank 0 produces the one
  canonical set of output files (`OUT3`, `OUT9`, `fort.*`, restart dumps,
  stdout), exactly as a serial run would. Stderr is left live on every rank so
  genuine runtime errors remain visible.
- **Manual termination**: `touch STOP` in the run directory works; rank 0
  probes the file and broadcasts the decision, so all ranks stop together.
- **Restart**: `fort.1`/`fort.2` are written by rank 0 only (all ranks hold
  identical state, so rank 0's dump is canonical); on restart every rank reads
  the same `fort.1`.

## Per-rank-directory mode (validation / debugging)

`NBODY_RANK0_IO=0` disables the guard: every rank performs the full original
I/O, as in the serial code. Use this **only** with one working directory per
rank (otherwise the ranks race on the same output files). This is the mode the
validation harnesses use to compare every rank's output bit-for-bit:

```bash
poc_validation/run_equiv.sh         # per-rank dirs, guard off: np 1/2/4 bit-identical
poc_validation/run_io_guard.sh      # shared dir, guard on: np4 outputs identical to np1
poc_validation/profile_scaling.sh   # scaling profile (per-rank dirs)
```

## Caveats

- **Bit-identity is a pure-MPI guarantee**: run with `OMP_NUM_THREADS=1`.
  Hybrid MPI+OpenMP (OMP > 1) is out of scope — the one-time FPOLY2 start-up
  loop is not thread-deterministic, so ranks can drift apart and deadlock
  (details in the project notes; serial OMP behaviour is unchanged).
- **`NBSTAT`** (neighbour-count histogram written by the GPUIRR library at
  close) is a per-process diagnostic: at np > 1 it covers rank 0's share of
  the irregular-force calls only.
- The COMMON dumps `fort.1`/`fort.2` embed CPU/wall-time accumulators, so they
  differ in a few bytes between *any* two runs (serial included) — compare
  them by size/structure, not byte-for-byte.
- Serial and AMUSE builds are completely unaffected: they link a no-op stub
  (`mpi_nbody_stub.f`) and execute the original I/O instruction stream.

Shipped at commit `841820d` (guard) on top of `8450d16`/`bf84510` (start-up
decomposition) and the force-decomposition increments; per-increment validation
logs live in `poc_validation/results_*.log`.
