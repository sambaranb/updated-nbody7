# Path B — Internal MPI parallelisation of the NBODY7 integrator
## Design document 01: `common6.h` data-distribution table + first proof-of-concept

Branch: `mpi-parallel`. Status: **design contract — no integrator code changed yet.**
Target file for the first PoC: `GPU2/intgrt.omp.f` (the OpenMP integrator driver actually
compiled into the cpu/sse/avx/gpu variants; it overrides `Block/intgrt.f`).
Scope of this document: classify every `common6.h` array by its role under MPI, pick the
decomposition strategy for the PoC, and name the exact code region + communication set.

This is the design contract referenced in the project memory (`project_path_b_mpi.md`).
Nothing here alters numerics — the whole point of the chosen PoC strategy is bit-for-bit
identical results to the serial run (see §4).

---

## 0. Sizing facts (from `Ncode/params.h`)

```
NMAX = 350010   ! max single bodies + 3*NBIN + NHI   -> per-particle array length
KMAX = 175005   ! max KS solutions                    -> per-pair array length
LMAX = 650      ! max neighbour-list length
```
Per-particle state is sized at `NMAX`; KS-pair state at `KMAX`; neighbour lists are
`LIST(LMAX,NMAX)`. No MPI scaffolding exists in the tree today (`grep -ri mpif.h` is empty),
so this is a greenfield port — this NBODY7 is the GPU/OpenMP lineage, *not* NBODY6++/MPI.

---

## 1. The integrator's hot loop (where the time goes)

`INTGRT` (in `GPU2/intgrt.omp.f`) per block step does:

1. Build the next block `NXTLST(1:NXTLEN)` (the particles due now).
2. Split it into irregular-only and regular-due lists: `IRR`, `IREG(1:NFR)` (lines 339–349).
3. **Irregular force** on the block: `GPUIRR_FIRR_VEC` + `NBINT`/`NBINTP` (lines 387–444).
4. **Regular force** on the `NFR` regular-due members (lines 447–552) — *the wallclock
   dominator*:
   - `CXVPRED(IFIRST,NTOT,…)` predicts **all** particles to `TIME` (line 453).
   - `GPUNB_SEND(NN, BODY(IFIRST), X(1,IFIRST), XDOT(1,IFIRST))` ships the full
     j-particle set to the force library (line 461).
   - Loop over `IREG` in chunks of `NIMAX=1024`:
       - `GPUNB_REGF(NI, H2I, DTR, XI, VI, GPUACC, GPUJRK, GPUPHI, LMAX, NBMAX, LISTGP)`
         — the O(N·N_nb) regular-force + neighbour-list evaluation (line 483). **This is
         the expensive call.**
       - per-i bookkeeping of the returned neighbour list (lines 486–524).
       - `GPUIRR_FIRR_VEC` then `GPUCOR(I, …)` correct each i-particle (lines 527–543).
5. Housekeeping: next block time, copy-back, `GPUIRR_SET_JP`, SE mass loss, etc.

The cost of step 4 is `NFR × (work per j-particle)`. Steps 3 and 5 are cheaper and already
OpenMP-threaded within a node.

### What `GPUCOR` writes per regular-due particle `I` (from `GPU2/gpucor.f`)
Reads: all j-particle predicted `X, XDOT, BODY`; globals `RDENS, RSCALE, RC22, BODYM, …`.
Writes (per-particle, index `I`): `X0, X0DOT, F, FDOT, FI, FR, FIDOT, FRDOT,
D0, D0R, D1, D1R, D2R, D3R, RS, T0R, STEPR, STEP, TNEW, LIST(:,I)`; and accumulates the
diagnostic counters `NBVOID, NBFULL, NBFLUX, NBSMIN, NICONV` and `ETIDE` (the last under
`!$omp critical`). This 20-array set is the recombination payload (§3).

---

## 2. Data-distribution table for `common6.h`

Role legend:
- **R — replicated**: every rank holds the full, identical copy. Control flags, physical
  constants, scalars, energy/diagnostic accumulators.
- **P — per-particle state (partition candidate)**: length-`NMAX` arrays indexed by particle.
  Replicated in the PoC (copy algorithm); the future full-decomposition target (§5).
- **K — KS-pair / regularisation state**: length-`KMAX` (or chain) arrays. Sequential by
  nature, few in count → stay on rank 0 / replicated; revisit last.
- **S — scratch**: rebuilt every block, never needs gathering.

| COMMON block | Members | Role | Notes |
|---|---|---|---|
| `/NAMES/`  | `NTOT,NPAIRS,N,NNBMAX,…,KZ(50),NNBOPT` | **R** | control + counts; KZ is the option array, must be identical on all ranks |
| `/PARAMS/` | `ETAI,ETAR,…,TIME,TADJ,ZKIN,POT,E(12),ENERGY accumulators,…` | **R** | integration scalars + energy bookkeeping; reductions needed for the accumulators (§3.3) |
| `/COUNTS/` | `NSTEPI,NSTEPR,NSTEPU,NBVOID,NBFULL,NBFLUX,…` | **R+reduce** | diagnostic counters; each rank counts its own slice → `MPI_Reduce` (sum) at output time |
| `/PLPOT/`  | `MP,AP2,VIR,MP0,MPDOT,…` | **R** | tidal-field scalars |
| `/BLOCKS/` | `TPREV,TBLOCK,DTK(40),KVEC(2*KMAX)` | **R** (KVEC: **K**) | block-time machinery replicated; `KVEC` is KS bookkeeping |
| `/STARS/`  | `EPOCH0,…,ZPARS(20),NBH,…,LISTR,LISTD,LISTV` | **R** | stellar-evolution scalars + special-particle lists; small, replicate |
| **`/NBODY/`** | `X,X0,X0DOT,F,FDOT,BODY,RS,XDOT, FI,D1,D2,D3, FR,D1R,D2R,D3R, STEP,T0,STEPR,T0R,TNEW, RADIUS,TEV,TEV0,BODY0,EPOCH,SPIN,XSTAR,ZLMSTY, FIDOT,D0,FRDOT,D0R,KSTAR` | **P** | **the partition candidates.** Force/derivative/step arrays are written by `GPUCOR`; the SE arrays (`RADIUS,TEV,…,KSTAR`) are touched only in `MDOT`/`MLOSS` (rank-0 for PoC) |
| `/PAIRS/`  | `U,U0,UDOT,FU,…,H,HDOT,…(KMAX)`; `KSLOW(KMAX)` | **K** | KS regularisation state; rank-0, not parallelised in PoC |
| `/PAIRS/`  | `NAME(NMAX)` | **P** | per-particle identity; travels with the particle |
| `/PAIRS/`  | `LIST(LMAX,NMAX)` | **P** | per-particle neighbour list; rewritten by `GPUCOR`, part of the gather set |
| `/LISTS/`  | `ILIST(NMAX),JLIST(NMAX),JPERT(5*LMAX)` | **S** | scratch index/perturber lists, rebuilt per use |
| `/SPIN2/`  | `SPN(NMAX),ASPN(NMAX)` | **P** | per-star spin (SE); rank-0 for PoC with the other SE arrays |
| BSE commons `/FLAGS*/,/VALUE*/` | flags + coefficients | **R** | set once in `data.f`, never change during integration |

**Reading of the table for the PoC:** under the copy algorithm (§3) *all* arrays stay
replicated. The table's value is forward-looking — it pre-classifies which arrays the
*future* full decomposition (§5) would partition (the **P** rows), which would live on
rank 0 (the **K**/SE rows), and which need reductions (the **R+reduce** counters and the
`/PARAMS/` energy accumulators).

---

## 3. PoC decomposition strategy — the "copy algorithm" (replicate data, partition work)

This is the NBODY6++/Spurzem approach and matches the memory's design note
("distribute the regular-force loop over the current block across ranks … `MPI_Allgatherv`
of the updated F/FDOT"). Chosen for the PoC because it is **minimally invasive and
provably identical to serial** (§4).

### 3.1 What is partitioned
Only the **work** of step 4's regular-force evaluation. The `NFR` regular-due i-particles
(`IREG(1:NFR)`) are split into `nranks` contiguous slices. Rank `r` computes
`GPUNB_REGF` + the neighbour-list bookkeeping for its slice only.

### 3.2 What is communicated, and when
Two clean options, in increasing redundancy-vs-traffic trade-off:

- **3.2(a) Gather forces, replicate the correction (recommended PoC).**
  Each rank fills its slice of `GPUACC, GPUJRK, GPUPHI` and the per-i neighbour lists
  `LISTGP`. One `MPI_Allgatherv` per block over those force arrays (`3+3+1 = 7` doubles
  per i-particle, plus the variable-length neighbour lists). Then **every rank runs
  `GPUCOR` for all `NFR` particles** (redundant but cheap — `GPUCOR` is correction
  arithmetic, not the O(N) force sum). Because the correction input is now identical on
  every rank, every rank writes identical `/NBODY/` updates → **no gather of the 20-array
  payload needed, and state stays bit-identical across ranks by construction.**
  *Pro:* smallest possible communication; trivially consistent; ~30-line change.
  *Con:* `GPUCOR` redundant compute (negligible vs force eval), and the neighbour-list
  Allgatherv needs care (variable length).

- **3.2(b) Gather the corrected state.** Each rank runs `GPUCOR` on its own slice and we
  `MPI_Allgatherv` the 20-array `/NBODY/` payload (§1) afterward. More traffic, no
  redundant compute. Keep as a fallback if `GPUCOR` redundancy ever shows up in a profile
  (it won't at PoC N).

### 3.3 Reductions
`/COUNTS/` counters and `/PARAMS/` energy accumulators (`ETIDE`, `EMDOT`, …) touched inside
the partitioned region are per-slice partial sums under 3.2(a)?  **No** — under 3.2(a) the
correction is replicated, so counters double-count. Two fixes: (i) only rank 0 commits
counter increments, or (ii) compute counters in the replicated `GPUCOR` and divide-by-none
because every rank has the same totals (they're diagnostics, read only at output on rank 0).
Simplest: let the replicated `GPUCOR` run everywhere, and at output time read counters from
rank 0 only. Document this so the diagnostics aren't misread as `nranks×` inflated.

### 3.4 Communicator hygiene (critical, from memory)
The AMUSE worker owns `MPI_COMM_WORLD` for IPC with the Python host. The integrator's
collectives **must not** collide with that. Bring up a private communicator early:
```
CALL MPI_COMM_DUP(MPI_COMM_WORLD, NBODY_COMM, IERR)   ! or COMM_SELF-derived under AMUSE
CALL MPI_COMM_RANK(NBODY_COMM, MYRANK, IERR)
CALL MPI_COMM_SIZE(NBODY_COMM, NRANKS, IERR)
```
All integrator `Allgatherv`/`Reduce` use `NBODY_COMM`. Under the AMUSE MPI-IPC worker the
spawn already consumes `COMM_WORLD`; the standalone-MPI build will `MPI_INIT` itself. The
PoC targets the **standalone** path first (simpler — no AMUSE in the loop); AMUSE-side MPI
is a later integration step.

---

## 4. Why the PoC is numerically identical to serial

Under 3.2(a) the *only* thing MPI changes is **who computes which i-particle's raw force**.
The force on i-particle `I` from the full j-set is independent of which rank evaluates it
(same inputs `X,XDOT,BODY`, same `GPUNB_REGF` code). After the Allgatherv every rank holds
the identical full `GPUACC/GPUJRK/GPUPHI` set, then runs the identical `GPUCOR`. There is no
reordering of the force summation across ranks (each i-particle's sum is still done whole on
one rank), so there is not even a floating-point reassociation difference. Expected result:
`dE/E` matches the serial cpu run to the last bit. This is the acceptance test for the PoC
(reuse the `nb7_equiv_test.py` harness: serial-cpu vs mpi-N-rank on the 128-body Plummer).

(Contrast: the AVX variant *does* show ~10× looser `dE/E` because its horizontal-add
reorders the per-i sum. The MPI PoC must **not** do that — it splits across i, never within
an i-particle's j-sum.)

---

## 5. Future direction (not PoC): full particle decomposition (Strategy B)

The table's **P** rows become genuinely partitioned: each rank owns `~N/nranks` particles
and their force/step/list arrays. Wins big on memory (the `(3,NMAX)` arrays dominate the
footprint) and lets the j-loop itself be distributed. Costs: predicted j-particles must be
`Allgather`ed every block before the force eval, the 20-array payload gathered after, and
KS/chain (rows **K**) need a home (rank 0 with broadcast of c.m. state). Defer until the
copy-algorithm PoC is proven and profiled — only move to B if memory (very large N) or the
j-loop cost forces it.

---

## 6. Concrete next coding steps (proposed, for review before I touch Fortran)

1. **Build plumbing:** add an `mpi` variant path (new `WITH_MPI` cpp gate + `MPIFC`/`mpif90`
   in `GPU2/Makefile`), so the serial build is untouched and MPI is opt-in. The AMUSE
   worker already builds with `MPIFC`; standalone needs the mpif90 wrapper.
2. **Communicator module:** a tiny `mpi_nbody.h` (COMMON `/MPICOMM/ NBODY_COMM, MYRANK,
   NRANKS, IS_PARALLEL`) + init in `START`/`NBODY6`, mirroring the `amuse.h` pattern.
3. **Decompose the §1 step-4 loop** in `GPU2/intgrt.omp.f` behind `IF (IS_PARALLEL)`:
   slice `IREG`, gate the `GPUNB_REGF` chunk loop to the local slice, `Allgatherv` the
   force arrays (3.2a), keep `GPUCOR` replicated. Serial path (`NRANKS==1`) is the exact
   current code — zero behaviour change when run on one rank.
4. **Acceptance:** `nb7_equiv_test.py`-style: serial-cpu vs `mpirun -np {1,2,4}` on the
   Mac (M4 Pro, MPICH/OpenMPI from conda), assert `dE/E` bit-identical and library
   identifiers unchanged.

**Open decisions for the user (expert) to steer before step 3:**
- 3.2(a) replicate-correction vs 3.2(b) gather-state — I recommend (a) for the PoC.
- PoC first on the **standalone** MPI build (recommended) vs straight into the AMUSE worker.
- Slice granularity: static contiguous (recommended PoC) vs step-level load balance (later;
  block membership varies, but PoC correctness doesn't need balancing).

---

## 7. Decisions taken (user sign-off, 2026-06-02)

1. **Recombination = 3.2(a)** — Allgather the raw force arrays, replicate `GPUCOR` on every
   rank. Smallest traffic; `/NBODY/` state identical across ranks by construction.
2. **First target = standalone MPI build** — prove the decomposition under `mpirun` on a
   plain NBODY7 run, no AMUSE in the loop. AMUSE-worker MPI is a later integration step.
3. **Slicing = static contiguous** — split `IREG(1:NFR)` into `NRANKS` equal contiguous
   chunks. No load balancing in the PoC.

## 8. Implementation architecture — wrapper + stub (no cpp in the hot file)

To keep the just-shipped serial cpu/sse/avx/gpu builds **byte-identical** and avoid
sprinkling `#ifdef WITH_MPI` through the 700-line `intgrt.omp.f`, MPI is confined to a small
wrapper API selected by the build variant — mirroring NBODY7's existing idiom of swapping a
`.f`/`.cpp`/`.cu` file per variant (gpunb.cpu / gpunb.velocity, etc.).

- **`Ncode/mpi_nbody.h`** (canonical, symlinked into `Block/`, `GPU2/` like `amuse.h`):
  declares `COMMON /MPICOMM/ NBODY_COMM, MYRANK, NRANKS, IS_PARALLEL`. Present in *every*
  build. In a serial build `NRANKS=1, MYRANK=0, IS_PARALLEL=.FALSE.`.
- **`GPU2/mpi_nbody.f`** (real, linked only in the `mpi` variant): `NBODY_MPI_INIT`
  (`MPI_INIT` if not already + `MPI_COMM_DUP` → `NBODY_COMM` + rank/size),
  `NBODY_MPI_FINALIZE`, and (increment 2) `NBODY_REGF_ALLGATHER` doing the `MPI_Allgatherv`
  of the force arrays.
- **`GPU2/mpi_nbody_stub.f`** (no-op, linked in serial cpu/sse/avx/gpu + the AMUSE build):
  `NBODY_MPI_INIT` just sets `NRANKS=1`; everything else returns immediately.

The integrator calls the wrapper API unconditionally and branches on `IF (NRANKS.GT.1)`.
When `NRANKS==1` (every serial build, and `mpirun -np 1`) the decomposition is the identity
and the gather is a no-op → the **exact original serial code path** executes. This is what
makes §4's bit-identical guarantee hold, and satisfies the project's never-alter-calculations
rule: the serial numerics are literally the same instructions.

**MPI lifecycle (standalone):** `NBODY_MPI_INIT` is called once at the top of the standalone
driver before `CALL NBODY6`; `NBODY_MPI_FINALIZE` is called at the genuine termination point
(`INTGRT` end: after `GPUNB_CLOSE`/`GPUIRR_CLOSE`, before the final `STOP`). `NBODY_MPI_INIT`
guards with `MPI_INITIALIZED` so it is safe if the AMUSE worker (later) has already
`MPI_INIT`'d the world.

**Build:** a new `mpi-cpu` variant in `GPU2/Makefile` compiles with `mpif90` (Open MPI in
the Mac `Amuse-env`: `/Users/sambaran/miniforge3/envs/Amuse-env/bin/mpif90`), links
`GPU2/mpi_nbody.f`, and passes `-DWITH_MPI` only if/where ever needed (not needed by the
wrapper approach). Serial variants link `GPU2/mpi_nbody_stub.f` and are otherwise unchanged.
The MPI variant family is named by force backend — `mpi-cpu` (this PoC), and later
`mpi-metal` (Apple GPU) and `mpi-gpu` (CUDA) — each producing `$(BINNAME).<variant>`.

**Staging:** increment 1 = scaffolding (header + wrappers + lifecycle wiring + build) proven
to compile and run as a 1-rank identity; increment 2 = the `IREG`-slice decomposition +
`NBODY_REGF_ALLGATHER` in `intgrt.omp.f`, validated with `nb7_equiv_test.py`.

---

## 9. Increment 1 status — DONE, run-verified on Mac (2026-06-03)

Standalone `mpi` variant builds and runs on the M4 Pro. `mpirun -np 1` →
`NRANKS = 1` + normal serial-identical startup; `mpirun -np 2` → `NRANKS = 2`.
Fully additive — to keep the shared `Block/nbody6_main.f` and all shipped builds
byte-identical, the `mpi` variant uses its own driver `GPU2/nbody6_main_mpi.f`
(= nbody6_main.f + `CALL NBODY_MPI_INIT`); the two Makefiles gained 58 inserted lines,
0 deletions.

Build: `conda activate Amuse-env; export SDKROOT=.../MacOSX15.4.sdk; cd GPU2;
make mpi-cpu CXX="$CXX"` → `run_versions/nbody7b.mpi-cpu`.

Mac toolchain notes (HPC x86 won't need these): pass `CXX="$CXX"` (conda clang with libomp;
`/usr/bin/g++` is Apple clang and rejects `-fopenmp`); link with the system linker via
`-B/usr/bin` (conda `ld64-956.6` asserts on NBODY7's huge static COMMON — the same
throwaway-link the AMUSE build avoids via `cpu-objects`); `-lc++` on macOS / `-lstdc++` on
Linux (`uname` switch in `LDMPI_OPT`).

Increment 2 (next) edits the shared `intgrt.omp.f` (decomposition behind `IF(NRANKS.GT.1)`)
and at that point wires `mpi_nbody_stub.o` into the serial + `*-objects` builds so the
shared file's wrapper calls resolve everywhere.

---

## 10. Increment 2 + 3 status — DONE, bit-identical (regular + irregular force)

**Increment 2 — regular force (commit `bfedf67`).** `NBODY_REGF_RANGE` (static
contiguous split of the `NI` regular-due block) + `NBODY_REGF_GATHER` (4×
`MPI_ALLGATHERV`: acc/jerk/phi + neighbour lists) in `mpi_nbody.f`; `intgrt.omp.f`
runs `GPUNB_REGF` per-rank slice → Allgather → replicated `GPUCOR`.

**Increment 3 — irregular force (2026-06-05).** The §1 step-3 block-wide
irregular force `GPUIRR_FIRR_VEC(NXTLEN,NXTLST,GF,GFD)` (line ~393) is now split
the same way: each rank evaluates its contiguous slice of the `NXTLEN` block,
`NBODY_IRRF_GATHER` Allgathers the two `(3,*)` arrays `GF`/`GFD`, then the
replicated `NBINT`/`NBINTP` correction runs on identical inputs. The range split
reuses the generic `NBODY_REGF_RANGE`; the gather mirrors `NBODY_REGF_GATHER` but
with only force + first derivative (no potential / neighbour list — the
irregular phase doesn't rebuild lists). Every rank holds the full, identical
GPUIRR j-particle state (the `GPUIRR_SET_JP`/`PRED_*` calls stay replicated), so
any rank can evaluate any i-particle's irregular force. The secondary
`GPUIRR_FIRR_VEC` inside the regular block (line ~543) stays replicated — it is
part of the replicated `GPUCOR` correction, not a separate decomposition target.

**Acceptance (both increments, pure-MPI `OMP_NUM_THREADS=1`).** On the 13049-body
BSE run (`standalone_version/test_cpu_with_bse`), `mpi_design/poc_validation/run_equiv.sh`
gives **serial cpu == mpi np=1 == np=2 == np=4**, bit-identical ADJUST + END RUN
diagnostics, all ranks reaching END RUN. Caveat unchanged: bit-identicality is a
pure-MPI (one OpenMP thread/rank) guarantee; hybrid MPI+OpenMP inherits the
`start.f` FPOLY2 external-field non-determinism documented at `GPU2/start.f:185`.
