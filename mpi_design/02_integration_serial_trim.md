# Path 2 — integration serial-trim: strategy sketch

*Status: Step 0 DONE (2026-06-10) — instrumentation shipped and the
attribution measured at N=5e4
(`poc_validation/results_phase_attribution_2026-06-10.log`); measured
re-ranking appended in §5. Next: C1. Builds on the force decomposition
(incr. 2+3), the start-up decomposition (FPOLY2 + FPOLY0), and the rank-0
I/O guard.*

## 1. The budget

Measured (N-scaling sweeps, `results_N_scaling_2026-06-06.log`): the
integration-only parallel fraction is **p ≈ 0.85, flat in N** (13049 → 1e5);
np8 integration speedup 3.9–4.1×; copy-algorithm ceiling at p = 0.85 is
**1/(1−p) ≈ 6.7×**. Path 2 = attack the ~15% replicated/serial residue.

| p (integ) | np8 | np16 | ceiling |
|-----------|-----|------|---------|
| 0.85 (now)| 3.9 | 5.2  | 6.7×    |
| 0.93      | 5.4 | 8.0  | 14×     |
| 0.95      | 5.9 | 9.5  | 20×     |

Two honesty notes on the 15%:
1. The sweeps used a **large DTADJ/DELTAT precisely to suppress the ADJUST
   diagnostics** (see `profile_scaling.sh` header). Production runs adjust far
   more often, so the episodic ADJUST cost (§3, C1) is *under-represented* in
   the measured 15%.
2. Part of the residue is not serial work but the **Allgather time itself**
   (grows with np, shrinks per-rank work). The instrumentation below separates
   the two — comm time is attacked with different tools (fewer/larger
   messages), not by decomposing more loops.

## 2. Step 0 — phase attribution (measure before cutting)

One instrumentation increment, no numerics change: WTIME accumulators around
each phase of the `intgrt.omp.f` cycle and the ADJUST chain, cumulative totals
printed at END RUN (rank 0; per-rank max/min for imbalance). Gate behind an
environment variable (`NBODY_PHASE_TIMERS=1`) so production output is
untouched.

Phases to time, per block step:

| # | phase | site (intgrt.omp.f) | today | scaling |
|---|-------|---------------------|-------|---------|
| P1 | INEXT scheduling + bookkeeping | ~239 | replicated | O(NQ) |
| P2 | prediction `GPUIRR_PRED_ACT/ALL` | 353–358 | replicated | O(NXTLEN) / O(N) |
| P3 | irregular force slice + gather | 403–409 | **decomposed** | O(NXTLEN·⟨NB⟩/np) + comm |
| P4 | irregular corrector `NBINT/NBINTP` | 413–460 | replicated | O(NXTLEN) |
| P5 | `CXVPRED` full predict (reg. step) | 471 | replicated | O(N) |
| P6 | `GPUNB_SEND` j-refresh (reg. step) | 479 | replicated | O(N) |
| P7 | regular force slice + gather | 508–516 | **decomposed** | O(NI·N/np) + comm |
| P8 | regular tail: list post + 2nd `FIRR(NI)` + `GPUCOR` | 530–571 | replicated | O(NI·⟨NB⟩) |
| P9 | `GPUIRR_SET_JP` block update | 622–624 | replicated | O(NXTLEN) |
| P10| `SUBINT` (KS/chain) + MDOT/BSE | 334, 674 | replicated (by contract) | — |
| A1 | ADJUST: `ENERGY2`→`GPUPOT` | adjust.f:18 | replicated | **O(N²) per DTADJ** |
| A2 | ADJUST: `LAGR` + rest + OUTPUT | adjust.f:167 | replicated | O(N log N) |

Run the table at N = 5e4 and 1e5, np = 1 vs 8, Config B, twice: once with the
sweep's large DTADJ and once with a production-like DTADJ. Rank the candidates
below by the np8 numbers. (~15 timer pairs; pure timing, results untouched.)

## 3. Candidates, ranked a priori

All decompositions follow the proven bit-identity recipe: **per-i split with
fully replicated j-state, each rank fills its slice in place, Allgatherv,
then any reduction runs replicated in fixed serial order on the gathered
arrays. Never MPI_Allreduce floating-point.**

### C1 — ENERGY2 / GPUPOT potential at ADJUST  ⟵ first increment
The only remaining **O(N²)** in the integration era, replicated on every rank
at every DTADJ. `lib/gpupot.cpp` is a plain i×j double loop → add a ranged
entry `gpupot_range(i0, ni, n, m, x, phi)` (old entry = full range, all
backends keep their unranged signature working); each rank computes φ for its
`NBODY_REGF_RANGE` i-slice, Allgatherv φ, and the existing PHICOR/POT
summation (energy2.f:45–56) runs replicated in original order → bit-identical.
- Yield: episodic but large; dominates production-DTADJ runs at large N.
- Effort/risk: small/low — exact FPOLY0 pattern.
- **GPU dividend**: the same ranged signature added to `gpupot.gpu.cu` gives
  multi-GPU potentials for free in the future mpi-gpu build.

### C2 — regular-phase tail (list post + 2nd FIRR + GPUCOR), P8
The regular block's replicated tail is force-eval-sized: for every regular
member it post-processes the new neighbour list, **re-evaluates the irregular
force (`GPUIRR_FIRR_VEC(NI)`)** and runs the `GPUCOR` corrector. Decompose
over the *same* i-slice the rank already used for `GPUNB_REGF`: move the
gather *after* the corrector and gather corrected state instead of raw kernel
output (X0, X0DOT, F, FDOT, FR/FRDOT, D-differences, T0R, STEPR, RS + LIST,
LMAX-stride — superset of the FPOLY2+FPOLY0 gathers); `GPUIRR_SET_LIST` /
`SET_JP` then run replicated on identical gathered state.
- Yield: likely the top non-episodic item (P8 ≈ another irregular-force pass
  over every regular block).
- Effort/risk: medium — biggest gather set yet; KS-trigger/neighbour-overflow
  side paths need care. Gate on the Step-0 number for P8.

### C3 — irregular corrector NBINT/NBINTP, P4
Same idea one level down: each rank corrects only its NXTLEN slice (it already
computed those forces), gather the corrected per-particle state (FPOLY2-family
arrays + `ISTAT`/IKS KS-trigger flags as an integer array).
- Yield: smaller flops than C2 (corrector is O(1)-per-neighbour-ish, not a
  force re-evaluation) and it **doubles the per-block-step gather count** —
  at small NXTLEN the added latency can eat the gain. Strictly Step-0-gated;
  attempt only if P4 ≫ P3's comm share.

### C4 — prediction (P2, P5): **do not decompose**
Predicted x,v of *all* particles are inputs to every rank's kernels, so a
split forces an O(N) Allgather per block step that costs more than the ~20
flops/particle prediction it saves. Replicate (NBODY6++GPU does the same).
Revisit only if Step 0 shows P2+P5 > 5% *and* we target single-node
shared-memory MPI only.

### C5 — GPUNB_SEND (P6): record for the GPU build
A replicated O(N) copy into the regular-force library per regular block. On
mpi-cpu it is memcpy-scale (cheap); on **mpi-gpu it becomes a full
host→device transfer of all N per rank per regular block** — the natural
fix there is an incremental j-update API in `gpunb` (as `gpuirr` already has
with SET_JP), not an MPI change. Not a Path-2 item; flagged for Phase G2.

### C6 — LAGR / OUTPUT diagnostics (A2), INEXT/SET_JP/KS (P1, P9, P10)
LAGR feeds RDENS back into physics, so it must stay consistent (decomposing a
sort is not worth it; rank-0 + Bcast gives no wall-time gain since the other
ranks just wait). KS/chain/BSE stay replicated by the design contract.
Deprioritised unless Step 0 says otherwise.

**Proposed sequence:** Step 0 (instrumentation + attribution table) → C1
(GPUPOT, also the first GPU-shared API) → C2 and/or C3 strictly as ranked by
the table, each validated with `run_equiv.sh` bit-identity np 1/2/4/8 before
the next.

## 4. The GPU endgame (why trim CPU serial now)

Ultimate goal: **MPI × GPU** — each rank drives its own GPU
(gpudyn1/gpudyn3-class nodes), regular force on the devices, everything else
on the CPUs.

1. **The decomposition API is already kernel-agnostic.** `NBODY_REGF_RANGE` /
   `NBODY_REGF_GATHER` wrap `GPUNB_REGF` identically whether the backend is
   `gpunb.cpp` (CPU) or `gpunb.velocity.cu` (CUDA). An `mpi-gpu` Makefile
   target is *wiring*, not new parallelisation: link the CUDA libs with
   `mpi_nbody.o`, bind rank → device (`MYRANK % nDevices`, or
   `CUDA_VISIBLE_DEVICES` per rank from the launcher), done. The same holds
   for `mpi-sse`/`mpi-avx`.
2. **Amdahl inverts on the GPU.** The GPU collapses the regular-force time
   (the bulk of today's parallel 85%) by an order of magnitude or more, so on
   an mpi-gpu build the *CPU-side replicated residue* — exactly the P-phases
   above — dominates the wall clock. Every Path-2 trim is therefore a direct
   investment in the GPU build; none of it is throwaway. Conversely, items
   that look minor on mpi-cpu (P6 = GPUNB_SEND, P8's second FIRR) grow in
   relative weight on GPU — another reason Step 0's table is run again on the
   GPU build before final ranking there.
3. **Validation semantics change on GPU.** CPU and GPU force kernels are not
   bit-identical *to each other* (summation order/precision), so the mpi-gpu
   acceptance is: bit-identity **np1-gpu vs npK-gpu** (the decomposition
   guarantee, same as today) + physics-level agreement gpu vs cpu (the
   existing four-host AMUSE sign-off methodology). Pure-MPI (one rank per
   GPU, OMP_NUM_THREADS=1) remains the bit-identical mode.

Suggested phase order across the two tracks:

- **G0 (cheap, anytime):** `mpi-gpu` build target + rank→device binding;
  validate np1-gpu ≡ npK-gpu with the existing harness on gpudyn1/3.
- **Path 2 on mpi-cpu:** Step 0 → C1 → C2/C3 (this document).
- **G1:** re-run the Step-0 attribution on mpi-gpu; re-rank; pull the next
  trim from §3 as indicated.
- **G2:** GPU-specific costs — incremental `gpunb` j-updates (C5), ranged
  `gpupot.gpu.cu` (falls out of C1), multi-GPU-per-node placement.

## 5. Step 0 results — measured re-ranking (2026-06-10, N=5e4)

Full table and analysis: `poc_validation/results_phase_attribution_2026-06-10.log`.
np8 phase totals (wall s), at the sweep ADJUST cadence (`asis`, 2 calls) and a
4× cadence (`adj4x`, 4 calls):

| phase (np8)              | asis | adj4x | note |
|--------------------------|------|-------|------|
| **A2 ADJUST (A1 ENERGY2)** | **8.9 (6.5)** | **14.8 (13.3)** | np-invariant, even ↑ vs np1 |
| P7+P7c regular force      | 6.6  | 4.7   | 7.4× vs np1, imbalance <0.5% |
| P3+P3c irregular force    | 0.9  | 0.7   | 4.7× — small-block latency |
| P2+P4+P9 replicated O(N·) | 2.3  | 1.8   | 1.5–2× *slower* than np1 (mem contention) |
| P12 OUTPUT                | 1.8  | 1.0   | file I/O, np-invariant |
| P8 regular tail           | 0.4  | 0.3   | **C2's guess refuted** |
| everything else           | <0.3 | <0.3  | |
| SUM                       | 21.1 | 23.3  | np1: 61.9 / 53.9 |

**Verdict: C1 is not a "15%" item — it is THE item.** ENERGY2/GPUPOT is 31%
of the np8 wall at the *suppressed* sweep cadence and 57% at the 4× cadence;
everything else on the trim list is ≤ 2.3 s combined. Post-C1 estimate:
phase-sum speedup 2.9→4.0× (asis), 2.3→4.5× (adj4x). C2 measured at 0.3–0.4 s
(the second FIRR covers only regular-block members — the "force-eval-sized"
guess was wrong) and is demoted until the GPU-build re-measurement. C3 stays
marginal (P4 = 1.0 s ≈ what its extra gather would cost). C4's
do-not-decompose call is vindicated (P5 = 0.08 s). P6 GPUNB_SEND = 0.03 s on
CPU — mpi-gpu concern only, as flagged.

**Next increment: C1** — ranged `gpupot` + Allgatherv + replicated
fixed-order summation (FPOLY0 pattern), then re-run this matrix.

## 6. C1 shipped — measured effect (2026-06-10)

Implementation and full numbers: `poc_validation/results_c1_gpupot_2026-06-10.log`.
`gpupot_range` (verbatim per-i kernel, full j-set) in `lib/gpupot.cpp` +
link-compatible fallbacks in the sse/avx/gpu/metal backends +
`NBODY_PHI_GATHER` + the `NRANKS>1` branch in `energy2.f`; serial path
untouched. Bit-identical np1/2/4 (run_equiv) and np8-vs-np1 at N=5e4 on both
cadences — and since np1 takes the *original* GPUPOT path, the equivalence
run directly proves `gpupot_range ≡ gpupot` per particle.

N=5e4 np8: **ENERGY2 6.5→0.83 s (7.9×) / 13.3→1.65 s (8.1×)**; phase-sum
speedup **2.94→4.07×** (sweep cadence) and **2.31→4.66×** (4× cadence) — the
§5 post-C1 estimates hit within 3%. The np8 residue is now ADJUST-rest
(2.5 s), OUTPUT file I/O (1.9 s), and the contended replicated O(N) phases
(2.3 s); the largest *compute* item is the already-decomposed regular force
slice. C2/C3 remain on hold pending the GPU-build re-measurement (§4 G1),
which is the natural next step together with G0 (mpi-gpu wiring).

## 7. Hybrid-mode irregular-force replication — kill the P3c Allgather (2026-06-30)

The Marvin A40 scaling campaign (OMP=8 hybrid runs) revealed that the **largest
single phase at high `np` is not compute but the irregular-force Allgather**
(`P3c`, `NBODY_IRRF_GATHER`). From the N=75k phase decomposition
(`~/nb7_mpi_scaling/phase_tables_N75k_eqmass.txt`, rank-0 cumulative wall s):

| phase (N=75k, OMP=8) | np1 | np2 | np4 | np8 | calls/run |
|---|---|---|---|---|---|
| **P3c irr-force Allgather** | 0.02 | 7.5 | 17.9 | **89.4** | ~348,670 |
| P7c reg-force Allgather | 0.00 | 4.3 | 11.4 | 50.3 | ~27,641 |
| P3 irr-force slice | 46.0 | 58.5 | 33.6 | 20.6 | |
| P7 reg-force slice | 10.6 | 9.3 | 8.4 | 8.1 | |

P3c is the largest phase in the whole run at np8 and exceeds P7c at every `np`,
because the irregular force fires **every block step** (~12.6× more often than
the regular-due steps that trigger P7c) — it is latency-bound, not bandwidth-
bound. This is the mechanism behind the campaign's U-shape scaling collapse.

**Correction to the scaling-paper prose.** `nb7_mpi_scaling/README.md` (and the
matching project memory) describe the exploding collective as *"the per-block-step
Allgather of all-N predicted positions."* That is wrong on two counts: prediction
is **not** decomposed (C4; `CXVPRED`/`GPUIRR_PRED_ALL` are replicated, never
gathered), and the collective that explodes is the irregular force + first
derivative `GF/GFD` of the **current block** (`NXTLEN` members, not all N). The
built-in timer label (`P3c irr force Allgather`) is correct; only the narrative
mislabeled it. Fix the README/paper text before publication.

**The fix (this increment).** When `OMP_NUM_THREADS>1` the irregular force should
not be MPI-decomposed at all: `GPUIRR_FIRR_VEC` is already OpenMP-threaded
(`irrlib/gpuirr.cpp`), so each rank can evaluate the **full block** with its
threads, and under the copy algorithm every rank then holds identical `GF/GFD`
with **zero communication** — `P3c` disappears. The MPI decomposition was only
ever a substitute for the OpenMP parallelism we disabled (`OMP=1`) to keep
bit-identity; once OMP is on in production, OpenMP is the right tool for this
intra-node, latency-bound phase and MPI should carry only the big *divisible*
regular force across nodes.

Implemented as a toggle (`IRR_REPLICATE`, common `/MPIIRR/` in `mpi_nbody.h`),
selected in `NBODY_MPI_INIT`: **auto** = replicate iff `OMP_NUM_THREADS>1`,
decompose at `OMP=1`; **override** `NBODY_IRR_MPI={0 force replicate, 1 force
decompose}`. A rank-0 banner reports the active mode. The integrator branch is
the single site `intgrt.omp.f` ~411–424 (force `MYL0=1,MYLEN=NXTLEN` and skip
`NBODY_IRRF_GATHER`); serial/AMUSE stubs set it `.FALSE.`. This is the **only**
per-step irregular gather — the `FPOLY0`/`FPOLY2` start-up gathers are one-time,
and the secondary `GPUIRR_FIRR_VEC` in the regular block is part of the
replicated `GPUCOR` (never gathered).

**Bit-identity.** `GPUIRR_FIRR_VEC(i)` depends only on the replicated identical
j-state, so at a fixed thread count the replicated full-block `GF/GFD` equals the
decomposed-then-gathered result bit-for-bit. The toggle is a pure compute-vs-comm
trade-off and cannot change numerics. **Validated** (`poc_validation/
results_irr_replicate_2026-06-30.log`, Slurm job 26403445, mpi-cpu, OMP=1):
np2/np4 in *both* forced modes are bit-identical to the np1 serial reference,
with the banner confirming each mode engaged.

**Companion facts (verified).**
- The regular force (`lib/gpunb.cpp`) threads over *i* only, each i's full j-sum
  on one thread → thread-deterministic. So keeping the regular-force MPI
  decomposition under OMP>1 is safe; this change is surgical to the irregular
  phase.
- Removing P3c does **not** by itself make mid-N scale: `P7c` (regular gather,
  50 s at np8) remains, and at N=75k the regular force is a sliver (8 s) not
  worth its gather, so single-GPU + OMP stays fastest until the divisible regular
  chunk grows large (~N=300k, P7=88.7 s). The fix removes a self-inflicted
  artifact and leaves the clean Amdahl story (divisible regular force vs its
  gather) — a *stronger* paper result, not a worse one.
- OMP>1 does **not** add a deadlock gate, and the campaign already proves it: the
  eqmass decks carry `KZ(14)=3` (MW tidal field) and ran at OMP=8 to END RUN,
  validated; the earlier OMP=4 hybrid run (job 26270705) confirmed it directly.
  The MPI Allgathers sit in serial funneled sections between closed `!$omp` loops,
  the per-block kernels/corrector are thread-deterministic (`start.f:195`), and the
  FPOLY2 start-up is gathered (`NBODY_FPOLY_GATHER`) so all ranks stay mutually
  consistent even though `XTRNLD` is not thread-safe with a field. The *only*
  consequence of OMP>1 is loss of bit-reproducibility **across thread counts**
  (`start.f:188–199`) — accepted by design in hybrid production. (An earlier
  caveat called this a "ranks drift and deadlock" risk; that was overcautious and
  is refuted by the runs above.) If strict cross-thread reproducibility *with* a
  field is ever wanted, serialise the FPOLY2 loop (drop its `!$omp`, one-time
  cost). Orthogonal to this increment either way.

**Next:** re-run the N=75k (and N=150k) hybrid sweep with replicate mode (the new
default at OMP=8) to quantify the wall-clock win (P3c→0) and regenerate the
scaling/phase figures; expect the U-shape to flatten markedly toward the clean
regular-force-vs-P7c story.
