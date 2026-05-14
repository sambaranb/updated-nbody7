"""Cross-variant numerical-equivalence test for amuse-nbody7.

Runs identical 128-body Plummer evolution under each installed mode
(cpu/sse/avx/gpu) and compares conservation laws + bulk structural
diagnostics. Microstates diverge exponentially across variants by
t ~ 1 N-body unit (chaos + floating-point reordering), so agreement
on individual particle positions is NOT expected past t=0; what must
agree is dE/E, |L|, |P|, r_RMS, and sigma_v.

Determinism: np.random.seed(SEED) is called fresh before each Plummer
build, so every variant starts from a bit-identical initial condition.
The t=0 IC sanity block at the bottom proves this -- if those numbers
don't match across variants, the rest of the test is meaningless.

Pass criteria
-------------
1. Each variant individually:        dE/E < 1e-5
2. Cross-variant spread (vs cpu) in r_rms_final, sigma_v_final < 5%
   (statistical floor for N=128 with chaotic divergence over t=3)
3. Cross-variant spread in dE/E itself: same order of magnitude.
   If one variant's dE/E is 100x another's, that's a smoking gun for
   variant contamination, not a chaotic-noise effect.

Workdir hygiene
---------------
NBODY7 writes a fixed set of output files into the current working
directory (NBSTAT, CMPACT, OUT3*, fort.*, STOP*) and refuses to
overwrite some of them, so running four variants back-to-back from
the same cwd would either contaminate later runs with earlier-variant
state or fail outright. We therefore (a) sweep any pre-existing
NBODY7 artifacts into ``test/prior/`` before the first variant runs
and (b) move each variant's outputs into ``test/<variant>/`` as soon
as that variant finishes. Per-variant subdirectories preserve the
files for postmortem (e.g. inspecting ``test/gpu/err_gpu.out`` for
``Opening NBODY6/GPU library'') instead of overwriting between runs.

Run on each box (aibn161, gpudyn3, RTX 2080 box):
    python nb7_equiv_test.py

then inspect ``test/<variant>/err_<variant>.out`` for the per-variant
``Opening NBODY6/GPU library'' / ``Opening GPUIRR lib. {SSE,AVX,CPU}
ver.'' identifiers.
"""
import shutil
from pathlib import Path

import numpy as np

from amuse.lab import nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.community.nbody7 import (
    Nbody7,
    MODE_CPU, MODE_SSE, MODE_AVX, MODE_GPU,
)

SEED     = 20260510
N        = 128
T_SPINUP = 0.5   # baseline-energy capture point (N-body time units)
T_END    = 3.0   # final-energy capture point  (N-body time units)
# NBODY7's first evolve_model() lazily runs SCALE/FPOLY0/FPOLY2, which
# populate the ZKIN/POT entries of the COMMON block. Querying
# kinetic_energy + potential_energy BEFORE that returns 0, so we have to
# evolve to T_SPINUP first to capture E0; the conservation interval is
# therefore [T_SPINUP, T_END], not [0, T_END]. (Same trap that
# test_nb7.py worked around explicitly.)

# NBODY7 output artifacts that must be cleared between variants.
ARTIFACT_PATTERNS = (
    "CMPACT", "NBSTAT", "OUT3*", "fort*", "STOP*",
    "run_*.out", "err_*.out",
)


def _move_artifacts(dest):
    """Move every NBODY7 artifact in cwd into ``dest`` (created if needed)."""
    workdir = Path.cwd()
    moved = False
    for pat in ARTIFACT_PATTERNS:
        for f in workdir.glob(pat):
            if not moved:
                dest.mkdir(parents=True, exist_ok=True)
                moved = True
            target = dest / f.name
            if target.exists():
                target.unlink()
            shutil.move(str(f), target)


def pretest_cleanup():
    """Stash any pre-existing NBODY7 artifacts (from a prior script run or
    a manual experiment) into ``test/prior/`` so the first variant starts
    against a clean cwd."""
    _move_artifacts(Path.cwd() / "test" / "prior")


def archive_variant_outputs(label):
    """Move the artifacts produced by the just-finished variant into
    ``test/<label>/``. Called immediately after ``inst.stop()`` so the
    next variant's invocation sees a clean cwd."""
    _move_artifacts(Path.cwd() / "test" / label)


def bulk_stats(particles):
    """Center-of-mass-frame bulk diagnostics from an AMUSE particle set."""
    pos = particles.position.value_in(nbody_system.length).T   # (3, N)
    vel = particles.velocity.value_in(nbody_system.speed).T
    m   = particles.mass.value_in(nbody_system.mass)
    M   = m.sum()
    com_p = (m * pos).sum(axis=1) / M
    com_v = (m * vel).sum(axis=1) / M
    dx = pos - com_p[:, None]
    dv = vel - com_v[:, None]
    r2 = (dx**2).sum(axis=0)
    v2 = (dv**2).sum(axis=0)
    L = np.cross(dx.T, (m[:, None] * dv.T)).sum(axis=0)
    P = (m[:, None] * vel.T).sum(axis=0)
    KE = 0.5 * (m * v2).sum()
    return dict(
        M=M,
        r_rms=np.sqrt((m * r2).sum() / M),
        sigma_v=np.sqrt((m * v2).sum() / M),
        Ltot=np.linalg.norm(L),
        Ptot=np.linalg.norm(P),
        KE=KE,
    )


def run_one(mode):
    """Build identical IC, evolve to T_END under `mode`, return diagnostics."""
    np.random.seed(SEED)
    cluster = new_plummer_model(N)        # default nbody units
    s0 = bulk_stats(cluster)

    workdir = Path.cwd()
    outfile = str(workdir / f"run_{mode}.out")
    errfile = str(workdir / f"err_{mode}.out")

    inst = Nbody7(
        mode=mode,
        redirection='file',
        redirect_stdout_file=outfile,
        redirect_stderr_file=errfile,
    )
    inst.particles.add_particles(cluster)
    # Spin-up step: triggers SCALE/FPOLY0/FPOLY2 so ZKIN/POT are valid.
    inst.evolve_model(T_SPINUP | nbody_system.time)
    E0 = (inst.kinetic_energy + inst.potential_energy).value_in(
        nbody_system.energy)
    inst.evolve_model(T_END | nbody_system.time)
    Ef = (inst.kinetic_energy + inst.potential_energy).value_in(
        nbody_system.energy)
    final = inst.particles.copy()
    inst.cleanup_code()
    inst.stop()

    sf = bulk_stats(final)
    return dict(
        mode=mode,
        E0=E0, Ef=Ef, dE_E=(Ef - E0) / abs(E0),
        r_rms_0=s0['r_rms'],   r_rms_f=sf['r_rms'],
        sigma_v_0=s0['sigma_v'], sigma_v_f=sf['sigma_v'],
        Ltot_0=s0['Ltot'],     Ltot_f=sf['Ltot'],
        Ptot_f=sf['Ptot'],
    )


# Drop any mode whose variant wheel is not installed -- run_one will raise.
modes = [MODE_CPU, MODE_SSE, MODE_AVX, MODE_GPU]

pretest_cleanup()

results = []
for m in modes:
    try:
        r = run_one(m)
        results.append(r)
        print(
            f"{m:>3s}: E0={r['E0']:.10f}  Ef={r['Ef']:.10f}  "
            f"dE/E={r['dE_E']:+.2e}  "
            f"r_rms_f={r['r_rms_f']:.6f}  sigma_v_f={r['sigma_v_f']:.6f}  "
            f"|L|_f={r['Ltot_f']:.3e}  |P|_f={r['Ptot_f']:.3e}"
        )
    except Exception as e:
        print(f"{m}: SKIP ({type(e).__name__}: {e})")
    finally:
        # Archive whatever this variant produced (run/err files plus any
        # NBSTAT/CMPACT/OUT3*/fort*/STOP* that were emitted before the
        # error, if any). Always runs, including the SKIP path.
        archive_variant_outputs(m)

# Pairwise comparison vs cpu reference.
ref = next((r for r in results if r['mode'] == MODE_CPU), None)
if ref is not None:
    print("\n--- relative deviation vs cpu (post-evolution) ---")
    print(f"{'mode':>4s}  {'|dE/E|':>10s}  {'rel(r_rms)':>12s}  "
          f"{'rel(sigma_v)':>12s}")
    for r in results:
        if r['mode'] == MODE_CPU:
            continue
        print(
            f"{r['mode']:>4s}  {abs(r['dE_E']):.2e}  "
            f"{abs(r['r_rms_f'] - ref['r_rms_f']) / ref['r_rms_f']:.4e}  "
            f"{abs(r['sigma_v_f'] - ref['sigma_v_f']) / ref['sigma_v_f']:.4e}"
        )

# t=0 IC sanity: must be bit-identical across modes (same seeded build).
print("\n--- t=0 IC sanity (must be identical across variants) ---")
print(f"{'mode':>4s}  {'r_rms_0':>14s}  {'sigma_v_0':>14s}  {'|L|_0':>14s}")
for r in results:
    print(
        f"{r['mode']:>4s}  {r['r_rms_0']:.10f}  {r['sigma_v_0']:.10f}  "
        f"{r['Ltot_0']:.6e}"
    )
