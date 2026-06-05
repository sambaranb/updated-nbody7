#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Path B force-decomposition LARGE-N scaling sweep.
#
# Times the standalone mpi-cpu binary in pure-MPI mode (OMP_NUM_THREADS=1) at a
# list of np values, REPS reps each (min wall reported), for one or more
# (label, input-file) physics configs that share an IC dir (Fort.10 + input_bse).
# Reports wall, speedup vs np=1, parallel efficiency, and an Amdahl serial-
# fraction estimate p (= parallel fraction) fit from the np>1 speedups.
# Each MPI rank runs in its own dir (per-rank fort.* would collide; the code
# also OPENs input_bse / Fort.10 from CWD). Success = rank0 reached "END RUN".
#
# Env overrides (all optional):
#   MPIRUN   path to mpirun                (default: Amuse-env mpirun)
#   BIN      mpi-cpu binary                (default: GPU2/run_versions/nbody7b.mpi-cpu)
#   IC_DIR   dir holding Fort.10+input_bse (default: standalone_version/test_cpu_with_bse)
#   NPLIST   space-sep np values           (default: "1 2 4 8")
#   REPS     reps per np, min taken        (default: 3)
#   CONFIGS  "label=path[,label=path...]"  (default: the two session configs)
#   WORK     scratch root                  (default: /tmp/nb7_sweep)
# -----------------------------------------------------------------------------
import os, sys, time, subprocess, shutil, stat

ROOT   = "/Users/sambaran/updated-nbody7"
MPIRUN = os.environ.get("MPIRUN", "/Users/sambaran/miniforge3/envs/Amuse-env/bin/mpirun")
BIN    = os.environ.get("BIN", f"{ROOT}/GPU2/run_versions/nbody7b.mpi-cpu")
IC_DIR = os.environ.get("IC_DIR", f"{ROOT}/standalone_version/test_cpu_with_bse")
NPLIST = [int(x) for x in os.environ.get("NPLIST", "1 2 4 8").split()]
REPS   = int(os.environ.get("REPS", "3"))
WORK   = os.environ.get("WORK", "/tmp/nb7_sweep")
PV     = f"{ROOT}/mpi_design/poc_validation"
DEFAULT_CONFIGS = f"A_bse_tidal={PV}/input_short,B_nobse_notidal={PV}/input_short_notidal"
CONFIGS = [c.split("=", 1) for c in os.environ.get("CONFIGS", DEFAULT_CONFIGS).split(",")]

for need in (MPIRUN, BIN, IC_DIR):
    if not os.path.exists(need):
        sys.exit(f"missing: {need}")

def run_one(label, inp, np_, base):
    """One mpirun -np np_ launch; per-rank dirs under base. Returns (wall_s, ok)."""
    shutil.rmtree(base, ignore_errors=True)
    os.makedirs(base, exist_ok=True)
    wrap = f"{base}/wrap.sh"
    with open(wrap, "w") as f:
        f.write(f"""#!/bin/bash
export OMP_NUM_THREADS=1
r=${{OMPI_COMM_WORLD_RANK:-0}}
d="{base}/rank$r"; rm -rf "$d"; mkdir -p "$d"; cd "$d"
ln -sf "{IC_DIR}/Fort.10" .
[ -e "{IC_DIR}/input_bse" ] && ln -sf "{IC_DIR}/input_bse" .
exec "{BIN}" < "{inp}" > run.out 2> err.out
""")
    os.chmod(wrap, os.stat(wrap).st_mode | stat.S_IEXEC)
    env = dict(os.environ); env["OMP_NUM_THREADS"] = "1"
    cmd = [MPIRUN, "--bind-to", "none", "-np", str(np_),
           "-x", "OMP_NUM_THREADS=1", wrap]
    t0 = time.perf_counter()
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
    wall = time.perf_counter() - t0
    out0 = f"{base}/rank0/run.out"
    ok = os.path.exists(out0) and "END RUN" in open(out0, errors="ignore").read()
    return wall, ok

def amdahl_p(np_, speedup):
    # S = 1/((1-p) + p/np)  ->  p = (1 - 1/S)/(1 - 1/np)
    if np_ <= 1 or speedup <= 0:
        return float("nan")
    return (1 - 1/speedup) / (1 - 1/np_)

print(f"binary : {BIN}")
print(f"IC dir : {IC_DIR}")
print(f"np list: {NPLIST}   reps: {REPS}   (min wall reported, OMP_NUM_THREADS=1)\n")

summary = {}
for label, inp in CONFIGS:
    if not os.path.exists(inp):
        print(f"!! skip {label}: missing input {inp}"); continue
    print(f"================ CONFIG {label}  ({os.path.basename(inp)}) ================")
    print(f"{'np':>3} {'wall_s':>9} {'speedup':>8} {'eff%':>6} {'Amdahl_p':>9} {'END RUN':>8}")
    ref = None; rows = []
    for np_ in NPLIST:
        best, ok_all = None, True
        for rep in range(REPS):
            base = f"{WORK}/{label}/np{np_}_r{rep}"
            wall, ok = run_one(label, inp, np_, base)
            ok_all = ok_all and ok
            if best is None or wall < best:
                best = wall
        if ref is None:
            ref = best
        sp = ref / best
        eff = 100 * sp / np_
        p = amdahl_p(np_, sp)
        rows.append((np_, best, sp, eff, p))
        print(f"{np_:>3} {best:>9.2f} {sp:>7.2f}x {eff:>5.0f} {p:>9.3f} {('yes' if ok_all else 'NO!'):>8}")
    summary[label] = rows
    print()

# side-by-side speedup comparison if >1 config
if len(summary) > 1:
    print("================ SPEEDUP COMPARISON (vs each config's np=1) ================")
    labels = list(summary)
    hdr = f"{'np':>3} " + " ".join(f"{l:>18}" for l in labels)
    print(hdr)
    npset = [r[0] for r in summary[labels[0]]]
    for i, np_ in enumerate(npset):
        cells = []
        for l in labels:
            sp = summary[l][i][2]; eff = summary[l][i][3]
            cells.append(f"{sp:>6.2f}x ({eff:>3.0f}%)")
        print(f"{np_:>3} " + " ".join(f"{c:>18}" for c in cells))
