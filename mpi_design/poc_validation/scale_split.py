#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Path B force-decomposition LARGE-N scaling with STARTUP/INTEGRATION split.
#
# For one input file, runs the mpi-cpu binary at np in NPLIST (REPS each, pure-
# MPI OMP=1, per-rank dirs). Parses rank0 run.out WTIME stamps (H M S) on the
# ADJUST lines:
#   startup    = WTIME at the FIRST ADJUST (TIME=0.00) -> the one-time START/
#                FPOLY2 setup, which is replicated-serial (NOT MPI-decomposed)
#                and therefore ~np-independent.
#   total      = WTIME at the LAST ADJUST (=TCRIT) -> whole in-code wall.
#   integration= total - startup -> the MPI-decomposed force-loop work, i.e.
#                the PRODUCTION-RELEVANT quantity (startup is amortised at large T).
# Reports integration-only speedup (+Amdahl parallel-fraction p) as the headline,
# raw total in-code speedup as a (startup-polluted) secondary, external wall for
# reference, and verifies bit-identity of the ADJUST physics across np.
#
# Env: INPUT (required), MPIRUN, BIN, IC_DIR, NPLIST (def "1 2 4 8"),
#      REPS (def 2), WORK (def /tmp/nb7_scale_split), LABEL.
# -----------------------------------------------------------------------------
import os, sys, time, subprocess, shutil, stat, re

ROOT   = "/Users/sambaran/updated-nbody7"
MPIRUN = os.environ.get("MPIRUN", "/Users/sambaran/miniforge3/envs/Amuse-env/bin/mpirun")
BIN    = os.environ.get("BIN", f"{ROOT}/GPU2/run_versions/nbody7b.mpi-cpu")
IC_DIR = os.environ.get("IC_DIR", f"{ROOT}/standalone_version/test_cpu_with_bse")
INPUT  = os.environ["INPUT"]
NPLIST = [int(x) for x in os.environ.get("NPLIST", "1 2 4 8").split()]
REPS   = int(os.environ.get("REPS", "2"))
WORK   = os.environ.get("WORK", "/tmp/nb7_scale_split")
LABEL  = os.environ.get("LABEL", os.path.basename(INPUT))

WT    = re.compile(r"WTIME\s*=\s*(\d+)\s+(\d+)\s+(\d+)")
STRIP = re.compile(r"WTIME.*$")

def parse(path):
    """Return (startup_s, total_s, adjust_physics_str, reached_end) or None."""
    if not os.path.exists(path):
        return None
    secs = []; phys = []; end = False
    for ln in open(path, errors="ignore"):
        if "ADJUST:" in ln:
            phys.append(STRIP.sub("", ln).rstrip())
            m = WT.search(ln)
            if m:
                secs.append(int(m.group(1)) * 3600 + int(m.group(2)) * 60 + int(m.group(3)))
        elif "END RUN" in ln:
            end = True
    if not secs:
        return None
    return secs[0], secs[-1], "\n".join(phys), end

def launch(np_, base):
    shutil.rmtree(base, ignore_errors=True); os.makedirs(base, exist_ok=True)
    wrap = f"{base}/wrap.sh"
    with open(wrap, "w") as f:
        f.write(f"""#!/bin/bash
export OMP_NUM_THREADS=1
r=${{OMPI_COMM_WORLD_RANK:-0}}
d="{base}/rank$r"; rm -rf "$d"; mkdir -p "$d"; cd "$d"
ln -sf "{IC_DIR}/Fort.10" .
[ -e "{IC_DIR}/input_bse" ] && ln -sf "{IC_DIR}/input_bse" .
exec "{BIN}" < "{INPUT}" > run.out 2> err.out
""")
    os.chmod(wrap, os.stat(wrap).st_mode | stat.S_IEXEC)
    env = dict(os.environ); env["OMP_NUM_THREADS"] = "1"
    t0 = time.perf_counter()
    subprocess.run([MPIRUN, "--bind-to", "none", "-np", str(np_),
                    "-x", "OMP_NUM_THREADS=1", wrap],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
    return time.perf_counter() - t0

def amdahl(np_, sp):
    if np_ <= 1 or sp <= 0:
        return float("nan")
    return (1 - 1 / sp) / (1 - 1 / np_)

print(f"LABEL={LABEL}  input={os.path.basename(INPUT)}  np={NPLIST}  reps={REPS}  OMP=1", flush=True)
shutil.rmtree(WORK, ignore_errors=True)
best = {}
for np_ in NPLIST:
    b = None
    for rep in range(REPS):
        base = f"{WORK}/np{np_}_r{rep}"
        ext = launch(np_, base)
        p = parse(f"{base}/rank0/run.out")
        if p is None:
            print(f"  np{np_} rep{rep}: NO ADJUST output (crash?)", flush=True); continue
        startup, total, phys, end = p
        integ = total - startup
        print(f"  np{np_} rep{rep}: ext={ext:.1f}s startup={startup}s total={total}s integ={integ}s end={end}", flush=True)
        if b is None or total < b[2]:
            b = (ext, startup, total, integ, phys, end)
    best[np_] = b

ref = best[NPLIST[0]]
print(f"\n==== {LABEL}: min-of-{REPS} ====")
print(f"{'np':>3} {'ext_s':>7} {'startup':>8} {'integ_s':>8} {'integ_sp':>9} {'integ_p':>8} {'tot_sp':>7} {'END':>4}")
for np_ in NPLIST:
    ext, startup, total, integ, phys, end = best[np_]
    isp = ref[3] / integ if integ > 0 else float("nan")
    tsp = ref[2] / total if total > 0 else float("nan")
    ip  = amdahl(np_, isp)
    print(f"{np_:>3} {ext:>7.1f} {startup:>8} {integ:>8} {isp:>8.2f}x {ip:>8.3f} {tsp:>6.2f}x {('yes' if end else 'NO'):>4}")

print("\nbit-identity (ADJUST physics vs np1 rank0):")
refphys = ref[4]
allok = True
for np_ in NPLIST:
    same = best[np_][4] == refphys
    allok = allok and same
    print(f"  np{np_}: {'OK bit-identical' if same else 'MISMATCH'}")
print(f"\n>>> bit-identity across np {NPLIST}: {'PASS' if allok else 'FAIL'}")
