#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Path B force-decomposition bit-identical EQUIVALENCE check (generalized).
#
# Runs the standalone mpi-cpu binary at np=1/2/4 for a given input file and a
# given OMP_NUM_THREADS, each rank in its own dir, then compares the
# deterministic physics diagnostics (ADJUST: / END RUN lines, WTIME column
# stripped) of every rank against np1_rank0. With the copy-algorithm force
# decomposition these MUST be bit-identical.
#
# Two regimes:
#   OMP=1  pure-MPI  -> bit-identical for ANY input (proven for Config A).
#   OMP>1  hybrid    -> bit-identical ONLY if start.f FPOLY2 is thread-
#                       deterministic, i.e. NO external field (KZ(14)=0). This
#                       is the extra property the no-tidal Config B should gain.
#
# Env overrides: MPIRUN, BIN, IC_DIR, INPUT (required-ish), OMP (default 1),
#                NPLIST (default "1 2 4"), WORK (default /tmp/nb7_equiv_gen).
# -----------------------------------------------------------------------------
import os, sys, subprocess, shutil, stat, re

ROOT   = "/Users/sambaran/updated-nbody7"
MPIRUN = os.environ.get("MPIRUN", "/Users/sambaran/miniforge3/envs/Amuse-env/bin/mpirun")
BIN    = os.environ.get("BIN", f"{ROOT}/GPU2/run_versions/nbody7b.mpi-cpu")
IC_DIR = os.environ.get("IC_DIR", f"{ROOT}/standalone_version/test_cpu_with_bse")
INPUT  = os.environ.get("INPUT", f"{ROOT}/mpi_design/poc_validation/input_short_notidal")
OMP    = os.environ.get("OMP", "1")
NPLIST = [int(x) for x in os.environ.get("NPLIST", "1 2 4").split()]
WORK   = os.environ.get("WORK", "/tmp/nb7_equiv_gen")
TIMEOUT = int(os.environ.get("TIMEOUT", "240"))   # per-launch guard (OMP>1 + tidal can deadlock)

WTIME = re.compile(r"WTIME.*$")
def extract(path):
    if not os.path.exists(path):
        return None
    lines = []
    for ln in open(path, errors="ignore"):
        if "ADJUST:" in ln or "END RUN" in ln:
            lines.append(WTIME.sub("", ln).rstrip())
    return "\n".join(lines)

def launch(np_, base):
    shutil.rmtree(base, ignore_errors=True); os.makedirs(base, exist_ok=True)
    wrap = f"{base}/wrap.sh"
    with open(wrap, "w") as f:
        f.write(f"""#!/bin/bash
export OMP_NUM_THREADS={OMP}
r=${{OMPI_COMM_WORLD_RANK:-0}}
d="{base}/rank$r"; rm -rf "$d"; mkdir -p "$d"; cd "$d"
ln -sf "{IC_DIR}/Fort.10" .
[ -e "{IC_DIR}/input_bse" ] && ln -sf "{IC_DIR}/input_bse" .
exec "{BIN}" < "{INPUT}" > run.out 2> err.out
""")
    os.chmod(wrap, os.stat(wrap).st_mode | stat.S_IEXEC)
    env = dict(os.environ); env["OMP_NUM_THREADS"] = OMP
    try:
        subprocess.run([MPIRUN, "--bind-to", "none", "-np", str(np_),
                        "-x", f"OMP_NUM_THREADS={OMP}", wrap],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       env=env, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        print(f"  (np{np_} exceeded {TIMEOUT}s -> killed; likely OMP drift deadlock)")
        subprocess.run(["pkill", "-9", "-f", "nbody7b.mpi-cpu"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

print(f"input: {os.path.basename(INPUT)}   OMP_NUM_THREADS={OMP}   np={NPLIST}")
shutil.rmtree(WORK, ignore_errors=True)
for np_ in NPLIST:
    launch(np_, f"{WORK}/np{np_}")

ref = extract(f"{WORK}/np{NPLIST[0]}/rank0/run.out")
if not ref:
    sys.exit(">>> FAIL: reference np1_rank0 produced no ADJUST/END RUN (crash?)")
ok = True
for np_ in NPLIST:
    for r in range(np_):
        d = f"{WORK}/np{np_}/rank{r}"
        cmp = extract(f"{d}/run.out")
        tag = f"np{np_}_rank{r}"
        if cmp is None or cmp == "":
            print(f"  EMPTY    {tag}  (hang/crash)"); ok = False
        elif cmp == ref:
            print(f"  OK       {tag}")
        else:
            print(f"  MISMATCH {tag}"); ok = False
end = "END RUN" in ref
print()
print(f">>> {'PASS' if ok else 'FAIL'}: bit-identical across np {NPLIST} at OMP={OMP}"
      f"   (reference reached END RUN: {end})")
sys.exit(0 if ok else 1)
