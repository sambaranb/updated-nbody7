#!/bin/bash
# -----------------------------------------------------------------------------
# Path B force-decomposition equivalence test (PoC, increments 2 + 3:
# regular-force AND irregular-force decomposition, both active here).
#
# Runs the standalone mpi-cpu binary under mpirun -np {1,2,4}. Each rank runs
# in its OWN working directory so (a) the per-rank fort.* outputs do not
# collide and (b) input_bse + Fort.10 are present in the CWD, which the
# standalone code REQUIRES (data.f does OPEN(UNIT=222,FILE='input_bse')).
# Each rank reads input_short from stdin and writes run.out / err.out locally.
#
# It then compares the ADJUST / END-RUN physics lines (the deterministic
# diagnostics, with the wall-clock WTIME column stripped) across all ranks of
# all np. With the copy-algorithm decomposition these MUST be bit-identical:
# np=1 executes the original serial path; np=2,4 must reproduce it exactly.
#
# Prereqs:
#   conda activate Amuse-env
#   (cd GPU2 && make mpi-cpu CXX="$CXX")     # -> run_versions/nbody7b.mpi-cpu
#
# Usage:  ./run_equiv.sh [workdir]           # default workdir: /tmp/nb7_mpi_equiv
#
# STATUS 2026-06-05: PASS. In pure-MPI mode (OMP_NUM_THREADS=1, set in wrap.sh)
# np=1/2/4 are bit-identical to each other AND to the serial cpu build, with
# both the regular- and irregular-force decompositions engaged (13049-body BSE
# run to END RUN). The earlier np>1 "hang" was cross-rank OpenMP drift, not an
# Allgatherv bug; forcing OMP_NUM_THREADS=1 removes it (see project memory).
# -----------------------------------------------------------------------------
set -u
ROOT=/Users/sambaran/updated-nbody7
BIN=$ROOT/GPU2/run_versions/nbody7b.mpi-cpu
EX=$ROOT/standalone_version/test_cpu_with_bse
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=${1:-/tmp/nb7_mpi_equiv}

[ -x "$BIN" ] || { echo "missing $BIN  (build: cd GPU2 && make mpi-cpu CXX=\"\$CXX\")"; exit 1; }
command -v mpirun >/dev/null || { echo "mpirun not on PATH -- conda activate Amuse-env"; exit 1; }

rm -rf "$WORK"; mkdir -p "$WORK/src"
cp "$EX/Fort.10" "$EX/input_bse" "$HERE/input_short" "$BIN" "$WORK/src/"

cat > "$WORK/wrap.sh" <<EOF
#!/bin/bash
export OMP_NUM_THREADS=1   # pure-MPI: avoid the start.f OMP non-determinism caveat
export NBODY_RANK0_IO=0    # per-rank-dir validation mode: every rank writes its own
                           # run.out / fort.* so they can be compared cross-rank
                           # (production default is the broad rank-0 I/O guard,
                           # which silences ranks > 0; see run_io_guard.sh)
r=\${OMPI_COMM_WORLD_RANK:-0}
d="$WORK/\${NPTAG}_rank\$r"; rm -rf "\$d"; mkdir -p "\$d"; cd "\$d"
ln -sf "$WORK/src/Fort.10" .; ln -sf "$WORK/src/input_bse" .
exec "$WORK/src/nbody7b.mpi-cpu" < "$WORK/src/input_short" > run.out 2> err.out
EOF
chmod +x "$WORK/wrap.sh"

for NP in 1 2 4; do
    echo "=== running np$NP (timeout 300s) ==="
    NPTAG=np$NP timeout 300 mpirun -x NPTAG=np$NP -np $NP "$WORK/wrap.sh" >/dev/null 2>&1
done

extract(){ grep -E 'ADJUST:|END RUN' "$1" 2>/dev/null | sed -E 's/WTIME.*$//'; }
extract "$WORK/np1_rank0/run.out" > "$WORK/ref.txt"
echo; echo "=== bit-identical check (reference = np1_rank0) ==="
ok=1
for d in "$WORK"/np2_rank* "$WORK"/np4_rank*; do
    [ -d "$d" ] || continue
    extract "$d/run.out" > "$d.cmp"
    if [ ! -s "$d.cmp" ]; then echo "  EMPTY   $(basename "$d")  (hang/crash)"; ok=0
    elif diff -q "$WORK/ref.txt" "$d.cmp" >/dev/null; then echo "  OK      $(basename "$d")"
    else echo "  MISMATCH $(basename "$d")"; ok=0; fi
done
echo
[ $ok -eq 1 ] && echo ">>> PASS: bit-identical across np 1/2/4" \
              || echo ">>> FAIL: see $WORK/np*_rank*/{run,err}.out"
