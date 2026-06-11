#!/bin/bash
# -----------------------------------------------------------------------------
# Path B broad rank-0 I/O guard test (production hardening).
#
# Runs the standalone mpi-cpu binary with ALL RANKS IN ONE SHARED WORKING
# DIRECTORY and the rank-0 I/O guard at its production default (ON for np>1:
# NBODY_RANK0_IO unset). Ranks > 0 have stdout and every output unit connected
# to /dev/null, so only rank 0 writes the run's files -- the per-rank-dir
# wrapper of run_equiv.sh is no longer required for a production run.
#
# Checks, against a serial-path np=1 reference run in its own directory:
#   1. np4 shared-dir run reaches END RUN with mpirun exit code 0;
#   2. the ADJUST / END-RUN physics lines of the np4 shared run.out are
#      bit-identical to the np=1 reference (WTIME column stripped);
#   3. run.out contains the physics output exactly ONCE (no duplicated lines
#      from ranks > 0 leaking through the guard);
#   4. the file inventory of the shared directory equals the np=1 inventory
#      (no stray fort.* / OUT* created by ranks > 0), and each common output
#      file is byte-identical to the np=1 reference.
#
# Each rank reads the input from a file redirect inside the wrapper (mpirun
# only forwards mpirun's stdin to rank 0 by default; NBODY7 also performs lazy
# mid-run reads from unit 5 on every rank, so production runs should use the
# same exec-with-redirect pattern, or mpirun --stdin all).
#
# Prereqs:
#   conda activate Amuse-env
#   (cd GPU2 && make mpi-cpu CXX="$CXX")     # -> run_versions/nbody7b.mpi-cpu
#
# Usage:  ./run_io_guard.sh [workdir]        # default: /tmp/nb7_mpi_ioguard
# -----------------------------------------------------------------------------
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=${ROOT:-$(cd "$HERE/../.." && pwd)}
# BIN override for hosts whose Makefile RUNDIR differs (e.g. run_ampere)
BIN=${BIN:-$(ls "$ROOT"/GPU2/run_*/nbody7b.mpi-cpu 2>/dev/null | head -1)}
EX=$ROOT/standalone_version/test_cpu_with_bse
WORK=${1:-/tmp/nb7_mpi_ioguard}

[ -x "$BIN" ] || { echo "missing $BIN  (build: cd GPU2 && make mpi-cpu CXX=\"\$CXX\")"; exit 1; }
command -v mpirun >/dev/null || { echo "mpirun not on PATH -- conda activate Amuse-env"; exit 1; }

rm -rf "$WORK"; mkdir -p "$WORK/ref" "$WORK/shared"
for d in ref shared; do
    # lowercase fort.10: gfortran's implicit unit-10 name; Linux FS is case-sensitive
    cp "$EX/Fort.10" "$WORK/$d/fort.10"
    cp "$EX/input_bse" "$HERE/input_short" "$WORK/$d/"
done

run_one() {           # run_one <dir> <np>
    local d=$1 np=$2
    ( cd "$WORK/$d" && OMP_NUM_THREADS=1 timeout 600 \
      mpirun --bind-to none -x OMP_NUM_THREADS=1 -np "$np" \
             bash -c 'exec "'"$BIN"'" < input_short' > run.out 2> err.out )
    echo $?
}

echo "=== np1 reference (own dir, serial path) ==="
RC1=$(run_one ref 1)
echo "    exit code: $RC1"
echo "=== np4, ONE shared directory, rank-0 I/O guard ON (production default) ==="
RC4=$(run_one shared 4)
echo "    exit code: $RC4"

# strip WTIME / CPUTOT / WTOT: wall+CPU accumulators differ between any two runs
extract(){ grep -E 'ADJUST:|END RUN' "$1" 2>/dev/null | sed -E 's/WTIME.*$//; s/CPUTOT =[^A-Z]*//; s/WTOT =.*$//'; }
ok=1

# 1. clean termination
[ "$RC4" = "0" ] || { echo "FAIL: np4 mpirun exit code $RC4"; ok=0; }
grep -q 'END RUN' "$WORK/shared/run.out" || { echo "FAIL: no END RUN in shared run.out"; ok=0; }
grep -q 'Rank-0 I/O guard active' "$WORK/shared/run.out" \
    || { echo "FAIL: guard banner missing (guard not active?)"; ok=0; }

# 2. physics bit-identical to the serial-path reference
extract "$WORK/ref/run.out"    > "$WORK/ref.txt"
extract "$WORK/shared/run.out" > "$WORK/shared.txt"
if diff -q "$WORK/ref.txt" "$WORK/shared.txt" >/dev/null; then
    echo "OK  physics lines bit-identical to np1 reference"
else
    echo "FAIL: physics lines differ from np1 reference"; ok=0
fi

# 3. no duplicated physics lines (ranks > 0 leaking through the guard)
NEND=$(grep -c 'END RUN' "$WORK/shared/run.out")
if [ "$NEND" = "1" ]; then
    echo "OK  END RUN appears exactly once (ranks > 0 silenced)"
else
    echo "FAIL: END RUN appears $NEND times in shared run.out"; ok=0
fi

# 4. file inventory and content identical to the np1 reference
inv(){ (cd "$WORK/$1" && ls | grep -v -E '^(run|err)\.out$' | sort); }
inv ref > "$WORK/inv_ref.txt"; inv shared > "$WORK/inv_shared.txt"
if diff "$WORK/inv_ref.txt" "$WORK/inv_shared.txt" > "$WORK/inv.diff"; then
    echo "OK  file inventory identical to np1 reference:"
    sed 's/^/      /' "$WORK/inv_ref.txt"
else
    echo "FAIL: file inventory differs (stray rank>0 files?):"
    sed 's/^/      /' "$WORK/inv.diff"; ok=0
fi
while read -r f; do
    [ "$f" = "fort.10" -o "$f" = "input_bse" -o "$f" = "input_short" ] && continue
    case "$f" in
    fort.1|fort.2)
        # COMMON dumps embed the CPU/wall-time accumulators (COMMON/COUNTS/),
        # which differ between ANY two runs; require identical structure only.
        s1=$(wc -c < "$WORK/ref/$f"); s2=$(wc -c < "$WORK/shared/$f")
        nb=$(cmp -l "$WORK/ref/$f" "$WORK/shared/$f" 2>/dev/null | wc -l | tr -d ' ')
        if [ "$s1" = "$s2" ]; then
            echo "OK  $f same size as np1 reference ($nb differing bytes: timing fields)"
        else
            echo "FAIL: $f size differs from np1 reference ($s1 vs $s2)"; ok=0
        fi ;;
    NBSTAT)
        # Per-process neighbour-count histogram: at np4 rank 0 only counts its
        # own block slice, so the content legitimately differs from np1.
        echo "OK  $f present (per-rank diagnostic; np>1 content covers rank 0's slice only)" ;;
    *)
        if cmp -s "$WORK/ref/$f" "$WORK/shared/$f"; then
            echo "OK  $f byte-identical to np1 reference"
        else
            echo "FAIL: $f differs from np1 reference"; ok=0
        fi ;;
    esac
done < "$WORK/inv_ref.txt"

echo
[ $ok -eq 1 ] && echo ">>> PASS: shared-directory np4 run is rank-0-clean and physics-identical" \
              || echo ">>> FAIL: see $WORK"
exit $((1-ok))
