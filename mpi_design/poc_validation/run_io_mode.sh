#!/bin/bash
# -----------------------------------------------------------------------------
# G2 decisive experiment: I/O-mode attribution of the np8 ADJUST-residue blow-up
# (mpi_design/02_integration_serial_trim.md; see project memory nbody7-path-b-mpi).
#
# The earlier intermediate-binary pass folded MYDUMP/ESCAPE/CHECK INTO the A2r
# residue (7-term formula) and saw A2r explode ~12x from np1->np8. The current
# binary (phase_timers.inc NPHT=24, fcce9d4) breaks A2d (MYDUMP COMMON save)
# out as its own phase. This harness runs the N=5e4 input at np1 and np8 under
# BOTH I/O modes so the blow-up can be attributed to a single sub-phase and the
# production guard tested:
#
#   io0 : NBODY_RANK0_IO=0, per-rank working dirs  -> EVERY rank writes the full
#         COMMON dump (8-way shared-disk contention).  Matches the recorded run.
#   io1 : NBODY_RANK0_IO=1, ONE shared working dir   -> only rank 0 writes
#         (broad rank-0 guard = production default for np>1).
#
# Hypothesis: under io0 the np8 blow-up lands in A2d (MYDUMP); under io1 A2d
# collapses to ~np1 because ranks>0 are silenced -> "nothing to cut in ADJUST,
# the residue is rank-0 restart-dump + OUTPUT I/O, np-invariant in production".
#
# Defaults to the mpi-gpu binary (P7 collapses on the A40 -> ADJUST is the wall
# AND runs are fast). Override:  BIN=.../nbody7b.mpi-cpu  NPLIST="1 8"  ./run_io_mode.sh
#
# Usage: ./run_io_mode.sh [workdir]            # default /tmp/nb7_g2_iomode
# -----------------------------------------------------------------------------
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=${ROOT:-$(cd "$HERE/../.." && pwd)}
# default to mpi-gpu; BIN env overrides (e.g. the mpi-cpu cross-check)
BIN=${BIN:-$(ls "$ROOT"/GPU2/run_*/nbody7b.mpi-gpu 2>/dev/null | head -1)}
EX=$ROOT/standalone_version/test_cpu_with_bse
IN=$ROOT/mpi_design/poc_validation/input_gen_N50000_eq
WORK=${1:-/tmp/nb7_g2_iomode}
NPLIST=${NPLIST:-"1 8"}
CFGLIST=${CFGLIST:-"asis adj4x"}

[ -x "$BIN" ] || { echo "missing $BIN"; exit 1; }
command -v mpirun >/dev/null || { echo "mpirun not on PATH"; exit 1; }
echo "BIN = $BIN"

rm -rf "$WORK"; mkdir -p "$WORK/src"
cp "$EX/Fort.10" "$EX/input_bse" "$WORK/src/"
cp "$IN" "$WORK/src/input_asis"
awk 'NR==3{$4=0.05} {print}' "$IN" > "$WORK/src/input_adj4x"

# io0: per-rank dirs, NBODY_RANK0_IO=0 -> every rank writes its own dump.
run_io0() {                       # run_io0 <tag> <inp> <np>
    local tag=$1 inp=$2 np=$3
    cat > "$WORK/wrap0.sh" <<EOF
#!/bin/bash
r=\${OMPI_COMM_WORLD_RANK:-0}
d="$WORK/${tag}_rank\$r"; rm -rf "\$d"; mkdir -p "\$d"; cd "\$d"
ln -sf "$WORK/src/Fort.10" ./fort.10; ln -sf "$WORK/src/input_bse" .
exec "$BIN" < "$WORK/src/$inp" > run.out 2> err.out
EOF
    chmod +x "$WORK/wrap0.sh"
    timeout 1800 mpirun --bind-to none -np "$np" \
        -x OMP_NUM_THREADS=1 -x NBODY_RANK0_IO=0 -x NBODY_PHASE_TIMERS=1 \
        "$WORK/wrap0.sh" >/dev/null 2>&1
    local rc=$?
    grep -q 'END RUN' "$WORK/${tag}_rank0/run.out" 2>/dev/null \
        && echo "    exit=$rc  (END RUN ok)" || echo "    exit=$rc  (NO END RUN)"
}

# io1: ONE shared dir, NBODY_RANK0_IO=1 -> only rank 0 writes (guard on).
run_io1() {                       # run_io1 <tag> <inp> <np>
    local tag=$1 inp=$2 np=$3
    local d="$WORK/${tag}"; rm -rf "$d"; mkdir -p "$d"
    cp "$WORK/src/Fort.10" "$d/fort.10"; cp "$WORK/src/input_bse" "$d/"
    ( cd "$d" && timeout 1800 mpirun --bind-to none -np "$np" \
        -x OMP_NUM_THREADS=1 -x NBODY_RANK0_IO=1 -x NBODY_PHASE_TIMERS=1 \
        bash -c 'exec "'"$BIN"'" < "'"$WORK/src/$inp"'"' > run.out 2> err.out )
    local rc=$?
    grep -q 'END RUN' "$d/run.out" 2>/dev/null \
        && echo "    exit=$rc  (END RUN ok)" || echo "    exit=$rc  (NO END RUN)"
}

for CFG in $CFGLIST; do
    for NP in $NPLIST; do
        echo "=== io0 ${CFG}_np${NP} ($(date +%H:%M:%S)) ==="; run_io0 "io0_${CFG}_np${NP}" "input_$CFG" "$NP"
        echo "=== io1 ${CFG}_np${NP} ($(date +%H:%M:%S)) ==="; run_io1 "io1_${CFG}_np${NP}" "input_$CFG" "$NP"
    done
done

# rank0 table lives in <tag>_rank0/run.out (io0) or <tag>/run.out (io1)
tbl(){ [ -f "$WORK/$1_rank0/run.out" ] && echo "$WORK/$1_rank0/run.out" || echo "$WORK/$1/run.out"; }

echo; echo "================ ADJUST sub-phase tables (rank 0) ================"
for CFG in $CFGLIST; do for NP in $NPLIST; do for IO in io0 io1; do
    tag="${IO}_${CFG}_np${NP}"; f=$(tbl "$tag")
    echo; echo "--- $tag ---"
    grep -A28 'PHASE TIMERS' "$f" 2>/dev/null \
      | grep -E 'A1 |A2 |A2o|A2e|A2d|A2k|A2r|SUM|P7 |P12' || echo "  (no table)"
done; done; done

echo; echo "================ KEY: A2d MYDUMP & A2r residue, io0 vs io1 ================"
printf "%-18s %10s %10s %10s\n" "cell" "A2d(MYDUMP)" "A2r(resid)" "A2(total)"
for CFG in $CFGLIST; do for NP in $NPLIST; do for IO in io0 io1; do
    tag="${IO}_${CFG}_np${NP}"; f=$(tbl "$tag")
    a2d=$(grep 'A2d MYDUMP' "$f" 2>/dev/null | awk '{print $5}')
    a2r=$(grep 'A2r ADJUST' "$f" 2>/dev/null | awk '{print $5}')
    a2=$(grep 'A2  ADJUST'  "$f" 2>/dev/null | awk '{print $6}')
    printf "%-18s %10s %10s %10s\n" "${CFG}_np${NP}_${IO}" "${a2d:-NA}" "${a2r:-NA}" "${a2:-NA}"
done; done; done
