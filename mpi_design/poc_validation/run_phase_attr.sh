#!/bin/bash
# -----------------------------------------------------------------------------
# Path 2 Step 0: phase-attribution matrix (mpi_design/02_integration_serial_trim.md).
#
# Runs the instrumented mpi-cpu binary (NBODY_PHASE_TIMERS=1) on the N=5e4
# generated-Plummer input at np in NPLIST, at two ADJUST cadences:
#   asis : input as committed   (DTADJ=0.2, DELTAT=0.2, TCRIT=0.2)
#   adj4x: DTADJ/4 = 0.05       (4x more frequent ENERGY2/LAGR diagnostics)
# Per-rank dirs + NBODY_RANK0_IO=0 so EVERY rank emits its own PHASE TIMERS
# table (rank imbalance visible). One rep per cell (attribution, not timing).
#
# Usage: ./run_phase_attr.sh [workdir]      # default /tmp/nb7_phase_attr
# -----------------------------------------------------------------------------
set -u
ROOT=/Users/sambaran/updated-nbody7
BIN=$ROOT/GPU2/run_versions/nbody7b.mpi-cpu
EX=$ROOT/standalone_version/test_cpu_with_bse
IN=$ROOT/mpi_design/poc_validation/input_gen_N50000_eq
WORK=${1:-/tmp/nb7_phase_attr}
NPLIST=${NPLIST:-"1 8"}

[ -x "$BIN" ] || { echo "missing $BIN"; exit 1; }
command -v mpirun >/dev/null || { echo "mpirun not on PATH"; exit 1; }

rm -rf "$WORK"; mkdir -p "$WORK/src"
cp "$EX/Fort.10" "$EX/input_bse" "$WORK/src/"
cp "$IN" "$WORK/src/input_asis"
awk 'NR==3{$4=0.05} {print}' "$IN" > "$WORK/src/input_adj4x"

cat > "$WORK/wrap.sh" <<EOF
#!/bin/bash
export OMP_NUM_THREADS=1
export NBODY_RANK0_IO=0
export NBODY_PHASE_TIMERS=1
r=\${OMPI_COMM_WORLD_RANK:-0}
d="$WORK/\${TAG}_rank\$r"; rm -rf "\$d"; mkdir -p "\$d"; cd "\$d"
ln -sf "$WORK/src/Fort.10" .; ln -sf "$WORK/src/input_bse" .
exec "$BIN" < "$WORK/src/\$INP" > run.out 2> err.out
EOF
chmod +x "$WORK/wrap.sh"

for CFG in asis adj4x; do
    for NP in $NPLIST; do
        TAG="${CFG}_np${NP}"
        echo "=== $TAG ($(date +%H:%M:%S)) ==="
        TAG=$TAG INP=input_$CFG timeout 3600 mpirun --bind-to none -np $NP \
            -x OMP_NUM_THREADS=1 -x NBODY_RANK0_IO=0 -x NBODY_PHASE_TIMERS=1 \
            -x TAG -x INP "$WORK/wrap.sh" >/dev/null 2>&1
        echo "    exit=$?"
        grep -q 'END RUN' "$WORK/${TAG}_rank0/run.out" \
            && echo "    END RUN ok" || echo "    MISSING END RUN"
    done
done

echo; echo "================ PHASE TABLES (rank 0) ================"
for CFG in asis adj4x; do
    for NP in $NPLIST; do
        TAG="${CFG}_np${NP}"
        echo; echo "--- $TAG ---"
        grep -A18 'PHASE TIMERS' "$WORK/${TAG}_rank0/run.out" 2>/dev/null
    done
done
echo; echo "rank imbalance (np8 asis, P3+P7 slice seconds per rank):"
for r in 0 1 2 3 4 5 6 7; do
    f="$WORK/asis_np8_rank$r/run.out"
    [ -f "$f" ] || continue
    p3=$(grep 'P3  irr force slice' "$f" | awk '{print $5}')
    p7=$(grep 'P7  reg force slice' "$f" | awk '{print $5}')
    echo "  rank$r: P3=$p3  P7=$p7"
done
