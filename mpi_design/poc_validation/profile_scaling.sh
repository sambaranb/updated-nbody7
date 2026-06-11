#!/bin/bash
# -----------------------------------------------------------------------------
# Path B force-decomposition SCALING profile (increments 2+3).
#
# Times the standalone mpi-cpu binary in pure-MPI mode (OMP_NUM_THREADS=1) at
# np=1/2/4/8 on the 13049-body BSE IC, REPS times each, reports the min wall
# clock and the speedup vs np=1. Each rank runs in its own dir (per-rank fort.*
# would otherwise collide; input_bse + Fort.10 must be in CWD).
#
# Args:  ./profile_scaling.sh [TCRIT] [DTADJ] [DELTAT] [REPS] [NPLIST]
#   defaults: TCRIT=1.0  DTADJ=0.2  DELTAT=1.0  REPS=2  NPLIST="1 2 4 8"
# Raising DTADJ/DELTAT cuts the (serial, unparallelised) ADJUST/energy/lagr2
# diagnostic overhead, exposing the force-eval scaling.
# -----------------------------------------------------------------------------
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=${ROOT:-$(cd "$HERE/../.." && pwd)}
# BIN override for hosts whose Makefile RUNDIR differs (e.g. run_ampere)
BIN=${BIN:-$(ls "$ROOT"/GPU2/run_*/nbody7b.mpi-cpu 2>/dev/null | head -1)}
EX=$ROOT/standalone_version/test_cpu_with_bse
IN=$ROOT/mpi_design/poc_validation/input_short
TCRIT=${1:-1.0}; DTADJ=${2:-0.2}; DELTAT=${3:-1.0}; REPS=${4:-2}; NPLIST=${5:-"1 2 4 8"}
WORK=/tmp/nb7_prof

[ -x "$BIN" ] || { echo "missing $BIN (build: cd GPU2 && make mpi-cpu CXX=\"\$CXX\")"; exit 1; }
command -v mpirun >/dev/null || { echo "mpirun not on PATH -- conda activate Amuse-env"; exit 1; }

rm -rf "$WORK"; mkdir -p "$WORK/src"
cp "$EX/Fort.10" "$EX/input_bse" "$BIN" "$WORK/src/"
# Patch line 3 fields 4,5,6 (DTADJ DELTAT TCRIT) of the input.
awk -v dtadj="$DTADJ" -v deltat="$DELTAT" -v tcrit="$TCRIT" \
    'NR==3{$4=dtadj;$5=deltat;$6=tcrit} {print}' "$IN" > "$WORK/src/input_run"

cat > "$WORK/wrap.sh" <<EOF
#!/bin/bash
export OMP_NUM_THREADS=1
export NBODY_RANK0_IO=0    # per-rank-dir mode (each rank keeps its own outputs)
r=\${OMPI_COMM_WORLD_RANK:-0}
d="$WORK/\${NPTAG}_rank\$r"; rm -rf "\$d"; mkdir -p "\$d"; cd "\$d"
# lowercase link: gfortran's implicit unit-10 name; Linux FS is case-sensitive
ln -sf "$WORK/src/Fort.10" ./fort.10; ln -sf "$WORK/src/input_bse" .
exec "$WORK/src/nbody7b.mpi-cpu" < "$WORK/src/input_run" > run.out 2> err.out
EOF
chmod +x "$WORK/wrap.sh"

echo "TCRIT=$TCRIT DTADJ=$DTADJ DELTAT=$DELTAT  reps=$REPS  (min wall reported)"
printf "%4s  %10s  %8s  %6s\n" "np" "wall_s" "speedup" "eff%"
REF=""
for NP in $NPLIST; do
    best=""
    for rep in $(seq 1 "$REPS"); do
        t=$(NPTAG=np$NP python3 -c "
import time,subprocess,os
env=dict(os.environ); env['NPTAG']='np$NP'
t=time.time()
subprocess.run(['mpirun','--bind-to','none','-np','$NP','-x','OMP_NUM_THREADS=1','-x','NPTAG','$WORK/wrap.sh'],
               stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,env=env)
print('%.3f'%(time.time()-t))")
        if [ -z "$best" ] || awk "BEGIN{exit !($t<$best)}"; then best=$t; fi
    done
    [ -z "$REF" ] && REF=$best
    sp=$(awk "BEGIN{printf \"%.2f\", $REF/$best}")
    ef=$(awk "BEGIN{printf \"%.0f\", 100*$REF/$best/$NP}")
    printf "%4s  %10s  %7sx  %5s\n" "$NP" "$best" "$sp" "$ef"
done
echo "(check $WORK/np1_rank0/run.out reached END RUN; per-rank dirs kept)"
