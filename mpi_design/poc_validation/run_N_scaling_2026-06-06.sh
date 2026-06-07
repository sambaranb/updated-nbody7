#!/bin/bash
# Path B equal-mass N-scaling: all three sweeps (13049 -> 50000 -> 100000),
# startup/integration split + bit-identity, pure-MPI OMP=1, np={1,2,4,8}, min-of-2.
# Decisive question: does integ-only parallel-fraction p climb at N=1e5 or stay ~0.85 flat?
set -u
PV=/Users/sambaran/updated-nbody7/mpi_design/poc_validation
OUT=/tmp/nb7_N_scaling
mkdir -p "$OUT"
export NPLIST="1 2 4 8"
export REPS=2

run() {  # label  input-file  workdir
  local label=$1 input=$2 work=$3
  {
    echo "######## START $label  $(date '+%F %T') ########"
    INPUT="$PV/$input" LABEL="$label" WORK="$work" REPS="$REPS" NPLIST="$NPLIST" \
      python3 "$PV/scale_split.py" 2>&1
    echo "######## DONE  $label  $(date '+%F %T') ########"
  } | tee "$OUT/$label.log"
  cat "$OUT/$label.log" >> "$OUT/master.log"
}

echo "==== N-scaling sweep started $(date '+%F %T') ====" > "$OUT/master.log"
run N13049  input_gen_N13049_eq  /tmp/nb7_scale_N13049
run N50000  input_gen_N50000_eq  /tmp/nb7_scale_N50000
run N100000 input_gen_N100000_eq /tmp/nb7_scale_N100000
echo "==== ALL THREE SWEEPS COMPLETE $(date '+%F %T') ====" | tee -a "$OUT/master.log"
