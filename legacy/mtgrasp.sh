#!/usr/bin/env bash
set -euo pipefail

############################################
# DEFAULTS
############################################

OUTDIR="mtgrasp_out"
THREADS=8
JOBS=1

TRANS_TABLE=""

RUN_QC=false
QC_THREADS=8
QC_MINLEN=30

############################################
# ARGUMENT PARSING
############################################

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o) OUTDIR="$2"; shift 2 ;;
    -t) THREADS="$2"; shift 2 ;;
    -j) JOBS="$2"; shift 2 ;;

    --mt-code) TRANS_TABLE="$2"; shift 2 ;;

    --run-qc) RUN_QC=true; shift ;;
    --qc-threads) QC_THREADS="$2"; shift 2 ;;
    --qc-set) QC_MINLEN="$2"; shift 2 ;;

    -h|--help)
      cat <<EOF
Usage:
  mtgrasp.sh --mt-code <N> [options]

Required:
  --mt-code        Translation table (e.g. 5, 9, 13)

Conventions:
  FASTQ: <sample>_R1_001.fastq.gz / _R2_
  REF:   <sample>.fasta

Options:
  -o               Output directory (default: mtgrasp_out)
  -t               Threads per mtGrasp job
  -j               Parallel jobs

QC options:
  --run-qc
  --qc-threads <n>
  --qc-set <minlen>

EOF
      exit 0
      ;;
    *) echo "ERROR: Unknown option $1"; exit 1 ;;
  esac
done

############################################
# PRECHECKS
############################################

[[ -n "$TRANS_TABLE" ]] || { echo "ERROR: --mt-code is required"; exit 1; }
mkdir -p "$OUTDIR"

############################################
# WORKER FUNCTION
############################################

run_one() {
  local R1="$1"
  local SAMPLE="${R1%_R1_001.fastq.gz}"
  local R2="${SAMPLE}_R2_001.fastq.gz"
  local REF="${SAMPLE}.fasta"

  echo "=============================="
  echo ">>> SAMPLE: $SAMPLE"

  [[ -f "$R2" ]] || { echo "SKIP: missing R2"; return; }
  [[ -f "$REF" ]] || { echo "ERROR: missing reference $REF"; return; }

  local WORK="${OUTDIR}/${SAMPLE}"
  mkdir -p "$WORK"

  local IN_R1="$R1"
  local IN_R2="$R2"

  if $RUN_QC; then
    echo ">>> Running QC"

    QC_reads.sh \
      "$QC_THREADS" \
      "$SAMPLE" \
      "$WORK" \
      "$R1" "$R2" \
      trimmomatic /home/ark/bin/Trimmomatic-0.39 flash2 pigz "$QC_MINLEN"

    IN_R1="${WORK}/${SAMPLE}.notCombined.final_1P.gz"
    IN_R2="${WORK}/${SAMPLE}.notCombined.final_2P.gz"

    [[ -f "$IN_R1" && -f "$IN_R2" ]] || {
      echo "ERROR: QC failed for $SAMPLE"
      return
    }
  fi

  echo ">>> Running mtGrasp"
  echo "R1=$IN_R1"
  echo "R2=$IN_R2"
  echo "REF=$REF"
  echo "MT_CODE=$TRANS_TABLE"

  eval "$(conda shell.bash hook)"
  conda activate mtgrasp

  mtgrasp.py \
    -r1 "$IN_R1" \
    -r2 "$IN_R2" \
    -r "$REF" \
    -m "$TRANS_TABLE" \
    -o "$WORK" \
    -t "$THREADS"

    conda deactivate
}

export -f run_one
export OUTDIR THREADS RUN_QC QC_THREADS QC_MINLEN TRANS_TABLE

############################################
# FIND FASTQS + PARALLEL
############################################

echo ">>> Launching mtGrasp jobs"

ls *_R1_001.fastq.gz | \
  parallel -j "$JOBS" \
    run_one {}

echo ">>> DONE"
