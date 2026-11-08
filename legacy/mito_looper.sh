#!/usr/bin/env bash
# mito_pipeline.sh
#
# End-to-end mitochondrial genome assembly and QC pipeline.
#
# For each sample (paired FASTQs in <fastq_dir>), the pipeline:
#   [1] QC_reads.sh      — adapter trim (Trimmomatic PE), read merging (FLASH2),
#                          quality trim, gzip
#   [2] GetOrganelle     — de novo animal_mt assembly from trimmed reads
#   [3] mtGrasp          — reference-guided mitogenome assembly from trimmed reads
#   [4] mito2summary.py  — per-sample run summary + gene completeness table
#                          (run independently for GetOrganelle and mtGrasp)
#   [5] aggregate_mito.sh— cohort-level concatenation of all per-sample CSVs
#
# Usage:
#   bash mito_pipeline.sh [options] <fastq_dir>
#
# FASTQ naming convention (Illumina):
#   <sample>_R1_001.fastq.gz  /  <sample>_R2_001.fastq.gz
#
# Reference FASTA convention for mtGrasp (must be in <fastq_dir> or --ref-dir):
#   <sample>.fasta
#
# See --help for full option list.

set -euo pipefail

##############################################################################
# DEFAULTS
##############################################################################

OUTDIR="mito_out"
THREADS=12
JOBS=1

# --- QC ---
RUN_QC=true
QC_MINLEN=30
TRIMMOMATIC_PATH="/home/ark/bin/Trimmomatic-0.39"   # directory containing adapters/

# --- GetOrganelle ---
RUN_GETORG=true
GETORG_ROUNDS=30
GETORG_DISENTANGLE=200000000
GETORG_MAX_READS="inf"

# --- mtGrasp ---
RUN_MTGRASP=true
MT_CODE=""          # REQUIRED if --mtgrasp is active; NCBI translation table number
REF_DIR=""          # directory containing <sample>.fasta references (default: FASTQ_DIR)

# --- mito2summary + aggregate ---
MITO2SUMMARY=""     # path to mito2summary.py (resolved below if in PATH)
AGGREGATE=""        # path to aggregate_mito.sh (resolved below if in PATH)
EXPECTED_CSV=""     # optional: mitochondrial_gene_average_lengths.csv

# --- Misc ---
SKIP_EXISTING=false
VERBOSE=false

##############################################################################
# USAGE
##############################################################################

usage() {
    cat <<'EOF'

Usage:
  bash mito_pipeline.sh [options] <fastq_dir>

Arguments:
  <fastq_dir>           Directory containing *_R1_001.fastq.gz / *_R2_001.fastq.gz pairs.
                        For mtGrasp, per-sample reference FASTAs (<sample>.fasta) must
                        also be present here unless --ref-dir is specified.

Core options:
  -o, --outdir DIR      Output directory (default: mito_out)
  -t, --threads N       Threads per tool per sample (default: 12)
  -j, --jobs N          Parallel samples via GNU parallel (default: 1)
  --skip-existing       Skip a stage for a sample if its output directory already exists
  -v, --verbose         Extra progress output
  -h, --help            Show this help and exit

QC options:
  --no-qc               Skip QC_reads.sh; raw reads fed directly to assemblers
  --qc-minlen N         SLIDINGWINDOW quality threshold for Trimmomatic (default: 30)
  --trimmomatic-path P  Path to Trimmomatic installation dir (default: /home/ark/bin/Trimmomatic-0.39)

GetOrganelle options:
  --no-getorg           Skip GetOrganelle assembly
  --getorg-rounds N     Maximum extension rounds (default: 30)

mtGrasp options:
  --no-mtgrasp          Skip mtGrasp assembly
  --mt-code N           NCBI translation table (required if mtGrasp is active)
                        Common values: 2=vertebrate, 5=invertebrate, 9=echinoderm, 13=ascidian
  --ref-dir DIR         Directory containing <sample>.fasta references for mtGrasp
                        (default: same as <fastq_dir>)

Post-processing options:
  --mito2summary PATH   Path to mito2summary.py (default: searches PATH)
  --aggregate PATH      Path to aggregate_mito.sh (default: searches PATH)
  --expected-csv FILE   CSV of expected gene lengths (Gene, Average_Length_bp);
                        passed to mito2summary.py (optional)

Examples:
  # Full run, all tools, vertebrate mitogenomes
  bash mito_pipeline.sh --mt-code 2 -o mito_results -t 16 -j 4 fastqs/

  # GetOrganelle only (no reference required)
  bash mito_pipeline.sh --no-mtgrasp -o getorg_results -t 12 fastqs/

  # mtGrasp only, references in a separate directory, skip QC
  bash mito_pipeline.sh --no-getorg --no-qc --mt-code 5 --ref-dir refs/ -o mtgrasp_results fastqs/

  # Skip samples already processed
  bash mito_pipeline.sh --mt-code 2 --skip-existing -o mito_results fastqs/

EOF
    exit 0
}

##############################################################################
# ARGUMENT PARSING
##############################################################################

if [[ $# -eq 0 ]]; then usage; fi

POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--outdir)           OUTDIR="$2";            shift 2 ;;
        -t|--threads)          THREADS="$2";           shift 2 ;;
        -j|--jobs)             JOBS="$2";              shift 2 ;;
        --skip-existing)       SKIP_EXISTING=true;     shift   ;;
        -v|--verbose)          VERBOSE=true;           shift   ;;
        -h|--help)             usage ;;
        --no-qc)               RUN_QC=false;           shift   ;;
        --qc-minlen)           QC_MINLEN="$2";         shift 2 ;;
        --trimmomatic-path)    TRIMMOMATIC_PATH="$2";  shift 2 ;;
        --no-getorg)           RUN_GETORG=false;       shift   ;;
        --getorg-rounds)       GETORG_ROUNDS="$2";     shift 2 ;;
        --no-mtgrasp)          RUN_MTGRASP=false;      shift   ;;
        --mt-code)             MT_CODE="$2";           shift 2 ;;
        --ref-dir)             REF_DIR="$2";           shift 2 ;;
        --mito2summary)        MITO2SUMMARY="$2";      shift 2 ;;
        --aggregate)           AGGREGATE="$2";         shift 2 ;;
        --expected-csv)        EXPECTED_CSV="$2";      shift 2 ;;
        -*) echo "ERROR: Unknown option: $1"; exit 1 ;;
        *)  POSITIONAL+=("$1"); shift ;;
    esac
done

if [[ ${#POSITIONAL[@]} -ne 1 ]]; then
    echo "ERROR: Exactly one positional argument required (<fastq_dir>)."
    echo "       Run with --help for usage."
    exit 1
fi

FASTQ_DIR="${POSITIONAL[0]%/}"
[[ -z "$REF_DIR" ]] && REF_DIR="$FASTQ_DIR"

##############################################################################
# PRECHECKS
##############################################################################

[[ -d "$FASTQ_DIR" ]] || { echo "ERROR: FASTQ directory not found: $FASTQ_DIR"; exit 1; }

# Both assemblers disabled?
if ! $RUN_GETORG && ! $RUN_MTGRASP; then
    echo "ERROR: Both --no-getorg and --no-mtgrasp specified — nothing to do."
    exit 1
fi

# mtGrasp requires translation table
if $RUN_MTGRASP && [[ -z "$MT_CODE" ]]; then
    echo "ERROR: --mt-code is required when running mtGrasp."
    echo "       Common values: 2=vertebrate, 5=invertebrate, 9=echinoderm, 13=ascidian"
    exit 1
fi

# Resolve helper scripts
resolve_script() {
    local name="$1" flag="$2" var="$3"
    if [[ -z "${!var}" ]]; then
        if command -v "$name" &>/dev/null; then
            eval "$var=$(command -v "$name")"
        else
            echo "ERROR: $name not found in PATH. Use $flag to specify the full path."
            exit 1
        fi
    elif [[ ! -f "${!var}" ]]; then
        echo "ERROR: $name not found at: ${!var}"
        exit 1
    fi
}

resolve_script QC_reads.sh      "--no-qc to skip"      _DUMMY_QC   2>/dev/null || true
resolve_script mito2summary.py  "--mito2summary"        MITO2SUMMARY
resolve_script aggregate_mito.sh "--aggregate"          AGGREGATE

# Tool availability (conda-based tools checked at activation time)
for tool in seqtk; do
    command -v "$tool" &>/dev/null || { echo "WARNING: $tool not found in PATH (needed for subsampling, not this pipeline)"; }
done

if [[ $JOBS -gt 1 ]] && ! command -v parallel &>/dev/null; then
    echo "ERROR: GNU parallel not found but -j $JOBS > 1 requested."
    echo "       Install with: conda install -c conda-forge parallel"
    exit 1
fi

[[ -n "$EXPECTED_CSV" && ! -f "$EXPECTED_CSV" ]] && {
    echo "ERROR: --expected-csv file not found: $EXPECTED_CSV"; exit 1
}

##############################################################################
# SETUP OUTPUT DIRS
##############################################################################

GETORG_ROOT="${OUTDIR}/getorganelle"
MTGRASP_ROOT="${OUTDIR}/mtgrasp"
SUMMARY_ROOT="${OUTDIR}/summaries"
LOG_DIR="${OUTDIR}/logs"

mkdir -p "$GETORG_ROOT" "$MTGRASP_ROOT" "$SUMMARY_ROOT" "$LOG_DIR"

PIPELINE_LOG="${OUTDIR}/mito_pipeline.log"
exec > >(tee -a "$PIPELINE_LOG") 2>&1

echo "============================================================"
echo "mito_pipeline.sh"
echo "Start : $(date '+%Y-%m-%d %H:%M:%S')"
echo "FASTQs: $FASTQ_DIR"
echo "Output: $OUTDIR"
echo "Stages: QC=${RUN_QC}  GetOrganelle=${RUN_GETORG}  mtGrasp=${RUN_MTGRASP}"
[[ $RUN_MTGRASP == true ]] && echo "MT code: $MT_CODE"
echo "Threads/job: $THREADS  |  Parallel jobs: $JOBS"
echo "============================================================"

##############################################################################
# HELPERS
##############################################################################

log_step() { echo ""; echo "--- $* ---"; }
log_ok()   { echo "[OK]  $*"; }
log_warn() { echo "[WARN] $*"; }
log_fail() { echo "[FAIL] $*"; }
vlog()     { $VERBOSE && echo "[INFO] $*" || true; }

##############################################################################
# STAGE 1 — QC (per sample, inline)
##############################################################################

run_qc() {
    local SAMPLE="$1" R1="$2" R2="$3" WORKDIR="$4"

    local QC_DONE="${WORKDIR}/${SAMPLE}.notCombined.final_1P.gz"
    if $SKIP_EXISTING && [[ -f "$QC_DONE" ]]; then
        log_ok "QC already done: $SAMPLE"
        return 0
    fi

    log_step "QC: $SAMPLE"
    mkdir -p "$WORKDIR"

    QC_reads.sh \
        "$THREADS" \
        "$SAMPLE" \
        "$WORKDIR" \
        "$R1" "$R2" \
        trimmomatic \
        "$TRIMMOMATIC_PATH" \
        flash2 \
        pigz \
        "$QC_MINLEN" \
        >> "${LOG_DIR}/${SAMPLE}.qc.log" 2>&1

    if [[ -f "${WORKDIR}/${SAMPLE}.notCombined.final_1P.gz" ]]; then
        log_ok "QC complete: $SAMPLE"
    else
        log_fail "QC failed: $SAMPLE (check ${LOG_DIR}/${SAMPLE}.qc.log)"
        return 1
    fi
}

##############################################################################
# STAGE 2 — GetOrganelle (per sample)
##############################################################################

run_getorganelle() {
    local SAMPLE="$1" WORKDIR="$2"

    local OUTD="${GETORG_ROOT}/${SAMPLE}-mitos"
    if $SKIP_EXISTING && [[ -d "$OUTD" ]]; then
        log_ok "GetOrganelle already done: $SAMPLE"
        return 0
    fi

    # Choose reads: QC output if available, else raw
    local R1 R2 U
    if [[ -f "${WORKDIR}/${SAMPLE}.notCombined.final_1P.gz" ]]; then
        R1="${WORKDIR}/${SAMPLE}.notCombined.final_1P.gz"
        R2="${WORKDIR}/${SAMPLE}.notCombined.final_2P.gz"
        U="${WORKDIR}/${SAMPLE}_reads.U.fq.gz"
    else
        R1="${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz"
        R2="${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz"
        U=""
    fi

    log_step "GetOrganelle: $SAMPLE"
    vlog "R1=$R1"
    vlog "R2=$R2"
    [[ -n "$U" ]] && vlog "U=$U"

    eval "$(conda shell.bash hook)"
    conda activate getorganelle

    local CMD=(
        get_organelle_from_reads.py
        -1 "$R1" -2 "$R2"
        -F animal_mt
        -o "$OUTD"
        -t "$THREADS"
        -R "$GETORG_ROUNDS"
        --max-reads "$GETORG_MAX_READS"
        --disentangle-time-limit "$GETORG_DISENTANGLE"
        --reduce-reads-for-coverage inf
    )
    [[ -n "$U" ]] && CMD+=(-u "$U")

    if "${CMD[@]}" >> "${LOG_DIR}/${SAMPLE}.getorg.log" 2>&1; then
        log_ok "GetOrganelle complete: $SAMPLE → $OUTD"
    else
        log_fail "GetOrganelle failed: $SAMPLE (check ${LOG_DIR}/${SAMPLE}.getorg.log)"
    fi

    conda deactivate
}

##############################################################################
# STAGE 3 — mtGrasp (per sample)
##############################################################################

run_mtgrasp() {
    local SAMPLE="$1" WORKDIR="$2"

    local OUTD="${MTGRASP_ROOT}/${SAMPLE}"
    if $SKIP_EXISTING && [[ -d "$OUTD" ]]; then
        log_ok "mtGrasp already done: $SAMPLE"
        return 0
    fi

    local REF="${REF_DIR}/${SAMPLE}.fasta"
    if [[ ! -f "$REF" ]]; then
        log_warn "mtGrasp skipped: no reference found at $REF"
        return 0
    fi

    local R1 R2
    if [[ -f "${WORKDIR}/${SAMPLE}.notCombined.final_1P.gz" ]]; then
        R1="${WORKDIR}/${SAMPLE}.notCombined.final_1P.gz"
        R2="${WORKDIR}/${SAMPLE}.notCombined.final_2P.gz"
    else
        R1="${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz"
        R2="${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz"
    fi

    log_step "mtGrasp: $SAMPLE"
    vlog "R1=$R1  R2=$R2  REF=$REF"

    mkdir -p "$OUTD"

    eval "$(conda shell.bash hook)"
    conda activate mtgrasp

    if mtgrasp.py \
            -r1 "$R1" -r2 "$R2" \
            -r  "$REF" \
            -m  "$MT_CODE" \
            -o  "$OUTD" \
            -t  "$THREADS" \
            >> "${LOG_DIR}/${SAMPLE}.mtgrasp.log" 2>&1; then
        log_ok "mtGrasp complete: $SAMPLE → $OUTD"
    else
        log_fail "mtGrasp failed: $SAMPLE (check ${LOG_DIR}/${SAMPLE}.mtgrasp.log)"
    fi

    conda deactivate
}

##############################################################################
# STAGE 4 — mito2summary (per sample, per assembler)
##############################################################################

run_mito2summary() {
    local SAMPLE="$1" ASSEM_DIR="$2" OUT_PREFIX="$3"

    # Locate required GetOrganelle outputs
    local FASTA HITS_CSV RUN_LOG

    # FASTA: prefer complete, fall back to any path_sequence
    FASTA=$(find "$ASSEM_DIR" -maxdepth 2 -name "*.path_sequence.fasta" \
                | sort | grep -i "complete" | head -1 || true)
    [[ -z "$FASTA" ]] && \
        FASTA=$(find "$ASSEM_DIR" -maxdepth 2 -name "*.path_sequence.fasta" \
                    | sort | head -1 || true)

    HITS_CSV=$(find "$ASSEM_DIR" -maxdepth 2 -name "*.hits.csv" | head -1 || true)
    [[ -z "$HITS_CSV" ]] && \
        HITS_CSV=$(find "$ASSEM_DIR" -maxdepth 2 -name "*.csv" \
                       | grep -v "_run_summary\|_gene_lengths" | head -1 || true)

    RUN_LOG=$(find "$ASSEM_DIR" -maxdepth 2 -name "get_org.log.txt" | head -1 || true)
    [[ -z "$RUN_LOG" ]] && \
        RUN_LOG=$(find "$ASSEM_DIR" -maxdepth 2 -name "*.log.txt" | head -1 || true)

    if [[ -z "$FASTA" || -z "$HITS_CSV" || -z "$RUN_LOG" ]]; then
        log_warn "mito2summary: incomplete outputs for $SAMPLE in $ASSEM_DIR — skipping"
        [[ -z "$FASTA"    ]] && log_warn "  missing: path_sequence.fasta"
        [[ -z "$HITS_CSV" ]] && log_warn "  missing: hits CSV"
        [[ -z "$RUN_LOG"  ]] && log_warn "  missing: log file"
        return 0
    fi

    local CMD=(python3 "$MITO2SUMMARY"
               -c "$HITS_CSV"
               -f "$FASTA"
               -l "$RUN_LOG"
               -o "$OUT_PREFIX")
    [[ -n "$EXPECTED_CSV" ]] && CMD+=(--expected-csv "$EXPECTED_CSV")

    if "${CMD[@]}" >> "${LOG_DIR}/${SAMPLE}.mito2summary.log" 2>&1; then
        log_ok "mito2summary: $SAMPLE ($(basename "$OUT_PREFIX"))"
    else
        log_fail "mito2summary failed: $SAMPLE (check ${LOG_DIR}/${SAMPLE}.mito2summary.log)"
    fi
}

##############################################################################
# PER-SAMPLE WORKER
##############################################################################

process_sample() {
    local R1="$1"
    local SAMPLE
    SAMPLE=$(basename "${R1%_R1_001.fastq.gz}")
    local WORKDIR="${OUTDIR}/qc/${SAMPLE}"

    echo ""
    echo "=============================="
    echo "SAMPLE: $SAMPLE"
    echo "=============================="

    # Stage 1 — QC
    if $RUN_QC; then
        run_qc "$SAMPLE" \
               "${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz" \
               "${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz" \
               "$WORKDIR" || return 0
    fi

    # Stage 2 — GetOrganelle
    if $RUN_GETORG; then
        run_getorganelle "$SAMPLE" "$WORKDIR"
    fi

    # Stage 3 — mtGrasp
    if $RUN_MTGRASP; then
        run_mtgrasp "$SAMPLE" "$WORKDIR"
    fi

    # Stage 4 — mito2summary (one call per assembler that ran)
    log_step "mito2summary: $SAMPLE"
    mkdir -p "${SUMMARY_ROOT}/getorganelle" "${SUMMARY_ROOT}/mtgrasp"

    if $RUN_GETORG && [[ -d "${GETORG_ROOT}/${SAMPLE}-mitos" ]]; then
        run_mito2summary "$SAMPLE" \
            "${GETORG_ROOT}/${SAMPLE}-mitos" \
            "${SUMMARY_ROOT}/getorganelle/${SAMPLE}"
    fi

    if $RUN_MTGRASP && [[ -d "${MTGRASP_ROOT}/${SAMPLE}" ]]; then
        run_mito2summary "$SAMPLE" \
            "${MTGRASP_ROOT}/${SAMPLE}" \
            "${SUMMARY_ROOT}/mtgrasp/${SAMPLE}"
    fi
}

export -f process_sample run_qc run_getorganelle run_mtgrasp run_mito2summary
export -f log_step log_ok log_warn log_fail vlog
export OUTDIR THREADS FASTQ_DIR REF_DIR
export RUN_QC QC_MINLEN TRIMMOMATIC_PATH
export RUN_GETORG GETORG_ROOT GETORG_ROUNDS GETORG_DISENTANGLE GETORG_MAX_READS
export RUN_MTGRASP MTGRASP_ROOT MT_CODE
export MITO2SUMMARY SUMMARY_ROOT EXPECTED_CSV
export SKIP_EXISTING VERBOSE LOG_DIR

##############################################################################
# DISCOVER SAMPLES
##############################################################################

mapfile -t R1_FILES < <(
    find "$FASTQ_DIR" -maxdepth 1 -name "*_R1_001.fastq.gz" | sort
)

N=${#R1_FILES[@]}
if [[ $N -eq 0 ]]; then
    echo "ERROR: No *_R1_001.fastq.gz files found in $FASTQ_DIR"
    exit 1
fi

echo "Found $N sample$([ $N -eq 1 ] && echo '' || echo s)."
for r1 in "${R1_FILES[@]}"; do
    echo "  $(basename "${r1%_R1_001.fastq.gz}")"
done
echo ""

##############################################################################
# DISPATCH (serial or parallel)
##############################################################################

if [[ $JOBS -eq 1 ]]; then
    for R1 in "${R1_FILES[@]}"; do
        process_sample "$R1"
    done
else
    printf '%s\n' "${R1_FILES[@]}" | parallel -j "$JOBS" process_sample {}
fi

##############################################################################
# STAGE 5 — COHORT AGGREGATION
##############################################################################

echo ""
echo "=============================="
echo "STAGE 5 — Cohort aggregation"
echo "=============================="

AGG_ARGS=(--mito2summary "$MITO2SUMMARY")
[[ -n "$EXPECTED_CSV" ]] && AGG_ARGS+=(--expected-csv "$EXPECTED_CSV")

if $RUN_GETORG && [[ -d "${SUMMARY_ROOT}/getorganelle" ]]; then
    bash "$AGGREGATE" \
        "${AGG_ARGS[@]}" \
        --suffix "" \
        -o "${OUTDIR}/cohort_getorganelle" \
        "${SUMMARY_ROOT}/getorganelle" \
        >> "${LOG_DIR}/aggregate_getorganelle.log" 2>&1
    log_ok "GetOrganelle cohort tables → ${OUTDIR}/cohort_getorganelle/"
fi

if $RUN_MTGRASP && [[ -d "${SUMMARY_ROOT}/mtgrasp" ]]; then
    bash "$AGGREGATE" \
        "${AGG_ARGS[@]}" \
        --suffix "" \
        -o "${OUTDIR}/cohort_mtgrasp" \
        "${SUMMARY_ROOT}/mtgrasp" \
        >> "${LOG_DIR}/aggregate_mtgrasp.log" 2>&1
    log_ok "mtGrasp cohort tables → ${OUTDIR}/cohort_mtgrasp/"
fi

##############################################################################
# DONE
##############################################################################

echo ""
echo "============================================================"
echo "Done : $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "Output layout:"
echo "  ${OUTDIR}/qc/                     per-sample QC reads"
echo "  ${OUTDIR}/getorganelle/           GetOrganelle assemblies"
echo "  ${OUTDIR}/mtgrasp/                mtGrasp assemblies"
echo "  ${OUTDIR}/summaries/              per-sample mito2summary CSVs"
echo "  ${OUTDIR}/cohort_getorganelle/    cohort run_summary + gene_lengths"
echo "  ${OUTDIR}/cohort_mtgrasp/         cohort run_summary + gene_lengths"
echo "  ${OUTDIR}/logs/                   per-sample and per-stage logs"
echo "  ${OUTDIR}/mito_pipeline.log       full pipeline log"
echo "============================================================"