#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ENTRY_SCRIPT="$1"
shift
[[ -f "$ENTRY_SCRIPT" ]] || { echo "ERROR: Missing controller: $ENTRY_SCRIPT" >&2; exit 1; }

# Start the controller in itsme automatically; individual stages select their own
# environments. Alternatives are documented in configure_belly.py beside SPECS:
# GetOrganelle: getorganelle first, mito-getorganelle also relevant.
# mtGrasp: mito-mtgrasp first, mtgrasp also relevant.
# SPAdes: itsme first, mito-spades also relevant.
if [[ -n "${MITO_PYTHON:-}" ]]; then
    exec "$MITO_PYTHON" "$ENTRY_SCRIPT" "$@"
fi

MITO_ENV_ROOT="${MITO_ENV_ROOT:-/home/ark/miniconda3/envs}"
MITO_CONDA_BIN="${MITO_CONDA_EXE:-${MITO_ENV_ROOT%/envs}/bin/conda}"
if [[ ! -x "$MITO_CONDA_BIN" && -n "${CONDA_EXE:-}" ]]; then
    MITO_CONDA_BIN="$CONDA_EXE"
fi
if [[ -x "$MITO_CONDA_BIN" && -x "$MITO_ENV_ROOT/itsme/bin/python" ]]; then
    exec "$MITO_CONDA_BIN" run --no-capture-output -p "$MITO_ENV_ROOT/itsme" \
        python "$ENTRY_SCRIPT" "$@"
fi
# For another host, an already activated Python environment can launch the
# controller; --env-root/--conda or MITO_ENV_ROOT/MITO_CONDA_EXE select tool roots.
exec python3 "$ENTRY_SCRIPT" "$@"
