#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Alternate Conda installations are noted beside SPECS in configure_belly.py.
exec bash "$SCRIPT_DIR/launch_python.sh" "$SCRIPT_DIR/mito_pipeline.py" "$@"
