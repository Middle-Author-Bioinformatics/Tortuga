#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Environment choices and alternatives are documented beside SPECS in configure_belly.py.
exec bash "$SCRIPT_DIR/launch_python.sh" "$SCRIPT_DIR/pipeline_controller.py" "$@"
