#!/bin/bash
# ==============================================================================
# Bash wrapper pro vyhodnocení výsledků EasyDock screeningu
#
# Použití:
#   ./evaluate_docking_results.sh
#   ./evaluate_docking_results.sh -r results_docking -o reports_docking
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Zamezení načítání starých nekompatibilních balíčků z ~/.local
export PYTHONNOUSERSITE=1

# Aktivace conda prostředí, pokud je k dispozici
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate TV_easydock 2>/dev/null || conda activate easydock 2>/dev/null || conda activate TV 2>/dev/null || true
fi

python3 "${SCRIPT_DIR}/evaluate_docking_results.py" "$@"
