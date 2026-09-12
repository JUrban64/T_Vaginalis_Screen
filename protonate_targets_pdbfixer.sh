#!/usr/bin/env bash
# ==============================================================================
# Spouštěcí skript pro dávkovou protonaci proteinů pomocí PDBFixer (OpenMM)
# ==============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONNOUSERSITE=1

# Aktivace conda prostředí s pdbfixer / openmm, pokud existuje
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook 2>/dev/null || true)"
    if conda env list | grep -q "pdbfixer"; then
        conda activate pdbfixer
    elif conda env list | grep -q "TV_easydock"; then
        conda activate TV_easydock
    elif conda env list | grep -q "cadd2026"; then
        conda activate cadd2026
    fi
fi

echo "======================================================================"
echo "Dávková protonace AlphaFold struktur (PDBFixer)"
echo "Aktivní Python: $(which python 2>/dev/null || which python3)"
echo "======================================================================"

python "${SCRIPT_DIR}/protonate_targets_pdbfixer.py" "$@"
