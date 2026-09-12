#!/usr/bin/env bash
# ==============================================================================
# Spouštěcí skript pro generování globálního Protein-Level Volcano Plotu
# ==============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONNOUSERSITE=1

# Detekce a aktivace conda prostředí
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook 2>/dev/null || true)"
    if conda env list | grep -q "TV_easydock"; then
        conda activate TV_easydock
    elif conda env list | grep -q "cadd2026"; then
        conda activate cadd2026
    fi
fi

INPUT_DIR="${1:-reports_docking}"

echo "======================================================================"
echo "Spouštím generování Protein-Level Volcano Plotu (Všech 49 proteinů)"
echo "Vstupní složka: ${INPUT_DIR}"
echo "======================================================================"

python "${SCRIPT_DIR}/plot_protein_volcano.py" \
    -i "${INPUT_DIR}"

echo "======================================================================"
echo "Hotovo! Graf a tabulku naleznete v: ${INPUT_DIR}/protein_volcano/"
echo "======================================================================"
