#!/usr/bin/env bash
# ==============================================================================
# Spouštěcí skript pro generování Target Selectivity Volcano Plotů (Varianta A)
# ==============================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONNOUSERSITE=1

# Detekce a aktivace conda prostředí
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook 2>/dev/null || true)"
    conda activate TV_easydock 2>/dev/null || \
    conda activate /storage/brno2/home/urbany/.conda/envs/TV_easydock 2>/dev/null || \
    source activate /storage/brno2/home/urbany/.conda/envs/TV_easydock 2>/dev/null || \
    conda activate cadd2026 2>/dev/null || true
fi

INPUT_DIR="${1:-reports_docking}"
TARGET_OPT="${2:-}"

echo "======================================================================"
echo "Spouštím generování Target Selectivity Volcano Plotů (Varianta A)"
echo "Vstupní složka: ${INPUT_DIR}"
echo "======================================================================"

if [ -n "${TARGET_OPT}" ]; then
    echo "Vybraný cíl: ${TARGET_OPT}"
    python "${SCRIPT_DIR}/plot_docking_volcano.py" \
        -i "${INPUT_DIR}" \
        -t "${TARGET_OPT}"
else
    echo "Nebyl zadán konkrétní cíl -> generuji Volcano ploty pro TOP 6 proteinů + sdružený grid..."
    python "${SCRIPT_DIR}/plot_docking_volcano.py" \
        -i "${INPUT_DIR}" \
        --top 6 \
        --grid
fi

echo "======================================================================"
echo "Hotovo! Grafy a tabulky naleznete v: ${INPUT_DIR}/volcano_plots/"
echo "======================================================================"
