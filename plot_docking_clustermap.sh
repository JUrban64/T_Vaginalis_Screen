#!/usr/bin/env bash
# ==============================================================================
# Spouštěcí skript pro generování Clustermapy dokovacích skóre
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

echo "======================================================================"
echo "Spouštím generování Clustermap (Proteiny x Ligandy)"
echo "Vstupní složka: ${INPUT_DIR}"
echo "======================================================================"

python "${SCRIPT_DIR}/plot_docking_clustermap.py" \
    -i "${INPUT_DIR}" \
    --both \
    --cmap viridis_r

echo ""
echo "Spouštím generování Z-score normalizované Clustermapy (pro selektivitu)..."
python "${SCRIPT_DIR}/plot_docking_clustermap.py" \
    -i "${INPUT_DIR}" \
    --standardize compound \
    --cmap coolwarm_r

echo "======================================================================"
echo "Hotovo! Grafy naleznete v: ${INPUT_DIR}/plots/"
echo "======================================================================"
