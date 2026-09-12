#!/bin/bash
# ==============================================================================
# Bash wrapper pro dávkový převod protein.pdb na protein.pdbqt pomocí Meeko
#
# Použití:
#   ./prepare_receptors_meeko.sh
#   ./prepare_receptors_meeko.sh -c 8 --force
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

python "${SCRIPT_DIR}/prepare_receptors_meeko.py" "$@"
