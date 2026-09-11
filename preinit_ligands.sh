#!/bin/bash
# ==============================================================================
# Skript pro lokální předpřípravu ligandů v EasyDocku (bez nutnosti PBS / qsub)
#
# Pro menší knihovny (~160 molekul) není potřeba zadávat výpočetní úlohu na PBS;
# výpočet trvá jen pár minut a lze jej spustit přímo lokálně nebo na frontendovém uzlu.
#
# Použití:
#   ./preinit_ligands.sh
#   ./preinit_ligands.sh -i moje_ligandy.smi -c 8
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Aktivace conda prostředí, pokud je k dispozici
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate easydock 2>/dev/null || true
fi

# Spuštění python skriptu se všemi předanými parametry
python3 "${SCRIPT_DIR}/preinit_ligands.py" "$@"
