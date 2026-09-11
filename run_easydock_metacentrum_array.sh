#!/bin/bash
# ==============================================================================
# OpenPBS / PBS Pro Array Job skript pro EasyDock na MetaCentru (CESNET)
#
# Účel:
#   Paralelní virtuální screening (docking) knihovny ligandů do všech proteinových
#   cílů pomocí PBS job array (1 podúloha = 1 proteinový cíl).
#
# Spuštění úlohy:
#   1. Připravte soubor se seznamem cílů (např. targets.txt, kde každý řádek je název cíle)
#   2. Zjistěte počet řádků/cílů: N=$(wc -l < targets.txt)
#   3. Odešlete úlohu: qsub -J 1-$N run_easydock_metacentrum_array.sh
# ==============================================================================

# ---------------------------- PBS Direktivity ---------------------------------
#SBATCH/PBS název úlohy
#PBS -N easydock_screen

# Zdroje na jeden uzel/cíl:
# - 16 CPU jader
# - 32 GB RAM
# - 30-50 GB rychlý lokální NVMe/SSD scratch (nutnost pro SQLite DB!)
#PBS -l select=1:ncpus=16:mem=32gb:scratch_local=40gb

# Maximální doba běhu (např. 24 hodin, upravte dle velikosti knihovny)
#PBS -l walltime=24:00:00

# Spojení STDOUT a STDERR do jednoho log souboru
#PBS -j oe

# Výchozí rozsah job array (lze přepsat při volání qsub -J 1-N ...)
#PBS -J 1-49

# ------------------------------------------------------------------------------
set -e
set -o pipefail

echo "========================================================================"
echo "EasyDock Screening na MetaCentru (OpenPBS Job Array)"
echo "Datum spuštění:      $(date)"
echo "Běží na uzlu:        $(hostname -f)"
echo "PBS Job ID:          ${PBS_JOBID}"
echo "PBS Array Index:     ${PBS_ARRAY_INDEX}"
echo "Pracovní adresář:    ${PBS_O_WORKDIR}"
echo "Lokální Scratch:     ${SCRATCHDIR}"
echo "========================================================================"

# Kontrola dostupnosti SCRATCHDIR (ochrana před zaplněním síťového disku)
if [ -z "$SCRATCHDIR" ]; then
    echo "[!] CHYBA: Proměnná \$SCRATCHDIR není definována!" >&2
    exit 1
fi

# -------------------- 1. Uživatelská konfigurace cest -------------------------
WORKDIR="${PBS_O_WORKDIR}"

# Seznam cílů (každý řádek odpovídá jednomu indexu v job array)
TARGETS_FILE="${WORKDIR}/targets.txt"

# Složky s daty
TARGETS_DIR="${WORKDIR}/targets"           # Zde jsou podadresáře pro jednotlivé proteiny
RESULTS_DIR="${WORKDIR}/results_docking"    # Kam se uloží výsledné .db a .sdf
LOGS_DIR="${WORKDIR}/logs_docking"          # Kam se zkopírují logy

# Knihovna ligandů (dvě možnosti - viz níže):
# Možnost A (DOPORUČENÁ): Šablona již předpřipravené databáze (obsahuje 3D struktury, protonaci, isomery)
LIGANDS_TEMPLATE_DB="${WORKDIR}/ligands_template.db"

# Možnost B: Surový SMILES soubor (pokud nemáte předpřipravenou DB)
LIGANDS_SMI="${WORKDIR}/ligands.smi"

# Docking parametry
DOCKING_PROGRAM="vina"      # vina, gnina, qvina, atd.
NCPUS=16                    # Celkový počet jader alokovaných PBS
PARALLEL_MOLS=8             # Kolik molekul dockovat paralelně (-c)
NCPU_PER_MOL=2              # Vlákna na jednu molekulu v konfigu (8 * 2 = 16 jader)

# ------------------- 2. Výběr konkrétního proteinového cíle -------------------
if [ ! -f "$TARGETS_FILE" ]; then
    echo "[!] CHYBA: Soubor se seznamem cílů '${TARGETS_FILE}' neexistuje!" >&2
    exit 1
fi

# Přečteme řádek odpovídající aktuálnímu PBS_ARRAY_INDEX
TARGET_ID=$(sed -n "${PBS_ARRAY_INDEX}p" "$TARGETS_FILE" | tr -d '\r')

if [ -z "$TARGET_ID" ]; then
    echo "[!] CHYBA: Řádek ${PBS_ARRAY_INDEX} v '${TARGETS_FILE}' je prázdný!" >&2
    exit 1
fi

echo "[*] Zpracovávám proteinový cíl č. ${PBS_ARRAY_INDEX}: '${TARGET_ID}'"

TARGET_SPECIFIC_DIR="${TARGETS_DIR}/${TARGET_ID}"
PROTEIN_PDBQT="${TARGET_SPECIFIC_DIR}/protein.pdbqt"
GRID_TXT="${TARGET_SPECIFIC_DIR}/grid.txt"
PROTEIN_PDB="${TARGET_SPECIFIC_DIR}/protein.pdb"   # Volitelné, pro PoseBusters / PLIF

if [ ! -f "$PROTEIN_PDBQT" ] || [ ! -f "$GRID_TXT" ]; then
    echo "[!] CHYBA: Chybí '${PROTEIN_PDBQT}' nebo '${GRID_TXT}' pro cíl '${TARGET_ID}'!" >&2
    exit 1
fi

# ----------------- 3. Aktivace prostředí / Moduly MetaCentra ------------------
echo "[*] Aktivuji prostředí pro EasyDock..."
# Zvolte podle vašeho nastavení na MetaCentru (Conda modul nebo vlastní miniconda):
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate easydock || {
        echo "[!] Nepodařilo se aktivovat conda env 'easydock', zkouším výchozí python"
    }
elif [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
    # module add conda-modules
    # conda activate easydock
fi

# Ověření nástroje
command -v easydock >/dev/null 2>&1 || {
    echo "[!] CHYBA: Příkaz 'easydock' nebyl nalezen v PATH!" >&2
    exit 1
}

# ----------------- 4. Příprava lokálního scratch disku ------------------------
# DŮLEŽITÉ: SQLite v WAL módu NIKDY nespouštějte přes síťové NFS/Ceph úložiště!
# Všechny výpočty a I/O operace musí běžet v $SCRATCHDIR a až na konci zkopírovat zpět.
mkdir -p "${RESULTS_DIR}" "${LOGS_DIR}"

# Nastavení úklidu scratch disku při jakémkoliv ukončení (úspěch i selhání)
cleanup() {
    EXIT_CODE=$?
    echo "[*] Probíhá úklid scratch disku: ${SCRATCHDIR} (Exit code: ${EXIT_CODE})..."
    clean_scratch
    exit ${EXIT_CODE}
}
trap cleanup TERM EXIT

cd "${SCRATCHDIR}"
echo "[*] Přesunuto do scratch: $(pwd)"

# Kopírování potřebných souborů receptoru do scratch
cp "$PROTEIN_PDBQT" "${SCRATCHDIR}/protein.pdbqt"
cp "$GRID_TXT" "${SCRATCHDIR}/grid.txt"
[ -f "$PROTEIN_PDB" ] && cp "$PROTEIN_PDB" "${SCRATCHDIR}/protein.pdb"

# Vytvoření cílového config.yml pro EasyDock
cat << EOF > "${SCRATCHDIR}/config.yml"
protein: ${SCRATCHDIR}/protein.pdbqt
protein_setup: ${SCRATCHDIR}/grid.txt
exhaustiveness: 8
seed: 42
n_poses: 5
ncpu: ${NCPU_PER_MOL}
EOF

# ----------------- 5. Příprava ligandové databáze -----------------------------
DB_NAME="${TARGET_ID}.db"
SDF_NAME="${TARGET_ID}.sdf"

if [ -f "$LIGANDS_TEMPLATE_DB" ]; then
    echo "[*] Nalezena šablona ligandové DB: ${LIGANDS_TEMPLATE_DB}"
    echo "[*] Vytvářím čistou kopii databáze pomocí 'make_clean_copy' pro ${TARGET_ID}..."
    make_clean_copy -i "$LIGANDS_TEMPLATE_DB" -o "${SCRATCHDIR}/${DB_NAME}"
    RUN_ARGS="-o ${SCRATCHDIR}/${DB_NAME} --program ${DOCKING_PROGRAM} --config ${SCRATCHDIR}/config.yml -c ${PARALLEL_MOLS} --sdf"
elif [ -f "$LIGANDS_SMI" ]; then
    echo "[*] Šablona DB nenalezena, použiji surové ligandy ze SMILES: ${LIGANDS_SMI}"
    cp "$LIGANDS_SMI" "${SCRATCHDIR}/ligands.smi"
    RUN_ARGS="-i ${SCRATCHDIR}/ligands.smi -o ${SCRATCHDIR}/${DB_NAME} --program ${DOCKING_PROGRAM} --config ${SCRATCHDIR}/config.yml --protonation molgpka -s 4 -c ${PARALLEL_MOLS} --sdf"
else
    echo "[!] CHYBA: Nebyla nalezena šablona ligandů '${LIGANDS_TEMPLATE_DB}' ani SMILES '${LIGANDS_SMI}'!" >&2
    exit 1
fi

# ----------------- 6. Spuštění EasyDockingu ----------------------------------
echo "------------------------------------------------------------------------"
echo "[*] Spouštím docking pro cíl '${TARGET_ID}'..."
echo "Příkaz: easydock ${RUN_ARGS}"
echo "------------------------------------------------------------------------"

easydock ${RUN_ARGS}

echo "[+] Docking pro '${TARGET_ID}' úspěšně dokončen!"

# ----------------- 7. Volitelná kontrola kvality a interakcí ------------------
# A) PoseBusters validace (pokud máme protein.pdb s vodíky)
if [ -f "${SCRATCHDIR}/protein.pdb" ] && command -v easydock_bust >/dev/null 2>&1; then
    echo "[*] Spouštím PoseBusters kontrolu kvality (easydock_bust)..."
    easydock_bust -i "${SCRATCHDIR}/${DB_NAME}" -p "${SCRATCHDIR}/protein.pdb" -c ${PARALLEL_MOLS} || {
        echo "[!] Varování: PoseBusters kontrola selhala nebo skončila s chybou."
    }
fi

# B) PLIF interakční fingerprinty (pokud máme protein.pdb s vodíky)
if [ -f "${SCRATCHDIR}/protein.pdb" ] && command -v easydock_plif >/dev/null 2>&1; then
    echo "[*] Počítám protein-ligand interakční fingerprinty (easydock_plif)..."
    easydock_plif -i "${SCRATCHDIR}/${DB_NAME}" -p "${SCRATCHDIR}/protein.pdb" -c ${PARALLEL_MOLS} || {
        echo "[!] Varování: PLIF výpočet selhal nebo skončil s chybou."
    }
fi

# ----------------- 8. Zkopírování výsledků zpět do úložiště --------------------
echo "[*] Kopíruji výsledky zpět do ${RESULTS_DIR}..."

# Zkopírujeme výslednou DB a vygenerovaný SDF
[ -f "${SCRATCHDIR}/${DB_NAME}" ] && cp -v "${SCRATCHDIR}/${DB_NAME}" "${RESULTS_DIR}/"
[ -f "${SCRATCHDIR}/${SDF_NAME}" ] && cp -v "${SCRATCHDIR}/${SDF_NAME}" "${RESULTS_DIR}/"

# Pokud vznikly dodatečné logy/reporty, uložíme je také
if compgen -G "${SCRATCHDIR}/*.log" > /dev/null; then
    cp -v "${SCRATCHDIR}"/*.log "${LOGS_DIR}/"
fi

echo "========================================================================"
echo "[+] Hotovo pro cíl '${TARGET_ID}' v: $(date)"
echo "Výsledky uloženy v:"
echo "  - DB:  ${RESULTS_DIR}/${DB_NAME}"
echo "  - SDF: ${RESULTS_DIR}/${SDF_NAME}"
echo "========================================================================"
