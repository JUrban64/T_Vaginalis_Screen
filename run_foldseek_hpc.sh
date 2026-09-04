#!/bin/bash
# ==============================================================================
# SLURM / Bash skript pro spuštění Foldseek TM-score klastrování na HPC
# ==============================================================================
#SBATCH --job-name=foldseek_cluster
#SBATCH --output=foldseek_cluster_%j.log
#SBATCH --error=foldseek_cluster_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -e

echo "=== Začátek úlohy: $(date) ==="
echo "Běží na uzlu: $(hostname)"
echo "Pracovní adresář: $(pwd)"

# 1. Aktivace prostředí s Foldseekem (upravte podle konfigurace vašeho clusteru)
if command -v conda &> /dev/null; then
    # conda activate foldseek_env
    echo "[*] Conda detekována"
elif [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
    # module load foldseek
    # module load python/3.10
    echo "[*] Environment modules detekovány"
fi

# Ověření dostupnosti nástrojů
command -v python3 >/dev/null 2>&1 || { echo "[!] Python3 nebyl nalezen!"; exit 1; }
command -v foldseek >/dev/null 2>&1 || { echo "[!] Varování: Foldseek nebyl nalezen v PATH. Ujistěte se, že je načten modul nebo conda prostředí."; }

# 2. Nastavení parametrů
INPUT_DIR="structures"
OUTPUT_DIR="structures_representative"
REPORTS_DIR="reports_clustering"
TM_THRESHOLD="0.75"      # 'Tak akorát přísný' práh pro redukci duplikátů v T. vaginalis
COV_THRESHOLD="0.75"     # Minimální pokrytí délky
THREADS=${SLURM_CPUS_PER_TASK:-8}

# 3. Spuštění klastrování s průzkumem citlivosti prahů
echo "[*] Spouštím cluster_methyltransferases.py s TM-score = ${TM_THRESHOLD} a ${THREADS} vlákny..."
python3 cluster_methyltransferases.py \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --reports-dir "${REPORTS_DIR}" \
    --tmscore-threshold "${TM_THRESHOLD}" \
    --coverage-threshold "${COV_THRESHOLD}" \
    --threads "${THREADS}" \
    --representative-metric "plddt" \
    --explore-thresholds

echo "=== Dokončeno: $(date) ==="
echo "Reprezentativní struktury jsou uloženy v: ${OUTPUT_DIR}/"
echo "Reporty naleznete v: ${REPORTS_DIR}/"
