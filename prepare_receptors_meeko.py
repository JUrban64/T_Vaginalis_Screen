#!/usr/bin/env python3
"""
prepare_receptors_meeko.py

Dávkový převod všech proteinových struktur PDB na PDBQT pomocí balíčku Meeko.
Tento skript projde všechny cíle v 'targets/' (definované v 'targets.txt')
a pro každý cíl vytvoří soubor 'protein.pdbqt' potřebný pro docking v EasyDocku / AutoDock Vina.

Použití:
    python3 prepare_receptors_meeko.py
    python3 prepare_receptors_meeko.py -c 8 --force
    ./prepare_receptors_meeko.sh
"""

import argparse
import concurrent.futures
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# Zamezení načítání starých nekompatibilních balíčků z ~/.local
os.environ["PYTHONNOUSERSITE"] = "1"


def convert_pdb_to_pdbqt(pdb_path: Path, pdbqt_path: Path, force: bool = False) -> tuple[str, bool, str]:
    """
    Převede jeden protein.pdb na protein.pdbqt pomocí Meeko (Python API nebo CLI).
    Vrací: (target_id, success, message)
    """
    target_id = pdb_path.parent.name

    if pdbqt_path.exists() and not force:
        return target_id, True, "Již existuje (přeskočeno)"

    if not pdb_path.exists():
        return target_id, False, f"Vstupní PDB neexistuje: {pdb_path}"

    errors = []

    # 1. Způsob: Přímé volání Meeko CLI modulu přes aktuální Python interpreter
    # (nejspolehlivější, nevyžaduje mít mk_prepare_receptor.py v PATH)
    python_bin = sys.executable
    out_base = pdbqt_path.with_suffix("")
    rigid_pdbqt = pdbqt_path.parent / f"{pdbqt_path.stem}_rigid.pdbqt"

    pdb_name = pdb_path.name

    def find_and_standardize_output() -> bool:
        """Zkontroluje a případně přejmenuje _rigid.pdbqt na požadovaný protein.pdbqt."""
        if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
            return True
        candidates = [
            rigid_pdbqt,
            pdbqt_path.parent / f"{pdbqt_path.name}_rigid.pdbqt",
            pdbqt_path.parent / f"{pdbqt_path.stem}_rigid.pdbqt",
            Path(str(out_base) + "_rigid.pdbqt"),
            Path(str(out_base) + ".pdbqt"),
        ]
        for cand in candidates:
            if cand.exists() and cand.stat().st_size > 0:
                shutil.move(cand, pdbqt_path)
                return True
        return False

    # 1. Způsob: Standardní ProDy čtení (-i) s explicitním --write_pdbqt <cesta>
    cmd_meeko_prody = [
        python_bin, "-m", "meeko.cli.mk_prepare_receptor",
        "-i", str(pdb_path),
        "-o", str(out_base),
        "--write_pdbqt", str(pdbqt_path)
    ]
    try:
        res = subprocess.run(cmd_meeko_prody, capture_output=True, text=True)
        if find_and_standardize_output():
            return target_id, True, f"Převedeno z {pdb_name} přes Meeko (ProDy -> pdbqt)"
        err = res.stderr.strip() or res.stdout.strip()
        errors.append(f"Meeko (ProDy) selhalo: {err[:250] if err else 'výstupní soubor nevznikl'}")
    except Exception as e:
        errors.append(f"Meeko (ProDy) výjimka: {e}")

    # 2. Způsob: Čtení bez ProDy přes --read_pdb s explicitním --write_pdbqt <cesta>
    cmd_meeko_direct = [
        python_bin, "-m", "meeko.cli.mk_prepare_receptor",
        "--read_pdb", str(pdb_path),
        "-o", str(out_base),
        "--write_pdbqt", str(pdbqt_path)
    ]
    try:
        res = subprocess.run(cmd_meeko_direct, capture_output=True, text=True)
        if find_and_standardize_output():
            return target_id, True, f"Převedeno z {pdb_name} přes Meeko (--read_pdb -> pdbqt)"
        err = res.stderr.strip() or res.stdout.strip()
        errors.append(f"Meeko (--read_pdb) selhalo: {err[:250] if err else 'výstupní soubor nevznikl'}")
    except Exception as e:
        errors.append(f"Meeko (--read_pdb) výjimka: {e}")

    # 3. Způsob: Samostatná binárka mk_prepare_receptor.py z PATH nebo sys.prefix/bin
    mk_bin = shutil.which("mk_prepare_receptor.py") or shutil.which("mk_prepare_receptor")
    if not mk_bin:
        conda_mk = Path(sys.prefix) / "bin" / "mk_prepare_receptor.py"
        if conda_mk.exists():
            mk_bin = str(conda_mk)

    if mk_bin:
        cmd_bin = [
            mk_bin,
            "-i", str(pdb_path),
            "-o", str(out_base),
            "--write_pdbqt", str(pdbqt_path)
        ]
        try:
            res = subprocess.run(cmd_bin, capture_output=True, text=True)
            if find_and_standardize_output():
                return target_id, True, f"Převedeno z {pdb_name} přes mk_prepare_receptor.py binárku"
            err = res.stderr.strip() or res.stdout.strip()
            errors.append(f"mk_prepare_receptor.py binárka selhala: {err[:250] if err else 'výstup nevznikl'}")
        except Exception as e:
            errors.append(f"mk_prepare_receptor.py binárka výjimka: {e}")
    else:
        errors.append("mk_prepare_receptor.py binárka nenalezena")

    summary_err = " | ".join(errors)
    return target_id, False, f"Chyba: {summary_err}"


def main():
    parser = argparse.ArgumentParser(
        description="Dávkový převod PDB struktur cílů na PDBQT pomocí Meeko."
    )
    parser.add_argument(
        "-t", "--targets-dir",
        default="targets",
        help="Složka s cíli (default: targets)"
    )
    parser.add_argument(
        "-l", "--targets-list",
        default="targets.txt",
        help="Seznam cílů (default: targets.txt)"
    )
    parser.add_argument(
        "-c", "--ncpu",
        type=int,
        default=min(os.cpu_count() or 4, 8),
        help=f"Počet paralelních vláken (default: {min(os.cpu_count() or 4, 8)})"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Přepsat již existující protein.pdbqt soubory"
    )
    args = parser.parse_args()

    targets_dir = Path(args.targets_dir)
    targets_list = Path(args.targets_list)

    if not targets_dir.exists():
        print(f"[!] CHYBA: Složka '{targets_dir}' neexistuje!", file=sys.stderr)
        sys.exit(1)

    # Načtení seznamu cílů
    if targets_list.exists():
        with open(targets_list, "r", encoding="utf-8") as f:
            target_ids = [line.strip() for line in f if line.strip() and not line.startswith("#")]
        print(f"[*] Načteno {len(target_ids)} cílů ze souboru '{targets_list}'.")
    else:
        # Jinak prohledáme přímo podsložky v targets/
        target_ids = [d.name for d in targets_dir.iterdir() if d.is_dir()]
        print(f"[*] Načteno {len(target_ids)} cílů z adresáře '{targets_dir}'.")

    if not target_ids:
        print("[!] Žádné cíle k převodu nebyly nalezeny!", file=sys.stderr)
        sys.exit(1)

    print("=" * 70)
    print("Dávková příprava receptorů PDBQT pomocí Meeko:")
    print(f"  Cílová složka:     {targets_dir}")
    print(f"  Celkem cílů:       {len(target_ids)}")
    print(f"  Paralelních jader: {args.ncpu}")
    print(f"  Přepsat existující:{args.force}")
    print("=" * 70)

    start_time = time.time()
    tasks = []
    for tid in target_ids:
        # Pokud existuje protonovaný protein_h.pdb (např. z PDBFixer), použijeme přednostně ten
        h_pdb = targets_dir / tid / "protein_h.pdb"
        std_pdb = targets_dir / tid / "protein.pdb"

        if h_pdb.exists():
            pdb_file = h_pdb
        else:
            pdb_file = std_pdb

        pdbqt_file = targets_dir / tid / "protein.pdbqt"
        tasks.append((pdb_file, pdbqt_file))

    success_count = 0
    fail_count = 0

    with concurrent.futures.ProcessPoolExecutor(max_workers=args.ncpu) as executor:
        futures = {
            executor.submit(convert_pdb_to_pdbqt, pdb, pdbqt, args.force): pdb.parent.name
            for pdb, pdbqt in tasks
        }
        for future in concurrent.futures.as_completed(futures):
            tid = futures[future]
            try:
                target_id, ok, msg = future.result()
                if ok:
                    success_count += 1
                    print(f"  [OK]  {target_id:15s} -> {msg}")
                else:
                    fail_count += 1
                    print(f"  [ERR] {target_id:15s} -> {msg}", file=sys.stderr)
            except Exception as e:
                fail_count += 1
                print(f"  [EXC] {tid:15s} -> Neočekávaná výjimka: {e}", file=sys.stderr)

    elapsed = time.time() - start_time
    print("=" * 70)
    print(f"Převod dokončen za {elapsed:.2f} s ({elapsed / 60:.2f} min).")
    print(f"Úspěšně připraveno: {success_count} / {len(target_ids)}")
    if fail_count > 0:
        print(f"Selhalo:            {fail_count} cílů!", file=sys.stderr)
    print("=" * 70)


if __name__ == "__main__":
    main()
