#!/usr/bin/env python3
# ==============================================================================
# protonate_targets_pdbfixer.py
#
# Dávková protonace a doplnění chybějících vodíků do AlphaFold PDB struktur
# pomocí PDBFixer (OpenMM) při fyziologickém pH (výchozí 7.4).
#
# Proč je tento krok nutný:
#   AlphaFold modely obsahují POUZE těžké atomy (C, N, O, S) a nemají ŽÁDNÉ vodíky.
#   Meeko (mk_prepare_receptor.py) vodíky samo negeneruje, ale vyžaduje je.
#   Bez polárních vodíků (HD) v PDBQT AutoDock Vina selže s chybovým kódem 1.
#
# Funkcionalita:
#   1. Projde všechny cíle v 'targets/' (nebo ze seznamu 'targets.txt').
#   2. Zálohuje původní strukturu jako 'protein_orig.pdb' (pokud ještě neexistuje).
#   3. Pomocí PDBFixer doplní chybějící těžké atomy a všechny vodíky při pH 7.4.
#   4. Uloží kompletní protonovanou strukturu jako 'protein.pdb' (a 'protein_h.pdb').
#   5. Běží paralelně přes zadaný počet CPU jader (ProcessPoolExecutor).
#
# Autor: Antigravity AI / Jáchym Urban
# ==============================================================================

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


def count_atoms_and_hydrogens(pdb_path: Path) -> tuple[int, int]:
    """Spočítá celkový počet atomů a počet vodíků v PDB souboru."""
    total_atoms = 0
    hydrogens = 0
    if not pdb_path.exists():
        return 0, 0

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                total_atoms += 1
                element = line[76:78].strip().upper() if len(line) >= 78 else ""
                atom_name = line[12:16].strip()
                if element == "H" or (not element and atom_name.startswith("H")):
                    hydrogens += 1

    return total_atoms, hydrogens


def fix_and_protonate_single_target(
    target_id: str,
    target_dir: Path,
    ph: float = 7.4,
    force: bool = False,
) -> tuple[str, bool, str]:
    """
    Protonuje jeden protein.pdb v daném adresáři pomocí PDBFixer.
    Vrací: (target_id, success, zpráva)
    """
    protein_pdb = target_dir / "protein.pdb"
    protein_orig = target_dir / "protein_orig.pdb"
    protein_h = target_dir / "protein_h.pdb"

    if not protein_pdb.exists() and not protein_orig.exists():
        return target_id, False, f"Chybí soubor protein.pdb v {target_dir}"

    # Pokud již existuje záloha orig.pdb, vstupem je orig.pdb (aby se nereprotonovalo dokola)
    source_pdb = protein_orig if protein_orig.exists() else protein_pdb

    # Kontrola, zda již máme protonováno
    if protein_orig.exists() and not force:
        _, h_count = count_atoms_and_hydrogens(protein_pdb)
        if h_count > 50:
            return target_id, True, f"Již protonováno ({h_count} vodíků, přeskočeno)"

    # Záloha původního PDB (pouze pokud ještě neexistuje)
    if not protein_orig.exists():
        shutil.copy2(protein_pdb, protein_orig)

    # 1. Způsob: Použití Python API PDBFixer / OpenMM
    api_success = False
    err_api = ""
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile

        fixer = PDBFixer(filename=str(source_pdb))

        # Kontrola a doplnění chybějících atomů
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        # Doplnění vodíků při zadaném pH
        fixer.addMissingHydrogens(ph)

        # Dočasný výstup pro ověření
        tmp_out = target_dir / "protein_fixed_tmp.pdb"
        with open(tmp_out, "w", encoding="utf-8") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

        if tmp_out.exists() and tmp_out.stat().st_size > 0:
            shutil.move(tmp_out, protein_pdb)
            shutil.copy2(protein_pdb, protein_h)
            api_success = True
    except ImportError:
        err_api = "PDBFixer Python API není dostupné v aktuálním prostředí"
    except Exception as e:
        err_api = f"PDBFixer API výjimka: {e}"

    if api_success:
        total, h_count = count_atoms_and_hydrogens(protein_pdb)
        return target_id, True, f"PDBFixer API -> {total} atomů ({h_count} vodíků, pH {ph})"

    # 2. Způsob: Použití CLI nástroje 'pdbfixer'
    pdbfixer_cli = shutil.which("pdbfixer")
    if pdbfixer_cli:
        cmd = [
            pdbfixer_cli,
            str(source_pdb),
            f"--output={protein_pdb}",
            f"--ph={ph}",
            "--add-atoms=heavy",
            "--add-hydrogens",
            "--keep-ids"
        ]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if protein_pdb.exists() and protein_pdb.stat().st_size > 0:
                shutil.copy2(protein_pdb, protein_h)
                total, h_count = count_atoms_and_hydrogens(protein_pdb)
                return target_id, True, f"PDBFixer CLI -> {total} atomů ({h_count} vodíků, pH {ph})"
            err_cli = res.stderr.strip() or res.stdout.strip()
            return target_id, False, f"Chyba CLI PDBFixer: {err_cli[:200]}"
        except Exception as e:
            return target_id, False, f"PDBFixer CLI výjimka: {e}"

    # 3. Způsob: Záložní protonace přes OpenBabel (výhradně PDB -> PDB, zachová formát PDB!)
    obabel_bin = shutil.which("obabel")
    if obabel_bin:
        cmd = [obabel_bin, "-ipdb", str(source_pdb), "-opdb", "-O", str(protein_pdb), "-p", str(ph), "-h"]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if protein_pdb.exists() and protein_pdb.stat().st_size > 0:
                shutil.copy2(protein_pdb, protein_h)
                total, h_count = count_atoms_and_hydrogens(protein_pdb)
                return target_id, True, f"OpenBabel fallback (PDB) -> {total} atomů ({h_count} vodíků, pH {ph})"
        except Exception as e:
            pass

    return target_id, False, f"Selhaly všechny metody: {err_api}"


def main():
    parser = argparse.ArgumentParser(
        description="Dávkové doplnění chybějících vodíků a atomů do receptorů pomocí PDBFixer.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Příklady použití:
  # Standardní protonace všech 49 cílů při pH 7.4 (8 jader):
  python protonate_targets_pdbfixer.py

  # Přepsání i dříve protonovaných cílů:
  python protonate_targets_pdbfixer.py --force -c 16

  # Spuštění přes bash wrapper:
  ./protonate_targets_pdbfixer.sh
        """
    )
    parser.add_argument(
        "-t", "--targets-dir",
        type=Path,
        default=Path("targets"),
        help="Složka s cíli (výchozí: 'targets')."
    )
    parser.add_argument(
        "-l", "--targets-list",
        type=Path,
        default=Path("targets.txt"),
        help="Seznam cílů (výchozí: 'targets.txt')."
    )
    parser.add_argument(
        "--ph",
        type=float,
        default=7.4,
        help="Hodnota pH pro výpočet protonačních stavů aminokyselin (výchozí: 7.4)."
    )
    parser.add_argument(
        "-c", "--ncpu",
        type=int,
        default=min(os.cpu_count() or 4, 8),
        help=f"Počet paralelních vláken (výchozí: {min(os.cpu_count() or 4, 8)})."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Vynutit reprotonaci i pro cíle, které již byly zpracovány."
    )

    args = parser.parse_args()

    if not args.targets_dir.exists():
        print(f"[!] CHYBA: Adresář '{args.targets_dir}' neexistuje!", file=sys.stderr)
        sys.exit(1)

    # Načtení seznamu cílů
    target_ids = []
    if args.targets_list.exists():
        with open(args.targets_list, "r", encoding="utf-8") as f:
            target_ids = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    if not target_ids:
        # Prohledáme podsložky
        target_ids = sorted([p.name for p in args.targets_dir.iterdir() if p.is_dir() and (p / "protein.pdb").exists()])

    if not target_ids:
        print(f"[!] CHYBA: V '{args.targets_dir}' nebyly nalezeny žádné cíle!", file=sys.stderr)
        sys.exit(1)

    print("=" * 70)
    print("Dávková protonace proteinových struktur (PDBFixer / OpenMM):")
    print(f"  Cílová složka:     {args.targets_dir}")
    print(f"  Celkem cílů:       {len(target_ids)}")
    print(f"  Hodnota pH:        {args.ph}")
    print(f"  Paralelních jader: {args.ncpu}")
    print(f"  Přepsat stávající: {args.force}")
    print("=" * 70)

    start_time = time.time()
    success_count = 0
    failed_targets = []

    tasks = []
    for tid in target_ids:
        tdir = args.targets_dir / tid
        tasks.append((tid, tdir, args.ph, args.force))

    if args.ncpu > 1 and len(tasks) > 1:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.ncpu) as executor:
            futures = [executor.submit(fix_and_protonate_single_target, *t) for t in tasks]
            for future in concurrent.futures.as_completed(futures):
                try:
                    tid, ok, msg = future.result()
                    if ok:
                        success_count += 1
                        print(f"  [OK]  {tid:<15} -> {msg}")
                    else:
                        failed_targets.append((tid, msg))
                        print(f"  [ERR] {tid:<15} -> {msg}", file=sys.stderr)
                except Exception as exc:
                    print(f"  [EXC] Úloha skončila výjimkou: {exc}", file=sys.stderr)
    else:
        for t in tasks:
            tid, ok, msg = fix_and_protonate_single_target(*t)
            if ok:
                success_count += 1
                print(f"  [OK]  {tid:<15} -> {msg}")
            else:
                failed_targets.append((tid, msg))
                print(f"  [ERR] {tid:<15} -> {msg}", file=sys.stderr)

    elapsed = time.time() - start_time
    print("=" * 70)
    print(f"Protonace dokončena za {elapsed:.2f} s:")
    print(f"  Úspěšně: {success_count}/{len(target_ids)}")

    if failed_targets:
        print(f"  Neúspěšně: {len(failed_targets)}")
        for tid, err in failed_targets:
            print(f"    - {tid}: {err}")
        sys.exit(1)
    else:
        print("\n[✓] Všechny receptory mají nyní kompletní sadu vodíků při pH 7.4!")
        print("    Nyní můžete spustit: ./prepare_receptors_meeko.sh --force")


if __name__ == "__main__":
    main()
