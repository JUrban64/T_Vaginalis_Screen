#!/usr/bin/env python3
"""
prepare_targets_setup.py

Příprava struktury adresářů cílů (targets/) a seznamu targets.txt
pro EasyDock screening na MetaCentru na základě ručně zadaných grid boxů (girds.csv).

Funkcionalita:
1. Načte 'structures_deduplicated/girds.csv'.
2. Vyřadí proteiny označené jako '///' (non-methyltransferázy).
3. Pro validní cíle extrahuje souřadnice středu [center_x, center_y, center_z]
   a vytvoří 'targets/<target_id>/grid.txt' se standardní velikostí vazebné kapsy (výchozí 20x20x20 Å).
4. Zkopíruje odpovídající PDB strukturu jako 'targets/<target_id>/protein.pdb'.
5. Vyčistí ze složky 'targets/' případné dříve vytvořené složky nevalidních cílů.
6. Vygeneruje aktualizovaný 'targets.txt' obsahující pouze validní methyltransferázy.

Použití:
    python3 prepare_targets_setup.py
    python3 prepare_targets_setup.py --box-size 22.0
"""

import argparse
import csv
import os
import re
import shutil
import sys
from pathlib import Path


def parse_accession(filename: str) -> str:
    """Extrahuje UniProt accession nebo unikátní ID z názvu souboru."""
    base = os.path.splitext(os.path.basename(filename))[0]
    match = re.search(r"AF-([A-Za-z0-9]+)-F\d+", base)
    if match:
        return match.group(1)
    return base


def parse_grids_csv(csv_path: Path):
    """
    Načte a naparsuje CSV soubor s grid boxy (např. girds.csv).
    Vrací:
        valid_targets: dict[accession] = (center_x, center_y, center_z)
        invalid_targets: list[accession]
    """
    if not csv_path.exists():
        raise FileNotFoundError(f"Soubor '{csv_path}' neexistuje!")

    with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    reader = csv.reader(content.splitlines(), delimiter=";")
    valid_targets = {}
    invalid_targets = []

    for row in reader:
        if not row or not row[0].strip():
            continue

        target_id = row[0].strip().replace("\n", "").replace("\r", "").strip('"').strip("'")
        if target_id.lower().startswith("protein"):
            continue

        raw_center = row[1].strip() if len(row) > 1 else ""

        if "///" in raw_center or not raw_center:
            invalid_targets.append(target_id)
            continue

        # Oprava známých překlepů z ručního zápisu
        clean_center = raw_center
        # Překlep u A2F7I7: 1-11 -> -11
        if "1-11" in clean_center:
            clean_center = clean_center.replace("1-11", "-11")
        # Čárka místo tečky u A2FYV2: -0,6 -> -0.6
        if "-0,6" in clean_center:
            clean_center = clean_center.replace("-0,6", "-0.6")

        # Rozdělení podle čárky
        parts = [p.strip() for p in clean_center.split(",") if p.strip()]
        if len(parts) != 3:
            print(f"[!] Varování: Cíl '{target_id}' má neplatný počet souřadnic '{raw_center}' -> {parts}", file=sys.stderr)
            invalid_targets.append(target_id)
            continue

        try:
            cx = float(parts[0])
            cy = float(parts[1])
            cz = float(parts[2])
            valid_targets[target_id] = (cx, cy, cz)
        except ValueError as e:
            print(f"[!] Varování: Nelze převést souřadnice u '{target_id}': '{raw_center}' ({e})", file=sys.stderr)
            invalid_targets.append(target_id)

    return valid_targets, invalid_targets


def main():
    parser = argparse.ArgumentParser(description="Příprava cílů a grid boxů pro EasyDock na základě grids.csv.")
    parser.add_argument(
        "-g", "--grids-csv",
        default="structures_deduplicated/girds.csv",
        help="Cesta k CSV souboru s grid boxy (default: structures_deduplicated/girds.csv)"
    )
    parser.add_argument(
        "-s", "--structures-dir",
        default="structures_deduplicated",
        help="Cesta ke složce s PDB strukturami (default: structures_deduplicated)"
    )
    parser.add_argument(
        "-t", "--targets-dir",
        default="targets",
        help="Cílová složka pro přípravu (default: targets)"
    )
    parser.add_argument(
        "-o", "--targets-list",
        default="targets.txt",
        help="Výstupní seznam validních cílů pro PBS job array (default: targets.txt)"
    )
    parser.add_argument(
        "--box-size",
        type=float,
        default=20.0,
        help="Velikost hrany grid boxu v Angströmech pro cílený docking (default: 20.0 Å)"
    )
    args = parser.parse_args()

    grids_csv = Path(args.grids_csv)
    struct_dir = Path(args.structures_dir)
    targets_dir = Path(args.targets_dir)

    print(f"[*] Načítám grid boxy z '{grids_csv}'...")
    valid_targets, invalid_targets = parse_grids_csv(grids_csv)

    print(f"[+] Nalezeno {len(valid_targets)} validních methyltransferáz s ručně definovaným středem.")
    print(f"[-] Vyřazeno {len(invalid_targets)} non-methyltransferáz (označeno jako '///'):")
    print(f"    {', '.join(invalid_targets)}")

    # Namapujeme PDB soubory podle accession
    pdb_map = {}
    for pdb in struct_dir.glob("*.pdb"):
        acc = parse_accession(pdb.name)
        pdb_map[acc] = pdb

    targets_dir.mkdir(parents=True, exist_ok=True)

    # 1. Odstraníme ze složky targets/ nevalidní proteiny, pokud tam dříve vznikly
    for inv in invalid_targets:
        inv_dir = targets_dir / inv
        if inv_dir.exists():
            print(f"[-] Odstraňuji nevalidní cíl ze složky targets/: '{inv}'")
            shutil.rmtree(inv_dir)

    # 2. Vytvoříme/aktualizujeme validní cíle
    sorted_valid_ids = sorted(list(valid_targets.keys()))
    prepared_ids = []

    for target_id in sorted_valid_ids:
        cx, cy, cz = valid_targets[target_id]
        target_folder = targets_dir / target_id
        target_folder.mkdir(exist_ok=True)

        # Zkopírování protein.pdb
        if target_id in pdb_map:
            dest_pdb = target_folder / "protein.pdb"
            shutil.copy(pdb_map[target_id], dest_pdb)
        else:
            print(f"[!] Varování: Pro validní cíl '{target_id}' nebyl nalezen PDB soubor v '{struct_dir}'!")

        # Zápis upraveného grid.txt s reálným středem kapsy a standardní velikostí
        grid_file = target_folder / "grid.txt"
        with open(grid_file, "w", encoding="utf-8") as gf:
            gf.write(f"center_x = {cx}\n")
            gf.write(f"center_y = {cy}\n")
            gf.write(f"center_z = {cz}\n")
            gf.write(f"size_x = {args.box_size}\n")
            gf.write(f"size_y = {args.box_size}\n")
            gf.write(f"size_z = {args.box_size}\n")

        prepared_ids.append(target_id)

    # 3. Zápis finálního targets.txt pro OpenPBS
    with open(args.targets_list, "w", encoding="utf-8") as f:
        for tid in prepared_ids:
            f.write(f"{tid}\n")

    print(f"\n[+] Hotovo! Vytvořen soubor '{args.targets_list}' s přesně {len(prepared_ids)} cíli.")
    print(f"[+] Každý cíl v '{targets_dir}/' má aktualizovaný 'grid.txt' se středem a rozměry {args.box_size}x{args.box_size}x{args.box_size} Å.")
    print(f"[*] PBS Job Array spusťte příkazem:")
    print(f"    qsub -J 1-{len(prepared_ids)} run_easydock_metacentrum_array.sh")


if __name__ == "__main__":
    main()
