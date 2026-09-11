#!/usr/bin/env python3
"""
preinit_ligands.py

Lokální skript pro jednorázovou předpřípravu knihovny ligandů v EasyDocku.
Pro menší knihovny (~100-500 molekul) není nutné zadávat PBS job na MetaCentru;
tento skript spustí EasyDock jako podproces přímo na lokálním stroji nebo na frontendovém uzlu.

Co provádí:
1. Validace a sanitizace molekul, odstranění solí.
2. Generování stereoisomerů pro nedefinovaná chirální centra (-s).
3. Generování 3D konformací.
4. Protonace při zadaném pH (--protonation molgpka / unipka).
5. Uložení do master šablony 'ligands_template.db'.

Použití:
    python3 preinit_ligands.py -i ligands.smi -o ligands_template.db
    python3 preinit_ligands.py -i ligands.smi -c 8 --protonation molgpka -s 4
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


def check_easydock():
    """Zkontroluje, zda je EasyDock dostupný v PATH."""
    easydock_bin = shutil.which("easydock")
    if not easydock_bin:
        print("[!] CHYBA: Nástroj 'easydock' nebyl nalezen v PATH!", file=sys.stderr)
        print("    Ujistěte se, že máte aktivované příslušné conda prostředí (např. 'conda activate easydock').", file=sys.stderr)
        sys.exit(1)
    return easydock_bin


def main():
    parser = argparse.ArgumentParser(
        description="Jednorázová příprava ligandů do SQLite šablony pomocí EasyDock."
    )
    parser.add_argument(
        "-i", "--input",
        default="ligands.smi",
        help="Vstupní soubor s ligandy (.smi nebo .sdf, default: ligands.smi)"
    )
    parser.add_argument(
        "-o", "--output",
        default="ligands_template.db",
        help="Výstupní SQLite databáze (default: ligands_template.db)"
    )
    parser.add_argument(
        "-c", "--ncpu",
        type=int,
        default=min(os.cpu_count() or 4, 8),
        help=f"Počet CPU jader (default: {min(os.cpu_count() or 4, 8)})"
    )
    parser.add_argument(
        "-s", "--stereoisomers",
        type=int,
        default=4,
        help="Maximální počet stereoisomerů na nedefinované centrum (default: 4)"
    )
    parser.add_argument(
        "-p", "--protonation",
        default="molgpka",
        choices=["molgpka", "molgpka_fix", "none", "chemaxon"],
        help="Metoda protonace (default: molgpka)"
    )
    parser.add_argument(
        "--pH",
        type=float,
        default=7.4,
        help="Cílové pH pro protonaci (default: 7.4)"
    )
    parser.add_argument(
        "--ring-sample",
        action="store_true",
        help="Zapnout vzorkování konformací nasycených kruhů (--ring_sample)"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Přepsat existující výstupní databázi bez ptaní"
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        print(f"[!] CHYBA: Vstupní soubor '{input_path}' neexistuje!", file=sys.stderr)
        sys.exit(1)

    if output_path.exists():
        if args.force:
            print(f"[*] Výstupní soubor '{output_path}' již existuje a bude smazán (--force).")
            output_path.unlink()
        else:
            print(f"[!] Varování: Výstupní databáze '{output_path}' již existuje.")
            print("    EasyDock neumožňuje přepsat existující databázi.")
            resp = input("    Chcete ji smazat a vytvořit znovu? [y/N]: ").strip().lower()
            if resp in ("y", "yes", "ano", "a"):
                output_path.unlink()
                print(f"[*] Smazáno '{output_path}'.")
            else:
                print("[*] Akce zrušena uživatelem.")
                sys.exit(0)

    easydock_cmd = check_easydock()

    # Sestavení příkazu pro EasyDock
    cmd = [
        easydock_cmd,
        "-i", str(input_path),
        "-o", str(output_path),
        "-c", str(args.ncpu),
    ]

    if args.stereoisomers > 1:
        cmd.extend(["-s", str(args.stereoisomers)])

    if args.protonation != "none":
        cmd.extend(["--protonation", args.protonation])
        cmd.extend(["--pH", str(args.pH)])

    if args.ring_sample:
        cmd.append("--ring_sample")

    print("=" * 70)
    print("Spouštím EasyDock předpřípravu ligandů:")
    print(f"  Vstupní soubor:   {input_path}")
    print(f"  Výstupní DB:      {output_path}")
    print(f"  CPU vláken:       {args.ncpu}")
    print(f"  Stereoisomery:    až {args.stereoisomers}")
    print(f"  Protonace:        {args.protonation} (pH {args.pH})")
    print(f"  Příkaz:           {' '.join(cmd)}")
    print("=" * 70)

    start_time = time.time()

    try:
        # Spuštění jako podproces se streamováním výstupu v reálném čase
        proc = subprocess.run(cmd, check=True)
        elapsed = time.time() - start_time
        print("\n" + "=" * 70)
        print(f"[+] Hotovo za {elapsed:.1f} sekund ({elapsed / 60:.2f} minut)!")
        print(f"[+] Ligandy jsou úspěšně uloženy v šabloně: '{output_path}'")
        print(f"[*] Nyní můžete spustit docking do všech cílů na MetaCentru přes:")
        print(f"    qsub run_easydock_metacentrum_array.sh")
        print("=" * 70)
    except subprocess.CalledProcessError as e:
        print(f"\n[!] EasyDock skončil s chybou (kód {e.returncode})!", file=sys.stderr)
        sys.exit(e.returncode)
    except KeyboardInterrupt:
        print("\n[!] Přerušeno uživatelem.", file=sys.stderr)
        sys.exit(130)


if __name__ == "__main__":
    main()
