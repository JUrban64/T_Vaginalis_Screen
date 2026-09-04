#!/usr/bin/env python3
"""
fetch_alphafold_structures.py

Skript pro automatické stahování predikovaných struktur (PDB) z AlphaFold Database (AFDB)
pro seznam UniProt identifikátorů (ze souboru Excel .xlsx, CSV, TSV nebo TXT).
Vytváří přehledný report o úspěšně stažených strukturách a podrobný seznam chybějících
proteinů, které v databázi AlphaFold DB nejsou dostupné.

Použití:
    python3 fetch_alphafold_structures.py
    python3 fetch_alphafold_structures.py -i uniprot_list.xlsx -o structures -r reports
    python3 fetch_alphafold_structures.py --format pdb --threads 4
"""

import argparse
import concurrent.futures
import csv
import datetime
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
import xml.etree.ElementTree as ET

AFDB_API_BASE = "https://alphafold.ebi.ac.uk/api/prediction"
DEFAULT_USER_AGENT = "AlphaFold-Fetcher/1.0 (Research Pipeline; mailto:user@example.com)"


def parse_xlsx_stdlib(filepath):
    """
    Načte Excel (.xlsx) soubor pouze pomocí standardní knihovny Pythonu (zipfile + xml),
    takže není nutné instalovat openpyxl ani pandas.
    Vrací seznam slovníků (řádek = dict {sloupec: hodnota}).
    """
    with zipfile.ZipFile(filepath, "r") as z:
        # 1. Načtení sdílených řetězců (sharedStrings.xml), pokud existují
        shared_strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            tree = ET.fromstring(z.read("xl/sharedStrings.xml"))
            ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
            for si in tree.findall(f"{ns}si"):
                t = si.find(f"{ns}t")
                if t is not None and t.text:
                    shared_strings.append(t.text)
                else:
                    parts = [node.text for node in si.findall(f".//{ns}t") if node.text]
                    shared_strings.append("".join(parts))

        # 2. Vyhledání prvního listu
        sheet_filename = "xl/worksheets/sheet1.xml"
        if sheet_filename not in z.namelist():
            sheets = [n for n in z.namelist() if n.startswith("xl/worksheets/sheet") and n.endswith(".xml")]
            if not sheets:
                raise ValueError("V souboru XLSX nebyl nalezen žádný platný list.")
            sheet_filename = sheets[0]

        sheet_tree = ET.fromstring(z.read(sheet_filename))
        ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
        rows = sheet_tree.findall(f"{ns}sheetData/{ns}row")
        if not rows:
            return []

        # 3. Zpracování záhlaví (první řádek)
        header_row = rows[0]
        col_names = {}
        for c in header_row.findall(f"{ns}c"):
            r_attr = c.get("r", "")
            col_letter = "".join(re.findall(r"[A-Z]+", r_attr))
            val = _extract_cell_value(c, ns, shared_strings)
            if val:
                col_names[col_letter] = str(val).strip()

        # 4. Načtení datových řádků
        data = []
        for r in rows[1:]:
            row_dict = {}
            for c in r.findall(f"{ns}c"):
                r_attr = c.get("r", "")
                col_letter = "".join(re.findall(r"[A-Z]+", r_attr))
                col_name = col_names.get(col_letter, col_letter)
                val = _extract_cell_value(c, ns, shared_strings)
                row_dict[col_name] = val
            if any(v is not None and str(v).strip() != "" for v in row_dict.values()):
                data.append(row_dict)

        return data


def _extract_cell_value(cell, ns, shared_strings):
    """Pomocná funkce pro extrakci textu buňky v XLSX."""
    t_type = cell.get("t")
    # inlineStr
    if t_type == "inlineStr":
        t_elem = cell.find(f"{ns}is/{ns}t")
        if t_elem is not None and t_elem.text:
            return t_elem.text.strip()
        parts = [node.text for node in cell.findall(f".//{ns}t") if node.text]
        return "".join(parts).strip() if parts else ""

    v_elem = cell.find(f"{ns}v")
    if v_elem is None or v_elem.text is None:
        return ""

    raw_val = v_elem.text.strip()
    if t_type == "s":
        try:
            idx = int(raw_val)
            return shared_strings[idx] if idx < len(shared_strings) else raw_val
        except (ValueError, IndexError):
            return raw_val

    return raw_val


def load_entries_from_file(input_path):
    """
    Načte záznamy ze zadaného souboru (.xlsx, .csv, .tsv, .txt).
    Vrací seznam slovníků reprezentujících jednotlivé proteiny.
    """
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Vstupní soubor nebyl nalezen: {input_path}")

    ext = os.path.splitext(input_path)[1].lower()
    entries = []

    if ext == ".xlsx":
        # Zkusit openpyxl/pandas pokud existují, jinak spolehlivý stdlib parser
        try:
            import pandas as pd
            df = pd.read_excel(input_path)
            entries = df.fillna("").to_dict(orient="records")
        except Exception:
            entries = parse_xlsx_stdlib(input_path)

    elif ext in [".csv", ".tsv", ".tab"]:
        delimiter = "\t" if ext in [".tsv", ".tab"] else ","
        with open(input_path, "r", encoding="utf-8-sig", errors="replace") as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            entries = [row for row in reader]

    elif ext in [".txt", ""]:
        # Seznam ID řádek po řádku
        with open(input_path, "r", encoding="utf-8-sig", errors="replace") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    parts = line.split()
                    entries.append({"Entry": parts[0]})
    else:
        raise ValueError(f"Nepodporovaný formát souboru: {ext}. Použijte .xlsx, .csv, .tsv nebo .txt.")

    # Normalizace klíče pro UniProt Accession
    normalized = []
    for item in entries:
        row = dict(item)
        accession = ""
        # Hledáme klíč odpovídající UniProt ID
        for candidate in ["Entry", "Accession", "UniProt", "uniprot_id", "ID", "entry"]:
            if candidate in row and row[candidate]:
                accession = str(row[candidate]).strip()
                break
        if not accession and len(row) > 0:
            # Vezmeme první hodnotu v řádku
            first_val = str(list(row.values())[0]).strip()
            # Pokud vypadá jako UniProt accession (6 nebo 10 znaků)
            if re.match(r"^[A-NR-Z0-9]{6,10}$", first_val, re.IGNORECASE):
                accession = first_val

        if accession:
            row["Entry"] = accession.upper()
            normalized.append(row)

    return normalized


def query_alphafold_api(accession, retries=3, delay_sec=0.1):
    """
    Dotáže se AlphaFold DB API na daný UniProt accession.
    Vrací tuple (status_code, data_or_none, error_message).
    """
    clean_acc = accession.strip().upper()
    url = f"{AFDB_API_BASE}/{clean_acc}"
    headers = {"User-Agent": DEFAULT_USER_AGENT, "Accept": "application/json"}
    req = urllib.request.Request(url, headers=headers)

    for attempt in range(retries):
        try:
            time.sleep(delay_sec)
            with urllib.request.urlopen(req, timeout=15) as response:
                if response.status == 200:
                    payload = json.loads(response.read().decode("utf-8"))
                    return 200, payload, None
                return response.status, None, f"Nečekaný HTTP kód {response.status}"
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return 404, None, "Struktura není v AlphaFold DB dostupná (404 Not Found)"
            elif e.code == 400:
                return 400, None, "Neplatný formát identifikátoru (400 Bad Request)"
            elif e.code in [429, 500, 502, 503, 504]:
                wait_time = (attempt + 1) * 2
                time.sleep(wait_time)
                if attempt == retries - 1:
                    return e.code, None, f"Chyba serveru {e.code} po {retries} pokusech"
            else:
                return e.code, None, f"HTTP chyba {e.code}: {e.reason}"
        except urllib.error.URLError as e:
            if attempt == retries - 1:
                return 0, None, f"Chyba sítě: {e.reason}"
            time.sleep((attempt + 1) * 1.5)
        except Exception as e:
            return -1, None, f"Neočekávaná chyba: {str(e)}"

    return -1, None, "Vypršel limit pokusů"


def download_file(url, target_path, retries=3):
    """Stáhne soubor z URL s podporou obnovení a opakování při chybě."""
    headers = {"User-Agent": DEFAULT_USER_AGENT}
    req = urllib.request.Request(url, headers=headers)

    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as response:
                if response.status == 200:
                    temp_path = target_path + ".tmp"
                    with open(temp_path, "wb") as out_f:
                        while True:
                            chunk = response.read(64 * 1024)
                            if not chunk:
                                break
                            out_f.write(chunk)
                    os.replace(temp_path, target_path)
                    return True, None
                return False, f"HTTP status {response.status}"
        except Exception as e:
            if attempt == retries - 1:
                return False, str(e)
            time.sleep((attempt + 1) * 1.5)
    return False, "Chyba stahování po vyčerpání pokusů"


def process_single_entry(row, output_dir, file_format="pdb", overwrite=False, dry_run=False, delay_sec=0.1):
    """
    Zpracuje jeden protein:
    1. Dotáže se API
    2. Stáhne požadované soubory (PDB / CIF / PAE)
    3. Vrátí souhrn o úspěchu nebo chybě
    """
    acc = row.get("Entry", "").strip().upper()
    result = {
        "Entry": acc,
        "Entry Name": row.get("Entry Name", ""),
        "Protein names": row.get("Protein names", ""),
        "Gene Names": row.get("Gene Names", ""),
        "Organism": row.get("Organism", ""),
        "Status": "Unknown",
        "HTTP_Code": 0,
        "pLDDT": None,
        "Sequence_Length": None,
        "PDB_File": "",
        "CIF_File": "",
        "Error_Message": "",
    }

    if not acc:
        result["Status"] = "Error"
        result["Error_Message"] = "Chybí UniProt accession"
        return result

    # 1. Dotaz na API
    status_code, data, err_msg = query_alphafold_api(acc, delay_sec=delay_sec)
    result["HTTP_Code"] = status_code

    if status_code != 200 or not data:
        result["Status"] = "Missing" if status_code in [400, 404] else "Error"
        result["Error_Message"] = err_msg or "Záznam nenalezen"
        return result

    # Výběr kanonického nebo nejdelšího záznamu
    entry_data = None
    for e in data:
        if e.get("uniprotAccession") == acc:
            entry_data = e
            break
    if entry_data is None:
        entry_data = max(data, key=lambda e: e.get("sequenceEnd", 0))

    result["pLDDT"] = entry_data.get("globalMetricValue")
    result["Sequence_Length"] = entry_data.get("sequenceEnd")
    if not result["Gene Names"] and entry_data.get("gene"):
        result["Gene Names"] = entry_data.get("gene")

    pdb_url = entry_data.get("pdbUrl")
    cif_url = entry_data.get("cifUrl")
    pae_url = entry_data.get("paeDocUrl")

    if dry_run:
        result["Status"] = "Found (Dry Run)"
        return result

    # 2. Stahování souborů podle formátu
    download_success = True
    files_to_get = []

    if file_format in ["pdb", "both", "all"] and pdb_url:
        pdb_filename = os.path.basename(pdb_url)
        pdb_path = os.path.join(output_dir, pdb_filename)
        files_to_get.append(("pdb", pdb_url, pdb_path))

    if file_format in ["cif", "both", "all"] and cif_url:
        cif_filename = os.path.basename(cif_url)
        cif_path = os.path.join(output_dir, cif_filename)
        files_to_get.append(("cif", cif_url, cif_path))

    if file_format == "all" and pae_url:
        pae_filename = os.path.basename(pae_url)
        pae_path = os.path.join(output_dir, pae_filename)
        files_to_get.append(("pae", pae_url, pae_path))

    for f_type, url, path in files_to_get:
        if os.path.exists(path) and os.path.getsize(path) > 0 and not overwrite:
            if f_type == "pdb":
                result["PDB_File"] = os.path.basename(path)
            elif f_type == "cif":
                result["CIF_File"] = os.path.basename(path)
            continue

        ok, dl_err = download_file(url, path)
        if ok:
            if f_type == "pdb":
                result["PDB_File"] = os.path.basename(path)
            elif f_type == "cif":
                result["CIF_File"] = os.path.basename(path)
        else:
            download_success = False
            result["Error_Message"] = f"Chyba při stahování {f_type}: {dl_err}"

    if download_success:
        result["Status"] = "Downloaded"
    else:
        result["Status"] = "Download Error"

    return result


def generate_reports(results, reports_dir, input_filename, file_format):
    """
    Vygeneruje souhrnné reporty:
    1. missing_structures.tsv a missing_structures.csv
    2. download_summary.tsv a download_summary.csv
    3. report.md (přehledné lidsky čitelné shrnutí s metrikami)
    """
    os.makedirs(reports_dir, exist_ok=True)

    missing_items = [r for r in results if r["Status"] in ["Missing", "Error", "Download Error"]]
    downloaded_items = [r for r in results if r["Status"] in ["Downloaded", "Found (Dry Run)"]]

    plddt_scores = [r["pLDDT"] for r in downloaded_items if r["pLDDT"] is not None]
    avg_plddt = sum(plddt_scores) / len(plddt_scores) if plddt_scores else 0.0

    # Kategorizace podle spolehlivosti pLDDT
    very_high = sum(1 for s in plddt_scores if s >= 90.0)
    confident = sum(1 for s in plddt_scores if 70.0 <= s < 90.0)
    low = sum(1 for s in plddt_scores if 50.0 <= s < 70.0)
    very_low = sum(1 for s in plddt_scores if s < 50.0)

    # 1. missing_structures.tsv & csv
    missing_fields = ["Entry", "Entry Name", "Protein names", "Gene Names", "Organism", "HTTP_Code", "Error_Message"]
    missing_tsv_path = os.path.join(reports_dir, "missing_structures.tsv")
    missing_csv_path = os.path.join(reports_dir, "missing_structures.csv")

    for path, sep in [(missing_tsv_path, "\t"), (missing_csv_path, ",")]:
        with open(path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=missing_fields, delimiter=sep, extrasaction="ignore")
            writer.writeheader()
            for item in missing_items:
                writer.writerow(item)

    # 2. download_summary.tsv & csv
    summary_fields = [
        "Entry",
        "Entry Name",
        "Gene Names",
        "Status",
        "HTTP_Code",
        "pLDDT",
        "Sequence_Length",
        "PDB_File",
        "Error_Message",
        "Protein names",
    ]
    summary_tsv_path = os.path.join(reports_dir, "download_summary.tsv")
    summary_csv_path = os.path.join(reports_dir, "download_summary.csv")

    for path, sep in [(summary_tsv_path, "\t"), (summary_csv_path, ",")]:
        with open(path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter=sep, extrasaction="ignore")
            writer.writeheader()
            for item in results:
                writer.writerow(item)

    # 3. report.md
    now_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    md_path = os.path.join(reports_dir, "report.md")

    with open(md_path, "w", encoding="utf-8") as md:
        md.write("# AlphaFold Database - Zpráva o stahování struktur\n\n")
        md.write(f"- **Datum spuštění:** `{now_str}`\n")
        md.write(f"- **Vstupní soubor:** `{os.path.basename(input_filename)}`\n")
        md.write(f"- **Zvolený formát:** `{file_format.upper()}`\n")
        md.write(f"- **Celkem zpracováno proteinů:** **{len(results)}**\n\n")

        md.write("## Souhrnná statistika\n\n")
        md.write("| Metrika | Počet | Podíl (%) |\n")
        md.write("| :--- | :---: | :---: |\n")
        total = len(results)
        dl_count = len(downloaded_items)
        miss_count = len(missing_items)
        dl_pct = (dl_count / total * 100) if total else 0
        miss_pct = (miss_count / total * 100) if total else 0

        md.write(f"| **Úspěšně dostupné / stažené** | **{dl_count}** | **{dl_pct:.1f} %** |\n")
        md.write(f"| **Chybějící v AlphaFold DB** | **{miss_count}** | **{miss_pct:.1f} %** |\n")
        md.write(f"| **Celkem v seznamu** | **{total}** | **100.0 %** |\n\n")

        if plddt_scores:
            md.write("## Distribuce kvality modelů (pLDDT)\n\n")
            md.write(f"- **Průměrné globální pLDDT skóre:** **`{avg_plddt:.2f}`**\n\n")
            md.write("| Interval spolehlivosti | Kategorie | Počet struktur | Podíl ze stažených |\n")
            md.write("| :--- | :--- | :---: | :---: |\n")
            md.write(f"| **pLDDT >= 90** | Velmi vysoká spolehlivost (Very high) | {very_high} | {very_high/dl_count*100:.1f} % |\n")
            md.write(f"| **70 <= pLDDT < 90** | Spolehlivá páteř (Confident) | {confident} | {confident/dl_count*100:.1f} % |\n")
            md.write(f"| **50 <= pLDDT < 70** | Nízká spolehlivost (Low) | {low} | {low/dl_count*100:.1f} % |\n")
            md.write(f"| **pLDDT < 50** | Velmi nízká / nestrukturovaná (Disordered) | {very_low} | {very_low/dl_count*100:.1f} % |\n\n")

        md.write("## Seznam chybějících struktur (nejsou v AlphaFold DB)\n\n")
        if not missing_items:
            md.write("Všechny struktury byly v AlphaFold DB úspěšně nalezeny.\n")
        else:
            md.write("Následující proteiny nemají v AlphaFold DB k dispozici 3D model:\n\n")
            md.write("| UniProt Entry | Entry Name | Gen | Důvod / Stav | Název proteinu |\n")
            md.write("| :--- | :--- | :--- | :--- | :--- |\n")
            for m in missing_items:
                entry = m.get("Entry", "")
                name = m.get("Entry Name", "-")
                gene = m.get("Gene Names", "-")
                reason = m.get("Error_Message", "404 Not Found")
                prot = m.get("Protein names", "-")
                if len(prot) > 60:
                    prot = prot[:57] + "..."
                md.write(f"| [`{entry}`](https://www.uniprot.org/uniprotkb/{entry}) | {name} | {gene} | {reason} | {prot} |\n")
            md.write("\n")

        md.write("## Vygenerované soubory\n\n")
        md.write(f"- **Adresář se strukturami:** [`{reports_dir}/../structures/`](file://{os.path.abspath(reports_dir)}/../structures)\n")
        md.write(f"- **Chybějící struktury (TSV):** [`{missing_tsv_path}`](file://{missing_tsv_path})\n")
        md.write(f"- **Chybějící struktury (CSV):** [`{missing_csv_path}`](file://{missing_csv_path})\n")
        md.write(f"- **Kompletní přehled (TSV):** [`{summary_tsv_path}`](file://{summary_tsv_path})\n")
        md.write(f"- **Kompletní přehled (CSV):** [`{summary_csv_path}`](file://{summary_csv_path})\n")

    return {
        "total": total,
        "downloaded": dl_count,
        "missing": miss_count,
        "missing_tsv": missing_tsv_path,
        "summary_tsv": summary_tsv_path,
        "report_md": md_path,
    }


def find_default_input_file():
    """Automaticky vyhledá vhodný vstupní soubor v pracovním adresáři."""
    candidates = [
        "uniprotkb_taxonomy_id_5722_AND_methyltr_2026_09_04.xlsx",
    ]
    for c in candidates:
        if os.path.exists(c):
            return c

    # Hledáme jakýkoliv .xlsx
    for f in os.listdir("."):
        if f.endswith(".xlsx") and not f.startswith("~$"):
            return f
    for f in os.listdir("."):
        if f.endswith(".tsv") or f.endswith(".csv"):
            return f
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Stahovač predikovaných struktur z AlphaFold DB s reportem chybějících záznamů."
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Vstupní soubor (.xlsx, .csv, .tsv, .txt). Výchozí: automaticky vyhledá .xlsx v aktuální složce.",
        default=None,
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help="Cílová složka pro stažené 3D struktury (výchozí: ./structures).",
        default="structures",
    )
    parser.add_argument(
        "-r",
        "--reports-dir",
        help="Cílová složka pro vygenerované reporty (výchozí: ./reports).",
        default="reports",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=["pdb", "cif", "both", "all"],
        default="pdb",
        help="Formát struktur ke stažení: 'pdb' (výchozí), 'cif', 'both' (PDB i CIF), 'all' (PDB, CIF i PAE).",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Počet paralelních vláken pro dotazování a stahování (výchozí: 4).",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.1,
        help="Pauza mezi požadavky v sekundách na vlákno (výchozí: 0.1s).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Zpracovat pouze prvních N proteinů (vhodné pro rychlý test).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Pouze ověří dostupnost v AlphaFold DB bez stahování velkých souborů.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Přepsat existující stažené soubory.",
    )

    args = parser.parse_args()

    input_file = args.input or find_default_input_file()
    if not input_file:
        print("[!] Chyba: Nebyl nalezen žádný vstupní soubor (.xlsx, .csv). Zadejte cestu pomocí -i <cesta>.")
        sys.exit(1)

    print("=" * 70)
    print("      AlphaFold DB - Stahovač struktur a generátor reportu")
    print("=" * 70)
    print(f"[*] Vstupní soubor : {input_file}")
    print(f"[*] Cílová složka  : {args.output_dir}")
    print(f"[*] Složka reportů : {args.reports_dir}")
    print(f"[*] Formát souborů : {args.format.upper()}")
    print(f"[*] Paralelní běh  : {args.threads} vláken (prodleva {args.delay} s)")
    if args.dry_run:
        print("[*] REŽIM DRY-RUN  : Struktury nebudou ukládány na disk.")
    print("-" * 70)

    # 1. Načtení vstupu
    print(f"[*] Načítám vstupní soubor '{input_file}'...")
    try:
        entries = load_entries_from_file(input_file)
    except Exception as e:
        print(f"[!] Chyba při načítání vstupního souboru: {e}")
        sys.exit(1)

    total_loaded = len(entries)
    if total_loaded == 0:
        print("[!] Ze vstupního souboru se nepodařilo načíst žádné UniProt záznamy.")
        sys.exit(1)

    if args.limit:
        entries = entries[: args.limit]
        print(f"[*] Omezeno na prvních {len(entries)} záznamů (--limit {args.limit}).")
    else:
        print(f"[*] Úspěšně načteno {total_loaded} proteinů.")

    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.reports_dir, exist_ok=True)

    # 2. Zpracování záznamů
    results = []
    completed_count = 0
    start_time = time.time()

    print(f"\n[*] Zahajuji dotazování AlphaFold DB API a stahování...")

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_entry = {
            executor.submit(
                process_single_entry,
                row,
                args.output_dir,
                file_format=args.format,
                overwrite=args.overwrite,
                dry_run=args.dry_run,
                delay_sec=args.delay,
            ): row
            for row in entries
        }

        for future in concurrent.futures.as_completed(future_to_entry):
            res = future.result()
            results.append(res)
            completed_count += 1

            status_symbol = "[+]" if res["Status"] in ["Downloaded", "Found (Dry Run)"] else "[!]"
            info_msg = ""
            if res["Status"] in ["Downloaded", "Found (Dry Run)"]:
                info_msg = f"pLDDT: {res['pLDDT']:.1f}" if res["pLDDT"] is not None else "OK"
            else:
                info_msg = f"{res['Status']} ({res['Error_Message']})"

            pct = (completed_count / len(entries)) * 100
            print(f"  [{completed_count:3d}/{len(entries):3d}] ({pct:5.1f}%) {status_symbol} {res['Entry']:<10} -> {info_msg}")

    # Seřadit výsledky podle původního pořadí v souboru
    entry_order = {row["Entry"]: idx for idx, row in enumerate(entries)}
    results.sort(key=lambda r: entry_order.get(r["Entry"], 999999))

    elapsed = time.time() - start_time
    print("-" * 70)
    print(f"[*] Dotazování a stahování dokončeno za {elapsed:.1f} s.")

    # 3. Vygenerování reportů
    print(f"[*] Generuji reporty v '{args.reports_dir}'...")
    rep_info = generate_reports(results, args.reports_dir, input_file, args.format)

    print("\n" + "=" * 70)
    print("                       VÝSLEDKY BĚHU")
    print("=" * 70)
    print(f"  Celkem zpracováno proteinů : {rep_info['total']}")
    print(f"  Úspěšně staženo            : {rep_info['downloaded']} ({(rep_info['downloaded']/rep_info['total']*100):.1f} %)")
    print(f"  Chybějící v AlphaFold DB   : {rep_info['missing']} ({(rep_info['missing']/rep_info['total']*100):.1f} %)")
    print("-" * 70)
    print(f"  Souhrnný Markdown report   : {rep_info['report_md']}")
    print(f"  Tabulka chybějících (TSV)  : {rep_info['missing_tsv']}")
    print(f"  Kompletní přehled (TSV)    : {rep_info['summary_tsv']}")
    print(f"  Složka se strukturami      : {os.path.abspath(args.output_dir)}")
    print("=" * 70)


if __name__ == "__main__":
    main()
