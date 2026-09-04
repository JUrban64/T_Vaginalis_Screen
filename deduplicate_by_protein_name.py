#!/usr/bin/env python3
"""
deduplicate_by_protein_name.py

Skript pro odstranění redundancí v seznamu methyltransferáz (Excel .xlsx, CSV, TSV)
čistě na základě anotace v sloupci 'Protein names'.

Vlastnosti:
1. Seskupí proteiny se shodným názvem (Protein names).
2. Jako reprezentanta skupiny vybere nejkvalitnější záznam:
   - Prioritně podle nejvyššího pLDDT (pokud existují stažené struktury nebo download_summary.tsv)
   - Sekundárně podle délky sekvence nebo prvního výskytu v tabulce.
3. Zachová všechny původní sloupce a přidá auditní sloupce:
   - Duplicate_Count (kolik duplikátů bylo sloučeno)
   - All_Merged_Entries (výčet všech UniProt ID)
   - All_Merged_Genes (výčet všech genů)
   - Representative_pLDDT (skóre kvality reprezentanta)
4. Vygeneruje nový deduplikovaný Excel soubor (.xlsx), TSV tabulku a přehledný Markdown report.

Použití:
    python3 deduplicate_by_protein_name.py
    python3 deduplicate_by_protein_name.py -i uniprotkb_taxonomy_id_5722_AND_methyltr_2026_09_04.xlsx -o deduplikovane.xlsx
    python3 deduplicate_by_protein_name.py --keep-uncharacterized-separate
"""

import argparse
import csv
import datetime
import os
import re
import sys
import xml.etree.ElementTree as ET
import xml.sax.saxutils
import zipfile


def parse_xlsx_stdlib(filepath):
    """Načte Excel (.xlsx) pomocí standardní knihovny Pythonu."""
    with zipfile.ZipFile(filepath, "r") as z:
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
            return [], []

        # Hlavička
        header_row = rows[0]
        col_names = {}
        headers_order = []
        for c in header_row.findall(f"{ns}c"):
            r_attr = c.get("r", "")
            col_letter = "".join(re.findall(r"[A-Z]+", r_attr))
            val = _extract_cell_value(c, ns, shared_strings)
            if val:
                val_str = str(val).strip()
                col_names[col_letter] = val_str
                headers_order.append(val_str)

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

        return data, headers_order


def _extract_cell_value(cell, ns, shared_strings):
    t_type = cell.get("t")
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


def load_input_data(filepath):
    """Načte vstupní data a vrátí seznam slovníků a původní pořadí sloupců."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Vstupní soubor nebyl nalezen: {filepath}")

    ext = os.path.splitext(filepath)[1].lower()
    if ext == ".xlsx":
        try:
            import pandas as pd
            df = pd.read_excel(filepath)
            headers = list(df.columns)
            records = df.fillna("").to_dict(orient="records")
            return records, headers
        except Exception:
            return parse_xlsx_stdlib(filepath)

    elif ext in [".csv", ".tsv", ".tab"]:
        delim = "\t" if ext in [".tsv", ".tab"] else ","
        with open(filepath, "r", encoding="utf-8-sig", errors="replace") as f:
            reader = csv.DictReader(f, delimiter=delim)
            headers = list(reader.fieldnames) if reader.fieldnames else []
            records = [row for row in reader]
            return records, headers
    else:
        raise ValueError(f"Nepodporovaný formát: {ext}")


def col_index_to_letter(idx):
    """Převede číslo sloupce (1-based) na písmeno Excelu (1 -> A, 27 -> AA)."""
    result = ""
    while idx > 0:
        idx, remainder = divmod(idx - 1, 26)
        result = chr(65 + remainder) + result
    return result


def write_xlsx_stdlib(filepath, headers, rows):
    """Zapíše tabulku do platného XLSX souboru pouze pomocí vestavěné knihovny zipfile."""
    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
</Types>"""

    rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>"""

    wb = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="Deduplicated_Proteins" sheetId="1" r:id="rId1"/>
  </sheets>
</workbook>"""

    wb_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
</Relationships>"""

    sheet_lines = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>',
        '<row r="1">',
    ]
    for col_idx, h in enumerate(headers, 1):
        col_letter = col_index_to_letter(col_idx)
        esc_val = xml.sax.saxutils.escape(str(h))
        sheet_lines.append(f'<c r="{col_letter}1" t="inlineStr"><is><t>{esc_val}</t></is></c>')
    sheet_lines.append("</row>")

    for row_idx, r in enumerate(rows, 2):
        sheet_lines.append(f'<row r="{row_idx}">')
        for col_idx, val in enumerate(r, 1):
            col_letter = col_index_to_letter(col_idx)
            val_str = "" if val is None else str(val)
            esc_val = xml.sax.saxutils.escape(val_str)
            sheet_lines.append(f'<c r="{col_letter}{row_idx}" t="inlineStr"><is><t>{esc_val}</t></is></c>')
        sheet_lines.append("</row>")

    sheet_lines.append("</sheetData></worksheet>")
    sheet_xml = "".join(sheet_lines)

    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    with zipfile.ZipFile(filepath, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", rels)
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        z.writestr("xl/worksheets/sheet1.xml", sheet_xml)


def load_plddt_metadata():
    """
    Pokusí se načíst pLDDT skóre z reports/download_summary.tsv nebo přímo ze souborů structures/*.pdb.
    Vrací dict {accession: plddt_float}.
    """
    plddt_map = {}
    summary_tsv = os.path.join("reports", "download_summary.tsv")
    if os.path.exists(summary_tsv):
        try:
            with open(summary_tsv, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for r in reader:
                    acc = r.get("Entry", "").strip()
                    val = r.get("pLDDT", "")
                    if acc and val:
                        try:
                            plddt_map[acc] = float(val)
                        except ValueError:
                            pass
        except Exception:
            pass

    # Pokud není v summary, zkusíme projít složku structures/
    if not plddt_map and os.path.exists("structures"):
        for fname in os.listdir("structures"):
            if fname.endswith(".pdb"):
                acc_match = re.search(r"AF-([A-NR-Z0-9]+)-F", fname)
                if acc_match:
                    acc = acc_match.group(1)
                    pdb_p = os.path.join("structures", fname)
                    bfactors = []
                    try:
                        with open(pdb_p, "r", encoding="utf-8", errors="replace") as f:
                            for line in f:
                                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                                    bfactors.append(float(line[60:66].strip()))
                        if bfactors:
                            plddt_map[acc] = round(sum(bfactors) / len(bfactors), 2)
                    except Exception:
                        pass

    return plddt_map


def normalize_protein_name(name, clean_primary=False):
    """Normalizuje text názvu proteinu pro porovnání."""
    if not name:
        return ""
    text = str(name).strip()
    # Nahradit vícenásobné mezery jednou
    text = re.sub(r"\s+", " ", text)

    if clean_primary:
        # Odstranit alternativní názvy v závorkách a EC čísla: "tRNA methyl... (EC 2.1.1.33) (...)" -> "tRNA methyl..."
        # Ponechat pouze první část před první otevírací závorkou, pokud je smysluplná
        parts = re.split(r"\s*\(", text, maxsplit=1)
        if parts and len(parts[0].strip()) > 3:
            text = parts[0].strip()

    return text.lower()


def deduplicate_entries(entries, name_col="Protein names", keep_uncharacterized_separate=False, clean_primary=False):
    """
    Seskupí položky podle normalizovaného Protein names.
    Vybere nejlepšího reprezentanta každé skupiny.
    """
    plddt_map = load_plddt_metadata()

    # Skupiny: norm_name -> list of row dicts
    groups = {}
    uncharacterized_counter = 0

    for idx, row in enumerate(entries):
        raw_name = str(row.get(name_col, "")).strip()
        entry_id = str(row.get("Entry", f"ID_{idx}")).strip()

        # Ošetření 'Uncharacterized protein'
        is_uncharacterized = "uncharacterized" in raw_name.lower()
        if is_uncharacterized and keep_uncharacterized_separate:
            uncharacterized_counter += 1
            group_key = f"__uncharacterized_{uncharacterized_counter}_{entry_id}__"
        else:
            group_key = normalize_protein_name(raw_name, clean_primary=clean_primary)
            if not group_key:
                group_key = f"__unnamed_{idx}_{entry_id}__"

        if group_key not in groups:
            groups[group_key] = []
        groups[group_key].append((idx, row))

    deduplicated = []
    group_stats = []

    for group_key, row_list in groups.items():
        # Výběr reprezentanta:
        # 1. nejvyšší pLDDT (pokud známo)
        # 2. nejnižší index (původní pořadí)
        def rep_key(item):
            orig_idx, r = item
            acc = str(r.get("Entry", "")).strip()
            plddt = plddt_map.get(acc, -1.0)
            return (plddt, -orig_idx)

        best_idx, best_row = max(row_list, key=rep_key)
        rep_acc = str(best_row.get("Entry", "")).strip()
        rep_plddt = plddt_map.get(rep_acc, None)

        # Seznam sloučených identifikátorů a genů
        all_entries = [str(r.get("Entry", "")).strip() for _, r in row_list if r.get("Entry")]
        all_genes = []
        for _, r in row_list:
            g = str(r.get("Gene Names", "")).strip()
            if g and g not in all_genes:
                all_genes.append(g)

        merged_row = dict(best_row)
        merged_row["Duplicate_Count"] = len(row_list)
        merged_row["Representative_pLDDT"] = rep_plddt if rep_plddt is not None else "N/A"
        merged_row["All_Merged_Entries"] = ", ".join(all_entries)
        merged_row["All_Merged_Genes"] = ", ".join(all_genes)

        deduplicated.append(merged_row)

        group_stats.append({
            "group_key": group_key,
            "display_name": best_row.get(name_col, ""),
            "representative_entry": rep_acc,
            "representative_gene": best_row.get("Gene Names", ""),
            "plddt": rep_plddt,
            "count": len(row_list),
            "members": all_entries,
        })

    # Seřadit zpět podle Duplicate_Count (největší rodiny nahoře) nebo podle původního pořadí
    group_stats.sort(key=lambda x: x["count"], reverse=True)

    return deduplicated, group_stats


def generate_markdown_report(group_stats, total_orig, total_dedup, output_path, excel_path):
    """Vygeneruje souhrnnou zprávu v Markdownu."""
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    now_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    removed = total_orig - total_dedup
    red_pct = (removed / total_orig * 100) if total_orig else 0

    with open(output_path, "w", encoding="utf-8") as md:
        md.write("# Redukce Redundance Methyltransferáz dle Názvu Proteinu (Protein Names)\n\n")
        md.write(f"- **Datum spuštění:** `{now_str}`\n")
        md.write(f"- **Původní počet záznamů v Excelu:** **{total_orig}**\n")
        md.write(f"- **Výsledný počet unikátních enzymů:** **{total_dedup}**\n")
        md.write(f"- **Odstraněno redundantních duplikátů:** **{removed}** (redukce o **{red_pct:.1f} %**)\n")
        md.write(f"- **Výstupní deduplikovaný Excel:** [`{os.path.basename(excel_path)}`](file://{os.path.abspath(excel_path)})\n\n")

        md.write("## 1. Souhrnná statistika\n\n")
        md.write("| Metrika | Hodnota |\n")
        md.write("| :--- | :---: |\n")
        md.write(f"| Celkem řádků před deduplikací | {total_orig} |\n")
        md.write(f"| Zachovaných reprezentantů | {total_dedup} |\n")
        md.write(f"| Odstraněných duplikátů | {removed} |\n")
        md.write(f"| Míra redukce | **{red_pct:.1f} %** |\n")
        multi_groups = [g for g in group_stats if g["count"] > 1]
        md.write(f"| Počet proteinových rodin s více duplikáty | {len(multi_groups)} |\n")
        md.write(f"| Počet unikátních proteinů (pouze 1 výskyt) | {len(group_stats) - len(multi_groups)} |\n\n")

        md.write("## 2. Největší rodiny duplikovaných methyltransferáz v *T. vaginalis*\n\n")
        md.write("| Počet duplikátů | Reprezentant | Gen | pLDDT | Název proteinu (Protein names) |\n")
        md.write("| :---: | :--- | :--- | :---: | :--- |\n")
        for g in multi_groups[:20]:
            plddt_str = f"{g['plddt']:.1f}" if g["plddt"] is not None else "-"
            name = g["display_name"]
            if len(name) > 65:
                name = name[:62] + "..."
            md.write(f"| **{g['count']}×** | [`{g['representative_entry']}`](https://www.uniprot.org/uniprotkb/{g['representative_entry']}) | {g['representative_gene']} | {plddt_str} | {name} |\n")
        md.write("\n")

        md.write("## 3. Výstupní soubory\n\n")
        md.write(f"- **Deduplikovaný Excel:** [`{excel_path}`](file://{os.path.abspath(excel_path)})\n")
        md.write(f"- **Tento report:** [`{output_path}`](file://{os.path.abspath(output_path)})\n")


def main():
    parser = argparse.ArgumentParser(
        description="Odstranění redundancí v seznamu methyltransferáz na základě sloupce 'Protein names'."
    )
    parser.add_argument(
        "-i",
        "--input",
        default="uniprotkb_taxonomy_id_5722_AND_methyltr_2026_09_04.xlsx",
        help="Cesta ke vstupnímu souboru (.xlsx, .csv, .tsv). Výchozí: uniprotkb_taxonomy_id_5722_AND_methyltr_2026_09_04.xlsx.",
    )
    parser.add_argument(
        "-o",
        "--output-excel",
        default="uniprot_methyltransferases_deduplicated.xlsx",
        help="Název výstupního deduplikovaného Excelu (výchozí: uniprot_methyltransferases_deduplicated.xlsx).",
    )
    parser.add_argument(
        "--output-tsv",
        default="reports/deduplicated_proteins.tsv",
        help="Cesta pro uložení TSV tabulky (výchozí: reports/deduplicated_proteins.tsv).",
    )
    parser.add_argument(
        "--name-column",
        default="Protein names",
        help="Název sloupce s popisem proteinu (výchozí: 'Protein names').",
    )
    parser.add_argument(
        "--keep-uncharacterized-separate",
        action="store_true",
        help="Ponechat proteiny s názvem 'Uncharacterized protein' jako samostatné unikátní položky (neslévat je do 1 řádku).",
    )
    parser.add_argument(
        "--clean-primary-name",
        action="store_true",
        help="Seskupovat podle primárního názvu (odstraní synonyma a EC čísla v závorkách před porovnáním).",
    )

    args = parser.parse_args()

    print("=" * 75)
    print("      Deduplikace Methyltransferáz dle 'Protein names'")
    print("=" * 75)
    print(f"[*] Vstupní soubor : {args.input}")
    print(f"[*] Výstupní Excel : {args.output_excel}")
    print(f"[*] Výstupní TSV   : {args.output_tsv}")
    print(f"[*] Sloupec názvu  : {args.name_column}")
    if args.keep_uncharacterized_separate:
        print("[*] Uncharacterized: Každý necharakterizovaný protein zůstane zachován samostatně.")
    else:
        print("[*] Uncharacterized: Slévají se pod jednotný název (výchozí).")
    print("-" * 75)

    # 1. Načtení dat
    print(f"[*] Načítám data z '{args.input}'...")
    entries, headers = load_input_data(args.input)
    print(f"[*] Načteno {len(entries)} záznamů s {len(headers)} sloupci.")

    if args.name_column not in headers:
        # Pokus najít alternativu
        candidates = [h for h in headers if "protein" in h.lower() and "name" in h.lower()]
        if candidates:
            args.name_column = candidates[0]
            print(f"[*] Sloupec nalezen jako: '{args.name_column}'")
        else:
            print(f"[!] Chyba: Sloupec '{args.name_column}' nebyl v tabulce nalezen!", file=sys.stderr)
            sys.exit(1)

    # 2. Deduplikace
    print(f"[*] Provádím deduplikaci podle sloupce '{args.name_column}'...")
    deduplicated, group_stats = deduplicate_entries(
        entries,
        name_col=args.name_column,
        keep_uncharacterized_separate=args.keep_uncharacterized_separate,
        clean_primary=args.clean_primary_name,
    )

    orig_count = len(entries)
    dedup_count = len(deduplicated)
    removed_count = orig_count - dedup_count
    red_pct = (removed_count / orig_count * 100) if orig_count else 0

    print(f"[+] Hotovo. Původní počet: {orig_count} -> Výsledný počet: {dedup_count}")
    print(f"[+] Odstraněno {removed_count} duplikátů (redukce o {red_pct:.1f} %).")

    # 3. Příprava výstupních sloupců
    audit_cols = ["Duplicate_Count", "Representative_pLDDT", "All_Merged_Entries", "All_Merged_Genes"]
    out_headers = [h for h in headers if h not in audit_cols] + audit_cols

    out_rows = []
    for r in deduplicated:
        row_vals = [r.get(h, "") for h in out_headers]
        out_rows.append(row_vals)

    # 4. Uložení XLSX
    print(f"[*] Ukládám deduplikovaný Excel: '{args.output_excel}'...")
    write_xlsx_stdlib(args.output_excel, out_headers, out_rows)

    # 5. Uložení TSV
    os.makedirs(os.path.dirname(os.path.abspath(args.output_tsv)), exist_ok=True)
    print(f"[*] Ukládám deduplikovaný TSV: '{args.output_tsv}'...")
    with open(args.output_tsv, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(out_headers)
        writer.writerows(out_rows)

    # 6. Generování Markdown reportu
    report_md_path = os.path.join("reports", "deduplication_by_name_report.md")
    print(f"[*] Generuji Markdown zprávu: '{report_md_path}'...")
    generate_markdown_report(group_stats, orig_count, dedup_count, report_md_path, args.output_excel)

    print("\n" + "=" * 75)
    print("                       VÝSLEDKY DEDUPLIKACE")
    print("=" * 75)
    print(f"  Původní počet proteinů      : {orig_count}")
    print(f"  Unikátních po deduplikaci   : {dedup_count}")
    print(f"  Odstraněných duplikátů      : {removed_count} ({red_pct:.1f} %)")
    print(f"  Výstupní Excel soubor       : {os.path.abspath(args.output_excel)}")
    print(f"  Výstupní TSV soubor         : {os.path.abspath(args.output_tsv)}")
    print(f"  Podrobný Markdown report    : {os.path.abspath(report_md_path)}")
    print("=" * 75)


if __name__ == "__main__":
    main()
