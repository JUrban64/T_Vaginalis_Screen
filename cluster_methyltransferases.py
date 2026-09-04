#!/usr/bin/env python3
"""
cluster_methyltransferases.py

Skript pro odstranění strukturně redundantních methyltransferáz pomocí Foldseeku
a metriky TM-score.

Klíčové vlastnosti:
1. Spustí Foldseek all-vs-all strukturní alignment (--alignment-type 1 pro TM-score).
2. Shlukuje struktury na základě zadaného TM-score a coverage thresholdů.
3. Pro každý klastr vybere nejkvalitnější reprezentativní strukturu (nejvyšší pLDDT skóre).
4. Vytvoří novou složku obsahující pouze reprezentativní PDB struktury.
5. Umožňuje otestovat citlivost různých TM-score prahů (--explore-thresholds).
6. Generuje detailní TSV a Markdown reporty o klastrech a redukci redundance.

Použití na HPC:
    python3 cluster_methyltransferases.py -i structures -o structures_representative
    python3 cluster_methyltransferases.py -i structures -o structures_representative --tmscore-threshold 0.75 --explore-thresholds
"""

import argparse
import csv
import datetime
import os
import re
import shutil
import subprocess
import sys


def extract_pdb_metadata(pdb_path):
    """
    Rychle extrahuje průměrné CA pLDDT skóre a délku sekvence přímo z PDB souboru.
    U AlphaFold modelů je per-reziduální pLDDT uloženo ve sloupci B-factor (pozice 61-66).
    """
    ca_bfactors = []
    seen_res = set()

    try:
        with open(pdb_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    res_id = (line[21:22].strip(), line[22:27].strip())  # chain, resSeq
                    if res_id not in seen_res:
                        seen_res.add(res_id)
                        try:
                            bf = float(line[60:66].strip())
                            ca_bfactors.append(bf)
                        except ValueError:
                            pass
    except Exception as e:
        print(f"[!] Varování: Nelze číst {pdb_path}: {e}", file=sys.stderr)

    mean_plddt = sum(ca_bfactors) / len(ca_bfactors) if ca_bfactors else 0.0
    seq_len = len(ca_bfactors)

    # Pokus o extrakci UniProt Accession z názvu souboru (např. AF-A2DE67-F1-model_v6.pdb -> A2DE67)
    basename = os.path.basename(pdb_path)
    acc_match = re.search(r"AF-([A-NR-Z0-9]+)-F", basename)
    accession = acc_match.group(1) if acc_match else os.path.splitext(basename)[0]

    return {
        "pdb_path": pdb_path,
        "filename": basename,
        "accession": accession,
        "plddt": round(mean_plddt, 2),
        "length": seq_len,
    }


def load_all_pdb_metadata(input_dir):
    """Načte metadata pro všechny PDB soubory v adresáři."""
    pdb_files = [
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith(".pdb") and not f.startswith(".")
    ]
    pdb_files.sort()

    if not pdb_files:
        raise FileNotFoundError(f"V adresáři '{input_dir}' nebyly nalezeny žádné .pdb soubory.")

    print(f"[*] Načítám metadata a pLDDT pro {len(pdb_files)} PDB souborů v '{input_dir}'...")
    metadata_map = {}
    for p in pdb_files:
        meta = extract_pdb_metadata(p)
        # Klíčem je jak celé jméno souboru, tak kmen (např. bez přípony) pro párování s Foldseek výstupem
        metadata_map[meta["filename"]] = meta
        metadata_map[os.path.splitext(meta["filename"])[0]] = meta
        metadata_map[meta["accession"]] = meta

    return metadata_map, pdb_files


def run_foldseek_all_vs_all(input_dir, output_aln, tmp_dir, foldseek_bin="foldseek", threads=8):
    """
    Spustí Foldseek easy-search (all-vs-all) s TM-align metrikou (--alignment-type 1).
    """
    os.makedirs(tmp_dir, exist_ok=True)
    out_dir = os.path.dirname(output_aln)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    format_fields = "query,target,alntmscore,qtmscore,ttmscore,qcov,tcov,fident,lddt,alnlen"
    cmd = [
        foldseek_bin,
        "easy-search",
        input_dir,
        input_dir,
        output_aln,
        tmp_dir,
        "--alignment-type", "1",  # TM-align výpočet TM-score
        "--format-output", format_fields,
        "--threads", str(threads),
        "-v", "2",
    ]

    print(f"[*] Spouštím Foldseek all-vs-all TM-align vyhledávání:")
    print(f"    {' '.join(cmd)}")

    try:
        res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(f"[+] Foldseek vyhledávání úspěšně dokončeno. Výstup uložen do: {output_aln}")
    except FileNotFoundError:
        print(f"\n[!] CHYBA: Příkaz '{foldseek_bin}' nebyl nalezen v systému!", file=sys.stderr)
        print("    Ujistěte se, že je Foldseek nainstalován a dostupný v PATH (např. conda activate foldseek nebo module load foldseek).", file=sys.stderr)
        print("    Pokud již máte předpočítaný výsledek, zadejte jej parametrem --existing-aln <cesta.m8>\n", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"\n[!] CHYBA při běhu Foldseeku (návratový kód {e.returncode}):", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        sys.exit(1)


def parse_foldseek_alignments(aln_file):
    """
    Načte Foldseek m8 tabulku do slovníku hran mezi proteiny s TM-score a pokrytím.
    """
    if not os.path.exists(aln_file):
        raise FileNotFoundError(f"Soubor s Foldseek zarovnáním nebyl nalezen: {aln_file}")

    print(f"[*] Načítám zarovnání z '{aln_file}'...")
    edges = {}  # query -> list of (target, qtmscore, ttmscore, qcov, tcov, alntmscore, fident)
    total_lines = 0

    with open(aln_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                parts = line.split()
            if len(parts) < 8:
                continue

            query = parts[0]
            target = parts[1]
            try:
                alntmscore = float(parts[2])
                qtmscore = float(parts[3])
                ttmscore = float(parts[4])
                qcov = float(parts[5])
                tcov = float(parts[6])
                fident = float(parts[7])
            except ValueError:
                continue

            total_lines += 1
            if query not in edges:
                edges[query] = []
            edges[query].append({
                "target": target,
                "alntmscore": alntmscore,
                "qtmscore": qtmscore,
                "ttmscore": ttmscore,
                "qcov": qcov,
                "tcov": tcov,
                "fident": fident,
            })

    print(f"[*] Načteno {total_lines} párů strukturních zarovnání pro {len(edges)} dotazů.")
    return edges


def cluster_by_tmscore(pdb_items, edges, tm_threshold=0.75, cov_threshold=0.75, representative_metric="plddt"):
    """
    Provádí greedy clustering na základě TM-score a pokrytí.

    Pravidlo ekvivalence (redundance):
    Dva proteiny A a B jsou v témže strukturním klastru, pokud:
      min(qtmscore, ttmscore) >= tm_threshold  A  qcov >= cov_threshold A tcov >= cov_threshold

    Výběr reprezentanta:
    1. Všechny proteiny jsou seřazeny podle zvoleného kritéria (nejvyšší pLDDT -> nejdelší sekvence).
    2. Nejlepší dostupný protein se stane reprezentantem klastru C.
    3. Všechny dosud nezařazené proteiny, které splňují strukturní podobnost s tímto reprezentantem,
       jsou přiřazeny do klastru C jako redundantní členové.
    4. Postup se opakuje, dokud nejsou zařazeny všechny proteiny.
    """
    # 1. Seřadit proteiny pro prioritizaci výběru reprezentanta
    if representative_metric == "plddt":
        sorted_items = sorted(pdb_items, key=lambda x: (x["plddt"], x["length"]), reverse=True)
    elif representative_metric == "length":
        sorted_items = sorted(pdb_items, key=lambda x: (x["length"], x["plddt"]), reverse=True)
    else:
        sorted_items = list(pdb_items)

    assigned = set()
    clusters = []  # list of dict: { "representative": meta, "members": [meta, ...], "edges_info": {target: edge} }

    # Rychlá indexace hran podle jmen souborů a kmenů
    # Vytvoříme normalizovanou mapu sousedů
    item_by_id = {item["filename"]: item for item in pdb_items}
    for item in pdb_items:
        item_by_id[os.path.splitext(item["filename"])[0]] = item
        item_by_id[item["accession"]] = item

    for candidate in sorted_items:
        cand_key = candidate["filename"]
        if cand_key in assigned:
            continue

        # Založíme nový klastr s tímto reprezentantem
        cluster_rep = candidate
        cluster_members = [cluster_rep]
        assigned.add(cand_key)
        rep_edges_map = {}

        # Vyhledat sousedy reprezentanta ve Foldseek zarovnání
        # Zkusíme různé tvary identifikátoru reprezentanta (plný název, bez .pdb, accession)
        cand_queries = [cand_key, os.path.splitext(cand_key)[0], candidate["accession"]]
        cand_edges = []
        for qk in cand_queries:
            if qk in edges:
                cand_edges = edges[qk]
                break

        for edge in cand_edges:
            target_raw = edge["target"]
            target_meta = item_by_id.get(target_raw)
            if target_meta is None:
                # Zkusit bez přípony
                target_meta = item_by_id.get(os.path.splitext(target_raw)[0])

            if target_meta is None:
                continue

            target_key = target_meta["filename"]
            if target_key in assigned:
                continue

            # Kontrola podmínek TM-score a pokrytí
            # qtmscore = TM-score normalizované na dotaz (reprezentant)
            # ttmscore = TM-score normalizované na cíl (člen)
            min_tm = min(edge["qtmscore"], edge["ttmscore"])
            min_cov = min(edge["qcov"], edge["tcov"])

            if min_tm >= tm_threshold and min_cov >= cov_threshold:
                cluster_members.append(target_meta)
                assigned.add(target_key)
                rep_edges_map[target_key] = edge

        clusters.append({
            "cluster_id": len(clusters) + 1,
            "representative": cluster_rep,
            "members": cluster_members,
            "edges_info": rep_edges_map,
            "size": len(cluster_members),
        })

    # Zkontrolujeme, zda nezbyl nějaký nezařazený protein (např. bez Foldseek hitu)
    for item in pdb_items:
        if item["filename"] not in assigned:
            clusters.append({
                "cluster_id": len(clusters) + 1,
                "representative": item,
                "members": [item],
                "edges_info": {},
                "size": 1,
            })
            assigned.add(item["filename"])

    return clusters


def explore_thresholds(pdb_items, edges, thresholds=None, cov_threshold=0.75, reports_dir="reports_clustering"):
    """
    Otestuje sadu různých TM-score thresholdů a vytvoří přehlednou tabulku citlivosti.
    Pomáhá uživateli zvolit ten 'tak akorát přísný' práh.
    """
    if thresholds is None:
        thresholds = [0.50, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]

    print("\n" + "=" * 75)
    print("        PRŮZKUM CITLIVOSTI TM-SCORE THRESHOLDŮ (REDUKCE REDUNDANCE)")
    print("=" * 75)
    print(f" Celkem vstupních struktur: {len(pdb_items)}")
    print(f" Minimální pokrytí (coverage): {cov_threshold * 100:.0f} %")
    print("-" * 75)
    print(f" {'TM-score':<10} | {'Počet klastrů':<14} | {'Redukce (%)':<12} | {'Biologický význam / interpretace'}")
    print("-" * 75)

    exploration_data = []

    for tm in thresholds:
        cls = cluster_by_tmscore(pdb_items, edges, tm_threshold=tm, cov_threshold=cov_threshold)
        num_clusters = len(cls)
        red_pct = (1.0 - (num_clusters / len(pdb_items))) * 100

        interp = ""
        if tm <= 0.50:
            interp = "Příliš volný: slévá různé enzymy se stejným záhybem (fold)"
        elif tm <= 0.60:
            interp = "Volný: slévá enzymy v rámci nadrodiny (superfamily)"
        elif tm <= 0.70:
            interp = "Střední: spojuje příbuzné rodiny se stejným uspořádáním domén"
        elif tm <= 0.75:
            interp = "DOPORUČENÝ: ideální pro odstranění duplikátů a paralogů v T. vaginalis"
        elif tm <= 0.80:
            interp = "DOPORUČENÝ: přísnější, zachová i jemné odchylky v konformaci smyček"
        elif tm <= 0.85:
            interp = "Přísný: rozpadá klastry s flexibilními N/C konci"
        else:
            interp = "Velmi přísný: téměř totožná páteř, ponechává většinu duplikátů"

        print(f" {tm:<10.2f} | {num_clusters:<14d} | {red_pct:<11.1f}% | {interp}")
        exploration_data.append({
            "tm_threshold": tm,
            "coverage_threshold": cov_threshold,
            "num_clusters": num_clusters,
            "total_input": len(pdb_items),
            "reduction_percent": round(red_pct, 1),
            "interpretation": interp,
        })

    print("=" * 75 + "\n")

    # Uložit do TSV
    os.makedirs(reports_dir, exist_ok=True)
    exp_tsv = os.path.join(reports_dir, "threshold_exploration.tsv")
    with open(exp_tsv, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(exploration_data[0].keys()), delimiter="\t")
        writer.writeheader()
        for row in exploration_data:
            writer.writerow(row)
    print(f"[*] Výsledky průzkumu prahů uloženy do: {exp_tsv}")

    return exploration_data


def copy_representative_structures(clusters, output_dir, use_symlink=False):
    """
    Zkopíruje (nebo nalinkuje) reprezentativní PDB soubory do nové výstupní složky.
    """
    os.makedirs(output_dir, exist_ok=True)
    copied_count = 0

    for c in clusters:
        rep = c["representative"]
        src = rep["pdb_path"]
        dst = os.path.join(output_dir, rep["filename"])

        if use_symlink:
            if os.path.exists(dst) or os.path.islink(dst):
                os.remove(dst)
            os.symlink(os.path.abspath(src), dst)
        else:
            shutil.copy2(src, dst)
        copied_count += 1

    action = "Nalinkováno" if use_symlink else "Zkopírováno"
    print(f"[+] {action} {copied_count} reprezentativních struktur do: '{output_dir}'")
    return copied_count


def generate_clustering_reports(clusters, pdb_items, tm_threshold, cov_threshold, output_dir, reports_dir, exploration_data=None):
    """
    Vygeneruje kompletní TSV a Markdown reporty o klastrování.
    """
    os.makedirs(reports_dir, exist_ok=True)

    # 1. Tabulka reprezentantů (cluster_representatives.tsv)
    rep_tsv = os.path.join(reports_dir, "cluster_representatives.tsv")
    with open(rep_tsv, "w", encoding="utf-8", newline="") as f:
        fields = ["Cluster_ID", "Representative_PDB", "UniProt_Entry", "pLDDT", "Length", "Cluster_Size", "Member_Entries"]
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for c in clusters:
            rep = c["representative"]
            members_accs = [m["accession"] for m in c["members"]]
            writer.writerow({
                "Cluster_ID": c["cluster_id"],
                "Representative_PDB": rep["filename"],
                "UniProt_Entry": rep["accession"],
                "pLDDT": rep["plddt"],
                "Length": rep["length"],
                "Cluster_Size": c["size"],
                "Member_Entries": ", ".join(members_accs),
            })

    # 2. Kompletní tabulka všech členů (cluster_all_members.tsv)
    members_tsv = os.path.join(reports_dir, "cluster_all_members.tsv")
    with open(members_tsv, "w", encoding="utf-8", newline="") as f:
        fields = ["UniProt_Entry", "PDB_File", "pLDDT", "Length", "Cluster_ID", "Role", "Representative_PDB", "TM_score_to_Rep"]
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for c in clusters:
            rep = c["representative"]
            for m in c["members"]:
                is_rep = (m["filename"] == rep["filename"])
                edge = c["edges_info"].get(m["filename"], {})
                tm_to_rep = 1.0 if is_rep else edge.get("alntmscore", edge.get("ttmscore", "N/A"))
                writer.writerow({
                    "UniProt_Entry": m["accession"],
                    "PDB_File": m["filename"],
                    "pLDDT": m["plddt"],
                    "Length": m["length"],
                    "Cluster_ID": c["cluster_id"],
                    "Role": "Representative" if is_rep else "Redundant_Member",
                    "Representative_PDB": rep["filename"],
                    "TM_score_to_Rep": tm_to_rep,
                })

    # 3. Markdown zpráva (clustering_report.md)
    md_path = os.path.join(reports_dir, "clustering_report.md")
    now_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    total_in = len(pdb_items)
    total_reps = len(clusters)
    total_redundant = total_in - total_reps
    red_pct = (total_redundant / total_in * 100) if total_in else 0

    multi_member_clusters = [c for c in clusters if c["size"] > 1]
    singletons = [c for c in clusters if c["size"] == 1]

    with open(md_path, "w", encoding="utf-8") as md:
        md.write("# Foldseek TM-score Strukturní Klastrování Methyltransferáz\n\n")
        md.write(f"- **Datum běhu:** `{now_str}`\n")
        md.write(f"- **Zvolený TM-score threshold:** **`{tm_threshold:.2f}`**\n")
        md.write(f"- **Minimální pokrytí (coverage):** **`{cov_threshold * 100:.0f} %`**\n")
        md.write(f"- **Složka nere redundantních reprezentantů:** [`{output_dir}`](file://{os.path.abspath(output_dir)})\n\n")

        md.write("## 1. Souhrn redukce redundance\n\n")
        md.write("| Metrika | Hodnota |\n")
        md.write("| :--- | :---: |\n")
        md.write(f"| **Celkem vstupních struktur** | **{total_in}** |\n")
        md.write(f"| **Zachovaných reprezentantů (nová složka)** | **{total_reps}** |\n")
        md.write(f"| **Odstraněných redundantních struktur** | **{total_redundant}** |\n")
        md.write(f"| **Celková redukce redundance** | **{red_pct:.1f} %** |\n")
        md.write(f"| Počet klastrů s více členy (paralogy/duplikáty) | {len(multi_member_clusters)} |\n")
        md.write(f"| Počet unikátních struktur (singletons) | {len(singletons)} |\n\n")

        md.write("## 2. Doporučení a zdůvodnění TM-score prahu\n\n")
        md.write("> [!TIP]\n")
        md.write(f"> **Proč je TM-score = {tm_threshold:.2f} 'tak akorát přísný'?**\n")
        md.write(">\n")
        md.write("> - **TM-score < 0.50:** Značí náhodnou shodu nebo pouze podobnost sekundárních struktur. Při takto nízkém prahu by došlo k chybnému sloučení zcela odlišných enzymů, které pouze sdílejí obecný Rossmannův záhyb.\n")
        md.write("> - **TM-score 0.60 – 0.70:** Sdružuje proteiny na úrovni enzymové nadrodiny. Může sloučit methyltransferázy s různou substrátovou specifitou (např. DNA vs tRNA methyltransferázy).\n")
        md.write(f"> - **TM-score 0.75 – 0.80 (Zvoleno: {tm_threshold:.2f}):** Zlatý střed. *Trichomonas vaginalis* prošel masivní expanzí genomu a obsahuje stovky duplikovaných genů. Tento práh spolehlivě zredukuje identické a téměř identické paralogy a izoformy, ale **ponechá oddělené různé funkční rodiny** (Rossmann Class I, SET domény, SPOUT domény, membránové ICMT).\n")
        md.write("> - **TM-score > 0.85:** Příliš přísné. Kvůli flexibilitě nestrukturovaných N/C konců u AlphaFold modelů by tento práh selhal v odstranění zjevně redundantních duplikátů.\n\n")

        if exploration_data:
            md.write("## 3. Průzkum citlivosti TM-score thresholdů\n\n")
            md.write("| TM-score práh | Počet klastrů | Redukce (%) | Interpretace |\n")
            md.write("| :---: | :---: | :---: | :--- |\n")
            for exp in exploration_data:
                mark = "**" if abs(exp["tm_threshold"] - tm_threshold) < 1e-4 else ""
                md.write(f"| {mark}{exp['tm_threshold']:.2f}{mark} | {mark}{exp['num_clusters']}{mark} | {mark}{exp['reduction_percent']:.1f} %{mark} | {exp['interpretation']} |\n")
            md.write("\n")

        md.write("## 4. Přehled klastrů s více členy (odhalené paralogy a rodiny)\n\n")
        md.write("| Klastr | Reprezentant (PDB) | pLDDT | Velikost | Členové klastru |\n")
        md.write("| :---: | :--- | :---: | :---: | :--- |\n")
        for c in sorted(multi_member_clusters, key=lambda x: x["size"], reverse=True):
            rep = c["representative"]
            mem_list = ", ".join([f"`{m['accession']}`" for m in c["members"]])
            md.write(f"| **Klastr {c['cluster_id']}** | [`{rep['filename']}`](file://{os.path.abspath(os.path.join(output_dir, rep['filename']))}) | {rep['plddt']:.1f} | **{c['size']}** | {mem_list} |\n")
        md.write("\n")

        md.write("## 5. Vygenerované soubory\n\n")
        md.write(f"- **Reprezentativní PDB soubory:** [`{output_dir}/`](file://{os.path.abspath(output_dir)})\n")
        md.write(f"- **Tabulka reprezentantů:** [`{rep_tsv}`](file://{os.path.abspath(rep_tsv)})\n")
        md.write(f"- **Detailní mapování všech členů:** [`{members_tsv}`](file://{os.path.abspath(members_tsv)})\n")
        if exploration_data:
            md.write(f"- **Průzkum prahů (TSV):** [`{os.path.join(reports_dir, 'threshold_exploration.tsv')}`](file://{os.path.abspath(os.path.join(reports_dir, 'threshold_exploration.tsv'))})\n")

    print(f"[+] Reporty byly uloženy do složky '{reports_dir}':")
    print(f"    - {md_path}")
    print(f"    - {rep_tsv}")
    print(f"    - {members_tsv}")


def main():
    parser = argparse.ArgumentParser(
        description="Odstranění redundantních methyltransferáz pomocí Foldseeku a metriky TM-score."
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        default="structures",
        help="Vstupní složka obsahující stažené PDB struktury (výchozí: structures).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="structures_representative",
        help="Nová cílová složka pro reprezentativní PDB struktury (výchozí: structures_representative).",
    )
    parser.add_argument(
        "-r",
        "--reports-dir",
        default="reports_clustering",
        help="Složka pro ukládání reportů o klastrování (výchozí: reports_clustering).",
    )
    parser.add_argument(
        "--tmscore-threshold",
        type=float,
        default=0.75,
        help="Práh TM-score pro redukci redundance (výchozí: 0.75 - 'tak akorát přísný').",
    )
    parser.add_argument(
        "--coverage-threshold",
        type=float,
        default=0.75,
        help="Minimální obousměrné pokrytí délky struktury (výchozí: 0.75).",
    )
    parser.add_argument(
        "--foldseek-bin",
        default="foldseek",
        help="Cesta k binárce Foldseek (výchozí: 'foldseek' z PATH).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Počet vláken pro Foldseek (výchozí: 8).",
    )
    parser.add_argument(
        "--tmp-dir",
        default="tmp_foldseek",
        help="Dočasná složka pro databáze Foldseeku (výchozí: tmp_foldseek).",
    )
    parser.add_argument(
        "--existing-aln",
        default=None,
        help="Cesta k již existujícímu Foldseek m8 výstupu (přeskočí běh Foldseeku a rovnou klastruje).",
    )
    parser.add_argument(
        "--explore-thresholds",
        action="store_true",
        help="Provede průzkum citlivosti pro sadu TM-score prahů (0.50 až 0.90) a vygeneruje srovnávací tabulku.",
    )
    parser.add_argument(
        "--symlink",
        action="store_true",
        help="Vytvořit ve výstupní složce pouze symbolické linky místo kopírování souborů.",
    )
    parser.add_argument(
        "--representative-metric",
        choices=["plddt", "length"],
        default="plddt",
        help="Kritérium výběru reprezentanta: 'plddt' (nejvyšší AlphaFold pLDDT, výchozí) nebo 'length' (nejdelší sekvence).",
    )

    args = parser.parse_args()

    print("=" * 75)
    print("      Foldseek TM-score Redukce Redundance Methyltransferáz")
    print("=" * 75)
    print(f"[*] Vstupní složka      : {args.input_dir}")
    print(f"[*] Cílová nová složka  : {args.output_dir}")
    print(f"[*] Složka reportů      : {args.reports_dir}")
    print(f"[*] TM-score threshold  : {args.tmscore_threshold}")
    print(f"[*] Coverage threshold  : {args.coverage_threshold * 100:.0f} %")
    print(f"[*] Výběr reprezentanta : {args.representative_metric.upper()} (nejvyšší kvalita v klastru)")
    print("-" * 75)

    # 1. Načtení PDB metadat
    metadata_map, pdb_files = load_all_pdb_metadata(args.input_dir)
    pdb_items = [metadata_map[os.path.basename(p)] for p in pdb_files]
    print(f"[*] Úspěšně analyzováno {len(pdb_items)} PDB struktur.")

    # 2. Foldseek zarovnání
    aln_file = args.existing_aln
    if not aln_file:
        aln_file = os.path.join(args.reports_dir, "foldseek_all_vs_all.m8")
        run_foldseek_all_vs_all(
            args.input_dir,
            aln_file,
            args.tmp_dir,
            foldseek_bin=args.foldseek_bin,
            threads=args.threads,
        )

    # 3. Načtení hran zarovnání
    edges = parse_foldseek_alignments(aln_file)

    # 4. Průzkum prahů, pokud je vyžádán
    exploration_data = None
    if args.explore_thresholds:
        exploration_data = explore_thresholds(
            pdb_items,
            edges,
            cov_threshold=args.coverage_threshold,
            reports_dir=args.reports_dir,
        )

    # 5. Klastrování podle zvoleného TM-score
    print(f"\n[*] Provádím greedy klastrování s TM-score >= {args.tmscore_threshold} a coverage >= {args.coverage_threshold}...")
    clusters = cluster_by_tmscore(
        pdb_items,
        edges,
        tm_threshold=args.tmscore_threshold,
        cov_threshold=args.coverage_threshold,
        representative_metric=args.representative_metric,
    )

    print(f"[+] Nalezeno {len(clusters)} nere redundantních klastrů z původních {len(pdb_items)} struktur.")
    reduction = len(pdb_items) - len(clusters)
    print(f"[+] Odstraněno {reduction} redundantních struktur (redukce o {(reduction / len(pdb_items) * 100):.1f} %).")

    # 6. Vytvoření nové složky a uložení reprezentantů
    print(f"\n[*] Ukládám reprezentativní PDB struktury do nové složky '{args.output_dir}'...")
    copy_representative_structures(clusters, args.output_dir, use_symlink=args.symlink)

    # 7. Generování reportů
    print(f"\n[*] Generuji podrobné reporty do '{args.reports_dir}'...")
    generate_clustering_reports(
        clusters,
        pdb_items,
        args.tmscore_threshold,
        args.coverage_threshold,
        args.output_dir,
        args.reports_dir,
        exploration_data=exploration_data,
    )

    print("\n" + "=" * 75)
    print("                     DOKONČENO!")
    print("=" * 75)
    print(f"  Původní počet struktur       : {len(pdb_items)}")
    print(f"  Výsledný počet reprezentantů : {len(clusters)}")
    print(f"  Nová složka s reprezentanty  : {os.path.abspath(args.output_dir)}")
    print(f"  Souhrnný Markdown report     : {os.path.abspath(os.path.join(args.reports_dir, 'clustering_report.md'))}")
    print("=" * 75)


if __name__ == "__main__":
    main()
