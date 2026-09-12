#!/usr/bin/env python3
"""
evaluate_docking_results.py

Skript pro komplexní vyhodnocení výsledků virtuálního screeningu EasyDock
napříč všemi proteinovými cíli.

Hlavní funkce:
1. Spojí výsledky všech proteinů (z SQLite databází .db nebo .sdf souborů) do jednoho master CSV.
2. Vypočítá metriky očištěné od vlivu molekulové hmotnosti (tzv. size-independent scores):
   - Ligand Efficiency (LE = -score / HAC)
   - Size-Independent Ligand Efficiency (SILE = -score / HAC^0.3)
   - Binding Efficiency Index (BEI = -score / (MW / 1000))
   - MW-Residual Score (reziduum z lineární regrese skóre vs. MW: zápornější = váže se lépe,
     než odpovídá pouhé velikosti molekuly).
3. Analyzuje proteiny a seřadí je podle toho:
   - Na které proteiny se ligandy vážou nejsilněji (průměr, medián, top skóre).
   - Které proteiny vážou největší počet potentních inhibitorů (skóre <= -8, -9, -10 kcal/mol).
4. Vytvoří matici skóre (ligandy x proteiny) pro rychlou analýzu selektivity / cross-reaktivity.
5. Generuje přehledný Markdown report se shrnutím výsledků.

Použití:
    python3 evaluate_docking_results.py
    python3 evaluate_docking_results.py --results-dir results_docking --output-dir reports_docking
"""

import argparse
import os
import re
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# RDKit pro výpočet deskriptorů
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


def calculate_mol_descriptors(smiles_dict: dict[str, str]) -> dict[str, dict]:
    """
    Vypočítá molekulovou hmotnost (MW), počet těžkých atomů (HAC) a LogP pro unikátní SMILES.
    """
    descriptors = {}
    if not RDKIT_AVAILABLE:
        print("[!] Varování: RDKit není k dispozici. Deskriptory MW a HAC budou aproximovány.", file=sys.stderr)

    for cid, smi in smiles_dict.items():
        mw = np.nan
        hac = np.nan
        logp = np.nan

        if RDKIT_AVAILABLE and smi and isinstance(smi, str):
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is not None:
                    mw = round(Descriptors.MolWt(mol), 2)
                    hac = int(mol.GetNumHeavyAtoms())
                    logp = round(Crippen.MolLogP(mol), 2)
            except Exception:
                pass

        # Záložní odhad, pokud RDKit selže nebo není k dispozici
        if np.isnan(hac) and smi and isinstance(smi, str):
            # Odhad počtu těžkých atomů ze SMILES (C, N, O, S, P, halogeny)
            hac = len(re.findall(r"[CNOFSPIBrcnlops]", smi))
            mw = hac * 13.0  # Hrubý odhad

        descriptors[cid] = {
            "MW": mw if not np.isnan(mw) else 300.0,
            "HAC": hac if not np.isnan(hac) and hac > 0 else 22,
            "LogP": logp
        }

    return descriptors


def load_from_sqlite_dbs(results_dir: Path) -> pd.DataFrame:
    """Načte dokovací výsledky ze všech .db souborů v zadané složce."""
    db_files = sorted(list(results_dir.glob("*.db")))
    if not db_files:
        return pd.DataFrame()

    print(f"[*] Nalezeno {len(db_files)} SQLite databází v '{results_dir}'.")
    records = []

    for db_path in db_files:
        target_id = db_path.stem
        try:
            conn = sqlite3.connect(str(db_path))
            cursor = conn.cursor()

            # Zjistíme dostupné sloupce v tabulce mols
            cursor.execute("PRAGMA table_info(mols)")
            cols = [info[1] for info in cursor.fetchall()]

            if "docking_score" not in cols:
                conn.close()
                continue

            smi_col = "smi" if "smi" in cols else ("smi_input" if "smi_input" in cols else "smi_protonated")

            query = f"""
                SELECT id, stereo_id, docking_score, {smi_col}
                FROM mols
                WHERE docking_score IS NOT NULL
            """
            cursor.execute(query)
            rows = cursor.fetchall()

            # Volitelná kontrola tabulky bust (PoseBusters)
            bust_dict = {}
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='bust'")
            if cursor.fetchone():
                try:
                    cursor.execute("SELECT id, stereo_id, result FROM bust WHERE pose=1")
                    for bid, bstereo, bres in cursor.fetchall():
                        bust_dict[(bid, bstereo)] = bool(bres)
                except Exception:
                    pass

            conn.close()

            for row in rows:
                cid, stereo_id, score, smi = row
                passed_bust = bust_dict.get((cid, stereo_id), np.nan)
                records.append({
                    "target_id": target_id,
                    "compound_id": cid,
                    "stereo_id": stereo_id,
                    "docking_score": float(score),
                    "smiles": smi,
                    "posebusters_pass": passed_bust
                })
        except Exception as e:
            print(f"[!] Chyba při čtení databáze '{db_path.name}': {e}", file=sys.stderr)

    return pd.DataFrame(records)


def load_from_sdf_files(results_dir: Path) -> pd.DataFrame:
    """Záložní načtení dokovacích výsledků ze souborů .sdf, pokud nejsou přítomny .db."""
    if not RDKIT_AVAILABLE:
        return pd.DataFrame()

    sdf_files = sorted(list(results_dir.glob("*.sdf")))
    if not sdf_files:
        return pd.DataFrame()

    print(f"[*] Načítám záložní data z {len(sdf_files)} .sdf souborů...")
    records = []

    for sdf_path in sdf_files:
        target_id = sdf_path.stem
        try:
            suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
            for mol in suppl:
                if mol is None:
                    continue
                cid = mol.GetProp("_Name") if mol.HasProp("_Name") else "unknown"
                score = np.nan
                for prop in ["docking_score", "affinity", "VINA_SCORE", "r_i_docking_score"]:
                    if mol.HasProp(prop):
                        try:
                            score = float(mol.GetProp(prop))
                            break
                        except ValueError:
                            pass

                if not np.isnan(score):
                    smi = Chem.MolToSmiles(mol)
                    records.append({
                        "target_id": target_id,
                        "compound_id": cid,
                        "stereo_id": 0,
                        "docking_score": score,
                        "smiles": smi,
                        "posebusters_pass": np.nan
                    })
        except Exception as e:
            print(f"[!] Chyba při čtení SDF '{sdf_path.name}': {e}", file=sys.stderr)

    return pd.DataFrame(records)


def compute_decorrelated_scores(df: pd.DataFrame) -> pd.DataFrame:
    """
    Vypočítá Ligand Efficiency (LE), SILE, BEI a lineární reziduální skóre očištěné od hmotnosti.
    """
    # 1. Extrakce unikátních molekul a výpočet deskriptorů
    unique_compounds = df.drop_duplicates(subset=["compound_id"])[["compound_id", "smiles"]].set_index("compound_id")["smiles"].to_dict()
    descriptors = calculate_mol_descriptors(unique_compounds)

    desc_df = pd.DataFrame.from_dict(descriptors, orient="index").reset_index()
    desc_df.rename(columns={"index": "compound_id"}, inplace=True)

    df = df.merge(desc_df, on="compound_id", how="left")

    # 2. Základní ligandová efektivita (kcal/mol na těžký atom - vyšší je lepší)
    # Vina skóre je záporné, proto LE = -score / HAC (kladné číslo, typicky 0.25 - 0.45)
    df["LE"] = round(-df["docking_score"] / df["HAC"], 3)

    # 3. Size-Independent Ligand Efficiency (SILE = -score / (HAC^0.3))
    df["SILE"] = round(-df["docking_score"] / (df["HAC"] ** 0.3), 3)

    # 4. Binding Efficiency Index (BEI = -score / (MW / 1000))
    df["BEI"] = round(-df["docking_score"] / (df["MW"] / 1000.0), 3)

    # 5. Skóre nekorelované s hmotností (MW-Residual Score):
    # Fitujeme lineární regresi: docking_score = a * MW + b
    # Očekávané skóre je čistě funkce hmotnosti.
    # Reziduum = docking_score - expected_score
    # Záporné reziduum znamená, že molekula se váže VÝRAZNĚ LÉPE, než by odpovídalo její velikosti!
    valid_mask = ~df["MW"].isna() & ~df["docking_score"].isna()
    if valid_mask.sum() > 5:
        mws = df.loc[valid_mask, "MW"].values
        scores = df.loc[valid_mask, "docking_score"].values

        # Lineární regrese OLS: y = a*x + b
        slope, intercept = np.polyfit(mws, scores, 1)
        r_corr = np.corrcoef(mws, scores)[0, 1]

        expected_scores = slope * df["MW"].values + intercept
        residuals = df["docking_score"].values - expected_scores

        df["MW_expected_score"] = np.round(expected_scores, 2)
        df["MW_residual_score"] = np.round(residuals, 2)

        # Standardizované Z-skóre reziduí
        std_res = np.std(residuals)
        if std_res > 1e-6:
            df["MW_residual_zscore"] = np.round((residuals - np.mean(residuals)) / std_res, 2)
        else:
            df["MW_residual_zscore"] = 0.0

        print(f"[*] Korelace mezi MW a Docking Score: r = {r_corr:.3f}")
        print(f"[*] Lineární model: Očekávané skóre = {slope:.4f} * MW + ({intercept:.2f})")
    else:
        df["MW_expected_score"] = np.nan
        df["MW_residual_score"] = np.nan
        df["MW_residual_zscore"] = np.nan

    return df


def generate_protein_summary(best_df: pd.DataFrame, potent_thresh: float = -9.0, mod_thresh: float = -8.0) -> pd.DataFrame:
    """
    Agreguje data podle proteinových cílů a hodnotí, na který protein se váže
    největší počet i nejsilnější inhibitory.
    """
    grouped = best_df.groupby("target_id")

    summary = pd.DataFrame({
        "target_id": grouped["target_id"].first(),
        "total_ligands": grouped["compound_id"].nunique(),
        "best_score": grouped["docking_score"].min(),
        "mean_score": grouped["docking_score"].mean().round(2),
        "median_score": grouped["docking_score"].median().round(2),
        "std_score": grouped["docking_score"].std().round(2),
        "top5_mean_score": grouped["docking_score"].apply(lambda s: s.nsmallest(5).mean() if len(s) >= 5 else s.min()).round(2),
        f"potent_binders (score<={potent_thresh})": grouped["docking_score"].apply(lambda s: (s <= potent_thresh).sum()),
        f"moderate_binders (score<={mod_thresh})": grouped["docking_score"].apply(lambda s: (s <= mod_thresh).sum()),
        "mean_LE": grouped["LE"].mean().round(3),
        "mean_SILE": grouped["SILE"].mean().round(3),
        "mean_MW_residual": grouped["MW_residual_score"].mean().round(2),
    }).reset_index(drop=True)

    # Procentuální podíl silných vazačů
    summary[f"potent_binders_%"] = round(100.0 * summary[f"potent_binders (score<={potent_thresh})"] / summary["total_ligands"], 1)

    # Celkové kompozitní skóre cíle (Target Druggability / Hit Enrichment Index)
    # Čím zápornější skóre a čím více silných vazačů, tím vyšší priorita.
    # Normalizujeme ranky:
    score_rank = summary["mean_score"].rank(ascending=True)  # Nejnižší skóre = rank 1
    count_rank = summary[f"potent_binders (score<={potent_thresh})"].rank(ascending=False) # Nejvíce vazačů = rank 1
    top5_rank = summary["top5_mean_score"].rank(ascending=True)
    res_rank = summary["mean_MW_residual"].rank(ascending=True)

    summary["overall_priority_rank"] = ((score_rank * 0.3) + (count_rank * 0.4) + (top5_rank * 0.2) + (res_rank * 0.1)).rank(ascending=True).astype(int)

    summary.sort_values(by="overall_priority_rank", inplace=True)
    return summary


def df_to_markdown(df: pd.DataFrame) -> str:
    """Převede DataFrame na Markdown tabulku bez nutnosti externí knihovny tabulate."""
    if df.empty:
        return ""
    headers = [str(c) for c in df.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |"
    ]
    for _, row in df.iterrows():
        row_str = [str(v) if pd.notna(v) else "" for v in row]
        lines.append("| " + " | ".join(row_str) + " |")
    return "\n".join(lines)


def write_markdown_report(summary_df: pd.DataFrame, top_hits_df: pd.DataFrame, report_path: Path, potent_thresh: float):
    """Vytvoří přehledný souhrnný report v Markdownu."""
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("# Souhrnný report virtuálního screeningu EasyDock\n\n")
        f.write(f"Vyhodnoceno celkem **{len(summary_df)}** proteinových cílů.\n\n")

        f.write("## 1. Nejlepší proteinové cíle (Seřazeno dle počtu a síly vazačů)\n\n")
        f.write("Kritéria: Kombinace počtu vazačů s vysokou afinitou, průměrného skóre a ligandové efektivity.\n\n")

        top_targets = summary_df.head(15).copy()
        cols_to_show = [
            "overall_priority_rank", "target_id", "total_ligands",
            f"potent_binders (score<={potent_thresh})", f"potent_binders_%",
            "best_score", "mean_score", "top5_mean_score", "mean_LE", "mean_MW_residual"
        ]
        cols_to_show = [c for c in cols_to_show if c in top_targets.columns]

        # Markdown tabulka
        f.write(df_to_markdown(top_targets[cols_to_show]))
        f.write("\n\n---\n\n")

        f.write("## 2. Top 20 nejlepších vazebných interakcí (Hit molekuly)\n\n")
        top_interactions = top_hits_df.head(20).copy()
        top_cols = ["target_id", "compound_id", "docking_score", "LE", "SILE", "MW_residual_score", "MW", "LogP"]
        top_cols = [c for c in top_cols if c in top_interactions.columns]
        f.write(df_to_markdown(top_interactions[top_cols]))
        f.write("\n\n---\n\n")

        f.write("## 3. Metodika a vysvětlení metrik\n\n")
        f.write("- **Docking Score (kcal/mol)**: Standardní AutoDock Vina afinitní skóre (nižší/zápornější = silnější vazba).\n")
        f.write("- **Ligand Efficiency (LE)**: $-\\frac{\\Delta G}{HAC}$ (skóre dělené počtem těžkých atomů). Hodnoty $> 0.30$ indikují vysoce efektivní vazače bez zbytečného balastu.\n")
        f.write("- **MW_residual_score (Skóre bez závislosti na hmotnosti)**: Reziduum z regrese mezi skóre a molekulovou hmotností. Záporná hodnota značí, že molekula má vyšší afinitu, než by odpovídalo pouze její molekulové velikosti.\n")
        f.write("- **Overall Priority Rank**: Celkové pořadí cílů kombinující počet potentních inhibitorů (40 %), průměrnou afinitu (30 %), top 5 průměr (20 %) a specifické reziduální skóre (10 %).\n")

    print(f"[+] Markdown report uložen do: '{report_path}'")


def main():
    parser = argparse.ArgumentParser(
        description="Vyhodnocení dokovacích výsledků EasyDock napříč všemi cíli."
    )
    parser.add_argument(
        "-r", "--results-dir",
        default="results_docking",
        help="Složka s výslednými .db a .sdf soubory (default: results_docking)"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="reports_docking",
        help="Výstupní složka pro CSV a reporty (default: reports_docking)"
    )
    parser.add_argument(
        "--potent-threshold",
        type=float,
        default=-9.0,
        help="Práh pro silné vazače v kcal/mol (default: -9.0)"
    )
    parser.add_argument(
        "--moderate-threshold",
        type=float,
        default=-8.0,
        help="Práh pro střední vazače v kcal/mol (default: -8.0)"
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)

    if not results_dir.exists():
        print(f"[!] CHYBA: Složka s výsledky '{results_dir}' neexistuje!", file=sys.stderr)
        print(f"    Ujistěte se, že dokování doběhlo a výsledky jsou v '{results_dir}'.", file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Vyhodnocení výsledků EasyDock screeningu")
    print(f"  Vstupní složka:     {results_dir}")
    print(f"  Výstupní složka:    {output_dir}")
    print(f"  Práh silných hitů:  <= {args.potent_threshold} kcal/mol")
    print("=" * 70)

    # 1. Načtení dat (nejprve SQLite DB, popř. záložní SDF)
    df = load_from_sqlite_dbs(results_dir)
    if df.empty:
        df = load_from_sdf_files(results_dir)

    if df.empty:
        print(f"[!] V '{results_dir}' nebyly nalezeny žádné výsledky (.db ani .sdf)!", file=sys.stderr)
        sys.exit(1)

    print(f"[+] Načteno celkem {len(df)} záznamů napříč {df['target_id'].nunique()} cíli a {df['compound_id'].nunique()} ligandy.")

    # 2. Výpočet metrik a očištění od hmotnosti
    print("[*] Počítám ligandovou efektivitu (LE, SILE, BEI) a skóre nekorelované s hmotností...")
    df = compute_decorrelated_scores(df)

    # Označení nejlepšího stereoisomeru pro každý pár (protein, sloučenina)
    df["is_best_stereoisomer"] = df.groupby(["target_id", "compound_id"])["docking_score"].transform(lambda s: s == s.min())

    # 3. Uložení Master CSV všech záznamů
    all_csv_path = output_dir / "docking_results_all.csv"
    df.to_csv(all_csv_path, index=False)
    print(f"[+] Všechny výsledky uloženy do: '{all_csv_path}'")

    # 4. Tabulka nejlepších stereoisomerů (hlavní pracovní soubor)
    best_df = df[df["is_best_stereoisomer"]].copy().drop_duplicates(subset=["target_id", "compound_id"])
    best_csv_path = output_dir / "docking_results_best_stereo.csv"
    best_df.to_csv(best_csv_path, index=False)
    print(f"[+] Nejlepší stereoisomery uloženy do: '{best_csv_path}'")

    # 5. Kontingenční matice skóre (Ligandy x Proteiny)
    matrix_scores = best_df.pivot_table(index="compound_id", columns="target_id", values="docking_score")
    matrix_scores.to_csv(output_dir / "docking_scores_matrix.csv")
    print(f"[+] Matice skóre (Ligandy x Proteiny) uložena do: '{output_dir / 'docking_scores_matrix.csv'}'")

    # Matice reziduálních skóre (bez vlivu hmotnosti)
    if "MW_residual_score" in best_df.columns and not best_df["MW_residual_score"].isna().all():
        matrix_res = best_df.pivot_table(index="compound_id", columns="target_id", values="MW_residual_score")
        matrix_res.to_csv(output_dir / "docking_mw_residuals_matrix.csv")
        print(f"[+] Matice reziduálních skóre uložena do: '{output_dir / 'docking_mw_residuals_matrix.csv'}'")

    # 6. Agregace cílů a ranking (hledáme proteiny s nejvíce vazači a nejsilnější afinitou)
    print("[*] Vytvářím ranking proteinových cílů...")
    protein_summary = generate_protein_summary(best_df, potent_thresh=args.potent_threshold, mod_thresh=args.moderate_threshold)
    summary_csv_path = output_dir / "protein_ranking_summary.csv"
    protein_summary.to_csv(summary_csv_path, index=False)
    print(f"[+] Souhrnný ranking proteinů uložen do: '{summary_csv_path}'")

    # 7. Top hity napříč celým screeningem
    top_hits = best_df.sort_values(by="docking_score", ascending=True).copy()
    top_hits_csv_path = output_dir / "top_overall_hits.csv"
    top_hits.to_csv(top_hits_csv_path, index=False)

    # 8. Markdown report
    md_report_path = output_dir / "docking_summary_report.md"
    write_markdown_report(protein_summary, top_hits, md_report_path, args.potent_threshold)

    print("\n" + "=" * 70)
    print("TOP 5 CÍLŮ S NEJVYŠŠÍ AFINITOU A POČTEM INHIBITORŮ:")
    print("=" * 70)
    cols_display = [
        "overall_priority_rank", "target_id",
        f"potent_binders (score<={args.potent_threshold})",
        "best_score", "mean_score", "mean_LE", "mean_MW_residual"
    ]
    cols_display = [c for c in cols_display if c in protein_summary.columns]
    print(protein_summary.head(5)[cols_display].to_string(index=False))
    print("=" * 70)
    print(f"[+] Kompletní analýza úspěšně dokončena! Všechny reporty naleznete v '{output_dir}/'.")


if __name__ == "__main__":
    main()
