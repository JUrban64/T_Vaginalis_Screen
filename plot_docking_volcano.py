#!/usr/bin/env python3
# ==============================================================================
# Skript pro vizualizaci dokovacích dat: Target Selectivity Volcano Plot
# (Varianta A: Selektivita konkrétního proteinu vůči zbytku panelu)
# ==============================================================================
# Pro vybraný protein (např. cíl z rodiny methyltransferáz T. vaginalis)
# srovná dokovací afinitu každého ligandu s jeho afinitou ke zbytku panelu (48 cílů).
#
# Osa X: Rozdíl afinity (Effect Size / Selektivita)
#        ΔScore = Score(cílový protein) - Medián(všechny ostatní proteiny)
#        [Záporné hodnoty = ligand se váže na cílový protein silněji než na zbytek]
#
# Osa Y: Statistická významnost
#        -log10(p-hodnota) z robustního Z-testu / t-testu vůči distribuci na panelu
#
# Výstup:
#   1. Publikační Volcano plot (PNG 300 DPI + PDF) s prahovými čarami a popisky hitů.
#   2. CSV tabulka selektivních hitů (seřazená dle selektivity a p-hodnoty).
#   3. Možnost multi-panel gridu pro TOP N proteinů v jednom přehledovém grafu.
#
# Autor: Antigravity AI / Jáchym Urban
# ==============================================================================

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Grafické a statistické knihovny
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats
except ImportError as e:
    print(f"[!] Chyba: Chybí potřebná knihovna: {e}", file=sys.stderr)
    print("    Nainstalujte: pip install matplotlib seaborn scipy pandas numpy", file=sys.stderr)
    sys.exit(1)

# Volitelná podpora adjustText pro zamezení překryvu popisků
try:
    from adjustText import adjust_text
    ADJUST_TEXT_AVAILABLE = True
except ImportError:
    ADJUST_TEXT_AVAILABLE = False


def load_data(input_path: Path) -> tuple[pd.DataFrame, pd.DataFrame | None, Path]:
    """
    Načte matici skóre a případný souhrnný ranking proteinů.
    """
    if input_path.is_dir():
        score_file = input_path / "docking_scores_matrix.csv"
        if not score_file.exists():
            candidates = list(input_path.glob("*scores_matrix*.csv"))
            if not candidates:
                raise FileNotFoundError(f"V adresáři '{input_path}' nebyl nalezen 'docking_scores_matrix.csv'!")
            score_file = candidates[0]
        matrix_df = pd.read_csv(score_file, index_col=0)

        ranking_file = input_path / "protein_ranking_summary.csv"
        ranking_df = pd.read_csv(ranking_file) if ranking_file.exists() else None
        output_dir = input_path
    else:
        if not input_path.exists():
            raise FileNotFoundError(f"Soubor '{input_path}' neexistuje!")
        matrix_df = pd.read_csv(input_path, index_col=0)
        ranking_file = input_path.parent / "protein_ranking_summary.csv"
        ranking_df = pd.read_csv(ranking_file) if ranking_file.exists() else None
        output_dir = input_path.parent

    # Převedeme na float
    matrix_df = matrix_df.apply(pd.to_numeric, errors="coerce")

    # Zkontrolujeme orientaci: chceme sloupce = proteiny, řádky = sloučeniny
    n_rows, n_cols = matrix_df.shape
    if n_cols > n_rows:
        # Původně: řádky = proteiny, sloupce = ligandy -> transponujeme
        matrix_df = matrix_df.T

    matrix_df.index.name = "compound_id"
    return matrix_df, ranking_df, output_dir


def calculate_target_selectivity(
    matrix_df: pd.DataFrame,
    target_id: str,
    potent_cutoff: float = -7.0,
    p_cutoff: float = 0.05,
    delta_cutoff: float = 1.0,
) -> pd.DataFrame:
    """
    Spočítá selektivitu a statistickou významnost pro každý ligand vůči zadanému proteinu.
    """
    if target_id not in matrix_df.columns:
        raise ValueError(f"Protein '{target_id}' nebyl nalezen v matici! Dostupné cíle: {list(matrix_df.columns[:5])}...")

    other_targets = [c for c in matrix_df.columns if c != target_id]
    if not other_targets:
        raise ValueError(f"Nelze provést analýzu selektivity: v matici není žádný jiný cíl kromě '{target_id}'!")

    results = []

    for compound_id, row in matrix_df.iterrows():
        target_score = row[target_id]
        if np.isnan(target_score):
            continue

        other_scores = row[other_targets].dropna().values
        if len(other_scores) < 3:
            continue

        median_others = float(np.median(other_scores))
        mean_others = float(np.mean(other_scores))
        std_others = float(np.std(other_scores, ddof=1))

        # Robustní směrodatná odchylka pomocí MAD (Median Absolute Deviation)
        mad = float(np.median(np.abs(other_scores - median_others)))
        robust_scale = max(1.4826 * mad, 0.08)  # Ochrana před nulovým rozptylem

        # Efekt: delta skóre (záporná hodnota = silnější na cíli)
        delta_score = target_score - median_others

        # Z-skóre a p-hodnota z robustního Studentova t-testu
        df_deg = len(other_scores) - 1
        z_score = (target_score - median_others) / robust_scale
        # Oboustranný p-value
        p_val = 2.0 * stats.t.sf(abs(z_score), df=df_deg)
        p_val = max(min(p_val, 1.0), 1e-15)  # Ochrana před 0.0 pro log10

        neg_log10_p = -np.log10(p_val)

        # Klasifikace
        # 1. Selektivní hit: silnější na cíli (delta <= -delta_cutoff) a statisticky signifikantní
        if delta_score <= -delta_cutoff and p_val <= p_cutoff:
            if target_score <= potent_cutoff:
                category = "Vysoce selektivní silný hit"
            else:
                category = "Selektivní (mírná afinita)"
        # 2. Preferuje jiné cíle: slabší na cíli než na zbytku
        elif delta_score >= delta_cutoff and p_val <= p_cutoff:
            category = "Preferuje ostatní proteiny"
        else:
            category = "Neselektivní / Nespecifický"

        results.append({
            "compound_id": str(compound_id),
            "target_id": target_id,
            "target_score": round(target_score, 2),
            "median_others": round(median_others, 2),
            "delta_score": round(delta_score, 2),
            "z_score": round(z_score, 2),
            "p_value": p_val,
            "neg_log10_p": round(neg_log10_p, 2),
            "category": category,
        })

    df_res = pd.DataFrame(results)
    if not df_res.empty:
        # Seřazení: nejvíce selektivní s nejvyšší statistickou významností nahoře
        df_res = df_res.sort_values(by=["delta_score", "neg_log10_p"], ascending=[True, False])

    return df_res


def plot_single_volcano(
    df_volcano: pd.DataFrame,
    target_id: str,
    out_png: Path,
    out_pdf: Path,
    p_cutoff: float = 0.05,
    delta_cutoff: float = 1.0,
    max_labels: int = 10,
    dpi: int = 300,
):
    """
    Vykreslí samostatný vysoce kvalitní Volcano Plot pro jeden cílový protein.
    """
    sns.set_theme(style="whitegrid", font="sans-serif")
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]

    fig, ax = plt.subplots(figsize=(9.0, 7.5), dpi=dpi)

    log10_p_thresh = -np.log10(p_cutoff)

    # Barevná mapa pro kategorie
    palette = {
        "Vysoce selektivní silný hit": "#b2182b",   # Sytá červená
        "Selektivní (mírná afinita)": "#ef8a62",    # Světle oranžová
        "Preferuje ostatní proteiny": "#2166ac",     # Modrá
        "Neselektivní / Nespecifický": "#999999",   # Neutrální šedá
    }

    # Vykreslení bodů podle kategorií (nejdříve šedé pozadí, pak významné body navrchu)
    for cat in ["Neselektivní / Nespecifický", "Preferuje ostatní proteiny", "Selektivní (mírná afinita)", "Vysoce selektivní silný hit"]:
        sub = df_volcano[df_volcano["category"] == cat]
        if sub.empty:
            continue
        alpha = 0.5 if cat.startswith("Neselektivní") else 0.9
        size = 35 if cat.startswith("Neselektivní") else 65
        zorder = 2 if cat.startswith("Neselektivní") else 4

        ax.scatter(
            sub["delta_score"],
            sub["neg_log10_p"],
            c=palette[cat],
            label=f"{cat} (n={len(sub)})",
            alpha=alpha,
            s=size,
            edgecolors="none" if cat.startswith("Neselektivní") else "black",
            linewidths=0.5,
            zorder=zorder,
        )

    # Prahové čáry
    ax.axhline(log10_p_thresh, color="black", linestyle="--", linewidth=1.0, alpha=0.7, zorder=1)
    ax.axvline(-delta_cutoff, color="crimson", linestyle=":", linewidth=1.0, alpha=0.8, zorder=1)
    ax.axvline(delta_cutoff, color="steelblue", linestyle=":", linewidth=1.0, alpha=0.8, zorder=1)

    # Textové anotace prahů
    x_min, x_max = ax.get_xlim()
    ax.text(x_max - 0.2, log10_p_thresh + 0.1, f"p = {p_cutoff} (-log10 = {log10_p_thresh:.1f})",
            fontsize=8, color="black", ha="right", va="bottom", fontstyle="italic")

    # Označení top hitů textovými popisky
    top_hits = df_volcano[
        (df_volcano["delta_score"] <= -delta_cutoff) &
        (df_volcano["neg_log10_p"] >= log10_p_thresh)
    ].head(max_labels)

    texts = []
    for _, r in top_hits.iterrows():
        t = ax.text(
            r["delta_score"],
            r["neg_log10_p"],
            f"{r['compound_id']} ({r['target_score']:.1f})",
            fontsize=8,
            fontweight="bold",
            color="#67001f",
        )
        texts.append(t)

    if texts:
        if ADJUST_TEXT_AVAILABLE:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="->", color="#67001f", lw=0.6))
        else:
            # Jednoduchý posun, pokud není adjustText
            for t in texts:
                t.set_position((t.get_position()[0], t.get_position()[1] + 0.15))

    # Šipky a popisky směrů selektivity v záhlaví
    ax.annotate("← Vyšší afinita k tomuto cíli (Selektivní)",
                xy=(0.03, 0.96), xycoords="axes fraction",
                fontsize=9.5, fontweight="bold", color="#b2182b")
    ax.annotate("Preferuje jiné methyltransferázy →",
                xy=(0.97, 0.96), xycoords="axes fraction",
                fontsize=9.5, fontweight="bold", color="#2166ac", ha="right")

    # Osy a titulky
    ax.set_xlabel("Rozdíl afinity: ΔScore = Score(cíl) - Medián(ostatní) [kcal/mol]", fontsize=11, fontweight="bold", labelpad=8)
    ax.set_ylabel("Statistická významnost: -log10(p-hodnota)", fontsize=11, fontweight="bold", labelpad=8)
    ax.set_title(f"Target Selectivity Volcano Plot: {target_id}\n(Srovnání vůči panelu ostatních proteinů)",
                 fontsize=13, fontweight="bold", pad=14)

    ax.legend(frameon=True, facecolor="white", edgecolor="lightgray", fontsize=8.5, loc="upper right")
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close("all")

    print(f"  [+] Volcano plot pro '{target_id}' uložen: {out_png}")


def plot_multi_volcano_grid(
    matrix_df: pd.DataFrame,
    targets: list[str],
    out_png: Path,
    out_pdf: Path,
    p_cutoff: float = 0.05,
    delta_cutoff: float = 1.0,
    dpi: int = 300,
):
    """
    Vykreslí přehledový multi-panel grid (např. 2x3 nebo 3x2) pro skupinu proteinů.
    """
    n_targets = len(targets)
    if n_targets == 0:
        return

    n_cols = 3 if n_targets >= 3 else n_targets
    n_rows = int(np.ceil(n_targets / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.0 * n_cols, 4.2 * n_rows), dpi=dpi, squeeze=False)
    log10_p_thresh = -np.log10(p_cutoff)

    for idx, target_id in enumerate(targets):
        row_idx = idx // n_cols
        col_idx = idx % n_cols
        ax = axes[row_idx][col_idx]

        try:
            df_v = calculate_target_selectivity(matrix_df, target_id, p_cutoff=p_cutoff, delta_cutoff=delta_cutoff)
        except Exception as e:
            ax.set_title(f"{target_id} (Chyba: {e})", fontsize=10)
            continue

        # Body
        non_sel = df_v[df_v["category"].str.startswith("Neselektivní")]
        sel = df_v[df_v["category"].str.contains("selektivní", case=False)]
        counter = df_v[df_v["category"].str.contains("Preferuje", case=False)]

        ax.scatter(non_sel["delta_score"], non_sel["neg_log10_p"], c="#aaaaaa", alpha=0.4, s=20, label="Nespecifické")
        ax.scatter(counter["delta_score"], counter["neg_log10_p"], c="#2166ac", alpha=0.8, s=35, label="Ostatní")
        ax.scatter(sel["delta_score"], sel["neg_log10_p"], c="#b2182b", alpha=0.9, s=45, edgecolors="black", linewidths=0.5, label="Selektivní hity")

        # Prahové linky
        ax.axhline(log10_p_thresh, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.axvline(-delta_cutoff, color="crimson", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.axvline(delta_cutoff, color="steelblue", linestyle=":", linewidth=0.8, alpha=0.7)

        # Anotace top 3 hitů
        top3 = sel.head(3)
        for _, r in top3.iterrows():
            ax.text(r["delta_score"], r["neg_log10_p"] + 0.1, r["compound_id"], fontsize=7, fontweight="bold", color="#67001f")

        ax.set_title(f"Cíl: {target_id} (Hity: {len(sel)})", fontsize=11, fontweight="bold")
        ax.set_xlabel("ΔScore [kcal/mol]", fontsize=9)
        ax.set_ylabel("-log10(p-val)", fontsize=9)
        ax.tick_params(labelsize=8)

    # Skrytí nepoužitých sub-plotů
    for idx in range(n_targets, n_rows * n_cols):
        r = idx // n_cols
        c = idx % n_cols
        axes[r][c].set_visible(False)

    fig.suptitle("Souhrnný přehled selektivity ligandů napříč vybranými methyltransferázami", fontsize=14, fontweight="bold", y=0.995)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close("all")

    print(f"  [+] Multi-panel Volcano grid uložen: {out_png}")


def main():
    parser = argparse.ArgumentParser(
        description="Generování Target Selectivity Volcano Plotů (Varianta A) pro virtuální screening.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Příklady použití:
  # 1. Volcano plot pro konkrétní cíl (např. A2DJA1):
  python plot_docking_volcano.py -i reports_docking -t A2DJA1

  # 2. Volcano ploty pro TOP 5 nejlepších proteinů ze souhrnného rankingu:
  python plot_docking_volcano.py -i reports_docking --top 5

  # 3. Vytvoření multi-panel gridu pro TOP 6 cílů v jednom obrázku:
  python plot_docking_volcano.py -i reports_docking --top 6 --grid

  # 4. Pro všechny proteiny:
  python plot_docking_volcano.py -i reports_docking --all
        """
    )
    parser.add_argument(
        "-i", "--input",
        type=Path,
        default=Path("reports_docking"),
        help="Cesta k 'docking_scores_matrix.csv' nebo adresáři s reporty (výchozí: 'reports_docking')."
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=None,
        help="Cílová složka pro grafy a tabulky (výchozí: <input_dir>/volcano_plots)."
    )
    parser.add_argument(
        "-t", "--target",
        type=str,
        default=None,
        help="Konkrétní ID cílového proteinu (např. A2DJA1). Lze zadat i více cílů oddělených čárkou."
    )
    parser.add_argument(
        "--top",
        type=int,
        default=None,
        help="Automaticky vybere N nejlepších proteinů z 'protein_ranking_summary.csv'."
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Vygeneruje volcano ploty pro všech 49 proteinů v panelu."
    )
    parser.add_argument(
        "--grid",
        action="store_true",
        help="Kromě samostatných grafů vygeneruje také sdružený multi-panel přehled (grid)."
    )
    parser.add_argument(
        "--p-cutoff",
        type=float,
        default=0.05,
        help="Prahová hodnota p-value pro statistickou významnost (výchozí: 0.05)."
    )
    parser.add_argument(
        "--delta-cutoff",
        type=float,
        default=1.0,
        help="Minimální rozdíl afinity v kcal/mol pro selektivitu (výchozí: 1.0 kcal/mol)."
    )
    parser.add_argument(
        "--potent-cutoff",
        type=float,
        default=-7.0,
        help="Prahová absolutní afinita v kcal/mol pro označení silného hitu (výchozí: -7.0 kcal/mol)."
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Rozlišení výstupních grafů (DPI, výchozí: 300)."
    )

    args = parser.parse_args()

    # 1. Načtení dat
    try:
        matrix_df, ranking_df, base_out_dir = load_data(args.input)
    except Exception as e:
        print(f"[!] Chyba při načítání dat: {e}", file=sys.stderr)
        sys.exit(1)

    out_dir = args.output_dir if args.output_dir else (base_out_dir / "volcano_plots")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 2. Určení seznamu cílů k analýze
    all_targets = list(matrix_df.columns)
    selected_targets = []

    if args.target:
        # Zadaný konkrétní cíl nebo seznam cílů
        for t in args.target.split(","):
            t_clean = t.strip()
            if t_clean in all_targets:
                selected_targets.append(t_clean)
            else:
                print(f"[!] Varování: Cíl '{t_clean}' nebyl nalezen v matici skóre!", file=sys.stderr)
    elif args.all:
        selected_targets = all_targets
    elif args.top:
        if ranking_df is not None and "target_id" in ranking_df.columns:
            ranked_list = [t for t in ranking_df["target_id"].tolist() if t in all_targets]
            selected_targets = ranked_list[:args.top]
        else:
            # Pokud nemáme ranking, vybereme prvních N
            selected_targets = all_targets[:args.top]
    else:
        # Výchozí chování: pokud máme ranking, vezmeme TOP 5 cílů, jinak první 3
        if ranking_df is not None and "target_id" in ranking_df.columns:
            selected_targets = [t for t in ranking_df["target_id"].tolist() if t in all_targets][:5]
            print(f"[*] Nebyl zadán parametr -t ani --top. Automaticky vybírám TOP 5 cílů z rankingu: {', '.join(selected_targets)}")
        else:
            selected_targets = all_targets[:3]
            print(f"[*] Nebyl zadán parametr -t ani --top. Vybírám první 3 cíle: {', '.join(selected_targets)}")

    if not selected_targets:
        print("[!] Nebyly vybrány žádné platné cíle k analýze!", file=sys.stderr)
        sys.exit(1)

    print("=" * 70)
    print("Generování Target Selectivity Volcano Plotů (Varianta A):")
    print(f"  Analyzované cíle ({len(selected_targets)}): {', '.join(selected_targets)}")
    print(f"  Prahová p-hodnota:     {args.p_cutoff} (-log10 = {-np.log10(args.p_cutoff):.2f})")
    print(f"  Prahová selektivita:   |ΔScore| >= {args.delta_cutoff} kcal/mol")
    print(f"  Prahová afinita hitu:  Score <= {args.potent_cutoff} kcal/mol")
    print(f"  Výstupní složka:       {out_dir}")
    print("=" * 70)

    # 3. Zpracování jednotlivých cílů
    summary_selective_rows = []

    for target_id in selected_targets:
        df_volcano = calculate_target_selectivity(
            matrix_df=matrix_df,
            target_id=target_id,
            potent_cutoff=args.potent_cutoff,
            p_cutoff=args.p_cutoff,
            delta_cutoff=args.delta_cutoff,
        )

        if df_volcano.empty:
            print(f"[!] Pro cíl '{target_id}' se nepodařilo spočítat žádné hodnoty.")
            continue

        # Uložení tabulky výsledků pro daný cíl
        csv_path = out_dir / f"selectivity_table_{target_id}.csv"
        df_volcano.to_csv(csv_path, index=False)

        # Uložení Volcano plotu
        png_path = out_dir / f"volcano_{target_id}.png"
        pdf_path = out_dir / f"volcano_{target_id}.pdf"
        plot_single_volcano(
            df_volcano=df_volcano,
            target_id=target_id,
            out_png=png_path,
            out_pdf=pdf_path,
            p_cutoff=args.p_cutoff,
            delta_cutoff=args.delta_cutoff,
            dpi=args.dpi,
        )

        # Agregace top selektivních hitů
        top_sel = df_volcano[df_volcano["category"].str.contains("selektivní", case=False)]
        if not top_sel.empty:
            summary_selective_rows.append(top_sel)

    # 4. Multi-panel přehledový grid (pokud je požadován nebo analyzujeme >= 4 cíle)
    if args.grid or len(selected_targets) in [4, 6, 8, 9, 12]:
        grid_png = out_dir / "volcano_grid_overview.png"
        grid_pdf = out_dir / "volcano_grid_overview.pdf"
        plot_multi_volcano_grid(
            matrix_df=matrix_df,
            targets=selected_targets[:12],  # Maximálně 12 panelů pro přehlednost
            out_png=grid_png,
            out_pdf=grid_pdf,
            p_cutoff=args.p_cutoff,
            delta_cutoff=args.delta_cutoff,
            dpi=args.dpi,
        )

    # 5. Společná souhrnná tabulka všech selektivních hitů napříč analyzovanými proteiny
    if summary_selective_rows:
        master_selective_df = pd.concat(summary_selective_rows, ignore_index=True)
        master_csv = out_dir / "all_selective_hits_summary.csv"
        master_selective_df.to_csv(master_csv, index=False)
        print(f"\n[+] Souhrnná tabulka všech selektivních hitů uložena: {master_csv}")

    print("\n[✓] Hotovo! Všechny Volcano ploty a tabulky selektivity byly vygenerovány.")


if __name__ == "__main__":
    main()
