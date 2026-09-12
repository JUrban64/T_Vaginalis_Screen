#!/usr/bin/env python3
# ==============================================================================
# Skript pro vizualizaci dokovacích dat: Protein-Level Volcano Plot
# (Globální srovnání všech 49 proteinů napříč celou knihovnou ligandů)
# ==============================================================================
# Každý bod v grafu představuje JEDEN PROTEIN (z rodiny methyltransferáz).
#
# Osa X: Posun afinity vůči celému panelu (Effect Size / Druggability Shift):
#        ΔAffinity = Medián(celý panel) - Průměr(daný protein) [kcal/mol]
#        [Kladné hodnoty = protein váže knihovnu signifikantně silněji než průměr]
#
# Osa Y: Statistická významnost:
#        -log10(p-hodnota) z párového t-testu napříč všemi ligandy (~160 měření)
#
# Velikost bodu: Počet silných vazačů (potent binders, např. skóre <= -8.0 kcal/mol)
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


def load_matrix(input_path: Path) -> tuple[pd.DataFrame, pd.DataFrame | None, Path]:
    """Načte matici skóre a případný souhrnný ranking proteinů."""
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

    matrix_df = matrix_df.apply(pd.to_numeric, errors="coerce")

    # Ujistíme se: řádky = sloučeniny, sloupce = proteiny
    n_rows, n_cols = matrix_df.shape
    if n_cols > n_rows:
        matrix_df = matrix_df.T

    matrix_df.index.name = "compound_id"
    return matrix_df, ranking_df, output_dir


def compute_protein_volcano_statistics(
    matrix_df: pd.DataFrame,
    potent_cutoff: float = -8.0,
    p_cutoff: float = 0.05,
    delta_cutoff: float = 0.5,
) -> pd.DataFrame:
    """
    Spočítá statistiky pro každý protein napříč všemi ligandy:
    - Osa X: Posun afinity vůči mediánu zbytku panelu
    - Osa Y: p-hodnota z párového t-testu napříč všemi ligandy
    """
    targets = list(matrix_df.columns)
    results = []

    # Pro každý ligand spočítáme medián napříč celým panelem jako referenci
    global_compound_medians = matrix_df.median(axis=1)

    for target_id in targets:
        target_series = matrix_df[target_id].dropna()
        if len(target_series) < 5:
            continue

        # Párová data: skóre tohoto proteinu vs medián všech ostatních proteinů pro stejné ligandy
        other_targets = [t for t in targets if t != target_id]
        if not other_targets:
            continue

        compounds_idx = target_series.index
        target_scores = target_series.values
        other_medians = matrix_df.loc[compounds_idx, other_targets].median(axis=1).values

        # Rozdíl pro každý ligand: (medián ostatních - cíl)
        # Protože nižší skóre = silnější vazba:
        # Pokud cíl má -9 a ostatní -7 -> rozdíl = -7 - (-9) = +2 kcal/mol (cíl je o 2 kcal/mol silnější!)
        diffs = other_medians - target_scores

        mean_diff = float(np.mean(diffs))
        median_diff = float(np.median(diffs))
        mean_score = float(np.mean(target_scores))
        best_score = float(np.min(target_scores))

        # Počet silných vazačů (<= potent_cutoff)
        potent_binders = int((target_scores <= potent_cutoff).sum())

        # Párový Studentův t-test napříč ligandy
        t_stat, p_val = stats.ttest_rel(target_scores, other_medians, nan_policy="omit")
        p_val = max(min(float(p_val), 1.0), 1e-15)
        neg_log10_p = -np.log10(p_val)

        # Kategorizace proteinu
        if mean_diff >= delta_cutoff and p_val <= p_cutoff:
            category = "Vysoce citlivý cíl (Super-binder)"
        elif mean_diff <= -delta_cutoff and p_val <= p_cutoff:
            category = "Rekalcitrantní cíl (Slabá afinita)"
        else:
            category = "Průměrný protein panelu"

        results.append({
            "target_id": target_id,
            "mean_docking_score": round(mean_score, 2),
            "best_docking_score": round(best_score, 2),
            "potent_binders_count": potent_binders,
            "delta_affinity": round(mean_diff, 2),
            "t_statistic": round(float(t_stat), 2),
            "p_value": p_val,
            "neg_log10_p": round(neg_log10_p, 2),
            "category": category,
        })

    df_res = pd.DataFrame(results)
    if not df_res.empty:
        # Seřazení: největší posun afinity a nejvyšší statistická významnost nahoře
        df_res = df_res.sort_values(by=["delta_affinity", "neg_log10_p"], ascending=[False, False])

    return df_res


def plot_protein_volcano_figure(
    df_stats: pd.DataFrame,
    title: str,
    out_png: Path,
    out_pdf: Path,
    p_cutoff: float = 0.05,
    delta_cutoff: float = 0.5,
    dpi: int = 300,
):
    """Vykreslí přehledný Protein-Level Volcano plot."""
    sns.set_theme(style="whitegrid", font="sans-serif")
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]

    fig, ax = plt.subplots(figsize=(10.5, 8.0), dpi=dpi)

    log10_p_thresh = -np.log10(p_cutoff)

    # Barevná paleta
    palette = {
        "Vysoce citlivý cíl (Super-binder)": "#d73027",  # Červená
        "Rekalcitrantní cíl (Slabá afinita)": "#4575b4",  # Modrá
        "Průměrný protein panelu": "#999999",            # Šedá
    }

    # Škálování velikosti bodů podle počtu silných vazačů (potent binders)
    min_size = 50
    max_size = 250
    max_binders = max(df_stats["potent_binders_count"].max(), 1)
    sizes = min_size + (df_stats["potent_binders_count"] / max_binders) * (max_size - min_size)

    # Vykreslení jednotlivých kategorií
    for cat in ["Průměrný protein panelu", "Rekalcitrantní cíl (Slabá afinita)", "Vysoce citlivý cíl (Super-binder)"]:
        sub = df_stats[df_stats["category"] == cat]
        if sub.empty:
            continue
        sub_sizes = sizes.loc[sub.index]
        alpha = 0.55 if cat.startswith("Průměrný") else 0.9
        zorder = 2 if cat.startswith("Průměrný") else 4

        ax.scatter(
            sub["delta_affinity"],
            sub["neg_log10_p"],
            c=palette[cat],
            label=f"{cat} (n={len(sub)})",
            s=sub_sizes,
            alpha=alpha,
            edgecolors="black" if not cat.startswith("Průměrný") else "none",
            linewidths=0.7,
            zorder=zorder,
        )

    # Prahové čáry
    ax.axhline(log10_p_thresh, color="black", linestyle="--", linewidth=1.0, alpha=0.7, zorder=1)
    ax.axvline(delta_cutoff, color="crimson", linestyle=":", linewidth=1.0, alpha=0.8, zorder=1)
    ax.axvline(-delta_cutoff, color="steelblue", linestyle=":", linewidth=1.0, alpha=0.8, zorder=1)

    # Označení všech signifikantních proteinů textovými popisky
    sig_targets = df_stats[
        (df_stats["neg_log10_p"] >= log10_p_thresh) &
        (df_stats["delta_affinity"].abs() >= delta_cutoff)
    ]

    texts = []
    for _, r in sig_targets.iterrows():
        t = ax.text(
            r["delta_affinity"],
            r["neg_log10_p"],
            f" {r['target_id']}",
            fontsize=8.5,
            fontweight="bold",
            color="#67001f" if r["delta_affinity"] > 0 else "#08306b",
        )
        texts.append(t)

    if texts:
        if ADJUST_TEXT_AVAILABLE:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="->", color="gray", lw=0.6))
        else:
            for t in texts:
                t.set_position((t.get_position()[0], t.get_position()[1] + 0.15))

    # Směrové šipky a popisy
    ax.annotate("Váže knihovnu inhibitorů významně SILNĚJI než průměr →\n(Primární cíle / Druggable Hotspots)",
                xy=(0.98, 0.95), xycoords="axes fraction",
                fontsize=9.5, fontweight="bold", color="#d73027", ha="right")
    ax.annotate("← Váže knihovnu SLABĚJI než průměr\n(Rekalcitrantní enzymy)",
                xy=(0.02, 0.95), xycoords="axes fraction",
                fontsize=9.5, fontweight="bold", color="#4575b4", ha="left")

    # Osy a popisky
    ax.set_xlabel("Posun afinity: ΔAffinity = Medián(panel) - Průměr(protein) [kcal/mol]\n(kladné = silnější vazba)",
                  fontsize=11, fontweight="bold", labelpad=8)
    ax.set_ylabel("Statistická významnost napříč ligandy: -log10(p-hodnota)\n(párový t-test přes ~160 sloučenin)",
                  fontsize=11, fontweight="bold", labelpad=8)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=14)

    # Legenda
    leg = ax.legend(frameon=True, facecolor="white", edgecolor="lightgray", fontsize=9, loc="upper center", bbox_to_anchor=(0.5, 0.90))
    leg.set_zorder(10)

    # Pomocná textová poznámka o velikosti bodů
    ax.text(0.02, 0.03, "Poznámka: Velikost bodu odpovídá počtu vysoce silných inhibitorů (skóre <= -8.0 kcal/mol)",
            transform=ax.transAxes, fontsize=8, fontstyle="italic", color="#555555")

    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close("all")

    print(f"  [+] Protein-level Volcano plot uložen (PNG): {out_png}")
    print(f"  [+] Protein-level Volcano plot uložen (PDF): {out_pdf}")


def main():
    parser = argparse.ArgumentParser(
        description="Generování Protein-Level Volcano Plotu (Globální srovnání všech proteinů napříč ligandy).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Příklady použití:
  python plot_protein_volcano.py -i reports_docking
  python plot_protein_volcano.py -i reports_docking --p-cutoff 0.01 --delta-cutoff 0.75
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
        help="Cílová složka pro graf a tabulku (výchozí: <input_dir>/protein_volcano)."
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
        default=0.5,
        help="Prahový posun afinity v kcal/mol (výchozí: 0.5 kcal/mol)."
    )
    parser.add_argument(
        "--potent-cutoff",
        type=float,
        default=-8.0,
        help="Hranice pro počítání silných vazačů (výchozí: -8.0 kcal/mol)."
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Rozlišení výstupních rastrů (DPI, výchozí: 300)."
    )

    args = parser.parse_args()

    # 1. Načtení dat
    try:
        matrix_df, ranking_df, base_out_dir = load_matrix(args.input)
    except Exception as e:
        print(f"[!] Chyba při načítání matice skóre: {e}", file=sys.stderr)
        sys.exit(1)

    out_dir = args.output_dir if args.output_dir else (base_out_dir / "protein_volcano")
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Generování Protein-Level Volcano Plotu (Globální srovnání proteinů):")
    print(f"  Celkem proteinů v panelu:  {matrix_df.shape[1]}")
    print(f"  Celkem ligandů v knihovně: {matrix_df.shape[0]}")
    print(f"  Prahová p-hodnota:          {args.p_cutoff} (-log10 = {-np.log10(args.p_cutoff):.2f})")
    print(f"  Prahový posun afinity:      |ΔAffinity| >= {args.delta_cutoff} kcal/mol")
    print(f"  Výstupní složka:            {out_dir}")
    print("=" * 70)

    # 2. Výpočet statistik pro každý protein
    df_stats = compute_protein_volcano_statistics(
        matrix_df=matrix_df,
        potent_cutoff=args.potent_cutoff,
        p_cutoff=args.p_cutoff,
        delta_cutoff=args.delta_cutoff,
    )

    if df_stats.empty:
        print("[!] Nepodařilo se spočítat statistiky pro žádný protein!", file=sys.stderr)
        sys.exit(1)

    # Uložení tabulky výsledků
    csv_path = out_dir / "protein_volcano_statistics.csv"
    df_stats.to_csv(csv_path, index=False)
    print(f"[+] Tabulka statistik proteinů uložena: {csv_path}")

    # 3. Vykreslení grafu
    png_path = out_dir / "protein_volcano_plot.png"
    pdf_path = out_dir / "protein_volcano_plot.pdf"

    title = f"Protein-Level Volcano Plot: Srovnání afinity 49 methyltransferáz T. vaginalis\n(Párová statistika přes {matrix_df.shape[0]} testovaných inhibitorů)"

    plot_protein_volcano_figure(
        df_stats=df_stats,
        title=title,
        out_png=png_path,
        out_pdf=pdf_path,
        p_cutoff=args.p_cutoff,
        delta_cutoff=args.delta_cutoff,
        dpi=args.dpi,
    )

    # Výpis TOP 5 nejlepších proteinů
    top_super = df_stats[df_stats["category"].str.contains("Super-binder")].head(5)
    if not top_super.empty:
        print("\n" + "=" * 70)
        print("TOP VYSOCE CITLIVÉ CÍLE KNIHOVNY (Super-binders):")
        print(top_super[["target_id", "delta_affinity", "mean_docking_score", "best_docking_score", "potent_binders_count", "neg_log10_p"]].to_string(index=False))
        print("=" * 70)

    print("\n[✓] Hotovo! Protein-Level Volcano plot byl úspěšně vygenerován.")


if __name__ == "__main__":
    main()
