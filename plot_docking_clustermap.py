#!/usr/bin/env python3
# ==============================================================================
# Skript pro vizualizaci dokovacích dat: Hierarchická Heatmapa (Clustermap)
# ==============================================================================
# Vytvoří oboustranně shlukovanou heatmapu (Clustermap) interakční matice:
#   Proteiny (49 methyltransferáz) x Ligandy (~160 inhibitorů)
#
# Funkce:
#   1. Shlukování proteinů (řádky) podle podobnosti vazebných profilů.
#   2. Shlukování ligandů (sloupce) podle podobnosti jejich selektivity.
#   3. Podpora surových skóre (kcal/mol), Z-skóre i reziduálních skóre bez vlivu MW.
#   4. Export do publikační kvality (PNG 300 DPI + vektorové PDF).
#
# Autor: Antigravity AI / Jáchym Urban
# ==============================================================================

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Kontrola dostupnosti grafických knihoven
try:
    import matplotlib
    matplotlib.use("Agg")  # Bezhlavý backend pro HPC / servery
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.cluster import hierarchy
    from scipy.spatial.distance import pdist
except ImportError as e:
    print(f"[!] Chyba: Chybí potřebná knihovna pro vykreslování: {e}", file=sys.stderr)
    print("    Nainstalujte: pip install matplotlib seaborn scipy pandas numpy", file=sys.stderr)
    print("    Nebo aktivujte prostředí s těmito balíčky (např. conda activate cadd2026).", file=sys.stderr)
    sys.exit(1)


def load_matrix(input_path: Path) -> tuple[pd.DataFrame, Path]:
    """
    Načte matici z CSV souboru nebo z adresáře s reporty.
    Vrací DataFrame a adresář pro uložení výstupů.
    """
    if input_path.is_dir():
        score_file = input_path / "docking_scores_matrix.csv"
        if not score_file.exists():
            # Zkusíme najít jakékoliv maticové CSV
            candidates = list(input_path.glob("*matrix*.csv"))
            if not candidates:
                raise FileNotFoundError(f"V adresáři '{input_path}' nebyl nalezen 'docking_scores_matrix.csv'!")
            score_file = candidates[0]
        matrix_df = pd.read_csv(score_file, index_col=0)
        output_dir = input_path
    else:
        if not input_path.exists():
            raise FileNotFoundError(f"Soubor '{input_path}' neexistuje!")
        matrix_df = pd.read_csv(input_path, index_col=0)
        output_dir = input_path.parent

    # Ujistíme se, že hodnoty jsou float
    matrix_df = matrix_df.apply(pd.to_numeric, errors="coerce")
    return matrix_df, output_dir


def prepare_matrix_orientation(df: pd.DataFrame) -> pd.DataFrame:
    """
    Upraví orientaci matice tak, aby řádky byly proteiny (typicky ~49)
    a sloupce byly ligandy (typicky ~160).
    Pokud je tomu naopak, matici transponuje.
    """
    n_rows, n_cols = df.shape
    # Předpoklad: proteinů je méně (49) než ligandů (160)
    if n_rows > n_cols:
        # Původně: řádky = sloučeniny, sloupce = proteiny -> transponujeme
        df_oriented = df.T
    else:
        df_oriented = df.copy()

    df_oriented.index.name = "Target Protein"
    df_oriented.columns.name = "Ligand / Compound"
    return df_oriented


def impute_missing(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ošetří případné chybějící hodnoty (např. selhání dokování pro daný pár).
    Chybějící buňky nahradí nejhorší pozorovanou hodnotou (nejméně stabilní vazba).
    """
    if df.isna().any().any():
        worst_score = df.max().max() + 0.5  # Mírně horší než nejhorší nalezené skóre
        print(f"[*] Varování: Matice obsahuje {df.isna().sum().sum()} prázdných hodnot (NaN).")
        print(f"    Chybějící hodnoty budou pro shlukování nahrazeny skóre: {worst_score:.2f} kcal/mol.")
        return df.fillna(worst_score)
    return df


def plot_clustermap_figure(
    df: pd.DataFrame,
    title: str,
    cbar_label: str,
    out_png: Path,
    out_pdf: Path,
    cmap: str = "viridis_r",
    standardize: str = "none",
    metric: str = "euclidean",
    method: str = "average",
    figsize: tuple[float, float] = (14.0, 10.0),
    dpi: int = 300,
    show_dendrograms: bool = True,
):
    """
    Vykreslí a uloží hierarchickou clustermapu pomocí Seaborn.
    """
    # Pokud se použije jiná metrika než euclidean, 'ward' nelze použít
    if metric != "euclidean" and method == "ward":
        method = "average"

    # Z-skóre standardizace
    z_score_val = None
    if standardize == "compound":
        # Standardizace napříč proteiny pro každý ligand (osa sloupců = 1)
        z_score_val = 1
    elif standardize == "target":
        # Standardizace napříč ligandy pro každý protein (osa řádků = 0)
        z_score_val = 0

    # Nastavení vizuálního stylu
    sns.set_theme(style="white", font="sans-serif")
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]

    dendro_ratio = (0.12, 0.12) if show_dendrograms else (0.001, 0.001)

    # Vytvoření clustermapy
    g = sns.clustermap(
        df,
        figsize=figsize,
        cmap=cmap,
        metric=metric,
        method=method,
        z_score=z_score_val,
        dendrogram_ratio=dendro_ratio,
        cbar_pos=(0.02, 0.82, 0.03, 0.14),  # Umístění colorbaru vlevo nahoře
        linewidths=0.2,
        linecolor="whitesmoke",
        row_cluster=show_dendrograms,
        col_cluster=show_dendrograms,
        xticklabels=True if df.shape[1] <= 80 else False,  # Popisky sloupců jen pokud se vejdou
        yticklabels=True,  # Proteiny vypisujeme vždy
    )

    # Úprava popisků os a fontů
    g.ax_heatmap.set_ylabel("Cílové proteiny (49 methyltransferáz)", fontsize=12, fontweight="bold", labelpad=12)
    g.ax_heatmap.set_xlabel(f"Knihovna inhibitorů ({df.shape[1]} sloučenin)", fontsize=12, fontweight="bold", labelpad=12)
    g.ax_heatmap.tick_params(axis="y", labelsize=8, rotation=0)

    if df.shape[1] <= 80:
        g.ax_heatmap.tick_params(axis="x", labelsize=7, rotation=90)

    # Titulek a popisek barevné škály
    g.ax_cbar.set_title(cbar_label, fontsize=9, pad=10, fontweight="bold")
    g.ax_cbar.tick_params(labelsize=8)

    # Hlavní nadpis celé figury
    g.fig.suptitle(title, fontsize=15, fontweight="bold", y=0.98)

    # Uložení do souborů
    out_png.parent.mkdir(parents=True, exist_ok=True)
    g.savefig(out_png, dpi=dpi, bbox_inches="tight")
    g.savefig(out_pdf, bbox_inches="tight")
    plt.close("all")

    print(f"  [+] Clustermapa uložena (PNG): {out_png}")
    print(f"  [+] Clustermapa uložena (PDF): {out_pdf}")


def main():
    parser = argparse.ArgumentParser(
        description="Generování hierarchické Clustermapy dokovacích skóre (Proteiny x Ligandy).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Příklady použití:
  python plot_docking_clustermap.py -i reports_docking
  python plot_docking_clustermap.py -i reports_docking/docking_scores_matrix.csv --both
  python plot_docking_clustermap.py -i reports_docking --standardize compound
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
        help="Cílová složka pro grafy (výchozí: <input_dir>/plots)."
    )
    parser.add_argument(
        "--both",
        action="store_true",
        help="Vykreslí jak surová skóre (kcal/mol), tak skóre očištěná od molekulové hmotnosti (MW-rezidua)."
    )
    parser.add_argument(
        "--standardize",
        choices=["none", "compound", "target"],
        default="none",
        help="Typ Z-skóre standardizace: 'none' (surová skóre), 'compound' (normalizace po sloučeninách pro odhalení selektivity), 'target' (normalizace po proteinech)."
    )
    parser.add_argument(
        "--cmap",
        type=str,
        default="viridis_r",
        help="Barevná paleta matplotlib/seaborn (výchozí: 'viridis_r', doporučeno také 'mako_r', 'coolwarm_r')."
    )
    parser.add_argument(
        "--metric",
        type=str,
        default="euclidean",
        choices=["euclidean", "correlation", "cosine", "cityblock"],
        help="Metrika vzdálenosti pro shlukování (výchozí: 'euclidean')."
    )
    parser.add_argument(
        "--method",
        type=str,
        default="ward",
        choices=["ward", "average", "complete", "single"],
        help="Metoda vazby shluků (výchozí: 'ward')."
    )
    parser.add_argument(
        "--figsize",
        type=str,
        default="14,10",
        help="Rozměry grafu v palcích jako 'šířka,výška' (výchozí: '14,10')."
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Rozlišení výstupního rastru (DPI, výchozí: 300)."
    )

    args = parser.parse_args()

    # Parsování rozměrů
    try:
        w, h = [float(x.strip()) for x in args.figsize.split(",")]
        figsize = (w, h)
    except Exception:
        figsize = (14.0, 10.0)

    # 1. Načtení základní matice dokovacích skóre
    try:
        matrix_df, base_out_dir = load_matrix(args.input)
    except Exception as e:
        print(f"[!] Chyba při načítání matice: {e}", file=sys.stderr)
        sys.exit(1)

    out_dir = args.output_dir if args.output_dir else (base_out_dir / "plots")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Úprava orientace a ošetření NaN
    df_oriented = prepare_matrix_orientation(matrix_df)
    df_imputed = impute_missing(df_oriented)

    print("=" * 70)
    print("Vytváření hierarchické Clustermapy dokovacích skóre:")
    print(f"  Vstupní matice:    {df_imputed.shape[0]} proteinů x {df_imputed.shape[1]} sloučenin")
    print(f"  Metrika/Metoda:    {args.metric} / {args.method}")
    print(f"  Standardizace:     {args.standardize}")
    print(f"  Výstupní adresář:  {out_dir}")
    print("=" * 70)

    # Definice popisků
    if args.standardize == "compound":
        cbar_label = "Z-skóre\n(selektivita\nligandu)"
        title_suffix = " [Standardizováno po ligandech (Z-skóre)]"
        cmap_to_use = "coolwarm_r" if args.cmap == "viridis_r" else args.cmap
    elif args.standardize == "target":
        cbar_label = "Z-skóre\n(relativní\nafinita k cíli)"
        title_suffix = " [Standardizováno po cílech (Z-skóre)]"
        cmap_to_use = "coolwarm_r" if args.cmap == "viridis_r" else args.cmap
    else:
        cbar_label = "Dokovací afinita\nΔG [kcal/mol]\n(nižší = silnější)"
        title_suffix = " [Dokovací afinita ΔG]"
        cmap_to_use = args.cmap

    # A) Vykreslení Clustermapy surových skóre
    title_main = f"Interakční Clustermapa: Methyltransferázy T. vaginalis vs. Inhibitory{title_suffix}"
    out_png_scores = out_dir / f"clustermap_docking_scores_{args.standardize}.png"
    out_pdf_scores = out_dir / f"clustermap_docking_scores_{args.standardize}.pdf"

    plot_clustermap_figure(
        df=df_imputed,
        title=title_main,
        cbar_label=cbar_label,
        out_png=out_png_scores,
        out_pdf=out_pdf_scores,
        cmap=cmap_to_use,
        standardize=args.standardize,
        metric=args.metric,
        method=args.method,
        figsize=figsize,
        dpi=args.dpi,
    )

    # B) Volitelné vykreslení matice reziduí očištěných od hmotnosti (--both)
    residuals_file = base_out_dir / "docking_mw_residuals_matrix.csv"
    if (args.both or "residuals" in str(args.input)) and residuals_file.exists():
        print("\n[*] Nalezena matice MW-reziduí. Generuji Clustermapu bez vlivu molekulové hmotnosti...")
        try:
            res_df = pd.read_csv(residuals_file, index_col=0).apply(pd.to_numeric, errors="coerce")
            res_oriented = prepare_matrix_orientation(res_df)
            res_imputed = impute_missing(res_oriented)

            out_png_res = out_dir / "clustermap_docking_mw_residuals.png"
            out_pdf_res = out_dir / "clustermap_docking_mw_residuals.pdf"

            plot_clustermap_figure(
                df=res_imputed,
                title="Clustermapa dokování očištěná od molekulové hmotnosti (MW Residuals)",
                cbar_label="Reziduum afinity\n[kcal/mol]\n(záporné = lepší než MW)",
                out_png=out_png_res,
                out_pdf=out_pdf_res,
                cmap="coolwarm_r",
                standardize="none",
                metric=args.metric,
                method=args.method,
                figsize=figsize,
                dpi=args.dpi,
            )
        except Exception as e:
            print(f"[!] Selhalo vykreslení clustermapy reziduí: {e}", file=sys.stderr)

    print("\n[✓] Hotovo! Všechny clustermapy byly úspěšně vygenerovány.")


if __name__ == "__main__":
    main()
