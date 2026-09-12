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
import re
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


def load_protein_annotations(search_dir: Path | None = None) -> dict[str, str]:
    """Načte anotace cílů pro výpis v titulku grafu (např. 'A2FE15 (ICMT)')."""
    candidates = [
        Path("reports/deduplicated_proteins.tsv"),
        Path("reports_deduplicated/download_summary.csv"),
        Path("targets.txt"),
    ]
    if search_dir:
        candidates.insert(0, search_dir / "reports/deduplicated_proteins.tsv")
        candidates.insert(1, search_dir.parent / "reports/deduplicated_proteins.tsv")

    annots = {}
    for c in candidates:
        if c.exists():
            try:
                df = pd.read_csv(c, sep="\t" if c.suffix == ".tsv" else ",")
                id_col = "Entry" if "Entry" in df.columns else df.columns[0]
                gene_col = "Gene Names" if "Gene Names" in df.columns else None
                entry_col = "Entry Name" if "Entry Name" in df.columns else None

                for _, row in df.iterrows():
                    tid = str(row[id_col]).strip()
                    label = ""
                    if gene_col and pd.notna(row[gene_col]):
                        label = str(row[gene_col]).strip().split()[0]
                    elif entry_col and pd.notna(row[entry_col]):
                        ename = str(row[entry_col]).strip().split("_")[0]
                        if ename != tid:
                            label = ename
                    annots[tid] = f"{tid} ({label})" if label else tid
                if annots:
                    break
            except Exception:
                pass
    return annots


def is_likely_target_id(val: str, known_targets: set[str] | None = None) -> bool:
    """Ověří, zda je daný řetězec skutečně protein (UniProt ID)."""
    if not isinstance(val, str):
        return False
    val_clean = val.strip()
    if known_targets and val_clean in known_targets:
        return True
    if any(c in val_clean for c in "=#()[]@/\\+"):
        return False
    if len(val_clean) > 20 or len(val_clean) < 4:
        return False
    return bool(re.match(r"^[A-Z0-9]{6,10}$", val_clean, re.IGNORECASE))


def detect_matrix_axes(df: pd.DataFrame, known_targets: set[str] | None = None) -> tuple[str, str]:
    """Deterministicky určí osu proteinů a osu sloučenin."""
    idx_name = str(df.index.name or "").lower()
    col_name = str(df.columns.name or "").lower()

    if "target" in idx_name or "protein" in idx_name:
        return "targets_in_rows", "compounds_in_cols"
    if "target" in col_name or "protein" in col_name:
        return "targets_in_cols", "compounds_in_rows"
    if "compound" in idx_name or "ligand" in idx_name or "smi" in idx_name:
        return "targets_in_cols", "compounds_in_rows"
    if "compound" in col_name or "ligand" in col_name or "smi" in col_name:
        return "targets_in_rows", "compounds_in_cols"

    idx_targets = sum(1 for x in df.index if is_likely_target_id(str(x), known_targets))
    col_targets = sum(1 for x in df.columns if is_likely_target_id(str(x), known_targets))

    if col_targets > idx_targets:
        return "targets_in_cols", "compounds_in_rows"
    elif idx_targets > col_targets:
        return "targets_in_rows", "compounds_in_cols"

    idx_has_smiles = any(any(c in str(x) for c in "=#()[]@") for x in df.index)
    col_has_smiles = any(any(c in str(x) for c in "=#()[]@") for x in df.columns)
    if idx_has_smiles and not col_has_smiles:
        return "targets_in_cols", "compounds_in_rows"
    if col_has_smiles and not idx_has_smiles:
        return "targets_in_rows", "compounds_in_cols"

    if abs(len(df.columns) - 49) < abs(len(df.index) - 49):
        return "targets_in_cols", "compounds_in_rows"
    return "targets_in_rows", "compounds_in_cols"


def sanitize_compound_label(label: str, max_len: int = 16) -> str:
    """Zkrátí příliš dlouhé popisky sloučenin (např. SMILES řetězce)."""
    s = str(label).strip()
    if len(s) > max_len:
        return s[:max_len - 3] + "..."
    return s


def load_data(input_path: Path, known_targets: set[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame | None, Path]:
    """
    Načte matici skóre a zajistí, že SLOUPCE jsou PROTEINY a ŘÁDKY jsou SLOUČENINY.
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

    matrix_df = matrix_df.apply(pd.to_numeric, errors="coerce")

    # DETERMINISTICKÁ KONTROLA ORIENTACE
    targets_axis, _ = detect_matrix_axes(matrix_df, known_targets)
    if targets_axis == "targets_in_rows":
        matrix_df = matrix_df.T

    # Odfiltrování neproteinových sloupců
    valid_cols = [c for c in matrix_df.columns if is_likely_target_id(str(c), known_targets)]
    if len(valid_cols) >= 3:
        matrix_df = matrix_df[valid_cols]

    matrix_df.columns.name = "target_id"
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
        if len(other_scores) < 1:
            continue

        median_others = float(np.median(other_scores))
        mean_others = float(np.mean(other_scores))

        # Robustní rozptyl (MAD)
        mad = float(np.median(np.abs(other_scores - median_others)))
        robust_scale = max(1.4826 * mad, 0.12)

        # Efekt: delta skóre (záporná hodnota = silnější na cíli)
        delta_score = target_score - median_others

        if len(other_scores) >= 3:
            df_deg = len(other_scores) - 1
            z_score = (target_score - median_others) / robust_scale
            p_val = 2.0 * stats.t.sf(abs(z_score), df=df_deg)
            p_val = max(min(p_val, 1.0), 1e-15)
        else:
            z_score = delta_score / 0.5
            p_val = 0.05 if abs(delta_score) >= delta_cutoff else 0.50

        neg_log10_p = -np.log10(p_val)

        # Klasifikace
        if delta_score <= -delta_cutoff and p_val <= p_cutoff:
            if target_score <= potent_cutoff:
                category = "Vysoce selektivní silný hit"
            else:
                category = "Selektivní (mírná afinita)"
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
        df_res = df_res.sort_values(by=["delta_score", "neg_log10_p"], ascending=[True, False])

    return df_res


def plot_single_volcano(
    df_volcano: pd.DataFrame,
    target_id: str,
    target_display: str,
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

    fig, ax = plt.subplots(figsize=(9.5, 7.5), dpi=dpi)

    log10_p_thresh = -np.log10(p_cutoff)

    palette = {
        "Vysoce selektivní silný hit": "#b2182b",   # Sytá červená
        "Selektivní (mírná afinita)": "#ef8a62",    # Světle oranžová
        "Preferuje ostatní proteiny": "#2166ac",     # Modrá
        "Neselektivní / Nespecifický": "#999999",   # Neutrální šedá
    }

    # Vykreslení bodů
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

    # Označení top hitů textovými popisky (sanitizovaný název, aby nebyl gigantický SMILES)
    top_hits = df_volcano[
        (df_volcano["delta_score"] <= -delta_cutoff) &
        (df_volcano["neg_log10_p"] >= log10_p_thresh)
    ].head(max_labels)

    texts = []
    for _, r in top_hits.iterrows():
        short_cid = sanitize_compound_label(r["compound_id"], max_len=14)
        t = ax.text(
            r["delta_score"],
            r["neg_log10_p"],
            f"{short_cid} ({r['target_score']:.1f})",
            fontsize=8,
            fontweight="bold",
            color="#67001f",
        )
        texts.append(t)

    if texts:
        if ADJUST_TEXT_AVAILABLE:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="->", color="#67001f", lw=0.6))
        else:
            for t in texts:
                t.set_position((t.get_position()[0], t.get_position()[1] + 0.15))

    # Šipky a popisky
    ax.annotate("← Vyšší afinita k tomuto cíli (Selektivní hity)",
                xy=(0.03, 0.96), xycoords="axes fraction",
                fontsize=9.5, fontweight="bold", color="#b2182b")
    ax.annotate("Preferuje jiné methyltransferázy →",
                xy=(0.97, 0.96), xycoords="axes fraction",
                fontsize=9.5, fontweight="bold", color="#2166ac", ha="right")

    # Osy a titulky
    ax.set_xlabel("Rozdíl afinity: ΔScore = Score(cíl) - Medián(ostatní) [kcal/mol]", fontsize=11, fontweight="bold", labelpad=8)
    ax.set_ylabel("Statistická významnost: -log10(p-hodnota)", fontsize=11, fontweight="bold", labelpad=8)
    ax.set_title(f"Target Selectivity Volcano Plot: {target_display}\n(Srovnání vůči panelu ostatních proteinů)",
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
    annotations: dict[str, str],
    out_png: Path,
    out_pdf: Path,
    p_cutoff: float = 0.05,
    delta_cutoff: float = 1.0,
    dpi: int = 300,
):
    """
    Vykreslí přehledový multi-panel grid pro skupinu proteinů.
    """
    n_targets = len(targets)
    if n_targets == 0:
        return

    n_cols = 3 if n_targets >= 3 else n_targets
    n_rows = int(np.ceil(n_targets / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.2 * n_cols, 4.4 * n_rows), dpi=dpi, squeeze=False)
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

        non_sel = df_v[df_v["category"].str.startswith("Neselektivní")]
        sel = df_v[df_v["category"].str.contains("selektivní", case=False)]
        counter = df_v[df_v["category"].str.contains("Preferuje", case=False)]

        ax.scatter(non_sel["delta_score"], non_sel["neg_log10_p"], c="#aaaaaa", alpha=0.4, s=20, label="Nespecifické")
        ax.scatter(counter["delta_score"], counter["neg_log10_p"], c="#2166ac", alpha=0.8, s=35, label="Ostatní")
        ax.scatter(sel["delta_score"], sel["neg_log10_p"], c="#b2182b", alpha=0.9, s=45, edgecolors="black", linewidths=0.5, label="Selektivní hity")

        ax.axhline(log10_p_thresh, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
        ax.axvline(-delta_cutoff, color="crimson", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.axvline(delta_cutoff, color="steelblue", linestyle=":", linewidth=0.8, alpha=0.7)

        # Anotace top 3 hitů (zkrácené)
        top3 = sel.head(3)
        for _, r in top3.iterrows():
            short_id = sanitize_compound_label(r["compound_id"], max_len=10)
            ax.text(r["delta_score"], r["neg_log10_p"] + 0.1, short_id, fontsize=7, fontweight="bold", color="#67001f")

        target_display = annotations.get(target_id, target_id)
        ax.set_title(f"{target_display} (Hity: {len(sel)})", fontsize=10.5, fontweight="bold")
        ax.set_xlabel("ΔScore [kcal/mol]", fontsize=8.5)
        ax.set_ylabel("-log10(p-val)", fontsize=8.5)
        ax.tick_params(labelsize=8)

    # Skrytí nepoužitých sub-plotů
    for idx in range(n_targets, n_rows * n_cols):
        r = idx // n_cols
        c = idx % n_cols
        axes[r][c].set_visible(False)

    fig.suptitle("Souhrnný přehled selektivity ligandů napříč vybranými methyltransferázami", fontsize=13.5, fontweight="bold", y=0.995)
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
  python plot_docking_volcano.py -i reports_docking -t A2DJA1
  python plot_docking_volcano.py -i reports_docking --top 5
  python plot_docking_volcano.py -i reports_docking --top 6 --grid
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

    # Načtení anotací
    annotations = load_protein_annotations()
    known_targets = set(annotations.keys())

    # Načtení dat
    try:
        matrix_df, ranking_df, base_out_dir = load_data(args.input, known_targets=known_targets)
    except Exception as e:
        print(f"[!] Chyba při načítání dat: {e}", file=sys.stderr)
        sys.exit(1)

    out_dir = args.output_dir if args.output_dir else (base_out_dir / "volcano_plots")
    out_dir.mkdir(parents=True, exist_ok=True)

    all_targets = list(matrix_df.columns)
    selected_targets = []

    if args.target:
        for t in args.target.split(","):
            t_clean = t.strip()
            if t_clean in all_targets:
                selected_targets.append(t_clean)
            else:
                print(f"[!] Varování: Cíl '{t_clean}' nebyl nalezen v matici!", file=sys.stderr)
    elif args.top:
        if ranking_df is not None and "target_id" in ranking_df.columns:
            top_ranked = ranking_df["target_id"].head(args.top).tolist()
            selected_targets = [t for t in top_ranked if t in all_targets]
        else:
            selected_targets = all_targets[:args.top]
    elif args.all:
        selected_targets = all_targets
    else:
        if ranking_df is not None and "target_id" in ranking_df.columns:
            selected_targets = ranking_df["target_id"].head(3).tolist()
        else:
            selected_targets = all_targets[:3]

    if not selected_targets:
        print("[!] Nebyly vybrány žádné platné proteiny pro analýzu!", file=sys.stderr)
        sys.exit(1)

    print("=" * 70)
    print(f"Generování Target Selectivity Volcano Plotů pro {len(selected_targets)} vybraných cílů:")
    print(f"  Celkem proteinů v panelu:  {len(all_targets)}")
    print(f"  Celkem testovaných sloučenin: {len(matrix_df.index)}")
    print(f"  Vybrané cíle:              {', '.join(selected_targets[:5])}{'...' if len(selected_targets)>5 else ''}")
    print(f"  Výstupní složka:           {out_dir}")
    print("=" * 70)

    # Samostatné grafy
    for target_id in selected_targets:
        df_volcano = calculate_target_selectivity(
            matrix_df=matrix_df,
            target_id=target_id,
            potent_cutoff=args.potent_cutoff,
            p_cutoff=args.p_cutoff,
            delta_cutoff=args.delta_cutoff,
        )

        csv_path = out_dir / f"volcano_selectivity_{target_id}.csv"
        df_volcano.to_csv(csv_path, index=False)

        png_path = out_dir / f"volcano_selectivity_{target_id}.png"
        pdf_path = out_dir / f"volcano_selectivity_{target_id}.pdf"

        target_display = annotations.get(target_id, target_id)
        plot_single_volcano(
            df_volcano=df_volcano,
            target_id=target_id,
            target_display=target_display,
            out_png=png_path,
            out_pdf=pdf_path,
            p_cutoff=args.p_cutoff,
            delta_cutoff=args.delta_cutoff,
            dpi=args.dpi,
        )

    # Multi-panel přehled
    if args.grid and len(selected_targets) > 1:
        grid_png = out_dir / "volcano_selectivity_grid.png"
        grid_pdf = out_dir / "volcano_selectivity_grid.pdf"
        plot_multi_volcano_grid(
            matrix_df=matrix_df,
            targets=selected_targets,
            annotations=annotations,
            out_png=grid_png,
            out_pdf=grid_pdf,
            p_cutoff=args.p_cutoff,
            delta_cutoff=args.delta_cutoff,
            dpi=args.dpi,
        )

    print(f"\n[✓] Hotovo! Všechny Volcano ploty selektivity uloženy do: '{out_dir}/'")


if __name__ == "__main__":
    main()
