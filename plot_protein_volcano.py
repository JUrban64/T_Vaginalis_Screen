#!/usr/bin/env python3
# ==============================================================================
# Skript pro vizualizaci dokovacích dat: Protein-Level Volcano Plot
# (Globální srovnání všech 49 proteinů napříč celou knihovnou inhibitorů)
# ==============================================================================
# Každý bod v grafu představuje JEDEN CÍLOVÝ PROTEIN (methyltransferázu T. vaginalis).
# ŽÁDNÉ LIGANDY ANI SMILES SE V GRAFU NEZOBRAZUJÍ JAKO BODY ANI POPISKY.
#
# Osa X: Posun afinity vůči celému panelu (Effect Size / Druggability Shift):
#        ΔAffinity = Medián(celý panel) - Průměr(daný protein) [kcal/mol]
#        [Kladné hodnoty = protein váže knihovnu signifikantně silněji než průměr]
#
# Osa Y: Statistická významnost:
#        -log10(FDR / q-hodnota) z Benjamini-Hochberg korekce (nebo Bonferroni / raw-p)
#        z párového t-testu napříč všemi testovanými ligandy
#
# Velikost bodu: Počet silných vazačů (potent binders, např. skóre <= -8.0 kcal/mol)
#
# Popisky bodů: UniProt ID + genový symbol / název z reports/deduplicated_proteins.tsv
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
    matplotlib.use("Agg")  # Bezhlavý backend pro HPC / servery
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


def multipletests_adjust(p_values: list[float] | np.ndarray, method: str = "fdr") -> np.ndarray:
    """
    Korekce p-hodnot na vícenásobné testování hypotéz (NumPy kompatibilní, bez nutnosti statsmodels).
    Podporované metody:
      - 'fdr' / 'fdr_bh': Benjamini-Hochberg False Discovery Rate (q-value)
      - 'bonferroni': Bonferroni Family-Wise Error Rate
      - 'none': bez korekce (původní p-hodnoty)
    """
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0 or method == "none":
        return p.copy()

    if method == "bonferroni":
        return np.clip(p * n, 0.0, 1.0)

    if method in ("fdr", "fdr_bh", "bh"):
        order = np.argsort(p)
        ranks = np.arange(1, n + 1)
        p_sorted = p[order]
        # Benjamini-Hochberg: p_adj = p * n / rank
        p_adj_sorted = p_sorted * (n / ranks)
        # Kumulativní minimum odzadu pro zachování monotónnosti
        p_adj_sorted = np.minimum.accumulate(p_adj_sorted[::-1])[::-1]
        p_adj_sorted = np.clip(p_adj_sorted, 0.0, 1.0)
        p_adj = np.empty_like(p_adj_sorted)
        p_adj[order] = p_adj_sorted
        return p_adj

    return p.copy()


def load_protein_annotations(search_dir: Path | None = None) -> dict[str, dict[str, str]]:
    """
    Načte biologické anotace pro 49 proteinů z deduplicated_proteins.tsv nebo jiných reportů.
    Vrací slovník: target_id -> {
        'gene': 'TVAG_...',
        'short_name': '...',
        'display_label': 'A2FE15 (ICMT)' nebo 'A2FE15 (TVAG_470080)'
    }
    """
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
                if c.suffix == ".tsv":
                    df = pd.read_csv(c, sep="\t")
                elif c.suffix == ".csv":
                    df = pd.read_csv(c)
                else:
                    continue

                id_col = "Entry" if "Entry" in df.columns else ("target_id" if "target_id" in df.columns else df.columns[0])
                gene_col = "Gene Names" if "Gene Names" in df.columns else None
                entry_name_col = "Entry Name" if "Entry Name" in df.columns else None
                name_col = "Protein names" if "Protein names" in df.columns else None

                for _, row in df.iterrows():
                    tid = str(row[id_col]).strip()
                    if not tid or tid == "nan":
                        continue

                    gene = ""
                    if gene_col and pd.notna(row[gene_col]):
                        first_gene = str(row[gene_col]).strip().split()[0]
                        if first_gene and first_gene.lower() != "nan":
                            gene = first_gene

                    short_name = ""
                    if entry_name_col and pd.notna(row[entry_name_col]):
                        ename = str(row[entry_name_col]).strip().split("_")[0]
                        if ename and ename != tid:
                            short_name = ename

                    if not short_name and name_col and pd.notna(row[name_col]):
                        pname = str(row[name_col]).strip().split("(")[0].strip()
                        if len(pname) > 22:
                            pname = pname[:20] + ".."
                        short_name = pname

                    best_label = short_name if (short_name and not short_name.startswith("TVAG")) else gene
                    if not best_label:
                        best_label = gene or short_name

                    display_label = f"{tid} ({best_label})" if best_label else tid

                    annots[tid] = {
                        "gene": gene,
                        "short_name": short_name,
                        "display_label": display_label
                    }

                if annots:
                    break
            except Exception:
                pass

    return annots


def is_likely_target_id(val: str, known_targets: set[str] | None = None) -> bool:
    """Ověří, zda je daný řetězec skutečně protein (UniProt ID) a ne ligand / SMILES."""
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


def load_matrix(input_path: Path, score_type: str = "adjusted", known_targets: set[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame | None, Path, str]:
    """
    Načte matici skóre a zajistí, že SLOUPCE jsou PROTEINY a ŘÁDKY jsou LIGANDY.
    Tím je 100% zaručeno, že Volcano plot bude vyhodnocovat proteiny, nikoli SMILES!

    Podporuje score_type:
      - 'adjusted': Prioritně načte 'docking_mw_adjusted_matrix.csv' (skóre očištěné od MW na škále kcal/mol).
      - 'raw':      Načte 'docking_scores_matrix.csv' (surová Vina skóre).
      - 'residual': Načte 'docking_mw_residuals_matrix.csv' (reziduální skóre kolem 0).
    """
    used_type = score_type
    if input_path.is_dir():
        score_file = None
        if score_type == "adjusted":
            cand = input_path / "docking_mw_adjusted_matrix.csv"
            if cand.exists():
                score_file = cand
                used_type = "MW-Adjusted (kcal/mol, bez vlivu hmotnosti)"
            else:
                score_file = input_path / "docking_scores_matrix.csv"
                used_type = "Raw (fallback)"
        elif score_type == "residual":
            score_file = input_path / "docking_mw_residuals_matrix.csv"
            used_type = "MW-Residuals (kolem 0.0)"
        else:
            score_file = input_path / "docking_scores_matrix.csv"
            used_type = "Raw (surová skóre)"

        if not score_file or not score_file.exists():
            candidates = list(input_path.glob("*matrix*.csv"))
            if not candidates:
                raise FileNotFoundError(f"V adresáři '{input_path}' nebyl nalezen 'docking_scores_matrix.csv'!")
            score_file = candidates[0]
            used_type = score_file.stem
        matrix_df = pd.read_csv(score_file, index_col=0)
        ranking_file = input_path / "protein_ranking_summary.csv"
        ranking_df = pd.read_csv(ranking_file) if ranking_file.exists() else None
        output_dir = input_path
    else:
        if not input_path.exists():
            raise FileNotFoundError(f"Soubor '{input_path}' neexistuje!")
        matrix_df = pd.read_csv(input_path, index_col=0)
        score_file = input_path
        used_type = input_path.stem
        ranking_file = input_path.parent / "protein_ranking_summary.csv"
        ranking_df = pd.read_csv(ranking_file) if ranking_file.exists() else None
        output_dir = input_path.parent

    matrix_df = matrix_df.apply(pd.to_numeric, errors="coerce")

    # DETERMINISTICKÁ ORIENTACE:
    targets_axis, _ = detect_matrix_axes(matrix_df, known_targets)
    if targets_axis == "targets_in_rows":
        matrix_df = matrix_df.T

    # Striktní filtrace neproteinových sloupců
    valid_protein_cols = [c for c in matrix_df.columns if is_likely_target_id(str(c), known_targets)]
    if len(valid_protein_cols) >= 3:
        matrix_df = matrix_df[valid_protein_cols]

    matrix_df.columns.name = "target_id"
    matrix_df.index.name = "compound_id"

    return matrix_df, ranking_df, output_dir, used_type


def compute_protein_volcano_statistics(
    matrix_df: pd.DataFrame,
    annotations: dict[str, dict[str, str]],
    p_adjust: str = "fdr",
    potent_cutoff: float = -8.0,
    p_cutoff: float = 0.05,
    delta_cutoff: float = 0.5,
) -> pd.DataFrame:
    """
    Spočítá statistiky pro každý PROTEIN napříč testovanými sloučeninami:
    - Osa X: ΔAffinity = Posun afinity vůči mediánu celého panelu [kcal/mol]
    - Osa Y: -log10(p-hodnota) s podporou FDR (Benjamini-Hochberg) i Bonferroni
    - Velikost bodu: počet silných vazačů (potent binders <= potent_cutoff)
    """
    targets = list(matrix_df.columns)
    raw_stats = []

    for target_id in targets:
        tid_str = str(target_id).strip()
        target_series = matrix_df[target_id].dropna()
        if len(target_series) < 1:
            continue

        other_targets = [t for t in targets if t != target_id]

        compounds_idx = target_series.index
        target_scores = target_series.values

        if len(other_targets) > 0 and len(compounds_idx) > 0:
            other_medians = matrix_df.loc[compounds_idx, other_targets].median(axis=1).values
        else:
            other_medians = target_scores

        diffs = other_medians - target_scores
        mean_diff = float(np.mean(diffs)) if len(diffs) > 0 else 0.0
        mean_score = float(np.mean(target_scores))
        best_score = float(np.min(target_scores))
        potent_binders = int((target_scores <= potent_cutoff).sum())

        # Párový test napříč ligandy
        if len(target_scores) >= 3 and np.std(diffs) > 1e-6:
            t_stat, p_val = stats.ttest_rel(target_scores, other_medians, nan_policy="omit")
            p_val = max(min(float(p_val), 1.0), 1e-15)
        else:
            t_stat = mean_diff / (np.std(diffs) + 1e-3) if len(diffs) > 1 else 0.0
            p_val = 0.05 if abs(mean_diff) >= delta_cutoff else 0.50

        annot_info = annotations.get(tid_str, {})
        gene_name = annot_info.get("gene", "")
        short_name = annot_info.get("short_name", "")
        display_label = annot_info.get("display_label", tid_str)

        raw_stats.append({
            "target_id": tid_str,
            "display_label": display_label,
            "gene_name": gene_name,
            "protein_name": short_name,
            "mean_docking_score": round(mean_score, 2),
            "best_docking_score": round(best_score, 2),
            "potent_binders_count": potent_binders,
            "delta_affinity": round(mean_diff, 2),
            "t_statistic": round(float(t_stat), 2),
            "raw_p_value": p_val,
        })

    if not raw_stats:
        return pd.DataFrame()

    df_res = pd.DataFrame(raw_stats)

    # VÝPOČET KOREKCÍ VÍCENÁSOBNÉHO TESTOVÁNÍ (FDR Benjamini-Hochberg i Bonferroni)
    p_values_array = df_res["raw_p_value"].values
    df_res["p_adj_fdr"] = np.round(multipletests_adjust(p_values_array, method="fdr"), 6)
    df_res["p_adj_bonferroni"] = np.round(multipletests_adjust(p_values_array, method="bonferroni"), 6)
    df_res["p_value"] = np.round(df_res["raw_p_value"], 6)

    # Zvolení aktivní p-hodnoty dle požadavku uživatele
    if p_adjust == "fdr":
        active_p = df_res["p_adj_fdr"].values
    elif p_adjust == "bonferroni":
        active_p = df_res["p_adj_bonferroni"].values
    else:  # 'none'
        active_p = df_res["p_value"].values

    # -log10 pro osu Y (ochrana před nulou)
    active_p_safe = np.clip(active_p, 1e-15, 1.0)
    df_res["neg_log10_p"] = np.round(-np.log10(active_p_safe), 2)

    # Kategorizace proteinů na základě korigované p-hodnoty
    categories = []
    for diff, p_act in zip(df_res["delta_affinity"], active_p):
        if diff >= delta_cutoff and p_act <= p_cutoff:
            categories.append("Vysoce citlivý cíl (Super-binder)")
        elif diff <= -delta_cutoff and p_act <= p_cutoff:
            categories.append("Rekalcitrantní cíl (Slabá afinita)")
        else:
            categories.append("Průměrný protein panelu")

    df_res["category"] = categories

    # Seřazení: největší posun afinity a nejvyšší statistická významnost nahoře
    cols_order = [
        "target_id", "display_label", "gene_name", "protein_name",
        "delta_affinity", "mean_docking_score", "best_docking_score", "potent_binders_count",
        "p_value", "p_adj_fdr", "p_adj_bonferroni", "neg_log10_p", "t_statistic", "category"
    ]
    df_res = df_res[cols_order].sort_values(by=["delta_affinity", "neg_log10_p"], ascending=[False, False])

    return df_res


def plot_protein_volcano_figure(
    df_stats: pd.DataFrame,
    title: str,
    n_compounds: int,
    out_png: Path,
    out_pdf: Path,
    p_adjust: str = "fdr",
    p_cutoff: float = 0.05,
    delta_cutoff: float = 0.5,
    dpi: int = 300,
):
    """
    Vykreslí přehledný, publikační Protein-Level Volcano plot.
    Každý bod = JEDEN PROTEIN s UniProt ID a genovým symbolem.
    """
    sns.set_theme(style="whitegrid", font="sans-serif")
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]

    fig, ax = plt.subplots(figsize=(11.0, 8.5), dpi=dpi)

    log10_p_thresh = -np.log10(p_cutoff)

    palette = {
        "Vysoce citlivý cíl (Super-binder)": "#d73027",  # Sytá červená
        "Rekalcitrantní cíl (Slabá afinita)": "#313695",  # Tmavě modrá
        "Průměrný protein panelu": "#999999",            # Neutrální šedá
    }

    # Škálování velikosti bodů podle počtu silných vazačů (potent binders)
    min_size = 55
    max_size = 280
    max_binders = max(df_stats["potent_binders_count"].max(), 1)
    sizes = min_size + (df_stats["potent_binders_count"] / max_binders) * (max_size - min_size)

    # Limity os pro barevné kvadranty
    x_min = min(df_stats["delta_affinity"].min() - 0.5, -delta_cutoff - 1.0)
    x_max = max(df_stats["delta_affinity"].max() + 0.5, delta_cutoff + 1.0)
    y_max = max(df_stats["neg_log10_p"].max() + 0.8, log10_p_thresh + 2.0)

    # Jemné podbarvení kvadrantů významnosti
    # 1. Pravý horní kvadrant: Vysoce citlivé cíle (Super-binders)
    ax.axvspan(delta_cutoff, x_max + 1.0, ymin=(log10_p_thresh / y_max), ymax=1.0,
               color="#fee8e7", alpha=0.45, zorder=0)
    # 2. Levý horní kvadrant: Rekalcitrantní enzymy
    ax.axvspan(x_min - 1.0, -delta_cutoff, ymin=(log10_p_thresh / y_max), ymax=1.0,
               color="#e7f0fa", alpha=0.45, zorder=0)

    # Vykreslení bodů
    for cat in ["Průměrný protein panelu", "Rekalcitrantní cíl (Slabá afinita)", "Vysoce citlivý cíl (Super-binder)"]:
        sub = df_stats[df_stats["category"] == cat]
        if sub.empty:
            continue
        sub_sizes = sizes.loc[sub.index]
        alpha = 0.55 if cat.startswith("Průměrný") else 0.92
        zorder = 2 if cat.startswith("Průměrný") else 4

        ax.scatter(
            sub["delta_affinity"],
            sub["neg_log10_p"],
            c=palette[cat],
            label=f"{cat} (n={len(sub)})",
            s=sub_sizes,
            alpha=alpha,
            edgecolors="black" if not cat.startswith("Průměrný") else "gray",
            linewidths=0.7,
            zorder=zorder,
        )

    # Prahové čáry
    ax.axhline(log10_p_thresh, color="black", linestyle="--", linewidth=1.1, alpha=0.75, zorder=1)
    ax.axvline(delta_cutoff, color="crimson", linestyle=":", linewidth=1.1, alpha=0.85, zorder=1)
    ax.axvline(-delta_cutoff, color="steelblue", linestyle=":", linewidth=1.1, alpha=0.85, zorder=1)

    # Textové označení prahů
    thresh_name = "FDR" if p_adjust == "fdr" else ("p-adj" if p_adjust == "bonferroni" else "p")
    ax.text(x_max - 0.1, log10_p_thresh + 0.08, f"{thresh_name} = {p_cutoff} (-log10 = {log10_p_thresh:.1f})",
            fontsize=8.5, color="black", ha="right", va="bottom", fontstyle="italic", zorder=3)

    # Textové popisky významných PROTEINŮ (zobrazuje se display_label, NIKDY SMILES!)
    sig_targets = df_stats[
        (df_stats["neg_log10_p"] >= log10_p_thresh) &
        (df_stats["delta_affinity"].abs() >= delta_cutoff)
    ]

    texts = []
    for _, r in sig_targets.iterrows():
        lbl = str(r["display_label"])
        if any(c in lbl for c in "=#()[]@") and len(lbl) > 25:
            lbl = str(r["target_id"])

        t = ax.text(
            r["delta_affinity"],
            r["neg_log10_p"],
            f" {lbl}",
            fontsize=8.5,
            fontweight="bold",
            color="#67001f" if r["delta_affinity"] > 0 else "#08306b",
            zorder=5,
        )
        texts.append(t)

    if texts:
        if ADJUST_TEXT_AVAILABLE:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="->", color="dimgray", lw=0.6))
        else:
            for i, t in enumerate(texts):
                offset_y = 0.12 if i % 2 == 0 else -0.12
                t.set_position((t.get_position()[0], t.get_position()[1] + offset_y))

    # Směrové šipky a vysvětlující texty v záhlaví
    ax.annotate("Váže testované ligandy významně SILNĚJI než průměr →\n(Primární cíle / Druggable Hotspots)",
                xy=(0.98, 0.96), xycoords="axes fraction",
                fontsize=9.5, fontweight="bold", color="#d73027", ha="right")
    ax.annotate("← Váže ligandy SLABĚJI než průměr\n(Rekalcitrantní enzymy)",
                xy=(0.02, 0.96), xycoords="axes fraction",
                fontsize=9.5, fontweight="bold", color="#313695", ha="left")

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-0.2, y_max)

    # Popisky os
    ax.set_xlabel("Posun afinity: ΔAffinity = Medián(celý panel) - Průměr(protein) [kcal/mol]\n(kladné hodnoty = silnější vazba proteinu vůči zbytku rodiny)",
                  fontsize=11, fontweight="bold", labelpad=8)

    if p_adjust == "fdr":
        y_axis_title = f"Statistická významnost: -log10(FDR / q-hodnota)\n(Benjamini-Hochberg korekce přes 49 proteinů, testováno přes {n_compounds} sloučenin)"
    elif p_adjust == "bonferroni":
        y_axis_title = f"Statistická významnost: -log10(p-adj)\n(Bonferroniho korekce přes 49 proteinů, testováno přes {n_compounds} sloučenin)"
    else:
        y_axis_title = f"Statistická významnost: -log10(p-hodnota)\n(párový t-test přes {n_compounds} testovaných sloučenin, bez korekce)"

    ax.set_ylabel(y_axis_title, fontsize=11, fontweight="bold", labelpad=8)
    ax.set_title(title, fontsize=13.5, fontweight="bold", pad=16)

    # Legenda
    leg = ax.legend(frameon=True, facecolor="white", edgecolor="lightgray", fontsize=9, loc="upper center", bbox_to_anchor=(0.5, 0.91))
    leg.set_zorder(10)

    # Informativní patička
    note_txt = f"Poznámka: Velikost bodu odpovídá počtu vysoce silných inhibitorů (dokovací skóre <= -8.0 kcal/mol) | Korekce: {p_adjust.upper()}"
    ax.text(0.02, 0.02, note_txt, transform=ax.transAxes, fontsize=8, fontstyle="italic", color="#555555")

    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close("all")

    print(f"  [+] Protein-Level Volcano plot uložen (PNG): {out_png}")
    print(f"  [+] Protein-Level Volcano plot uložen (PDF): {out_pdf}")


def main():
    parser = argparse.ArgumentParser(
        description="Generování Protein-Level Volcano Plotu (Srovnání 49 proteinů napříč ligandy).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Příklady použití:
  python plot_protein_volcano.py -i reports_docking
  python plot_protein_volcano.py -i reports_docking --p-adjust fdr --p-cutoff 0.05
  python plot_protein_volcano.py -i reports_docking --p-adjust bonferroni
  python plot_protein_volcano.py -i reports_docking --p-adjust none
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
        "--score-type",
        choices=["adjusted", "raw", "residual"],
        default="adjusted",
        help="Typ dokovacího skóre: 'adjusted' (výchozí: očištěné od hmotnosti na škále kcal/mol), 'raw' (surová skóre), nebo 'residual' (čistá rezidua kolem 0)."
    )
    parser.add_argument(
        "--p-adjust",
        choices=["fdr", "bonferroni", "none"],
        default="fdr",
        help="Metoda korekce na vícenásobné testování: 'fdr' (Benjamini-Hochberg, výchozí), 'bonferroni', nebo 'none' (surové p-hodnoty)."
    )
    parser.add_argument(
        "--p-cutoff",
        type=float,
        default=0.05,
        help="Prahová hodnota pro statistickou významnost (výchozí: 0.05 pro FDR/p-value)."
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

    # 1. Načtení biologických anotací proteinů (49 methyltransferáz)
    annotations = load_protein_annotations()
    known_targets = set(annotations.keys())

    # 2. Načtení matice skóre (garantovaně: sloupce = proteiny, řádky = sloučeniny)
    try:
        matrix_df, ranking_df, base_out_dir, used_score_type = load_matrix(
            args.input, score_type=args.score_type, known_targets=known_targets
        )
    except Exception as e:
        print(f"[!] Chyba při načítání matice skóre: {e}", file=sys.stderr)
        sys.exit(1)

    out_dir = args.output_dir if args.output_dir else (base_out_dir / "protein_volcano")
    out_dir.mkdir(parents=True, exist_ok=True)

    n_proteins = matrix_df.shape[1]
    n_compounds = matrix_df.shape[0]

    method_desc = {
        "fdr": "FDR (Benjamini-Hochberg / q-hodnota)",
        "bonferroni": "Bonferroni (p-adj)",
        "none": "Nekorigovaná surová p-hodnota"
    }[args.p_adjust]

    print("=" * 70)
    print("Generování Protein-Level Volcano Plotu (Srovnání proteinových cílů):")
    print(f"  Celkem vyhodnocených PROTEINŮ: {n_proteins}")
    print(f"  Celkem testovaných sloučenin:  {n_compounds}")
    print(f"  Použitá matice skóre:          {used_score_type}")
    print(f"  Korekce vícenásobného testu:   {method_desc}")
    print(f"  Prahová významnost ({args.p_adjust.upper()}):     <= {args.p_cutoff} (-log10 = {-np.log10(args.p_cutoff):.2f})")
    print(f"  Prahový posun afinity:         |ΔAffinity| >= {args.delta_cutoff} kcal/mol")
    print(f"  Výstupní složka:               {out_dir}")
    print("=" * 70)

    # 3. Výpočet statistik pro každý PROTEIN s FDR korekcí
    df_stats = compute_protein_volcano_statistics(
        matrix_df=matrix_df,
        annotations=annotations,
        p_adjust=args.p_adjust,
        potent_cutoff=args.potent_cutoff,
        p_cutoff=args.p_cutoff,
        delta_cutoff=args.delta_cutoff,
    )

    if df_stats.empty:
        print("[!] Nepodařilo se spočítat statistiky pro žádný protein!", file=sys.stderr)
        sys.exit(1)

    # Uložení kompletní tabulky statistik proteinů
    csv_path = out_dir / "protein_volcano_statistics.csv"
    df_stats.to_csv(csv_path, index=False)
    print(f"[+] Tabulka statistik proteinů uložena: {csv_path}")

    # 4. Vykreslení Volcano grafu
    png_path = out_dir / "protein_volcano_plot.png"
    pdf_path = out_dir / "protein_volcano_plot.pdf"

    title = f"Protein-Level Volcano Plot: Srovnání afinity {n_proteins} methyltransferáz T. vaginalis\n(Párová statistika přes {n_compounds} sloučenin, {used_score_type}, {method_desc})"

    plot_protein_volcano_figure(
        df_stats=df_stats,
        title=title,
        n_compounds=n_compounds,
        out_png=png_path,
        out_pdf=pdf_path,
        p_adjust=args.p_adjust,
        p_cutoff=args.p_cutoff,
        delta_cutoff=args.delta_cutoff,
        dpi=args.dpi,
    )

    # 5. Výpis TOP cílů v terminálu
    top_super = df_stats[df_stats["category"].str.contains("Super-binder")].head(5)
    if not top_super.empty:
        print("\n" + "=" * 70)
        print(f"TOP VYSOCE CITLIVÉ CÍLE KNIHOVNY (Super-binders / Druggable Hotspots, {args.p_adjust.upper()} <= {args.p_cutoff}):")
        cols_print = ["target_id", "display_label", "delta_affinity", "mean_docking_score", "p_value", "p_adj_fdr", "potent_binders_count"]
        cols_print = [c for c in cols_print if c in top_super.columns]
        print(top_super[cols_print].to_string(index=False))
        print("=" * 70)

    print("\n[✓] Hotovo! Protein-Level Volcano plot byl úspěšně vygenerován.")


if __name__ == "__main__":
    main()
