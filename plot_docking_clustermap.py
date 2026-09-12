#!/usr/bin/env python3
# ==============================================================================
# Skript pro vizualizaci dokovacích dat: Hierarchická Heatmapa (Clustermap)
# ==============================================================================
# Vytvoří oboustranně shlukovanou heatmapu (Clustermap) interakční matice:
#   Proteiny (49 methyltransferáz) x Ligandy (knihovna inhibitorů)
#
# Klíčové vlastnosti:
#   1. Robustní automatická detekce os: spolehlivě pozná proteiny (UniProt ID)
#      a ligandy bez ohledu na velikost matice (např. i při testu s 1 ligandem).
#   2. Možnost volby orientace: --rows {proteins,ligands} (výchozí: proteiny na ose Y).
#   3. Správné a konzistentní popisky os (nikdy neobrácené).
#   4. Anotace proteinů genovými symboly z reports/deduplicated_proteins.tsv.
#   5. Zkrácení dlouhých SMILES popisků, aby nerozbíjely layout.
#   6. Podpora surových skóre (kcal/mol), Z-skóre i MW-reziduálních skóre.
#   7. Export v publikační kvalitě (PNG 300 DPI + vektorové PDF).
#
# Autor: Antigravity AI / Jáchym Urban
# ==============================================================================

import argparse
import re
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


def load_target_annotations(search_dir: Path | None = None) -> dict[str, str]:
    """
    Načte anotace proteinových cílů z 'deduplicated_proteins.tsv' nebo 'targets.txt'.
    Vrací slovník: target_id -> 'target_id (Gene / Name)' nebo 'target_id'.
    """
    candidates = [
        Path("reports/deduplicated_proteins.tsv"),
        Path("reports_deduplicated/download_summary.csv"),
        Path("targets.txt"),
    ]
    if search_dir:
        candidates.insert(0, search_dir / "reports/deduplicated_proteins.tsv")
        candidates.insert(1, search_dir.parent / "reports/deduplicated_proteins.tsv")

    mapping = {}

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

                    label = ""
                    # 1. Zkusíme genový symbol (např. TVAG_069970 nebo TRM5)
                    if gene_col and pd.notna(row[gene_col]):
                        first_gene = str(row[gene_col]).strip().split()[0]
                        if first_gene and first_gene.lower() != "nan":
                            label = first_gene

                    # 2. Pokud není gen, zkusíme Entry Name (např. TRM5_TRIV3 -> TRM5)
                    if not label and entry_name_col and pd.notna(row[entry_name_col]):
                        ename = str(row[entry_name_col]).strip().split("_")[0]
                        if ename and ename != tid:
                            label = ename

                    # 3. Zkusíme zkrácený název proteinu
                    if not label and name_col and pd.notna(row[name_col]):
                        pname = str(row[name_col]).strip().split("(")[0].strip()
                        if len(pname) > 20:
                            pname = pname[:18] + ".."
                        label = pname

                    if label:
                        mapping[tid] = f"{tid} ({label})"
                    else:
                        mapping[tid] = tid

                if mapping:
                    break
            except Exception:
                pass

    return mapping


def is_likely_target_id(val: str, known_targets: set[str] | None = None) -> bool:
    """Ověří, zda daný řetězec odpovídá proteinovému cíli (UniProt formát)."""
    if not isinstance(val, str):
        return False
    val_clean = val.strip()
    if known_targets and val_clean in known_targets:
        return True
    # Pokud obsahuje chemické znaky ze SMILES, není to protein
    if any(c in val_clean for c in "=#()[]@/\\+"):
        return False
    if len(val_clean) > 20 or len(val_clean) < 4:
        return False
    # UniProt formát (např. A2FE15, A0A7D6J3B8)
    return bool(re.match(r"^[A-Z0-9]{6,10}$", val_clean, re.IGNORECASE))


def detect_matrix_axes(df: pd.DataFrame, known_targets: set[str] | None = None) -> tuple[str, str]:
    """
    Deterministicky určí, která osa obsahuje proteiny a která sloučeniny.
    Vrací ('targets_in_rows', 'compounds_in_cols') nebo ('targets_in_cols', 'compounds_in_rows').
    """
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

    # Kontrola hodnot
    idx_target_matches = sum(1 for x in df.index if is_likely_target_id(str(x), known_targets))
    col_target_matches = sum(1 for x in df.columns if is_likely_target_id(str(x), known_targets))

    idx_ratio = idx_target_matches / max(len(df.index), 1)
    col_ratio = col_target_matches / max(len(df.columns), 1)

    if col_ratio > idx_ratio and col_target_matches > 0:
        return "targets_in_cols", "compounds_in_rows"
    elif idx_ratio > col_ratio and idx_target_matches > 0:
        return "targets_in_rows", "compounds_in_cols"

    # Kontrola SMILES znaků
    idx_has_smiles = any(any(c in str(x) for c in "=#()[]@") for x in df.index)
    col_has_smiles = any(any(c in str(x) for c in "=#()[]@") for x in df.columns)
    if idx_has_smiles and not col_has_smiles:
        return "targets_in_cols", "compounds_in_rows"
    if col_has_smiles and not idx_has_smiles:
        return "targets_in_rows", "compounds_in_cols"

    # Default podle 49 známých proteinů
    if abs(len(df.columns) - 49) < abs(len(df.index) - 49):
        return "targets_in_cols", "compounds_in_rows"
    return "targets_in_rows", "compounds_in_cols"


def sanitize_compound_label(label: str, max_len: int = 16) -> str:
    """Zkrátí příliš dlouhé popisky sloučenin (např. SMILES řetězce)."""
    s = str(label).strip()
    if len(s) > max_len:
        return s[:max_len - 3] + "..."
    return s


def load_matrix(input_path: Path, score_type: str = "adjusted") -> tuple[pd.DataFrame, Path, str]:
    """Načte matici skóre z CSV nebo adresáře reportů s podporou MW-adjusted matice."""
    used_type = score_type
    if input_path.is_dir():
        score_file = None
        if score_type == "adjusted":
            cand = input_path / "docking_mw_adjusted_matrix.csv"
            if cand.exists():
                score_file = cand
                used_type = "MW-Adjusted (kcal/mol)"
            else:
                score_file = input_path / "docking_scores_matrix.csv"
                used_type = "Raw (fallback)"
        elif score_type == "residual":
            score_file = input_path / "docking_mw_residuals_matrix.csv"
            used_type = "MW-Residuals"
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
        output_dir = input_path
    else:
        if not input_path.exists():
            raise FileNotFoundError(f"Soubor '{input_path}' neexistuje!")
        matrix_df = pd.read_csv(input_path, index_col=0)
        output_dir = input_path.parent
        used_type = input_path.stem

    matrix_df = matrix_df.apply(pd.to_numeric, errors="coerce")
    return matrix_df, output_dir, used_type


def prepare_matrix_orientation(
    df: pd.DataFrame,
    desired_rows: str = "proteins",
    annotations: dict[str, str] | None = None
) -> tuple[pd.DataFrame, str, str]:
    """
    Zajistí přesnou orientaci matice dle volby uživatele:
    desired_rows = 'proteins' -> řádky = proteiny, sloupce = ligandy
    desired_rows = 'ligands'  -> řádky = ligandy, sloupce = proteiny

    Vrací (df_oriented, row_entity, col_entity).
    """
    known_targets = set(annotations.keys()) if annotations else None
    current_targets_axis, _ = detect_matrix_axes(df, known_targets)

    if desired_rows == "proteins":
        if current_targets_axis == "targets_in_cols":
            # Proteiny jsou ve sloupcích -> transponujeme do řádků
            df_oriented = df.T
        else:
            df_oriented = df.copy()
        row_entity = "proteins"
        col_entity = "ligands"
    else:  # desired_rows == "ligands"
        if current_targets_axis == "targets_in_rows":
            # Proteiny jsou v řádcích -> transponujeme do sloupců
            df_oriented = df.T
        else:
            df_oriented = df.copy()
        row_entity = "ligands"
        col_entity = "proteins"

    # Aplikace anotací a sanitizace popisků
    if annotations:
        if row_entity == "proteins":
            df_oriented.index = [annotations.get(str(x), str(x)) for x in df_oriented.index]
            df_oriented.columns = [sanitize_compound_label(str(c)) for c in df_oriented.columns]
        else:
            df_oriented.columns = [annotations.get(str(x), str(x)) for x in df_oriented.columns]
            df_oriented.index = [sanitize_compound_label(str(c)) for c in df_oriented.index]
    else:
        if row_entity == "proteins":
            df_oriented.columns = [sanitize_compound_label(str(c)) for c in df_oriented.columns]
        else:
            df_oriented.index = [sanitize_compound_label(str(c)) for c in df_oriented.index]

    df_oriented.index.name = "Target Protein" if row_entity == "proteins" else "Ligand / Compound"
    df_oriented.columns.name = "Ligand / Compound" if col_entity == "ligands" else "Target Protein"

    return df_oriented, row_entity, col_entity


def impute_missing(df: pd.DataFrame) -> pd.DataFrame:
    """Ošetří prázdné buňky (NaN) nahrazením nejméně stabilním skóre."""
    if df.isna().any().any():
        worst_score = df.max().max() + 0.5
        print(f"[*] Varování: Matice obsahuje {df.isna().sum().sum()} prázdných hodnot (NaN).")
        print(f"    Chybějící buňky budou nahrazeny skóre: {worst_score:.2f} kcal/mol.")
        return df.fillna(worst_score)
    return df


def plot_clustermap_figure(
    df: pd.DataFrame,
    title: str,
    cbar_label: str,
    out_png: Path,
    out_pdf: Path,
    row_entity: str = "proteins",
    col_entity: str = "ligands",
    cmap: str = "viridis_r",
    standardize: str = "none",
    metric: str = "euclidean",
    method: str = "average",
    figsize: tuple[float, float] = (14.0, 10.5),
    dpi: int = 300,
    show_dendrograms: bool = True,
):
    """
    Vykreslí a uloží hierarchickou clustermapu pomocí Seaborn s garantovaně správnými osami.
    """
    if metric != "euclidean" and method == "ward":
        method = "average"

    # Z-skóre normalizace:
    # 'compound': normalizace po sloučeninách (ať už jsou v řádcích nebo sloupcích)
    # 'target': normalizace po cílech
    z_score_val = None
    if standardize == "compound":
        z_score_val = 0 if row_entity == "ligands" else 1
    elif standardize == "target":
        z_score_val = 0 if row_entity == "proteins" else 1

    # Seaborn clustermap vyžaduje alespoň 2 prvky v každé dimenzi pro shlukování
    can_cluster_rows = show_dendrograms and (df.shape[0] >= 2)
    can_cluster_cols = show_dendrograms and (df.shape[1] >= 2)

    sns.set_theme(style="white", font="sans-serif")
    plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]

    dendro_ratio = (
        0.12 if can_cluster_rows else 0.001,
        0.12 if can_cluster_cols else 0.001
    )

    g = sns.clustermap(
        df,
        figsize=figsize,
        cmap=cmap,
        metric=metric,
        method=method,
        z_score=z_score_val,
        dendrogram_ratio=dendro_ratio,
        cbar_pos=(0.02, 0.80, 0.03, 0.15),  # Čisté umístění colorbaru vlevo nahoře
        linewidths=0.2,
        linecolor="whitesmoke",
        row_cluster=can_cluster_rows,
        col_cluster=can_cluster_cols,
        xticklabels=True if df.shape[1] <= 80 else False,
        yticklabels=True if df.shape[0] <= 80 else False,
    )

    # Popisky os DYNAMICKY a PŘESNĚ podle toho, co je na které ose!
    n_rows, n_cols = df.shape
    if row_entity == "proteins":
        y_label = f"Cílové proteiny ({n_rows} methyltransferáz)"
        x_label = f"Knihovna inhibitorů ({n_cols} sloučenin)"
        g.ax_heatmap.tick_params(axis="y", labelsize=8, rotation=0)
        if n_cols <= 80:
            g.ax_heatmap.tick_params(axis="x", labelsize=7, rotation=90)
    else:
        y_label = f"Knihovna inhibitorů ({n_rows} sloučenin)"
        x_label = f"Cílové proteiny ({n_cols} methyltransferáz)"
        g.ax_heatmap.tick_params(axis="y", labelsize=7, rotation=0)
        if n_cols <= 80:
            g.ax_heatmap.tick_params(axis="x", labelsize=8, rotation=90)

    g.ax_heatmap.set_ylabel(y_label, fontsize=12, fontweight="bold", labelpad=12)
    g.ax_heatmap.set_xlabel(x_label, fontsize=12, fontweight="bold", labelpad=12)

    # Colorbar
    g.ax_cbar.set_title(cbar_label, fontsize=9, pad=10, fontweight="bold")
    g.ax_cbar.tick_params(labelsize=8)

    # Titulek
    g.fig.suptitle(title, fontsize=14, fontweight="bold", y=0.99)

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
  python plot_docking_clustermap.py -i reports_docking --both
  python plot_docking_clustermap.py -i reports_docking --standardize compound
  python plot_docking_clustermap.py -i reports_docking --rows ligands  # otočení os
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
        "--score-type",
        choices=["adjusted", "raw", "residual"],
        default="adjusted",
        help="Typ dokovacího skóre: 'adjusted' (výchozí: MW-očištěné na škále kcal/mol), 'raw' (surové), nebo 'residual' (kolem 0)."
    )
    parser.add_argument(
        "--rows",
        choices=["proteins", "ligands"],
        default="proteins",
        help="Která entita má tvořit řádky (osa Y) grafu: 'proteins' (výchozí) nebo 'ligands'."
    )
    parser.add_argument(
        "--no-annotations",
        action="store_true",
        help="Nepřidávat genové symboly k UniProt ID cílů (zobrazit pouze čistá ID)."
    )
    parser.add_argument(
        "--both",
        action="store_true",
        help="Vykreslí jak surová skóre (kcal/mol), tak MW-reziduální/očištěná skóre."
    )
    parser.add_argument(
        "--standardize",
        choices=["none", "compound", "target"],
        default="none",
        help="Typ Z-skóre: 'none' (surová skóre), 'compound' (normalizace po sloučeninách), 'target' (po proteinech)."
    )
    parser.add_argument(
        "--cmap",
        type=str,
        default="viridis_r",
        help="Barevná paleta (výchozí: 'viridis_r', doporučeno 'mako_r', 'coolwarm_r')."
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
        default="14,10.5",
        help="Rozměry grafu v palcích 'šířka,výška' (výchozí: '14,10.5')."
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Rozlišení výstupního rastru (DPI, výchozí: 300)."
    )

    args = parser.parse_args()

    try:
        w, h = [float(x.strip()) for x in args.figsize.split(",")]
        figsize = (w, h)
    except Exception:
        figsize = (14.0, 10.5)

    # 1. Načtení matice
    try:
        matrix_df, base_out_dir, used_score_type = load_matrix(args.input, score_type=args.score_type)
    except Exception as e:
        print(f"[!] Chyba při načítání matice: {e}", file=sys.stderr)
        sys.exit(1)

    out_dir = args.output_dir if args.output_dir else (base_out_dir / "plots")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 2. Načtení anotací proteinů (pokud není zakázáno)
    annotations = None
    if not args.no_annotations:
        annotations = load_target_annotations(base_out_dir)

    # 3. Orientace matice dle zadání
    df_oriented, row_entity, col_entity = prepare_matrix_orientation(
        matrix_df, desired_rows=args.rows, annotations=annotations
    )
    df_imputed = impute_missing(df_oriented)

    print("=" * 70)
    print("Vytváření hierarchické Clustermapy dokovacích skóre:")
    print(f"  Typ skóre:         {used_score_type}")
    print(f"  Orientace os:      Řádky (Y) = {row_entity.upper()} ({df_imputed.shape[0]}), Sloupce (X) = {col_entity.upper()} ({df_imputed.shape[1]})")
    print(f"  Metrika/Metoda:    {args.metric} / {args.method}")
    print(f"  Standardizace:     {args.standardize}")
    print(f"  Výstupní složka:   {out_dir}")
    print("=" * 70)

    # Definice popisků
    if args.standardize == "compound":
        cbar_label = "Z-skóre\n(selektivita\nligandu)"
        title_suffix = " [Standardizováno po sloučeninách (Z-skóre)]"
        cmap_to_use = "coolwarm_r" if args.cmap == "viridis_r" else args.cmap
    elif args.standardize == "target":
        cbar_label = "Z-skóre\n(relativní\nafinita k cíli)"
        title_suffix = " [Standardizováno po cílech (Z-skóre)]"
        cmap_to_use = "coolwarm_r" if args.cmap == "viridis_r" else args.cmap
    else:
        cbar_label = "Dokovací afinita\nΔG [kcal/mol]\n(nižší = silnější)"
        title_suffix = " [Dokovací afinita ΔG]"
        cmap_to_use = args.cmap

    # A) Vykreslení Clustermapy dokovacích skóre
    title_main = f"Interakční Clustermapa: Methyltransferázy T. vaginalis vs. Inhibitory{title_suffix}"
    out_png_scores = out_dir / f"clustermap_docking_scores_{args.standardize}.png"
    out_pdf_scores = out_dir / f"clustermap_docking_scores_{args.standardize}.pdf"

    plot_clustermap_figure(
        df=df_imputed,
        title=title_main,
        cbar_label=cbar_label,
        out_png=out_png_scores,
        out_pdf=out_pdf_scores,
        row_entity=row_entity,
        col_entity=col_entity,
        cmap=cmap_to_use,
        standardize=args.standardize,
        metric=args.metric,
        method=args.method,
        figsize=figsize,
        dpi=args.dpi,
    )

    # B) Volitelné vykreslení MW-očištěné matice / reziduí (--both)
    adj_file = base_out_dir / "docking_mw_adjusted_matrix.csv"
    res_file = base_out_dir / "docking_mw_residuals_matrix.csv"

    target_second_file = adj_file if adj_file.exists() else (res_file if res_file.exists() else None)

    if (args.both or "residuals" in str(args.input)) and target_second_file and target_second_file.exists():
        is_adj = (target_second_file == adj_file)
        lbl_type = "MW-Adjusted (kcal/mol)" if is_adj else "MW-Residuals"
        print(f"\n[*] Nalezena matice {lbl_type}. Generuji srovnávací Clustermapu bez vlivu molekulové hmotnosti...")
        try:
            sec_df = pd.read_csv(target_second_file, index_col=0).apply(pd.to_numeric, errors="coerce")
            sec_oriented, r_ent, c_ent = prepare_matrix_orientation(
                sec_df, desired_rows=args.rows, annotations=annotations
            )
            sec_imputed = impute_missing(sec_oriented)

            file_stem = "clustermap_docking_mw_adjusted" if is_adj else "clustermap_docking_mw_residuals"
            out_png_sec = out_dir / f"{file_stem}.png"
            out_pdf_sec = out_dir / f"{file_stem}.pdf"

            cbar_sec = "MW-Očištěná afinita\n[kcal/mol]\n(nižší = silnější)" if is_adj else "Reziduum afinity\n[kcal/mol]\n(záporné = lepší než MW)"
            title_sec = "Clustermapa dokování: Skóre očištěné od vlivu molekulové hmotnosti (MW-Adjusted)" if is_adj else "Clustermapa dokování: Rezidua od modelu hmotnosti (MW Residuals)"

            plot_clustermap_figure(
                df=sec_imputed,
                title=title_sec,
                cbar_label=cbar_sec,
                out_png=out_png_sec,
                out_pdf=out_pdf_sec,
                row_entity=r_ent,
                col_entity=c_ent,
                cmap="viridis_r" if is_adj else "coolwarm_r",
                standardize="none",
                metric=args.metric,
                method=args.method,
                figsize=figsize,
                dpi=args.dpi,
            )
        except Exception as e:
            print(f"[!] Selhalo vykreslení srovnávací clustermapy: {e}", file=sys.stderr)

    print("\n[✓] Hotovo! Všechny clustermapy byly úspěšně vygenerovány.")


if __name__ == "__main__":
    main()
