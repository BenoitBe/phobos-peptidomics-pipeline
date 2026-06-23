# ==============================================================================
# ptm_analysis.py — Targeted PTM sub-analyses for Phobos
# ==============================================================================
# In addition to the global peptide-level analysis, Phobos can run the SAME
# pipeline (volcano / ANOVA / heatmap / UpSet / scatter) on PTM-filtered
# subsets: each selected modification yields an independent differential
# analysis on the peptides carrying that PTM.
#
# Detection is COMBINED: a peptide carries a PTM if it matches EITHER
#   (a) the text annotation in the PTM column (e.g. "Oxidation (M)"), OR
#   (b) an inline delta mass in the sequence (e.g. "M(+15.99)")
# within a tight mass tolerance.
#
# Trimethylation (+42.047) vs Acetylation (+42.011) differ by only 0.036 Da:
# the tolerance is set tight (±0.01 Da) and the text annotation takes
# precedence to disambiguate when present.
# ==============================================================================

import re
import numpy as np
import pandas as pd


# ------------------------------------------------------------------------------
# PTM registry — delta masses (monoisotopic) + text patterns
# ------------------------------------------------------------------------------
# Each entry:
#   key        : short id (used for sheet names / file prefixes)
#   label      : human-readable name
#   deltas     : list of (mass, tol) accepted for inline detection
#   text_regex : regex matched (case-insensitive) against the PTM column
#   priority_text_excludes : optional list of regexes that, if matched, EXCLUDE
#                            the peptide (disambiguation, e.g. trimethyl vs acetyl)
PTM_REGISTRY = {
    "Phospho": {
        "label": "Phosphorylation",
        "deltas": [(79.9663, 0.01)],
        "text_regex": r"phospho",
    },
    "Oxidation": {
        "label": "Oxidation",
        "deltas": [(15.9949, 0.01)],
        "text_regex": r"oxidation|oxidative",
    },
    "Methyl": {
        "label": "Methylation (mono)",
        "deltas": [(14.0157, 0.008)],
        "text_regex": r"(?<!di)(?<!tri)methyl(?!.*di)(?!.*tri)",
    },
    "Dimethyl": {
        "label": "Di-methylation",
        "deltas": [(28.0313, 0.008)],
        "text_regex": r"dimethyl|di-methyl",
    },
    "Trimethyl": {
        "label": "Tri-methylation",
        "deltas": [(42.0470, 0.008)],
        "text_regex": r"trimethyl|tri-methyl",
        # +42.047 (trimethyl) is 0.036 Da from +42.011 (acetyl): if the text
        # says acetyl, this is NOT a trimethyl.
        "exclude_text": r"acetyl",
    },
    "Methyl_all": {
        "label": "Methylation (mono+di+tri combined)",
        "deltas": [(14.0157, 0.008), (28.0313, 0.008), (42.0470, 0.008)],
        "text_regex": r"methyl",
        "exclude_text": r"acetyl",
    },
    "Acetyl": {
        "label": "Acetylation",
        "deltas": [(42.0106, 0.01)],
        "text_regex": r"acetyl",
        # +42.011 (acetyl) vs +42.047 (trimethyl): if text says (tri)methyl,
        # not an acetyl.
        "exclude_text": r"methyl",
    },
}

# Default selectable PTMs offered interactively (order matters for the menu)
DEFAULT_PTM_MENU = ["Phospho", "Oxidation", "Methyl_all", "Acetyl"]


# ------------------------------------------------------------------------------
# Detection
# ------------------------------------------------------------------------------

def _inline_masses(seq: str) -> list:
    """Extract inline delta masses from a sequence string.

    Handles PEAKS' occasional double sign (e.g. '(++0.98)' -> +0.98,
    '(+-0.98)' -> -0.98).
    """
    out = []
    for raw in re.findall(r"\(([+-]{1,2}\d+\.?\d*)\)", str(seq)):
        s = raw.replace("++", "+").replace("+-", "-").replace("-+", "-")
        try:
            out.append(float(s))
        except ValueError:
            continue
    return out


def _matches_ptm(seq: str, ptm_text: str, spec: dict) -> bool:
    """True if the peptide carries the PTM, by text OR by inline mass."""
    txt = (ptm_text or "").lower()

    # Exclusion first (disambiguation)
    excl = spec.get("exclude_text")
    text_says_excluded = bool(excl and re.search(excl, txt))

    # (a) text annotation
    if spec.get("text_regex") and re.search(spec["text_regex"], txt):
        if not text_says_excluded:
            return True

    # (b) inline delta mass — only if text does not explicitly exclude it
    if not text_says_excluded:
        masses = _inline_masses(seq)
        for m in masses:
            for (target, tol) in spec["deltas"]:
                if abs(m - target) <= tol:
                    return True
    return False


def detect_ptm_mask(df_meta: pd.DataFrame, ptm_key: str,
                    seq_col: str = "Sequence",
                    ptm_col: str = "Modifications") -> np.ndarray:
    """
    Boolean mask over rows of df_meta for peptides carrying ptm_key.

    df_meta must contain the sequence (with inline masses) and ideally the
    PTM text column. Missing columns are tolerated (detection falls back to
    whichever is available).
    """
    spec = PTM_REGISTRY[ptm_key]
    seqs = (df_meta[seq_col].astype(str) if seq_col in df_meta.columns
            else pd.Series([""] * len(df_meta)))
    txts = (df_meta[ptm_col].astype(str) if ptm_col in df_meta.columns
            else pd.Series([""] * len(df_meta)))
    mask = np.array([_matches_ptm(s, t, spec)
                     for s, t in zip(seqs, txts)], dtype=bool)
    return mask


def inventory_ptms(df_meta: pd.DataFrame, seq_col: str = "Sequence",
                   ptm_col: str = "Modifications") -> dict:
    """
    Count how many peptides match each registry PTM (for the interactive menu
    and for logging). Returns {ptm_key: n_peptides}.
    """
    counts = {}
    for key in PTM_REGISTRY:
        counts[key] = int(detect_ptm_mask(df_meta, key, seq_col, ptm_col).sum())
    return counts


# ------------------------------------------------------------------------------
# Interactive selection
# ------------------------------------------------------------------------------

def ask_ptm_selection(df_meta: pd.DataFrame,
                      seq_col: str = "Sequence",
                      ptm_col: str = "Modifications") -> list:
    """
    Interactive menu: show available PTMs (with peptide counts) and let the
    user select which ones to analyse. Returns a list of ptm_keys.
    Returns [] if the user selects none.
    """
    counts = inventory_ptms(df_meta, seq_col, ptm_col)

    print("\n" + "=" * 60)
    print("  PTM-TARGETED SUB-ANALYSES")
    print("=" * 60)
    print("  Run the full pipeline (volcano/ANOVA/heatmap) on PTM subsets.")
    print("  Detected peptides per modification (text OR inline mass):\n")

    menu = [k for k in DEFAULT_PTM_MENU]
    # Append any other registry PTM that is present but not in the default menu
    for k in PTM_REGISTRY:
        if k not in menu and counts.get(k, 0) > 0 and k not in ("Methyl", "Dimethyl", "Trimethyl"):
            menu.append(k)

    for i, k in enumerate(menu, 1):
        n = counts.get(k, 0)
        flag = "" if n > 0 else "  (none detected — will be skipped)"
        print(f"    [{i}] {PTM_REGISTRY[k]['label']:<34} {n:>6} peptides{flag}")
    print(f"    [0] None — skip PTM sub-analyses")

    raw = input("\n  Select PTMs to analyse (comma-separated, e.g. 1,2,4) -> ").strip()
    if raw in ("", "0"):
        print("  -> No PTM sub-analysis selected.")
        return []

    chosen = []
    for tok in raw.split(","):
        tok = tok.strip()
        if not tok.isdigit():
            continue
        idx = int(tok)
        if 1 <= idx <= len(menu):
            key = menu[idx - 1]
            if counts.get(key, 0) == 0:
                print(f"  [WARN] '{PTM_REGISTRY[key]['label']}' has no peptide "
                      f"— skipped.")
                continue
            if key not in chosen:
                chosen.append(key)

    if chosen:
        print(f"\n  [OK] PTM analyses selected: "
              f"{', '.join(PTM_REGISTRY[k]['label'] for k in chosen)}")
    else:
        print("  -> No valid PTM selected.")
    return chosen


# ------------------------------------------------------------------------------
# Sub-analysis driver
# ------------------------------------------------------------------------------

def run_ptm_subanalyses(ptm_keys: list,
                        mat_filt: pd.DataFrame, mat_imp: pd.DataFrame,
                        df_pep_filt: pd.DataFrame, design: pd.DataFrame,
                        params: dict, out_dir: str,
                        analysis_callbacks: dict,
                        n_iter: int = 100,
                        min_peptides: int = 12) -> dict:
    """
    For each selected PTM, filter the imputed/filtered matrices to the carrying
    peptides and re-run the differential + ANOVA pipeline on that subset.

    Parameters
    ----------
    ptm_keys          : list of registry keys selected by the user
    mat_filt/mat_imp  : full filtered / imputed log2 matrices (peptides × samples)
    df_pep_filt       : peptide metadata aligned with the matrices (same order)
    design            : experimental design
    params            : pipeline params (thresholds)
    out_dir           : output directory (PTM figures go in a subfolder)
    analysis_callbacks: dict of the phobos.py functions to reuse, with keys:
                        'differential', 'volcanoes', 'anova', 'upset'
                        (so this module stays decoupled from phobos internals)
    n_iter            : robustness iterations
    min_peptides      : minimum peptides required to run a subset (else skipped)

    Returns
    -------
    dict {ptm_key: result_bundle} where result_bundle holds the dataframes and
    figure paths needed by the Excel exporter. Empty subsets are skipped.
    """
    import os

    results = {}
    seq_col = "Sequence" if "Sequence" in df_pep_filt.columns else "peptide_id"
    ptm_col = "Modifications" if "Modifications" in df_pep_filt.columns else None

    for key in ptm_keys:
        spec = PTM_REGISTRY[key]
        label = spec["label"]
        mask = detect_ptm_mask(df_pep_filt, key, seq_col,
                               ptm_col or "Modifications")
        n_sub = int(mask.sum())

        print(f"\n{'─'*60}")
        print(f"  PTM SUB-ANALYSIS — {label}  ({n_sub} peptides)")
        print(f"{'─'*60}")

        if n_sub < min_peptides:
            print(f"  [SKIP] Only {n_sub} peptides (< {min_peptides} required) "
                  f"— sub-analysis skipped.")
            continue

        # Dedicated subfolder for this PTM's figures
        ptm_dir = os.path.join(out_dir, f"ptm_{key}")
        os.makedirs(ptm_dir, exist_ok=True)

        # Subset matrices + metadata (preserve order, reset index)
        sub_filt = mat_filt.loc[mask].reset_index(drop=True)
        sub_imp  = mat_imp.loc[mask].reset_index(drop=True)
        sub_meta = df_pep_filt.loc[mask].reset_index(drop=True)
        sub_filt.columns = mat_filt.columns
        sub_imp.columns  = mat_imp.columns

        bundle = {"label": label, "key": key, "n_peptides": n_sub,
                  "dir": ptm_dir}

        # 1) Differential analysis (limma eBayes) + robustness + scatter
        try:
            df_res, contrasts, scatter_files = analysis_callbacks["differential"](
                sub_imp, sub_filt, sub_meta, design, params, ptm_dir,
                n_iter=n_iter)
            bundle["df_results"] = df_res
            bundle["contrasts"] = contrasts
            bundle["scatter_files"] = scatter_files
        except Exception as e:
            print(f"  [WARN] Differential failed for {label}: "
                  f"{type(e).__name__}: {e}")
            continue

        # 2) Volcanoes
        try:
            volc_files, facet = analysis_callbacks["volcanoes"](
                df_res, contrasts, params, ptm_dir)
            bundle["volc_files"] = volc_files
            bundle["facet_volc"] = facet
        except Exception as e:
            print(f"  [WARN] Volcanoes failed for {label}: {e}")
            bundle["volc_files"], bundle["facet_volc"] = [], None

        # 3) ANOVA + heatmaps
        try:
            (df_anova, hm_classic, hm_clusters, hm_violin,
             mat_z, cluster_map, sig_ids) = analysis_callbacks["anova"](
                sub_imp, sub_meta, design, params, ptm_dir)
            bundle.update({
                "df_anova": df_anova, "hm_classic": hm_classic,
                "hm_clusters": hm_clusters, "hm_violin": hm_violin,
                "mat_zscore": mat_z, "sig_ids": sig_ids,
            })
        except Exception as e:
            print(f"  [WARN] ANOVA failed for {label}: {e}")

        # 4) UpSet
        try:
            upset_file, df_inter = analysis_callbacks["upset"](
                df_res, contrasts, params, ptm_dir)
            bundle["upset_file"] = upset_file
            bundle["df_intersections"] = df_inter
        except Exception as e:
            print(f"  [WARN] UpSet failed for {label}: {e}")
            bundle["upset_file"], bundle["df_intersections"] = None, None

        results[key] = bundle
        print(f"  [OK] {label} sub-analysis complete "
              f"({len(contrasts)} contrasts, {n_sub} peptides).")

    return results
