# ==============================================================================
# make_design_template.py — Generate an ExperimentalDesign template from a
# Peaks peptide export (Area columns), for the Phobos pipeline.
# ==============================================================================
# Reads the Peaks file (.xlsx 'Peptides' sheet or .csv), extracts every
# 'Area <sample>' column, and writes an ExperimentalDesign_template.csv with:
#     label;condition;replicate
# where 'label' is the suffix after 'Area ' (e.g. 'Area MA1' -> 'MA1').
#
# 'condition' and 'replicate' are PRE-FILLED with a naive guess (alphabetic
# prefix = condition, trailing digits = replicate) that the user MUST review.
#
# Usage:
#   python make_design_template.py Classeur1.xlsx
#   python make_design_template.py protein-peptides.csv -o my_design.csv
# ==============================================================================

import os
import re
import argparse
import pandas as pd


def extract_area_labels(peaks_path: str) -> list:
    """Return the list of sample labels (suffix after 'Area ') in file order."""
    ext = os.path.splitext(peaks_path)[1].lower()
    if ext in (".xlsx", ".xlsm", ".xls"):
        xl = pd.ExcelFile(peaks_path)
        sheet = "Peptides" if "Peptides" in xl.sheet_names else xl.sheet_names[0]
        header = pd.read_excel(peaks_path, sheet_name=sheet, nrows=0)
    else:
        header = pd.read_csv(peaks_path, sep=",", nrows=0)
    cols = [str(c).strip() for c in header.columns]

    area_cols = [c for c in cols if re.match(r"^area\b", c, re.IGNORECASE)]
    if not area_cols:
        area_cols = [c for c in cols if "area" in c.lower()]
    if not area_cols:
        raise ValueError("No 'Area' column found in the file.")

    labels = [re.sub(r"^area\s+", "", c, flags=re.IGNORECASE).strip()
              for c in area_cols]
    return labels


def _guess_condition_replicate(label: str) -> tuple:
    """
    Naive split: alphabetic (or alphanumeric) prefix = condition,
    trailing digits = replicate. Always review manually.
    E.g. 'MA1' -> ('MA', '1'); 'Ctrl_3' -> ('Ctrl', '3'); 'WT-2' -> ('WT', '2')
    """
    m = re.match(r"^(.*?)[\s_\-]?(\d+)$", label)
    if m:
        cond = m.group(1).strip("_- ")
        rep  = m.group(2)
        return (cond if cond else label, rep)
    return (label, "1")


def build_template(peaks_path: str, out_path: str = None) -> str:
    """Write ExperimentalDesign_template.csv next to the input (or out_path)."""
    labels = extract_area_labels(peaks_path)
    rows = []
    for lbl in labels:
        cond, rep = _guess_condition_replicate(lbl)
        rows.append({"label": lbl, "condition": cond, "replicate": rep})
    design = pd.DataFrame(rows, columns=["label", "condition", "replicate"])

    if out_path is None:
        base = os.path.dirname(os.path.abspath(peaks_path))
        out_path = os.path.join(base, "ExperimentalDesign_template.csv")

    # IMPORTANT: separator ';' (Phobos/Deimos convention)
    design.to_csv(out_path, sep=";", index=False)

    print(f"[OK] {len(labels)} samples detected.")
    print(f"     Template written: {out_path}")
    print(f"     Separator: ';'  |  columns: label;condition;replicate")
    print("\n  Detected labels (label -> guessed condition / replicate):")
    for r in rows:
        print(f"     {r['label']:<14} -> {r['condition']} / {r['replicate']}")
    print("\n  [!] REVIEW the 'condition' and 'replicate' columns before running "
          "Phobos.\n      The guess is naive (prefix=condition, digits=replicate).")
    return out_path


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Generate an ExperimentalDesign template from a Peaks export.")
    ap.add_argument("peaks_file", help="Peaks .xlsx ('Peptides' sheet) or .csv")
    ap.add_argument("-o", "--out", default=None,
                    help="Output CSV path (default: ExperimentalDesign_template.csv)")
    args = ap.parse_args()
    build_template(args.peaks_file, args.out)
