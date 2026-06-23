# ==============================================================================
# fasta_descriptions.py — Backfill protein descriptions from a FASTA file
# ==============================================================================
# PEAKS peptide exports do not always include a protein Description column.
# If a UniProt-style FASTA is present in the working directory, Phobos uses it
# to recover the description (and gene name) for each peptide's Accession.
#
# Behaviour (configured for Phobos):
#   • Auto-detection: the FIRST *.fasta / *.fa file found in the directory is
#     used automatically (filename does not matter — it can change per project).
#   • Join key: the FIRST UniProt ID, before the first '|', of the first protein
#     in the Accession group (e.g. "Q503D3|..:A0A..|.." -> "q503d3").
#   • Multi-protein groups (Accession separated by ':'): only the FIRST
#     protein's description is kept.
#
# Pure-Python (no Biopython dependency): a minimal FASTA header parser handles
# UniProt headers like:
#   >sp|A0A8M9QFQ9|A0A8M9QFQ9_DANRE Klhl13 protein OS=Danio rerio OX=7955 GN=klhl13 ...
#   >tr|Q503D3|Q503D3_DANRE Some description OS=...
#   >Q503D3|Q503D3_DANRE Some description OS=...    (no sp/tr prefix)
# ==============================================================================

import os
import re
import glob


# ------------------------------------------------------------------------------
# FASTA discovery
# ------------------------------------------------------------------------------

def find_fasta(search_dir: str) -> "str | None":
    """
    Return the path of the first FASTA file found in search_dir, or None.
    Extensions tried (case-insensitive): .fasta, .fa, .faa
    """
    if not search_dir or not os.path.isdir(search_dir):
        return None
    patterns = ["*.fasta", "*.fa", "*.faa", "*.FASTA", "*.FA", "*.FAA"]
    hits = []
    for pat in patterns:
        hits.extend(glob.glob(os.path.join(search_dir, pat)))
    # Deduplicate (case-insensitive filesystems can return dupes) and sort
    hits = sorted(set(os.path.abspath(h) for h in hits))
    return hits[0] if hits else None


# ------------------------------------------------------------------------------
# FASTA header parsing
# ------------------------------------------------------------------------------

def _parse_header(header: str) -> dict:
    """
    Parse a single FASTA header line (without the leading '>').
    Returns {accession, gene, description}.

    Accession: the ID between the first two '|' if present (UniProt sp/tr),
    otherwise the first '|'-delimited field, otherwise the first whitespace
    token.
    Description: text between the accession block and the first ' OS=' (or end).
    Gene: value after 'GN=' up to the next space, if present.
    """
    h = header.strip()
    if h.startswith(">"):
        h = h[1:].strip()

    # --- Accession ---
    acc = None
    # UniProt style: prefix|ACC|ENTRY  -> take the middle field
    m_mid = re.search(r"^[^|]*\|([^|]+)\|", h)
    if m_mid:
        acc = m_mid.group(1)
    else:
        # Fallback: first '|'-field, else first token
        first_field = h.split()[0] if h.split() else h
        acc = first_field.split("|")[0]

    # --- Gene (GN=) ---
    gene = None
    m_gn = re.search(r"GN=([^\s]+)", h)
    if m_gn:
        gene = m_gn.group(1)

    # --- Description ---
    # Everything before ' OS=' (UniProt) ...
    before_os = re.split(r"\sOS=", h, maxsplit=1)[0]
    # ... minus the leading "prefix|ACC|ENTRY " or "ACC|ENTRY " block.
    desc = re.sub(r"^>?(?:[^|\s]*\|)?[^|\s]+\|[^\s]+\s+", "", before_os)
    if desc == before_os:
        # Header had no '|' block (e.g. ">ACC description ...")
        desc = re.sub(r"^>?[^\s]+\s+", "", before_os)
    desc = desc.strip()

    return {"accession": acc, "gene": gene, "description": desc}


def load_fasta_index(fasta_path: str) -> dict:
    """
    Build a lookup {join_key: {description, gene}} from the FASTA.
    join_key = lowercased, trimmed UniProt accession (middle '|' field).
    Only header lines are read (sequences are skipped) — fast and memory-light.
    """
    index = {}
    n_entries = 0
    try:
        with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                n_entries += 1
                info = _parse_header(line)
                acc = info.get("accession")
                if not acc:
                    continue
                key = acc.strip().lower()
                # First occurrence wins (stable, like a left join on first hit)
                if key not in index:
                    index[key] = {"description": info.get("description") or "",
                                  "gene": info.get("gene") or ""}
    except Exception as e:
        print(f"  [WARN] FASTA read error ({type(e).__name__}: {str(e)[:80]}) "
              f"— descriptions not recovered.")
        return {}
    print(f"  -> FASTA parsed: {n_entries} headers, {len(index)} unique accessions indexed")
    return index


# ------------------------------------------------------------------------------
# Accession → join key
# ------------------------------------------------------------------------------

def _first_accession_key(accession: str) -> str:
    """
    From a PEAKS Accession field, derive the join key = first UniProt ID
    (before the first '|') of the FIRST protein (before the first ':').

    Examples
    --------
    "Q503D3|Q503D3_DANRE:A0A8M9QFQ9|A0A8M9QFQ9_DANRE" -> "q503d3"
    "Q9I8V0|PRV2_DANRE"                               -> "q9i8v0"
    "sp|P12345|NAME_HUMAN"                            -> "sp"  (see note)

    Note: PEAKS accessions are typically "ID|ENTRY", so the first field IS the
    ID. If the export uses the "sp|ID|ENTRY" UniProt form, the first field is
    the prefix; we then fall back to the middle field.
    """
    if accession is None:
        return ""
    first_protein = str(accession).split(":")[0].strip()
    fields = first_protein.split("|")
    if len(fields) >= 3 and fields[0].lower() in ("sp", "tr"):
        # sp|ID|ENTRY form -> ID is the middle field
        key = fields[1]
    else:
        # ID|ENTRY form (PEAKS default) -> ID is the first field
        key = fields[0]
    return key.strip().lower()


# ------------------------------------------------------------------------------
# Public entry point
# ------------------------------------------------------------------------------

def add_descriptions(df, fasta_path: str = None, search_dir: str = None,
                     accession_col: str = "Accession") -> "tuple":
    """
    Add a 'Description' column (and fill an empty 'Gene' if needed) to df by
    joining on the protein accession against a FASTA index.

    Parameters
    ----------
    df            : DataFrame containing an Accession column
    fasta_path    : explicit FASTA path; if None, auto-detect in search_dir
    search_dir    : directory to scan for a FASTA (defaults to df's location
                    is unknown here, so caller should pass it)
    accession_col : name of the accession column (default 'Accession')

    Returns
    -------
    (df, info) where df has a 'Description' column (added or filled) and info is
    a dict with {used, fasta, n_total, n_matched, n_missing}. If no FASTA is
    found or the column is absent, df is returned unchanged with used=False.
    """
    import pandas as pd

    info = {"used": False, "fasta": None, "n_total": len(df),
            "n_matched": 0, "n_missing": 0}

    if accession_col not in df.columns:
        return df, info

    # Resolve FASTA
    if fasta_path is None:
        fasta_path = find_fasta(search_dir)
    if not fasta_path or not os.path.exists(fasta_path):
        return df, info

    print(f"\n[FASTA] Recovering protein descriptions from: "
          f"{os.path.basename(fasta_path)}")
    index = load_fasta_index(fasta_path)
    if not index:
        return df, info

    keys = df[accession_col].map(_first_accession_key)
    descs = keys.map(lambda k: index.get(k, {}).get("description", ""))
    genes_fa = keys.map(lambda k: index.get(k, {}).get("gene", ""))

    n_matched = int((descs.str.len() > 0).sum())
    n_missing = len(df) - n_matched

    out = df.copy()
    out["Description"] = descs.values

    # Fill Gene only where it is missing/empty in the export
    if "Gene" in out.columns:
        gene_empty = out["Gene"].isna() | (out["Gene"].astype(str).str.strip() == "")
        out.loc[gene_empty, "Gene"] = genes_fa[gene_empty].values
    else:
        out["Gene"] = genes_fa.values

    info.update({"used": True, "fasta": fasta_path,
                 "n_matched": n_matched, "n_missing": n_missing})
    print(f"  -> {n_matched}/{len(df)} peptides matched a description "
          f"({n_missing} without match)")
    if n_missing > 0:
        print(f"  [INFO] Unmatched accessions keep an empty description "
              f"(ID not found in this FASTA).")
    return out, info


if __name__ == "__main__":
    import argparse
    import pandas as pd
    ap = argparse.ArgumentParser(
        description="Backfill protein descriptions in a PEAKS export from a FASTA.")
    ap.add_argument("peaks_csv", help="PEAKS peptide CSV/XLSX")
    ap.add_argument("fasta", nargs="?", default=None,
                    help="FASTA file (optional; auto-detected if omitted)")
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    ext = os.path.splitext(args.peaks_csv)[1].lower()
    if ext in (".xlsx", ".xlsm", ".xls"):
        df = pd.read_excel(args.peaks_csv)
    else:
        df = pd.read_csv(args.peaks_csv, sep=None, engine="python")

    search = os.path.dirname(os.path.abspath(args.peaks_csv))
    df2, info = add_descriptions(df, fasta_path=args.fasta, search_dir=search)

    out = args.out or os.path.splitext(args.peaks_csv)[0] + "_with_desc.csv"
    df2.to_csv(out, index=False)
    print(f"\n[OK] Written: {out}")
    print(f"     {info['n_matched']}/{info['n_total']} descriptions recovered.")
