# ==============================================================================
# ptm_validation.py — FASTA-based validation flags for Phobos
# ==============================================================================
# Two sequence-based annotations, computed only when a FASTA is available:
#
#  1. Amidation validation (C-terminal +1 glycine)
#     C-terminal alpha-amidation requires a glycine immediately C-terminal to
#     the amidation site in the PRECURSOR: ...X-G(-K/R)... PAM consumes the
#     glycine, transferring its nitrogen onto X -> mature peptide ends in X-NH2,
#     and the glycine is ABSENT from the observed peptide.
#     => To validate an amidated peptide we map it onto its protein sequence and
#        check that the residue at position +1 (just after the peptide's last
#        residue, in the protein) is a Glycine. No glycine downstream = likely
#        false positive (chemical/in-source artefact, mis-assignment).
#
#  2. Signal peptide flag (yes/no)
#     Lightweight von Heijne-style heuristic on the protein N-terminus
#     (n-region positive charge, h-region hydrophobic stretch, polar c-region).
#     This is an APPROXIMATION, not SignalP. It returns a boolean per protein,
#     no cleavage position. Pure Python, no external dependency.
#
# Both are best-effort: any peptide whose protein is not found in the FASTA gets
# an empty / NA flag rather than an error.
# ==============================================================================

import os
import re


# ------------------------------------------------------------------------------
# FASTA sequences (accession -> sequence)
# ------------------------------------------------------------------------------

def load_fasta_sequences(fasta_path: str) -> dict:
    """
    Parse a FASTA into {join_key: sequence}, join_key = lowercased UniProt
    accession (middle '|' field, or first field, mirroring fasta_descriptions).
    Reads sequences (needed for mapping), still light: one pass, str concat.
    """
    if not fasta_path or not os.path.exists(fasta_path):
        return {}
    seqs = {}
    key = None
    buf = []
    def _flush():
        if key and buf:
            seqs.setdefault(key, "".join(buf))
    try:
        with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    _flush()
                    buf = []
                    h = line[1:].strip()
                    m = re.search(r"^[^|]*\|([^|]+)\|", h)
                    if m:
                        acc = m.group(1)
                    else:
                        first = h.split()[0] if h.split() else h
                        acc = first.split("|")[0]
                    key = acc.strip().lower()
                else:
                    buf.append(line.strip())
            _flush()
    except Exception as e:
        print(f"  [WARN] FASTA sequence read error "
              f"({type(e).__name__}: {str(e)[:80]}) — PTM validation skipped.")
        return {}
    print(f"  -> FASTA sequences loaded: {len(seqs)} proteins")
    return seqs


# ------------------------------------------------------------------------------
# Accession → join key (same rule as fasta_descriptions)
# ------------------------------------------------------------------------------

def _acc_key(accession: str) -> str:
    if accession is None:
        return ""
    first_protein = str(accession).split(":")[0].strip()
    fields = first_protein.split("|")
    if len(fields) >= 3 and fields[0].lower() in ("sp", "tr"):
        key = fields[1]
    else:
        key = fields[0]
    return key.strip().lower()


def _strip_sequence(seq: str) -> str:
    """Remove inline modifications/charges to get the bare amino-acid string.

    Handles PEAKS inline forms:  K(+42.01) -> K, (sub A) substitutions, |z=2 …
    Substitutions '(sub X)' indicate the OBSERVED residue differs from the
    reference; for mapping we keep the reference letter that precedes it.
    """
    s = str(seq).split("|z=")[0]            # drop charge suffix if present
    s = re.sub(r"\(sub [A-Z]\)", "", s)     # drop substitution annotations
    s = re.sub(r"\([^)]*\)", "", s)         # drop inline mod masses / names
    s = re.sub(r"[^A-Za-z]", "", s)         # keep letters only
    return s.upper()


# ------------------------------------------------------------------------------
# 1. Amidation validation — C-terminal +1 glycine in the precursor
# ------------------------------------------------------------------------------

def validate_amidation(df_meta, sequences: dict,
                       seq_col: str = "Sequence",
                       acc_col: str = "Accession") -> "list":
    """
    For each peptide, return one of:
        True   -> protein found AND residue at +1 (after the peptide) is Glycine
        False  -> protein found but no downstream glycine (likely false positive)
        None   -> protein not found / peptide not located (cannot validate)

    The check is applied to ALL peptides (caller decides to use it only on the
    amidated subset). Mapping uses the bare (stripped) peptide sequence.
    """
    flags = []
    for _, row in df_meta.iterrows():
        pep = _strip_sequence(row.get(seq_col, ""))
        key = _acc_key(row.get(acc_col, ""))
        prot = sequences.get(key)
        if not pep or not prot:
            flags.append(None)
            continue
        pos = prot.find(pep)
        if pos < 0:
            # try a looser match: Leu/Ile ambiguity (I<->L are isobaric)
            prot_li = prot.replace("I", "L")
            pep_li = pep.replace("I", "L")
            pos = prot_li.find(pep_li)
        if pos < 0:
            flags.append(None)
            continue
        end = pos + len(pep)                 # index of residue just AFTER peptide
        nxt = prot[end] if end < len(prot) else ""
        flags.append(nxt == "G")
    return flags


# ------------------------------------------------------------------------------
# 2. Signal peptide flag — lightweight von Heijne-style heuristic
# ------------------------------------------------------------------------------

# Kyte-Doolittle hydrophobicity (for the h-region detection)
_KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,
       "H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,
       "T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}


def _has_signal_peptide(prot: str) -> bool:
    """
    Heuristic yes/no signal-peptide call on a protein N-terminus.
    Tripartite structure (von Heijne):
      • n-region : first ~1-5 residues, short, NOT strongly acidic
      • h-region : hydrophobic stretch (>=7 residues, high mean hydrophobicity)
                   beginning within the first ~12 residues
      • cleavage : typically within residues 15-30
    This is an APPROXIMATION (no neural net): tuned to flag the canonical
    hydrophobic N-terminal stretch while rejecting acidic cytoplasmic N-termini.
    """
    if not prot or len(prot) < 15:
        return False
    n = prot[:35].upper()
    n = "".join(c for c in n if c in _KD)
    if len(n) < 15:
        return False

    # Disqualifier : N-terminus très acide (ex. actine MDDDIAA…) — typique
    # cytoplasmique, jamais un signal peptide.
    if sum(1 for c in n[:6] if c in "DE") >= 2:
        return False

    # h-region : fenêtre de 7-9 résidus, démarrant tôt (<=10), très hydrophobe.
    best_h = -10.0
    best_start = -1
    best_win = 7
    for win in (7, 8, 9):
        for i in range(0, min(len(n) - win, 11)):
            seg = n[i:i + win]
            n_charged = sum(1 for c in seg if c in "DEKR")
            if n_charged > 1:           # cœur hydrophobe = quasi aucun chargé
                continue
            h = sum(_KD[c] for c in seg) / win
            if h > best_h:
                best_h = h; best_start = i; best_win = win
    if best_start < 0:
        return False

    core = n[best_start:best_start + best_win]
    has_pro = "P" in core[1:-1]         # proline casse l'hélice -> pas SP
    # Compter les résidus fortement hydrophobes dans le cœur (L/I/V/F/M/A/W/C)
    strong = sum(1 for c in core if c in "LIVFMAWC")

    # Critères stricts : cœur très hydrophobe (mean KD >= 2.0), au moins 5/7
    # résidus fortement hydrophobes, démarrage précoce, pas de proline.
    return bool(best_h >= 2.0 and strong >= 5 and not has_pro
                and best_start <= 10)


def signal_peptide_flags(df_meta, sequences: dict,
                         acc_col: str = "Accession") -> "list":
    """
    Per-peptide signal-peptide flag (yes/no), based on the PARENT protein.
    All peptides of a protein get the same flag. None if protein not found.
    Computed once per protein and cached.
    """
    cache = {}
    flags = []
    for _, row in df_meta.iterrows():
        key = _acc_key(row.get(acc_col, ""))
        prot = sequences.get(key)
        if not prot:
            flags.append(None)
            continue
        if key not in cache:
            cache[key] = _has_signal_peptide(prot)
        flags.append(cache[key])
    return flags


# ------------------------------------------------------------------------------
# Public entry point
# ------------------------------------------------------------------------------

def annotate(df_meta, fasta_path: str = None, search_dir: str = None,
             seq_col: str = "Sequence", acc_col: str = "Accession",
             do_amidation: bool = True, do_signalp: bool = True) -> "tuple":
    """
    Add validation columns to df_meta when a FASTA is available:
      • amidation_Gflank : True/False/None  (G at C-terminal +1 in the protein)
      • signal_peptide   : True/False/None  (heuristic SP call on parent protein)

    Returns (df_meta, info). If no FASTA / no Accession column, df is unchanged
    and info['used'] = False.
    """
    info = {"used": False, "fasta": None}
    if acc_col not in df_meta.columns:
        return df_meta, info

    if fasta_path is None and search_dir:
        try:
            from fasta_descriptions import find_fasta
            fasta_path = find_fasta(search_dir)
        except ImportError:
            fasta_path = None
    if not fasta_path or not os.path.exists(fasta_path):
        return df_meta, info

    sequences = load_fasta_sequences(fasta_path)
    if not sequences:
        return df_meta, info

    out = df_meta.copy()
    def _fmt(v):
        return "yes" if v is True else ("no" if v is False else "")
    if do_amidation:
        raw = validate_amidation(out, sequences, seq_col, acc_col)
        n_true = sum(1 for v in raw if v is True)
        n_false = sum(1 for v in raw if v is False)
        out["amidation_Gflank"] = [_fmt(v) for v in raw]
        print(f"  -> Amidation C-term+1 glycine: {n_true} with G (validated), "
              f"{n_false} without G (suspect)")
    if do_signalp:
        raw = signal_peptide_flags(out, sequences, acc_col)
        n_sp = sum(1 for v in raw if v is True)
        out["signal_peptide"] = [_fmt(v) for v in raw]
        print(f"  -> Signal peptide (heuristic): {n_sp} peptides from SP+ proteins")

    info.update({"used": True, "fasta": fasta_path})
    return out, info
