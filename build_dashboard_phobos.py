# ==============================================================================
# build_dashboard_phobos.py — Interactive HTML dashboard for Phobos (DDA peptide)
# ==============================================================================
# Companion to deimos build_dashboardv7.py, specialised for PEPTIDE-level data.
#
# Key differences vs the Deimos dashboard:
#   • Peptide/Protein toggle (rollup view from Protein_Aggregation sheet)
#   • Peptide-specific facets: sequence, charge state, modifications
#   • Volcano coloured by charge / modification status (optional)
#   • Reads the Phobos Excel workbook (Differential_Expression, Log2_Impute,
#     Protein_Aggregation, UMAP, ANOVA_Results)
#
# Output: a SINGLE self-contained .html file (data embedded as JSON, Chart.js
# loaded from local file if present, otherwise CDN). Works fully offline.
#
# Usage:
#   from build_dashboard_phobos import build_dashboard
#   build_dashboard("Phobos_output/PeptideAnalysis_Results.xlsx",
#                   design_df, "phobos_dashboard.html", params=params)
# ==============================================================================

import os
import re
import json
import base64
import numpy as np
import pandas as pd


# ------------------------------------------------------------------------------
# Local Chart.js (offline) — falls back to CDN if the file is missing
# ------------------------------------------------------------------------------
_CHARTJS_LOCAL = "chart.umd.min.js"
_CHARTJS_CDN   = "https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js"


def _chartjs_tag(base_dir: str) -> str:
    """Inline the local Chart.js for a fully offline dashboard, else CDN."""
    local = os.path.join(base_dir, _CHARTJS_LOCAL)
    if os.path.exists(local):
        with open(local, "r", encoding="utf-8") as f:
            return f"<script>{f.read()}</script>"
    return f'<script src="{_CHARTJS_CDN}"></script>'


# ------------------------------------------------------------------------------
# Data extraction from the Phobos workbook
# ------------------------------------------------------------------------------

def _safe_read(xl_path: str, sheet: str) -> "pd.DataFrame | None":
    try:
        return pd.read_excel(xl_path, sheet_name=sheet)
    except Exception:
        return None


def _detect_contrasts(diff_cols: list) -> list:
    """Extract contrast names from '<contrast>_diff' columns."""
    return [c[:-5] for c in diff_cols if c.endswith("_diff")]


def _df_to_records(df: pd.DataFrame, cols: list) -> list:
    """Subset + NaN-safe conversion to list of dicts (JSON-ready)."""
    sub = df[[c for c in cols if c in df.columns]].copy()
    sub = sub.replace({np.nan: None})
    return sub.to_dict(orient="records")


def _extract_payload(xlsx_path: str, design: pd.DataFrame,
                     params: dict) -> dict:
    """Build the JSON payload consumed by the dashboard JS."""
    payload = {"meta": {}, "contrasts": [], "peptides": [],
               "proteins": [], "samples": [], "umap": [],
               "thresholds": {}, "facets": {}}

    # --- thresholds ---
    payload["thresholds"] = {
        "use_padj":   bool(params.get("volcano_use_padj", False)),
        "p_thresh":   float(params.get("volcano_p_thresh", 0.05)),
        "lfc_min":    float(params.get("volcano_lfc_min", np.log2(1.5))),
        "ratio_min":  float(params.get("volcano_ratio_min", 1.5)),
    }

    # --- samples / design ---
    cond_map = {}
    for _, r in design.iterrows():
        cond_map[str(r["label"])] = str(r["condition"])
    payload["samples"] = [{"label": str(r["label"]),
                           "condition": str(r["condition"])}
                          for _, r in design.iterrows()]
    payload["meta"]["conditions"] = sorted(set(cond_map.values()))

    # --- differential expression (peptide) ---
    diff = _safe_read(xlsx_path, "Differential_Expression")
    contrasts = []
    if diff is not None:
        contrasts = _detect_contrasts(diff.columns.tolist())
        payload["contrasts"] = contrasts

        keep = ["peptide_id", "Sequence", "Modifications", "Charge",
                "Accession", "Gene", "Description", "Score_10lgP",
                "imputed", "num_imputed"]
        keep = [k for k in keep if k in diff.columns]
        for c in contrasts:
            keep += [f"{c}_diff", f"{c}_p.val", f"{c}_p.adj",
                     f"Pi_Score_{c}", f"Robustness_Score_{c}"]
        payload["peptides"] = _df_to_records(diff, keep)

        # facet inventories
        if "Charge" in diff.columns:
            charges = sorted(pd.Series(diff["Charge"].dropna().unique())
                             .astype(str).tolist())
            payload["facets"]["charges"] = charges
        if "Modifications" in diff.columns:
            mods = diff["Modifications"].fillna("").astype(str)
            payload["facets"]["has_mod_count"] = int((mods.str.len() > 0).sum())
            payload["facets"]["no_mod_count"]  = int((mods.str.len() == 0).sum())

    # --- per-sample intensities (peptide) ---
    log2 = _safe_read(xlsx_path, "Log2_Impute")
    if log2 is not None:
        sample_cols = [c for c in log2.columns
                       if c in cond_map]   # only real sample columns
        payload["meta"]["sample_cols"] = sample_cols
        # attach intensities keyed by peptide_id (for boxplots / PCA)
        intens = {}
        for _, row in log2.iterrows():
            pid = row.get("peptide_id")
            intens[str(pid)] = [None if pd.isna(row[c]) else round(float(row[c]), 3)
                                for c in sample_cols]
        payload["peptide_intensities"] = intens

    # --- per-sample PRE-imputation matrix (real missing values) for QC ---
    pre = _safe_read(xlsx_path, "Log2_PreImpute")
    if pre is not None:
        pre_cols = [c for c in pre.columns if c in cond_map]
        # Comptage des vrais trous par échantillon (NaN = non détecté avant imputation)
        miss_by = {c: int(pre[c].isna().sum()) for c in pre_cols}
        det_by  = {c: int(pre[c].notna().sum()) for c in pre_cols}
        payload["qc_preimpute"] = {
            "sample_cols": pre_cols,
            "missing": miss_by,
            "detected": det_by,
            "n_peptides": int(len(pre)),
        }

    # --- protein aggregation (rollup) ---
    prot = _safe_read(xlsx_path, "Protein_Aggregation")
    if prot is not None and "Accession" in prot.columns:
        sample_cols_p = [c for c in prot.columns if c in cond_map]
        payload["meta"]["protein_sample_cols"] = sample_cols_p
        prot_recs = []
        for _, row in prot.iterrows():
            prot_recs.append({
                "Accession": str(row["Accession"]),
                "intensities": [None if pd.isna(row[c]) else round(float(row[c]), 3)
                                for c in sample_cols_p],
            })
        payload["proteins"] = prot_recs

        # peptide → protein count (for the rollup view)
        if diff is not None and "Accession" in diff.columns:
            pep_per_prot = (diff.groupby("Accession")["peptide_id"]
                            .count().to_dict())
            payload["meta"]["pep_per_prot"] = {str(k): int(v)
                                               for k, v in pep_per_prot.items()}

    # --- UMAP coords ---
    umap = _safe_read(xlsx_path, "UMAP")
    if umap is not None:
        payload["umap"] = _df_to_records(
            umap, ["label", "condition", "UMAP1", "UMAP2"])

    # --- ANOVA significant (for the heatmap tab) ---
    anova = _safe_read(xlsx_path, "ANOVA_Results")
    if anova is not None:
        sig = anova[anova.get("significant", False) == True] if "significant" in anova.columns else anova.head(0)
        payload["meta"]["n_sig_anova"] = int(len(sig))
        payload["meta"]["n_peptides"]  = int(len(anova))

    # --- summary counts ---
    payload["meta"]["n_samples"]   = len(design)
    payload["meta"]["n_contrasts"] = len(contrasts)
    if diff is not None:
        payload["meta"]["n_peptides_total"] = int(len(diff))
        if "Accession" in diff.columns:
            payload["meta"]["n_proteins_total"] = int(diff["Accession"].nunique())

    return payload


# ------------------------------------------------------------------------------
# HTML template
# ------------------------------------------------------------------------------

def _html_template(payload_json: str, chartjs_tag: str,
                   title: str = "Phobos — Peptide-level Heuristics for Bottom-up Omics Suite") -> str:
    """Returns the full self-contained HTML string."""
    # NOTE: kept as one f-string with doubled braces for the CSS/JS blocks.
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{title}</title>
{chartjs_tag}
<style>
  :root {{
    --bg:#0f1419; --panel:#1a2129; --panel2:#222c37; --ink:#e6edf3;
    --muted:#8b97a6; --accent:#5B8FF9; --up:#E74C3C; --down:#3498DB;
    --ns:#4a5568; --line:#2d3742; --good:#5AD8A6;
  }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; font-family:'Segoe UI',system-ui,sans-serif;
         background:var(--bg); color:var(--ink); font-size:14px; }}
  header {{ background:linear-gradient(135deg,#1a2129,#222c37);
           padding:16px 24px; border-bottom:1px solid var(--line);
           display:flex; align-items:center; gap:16px; flex-wrap:wrap; }}
  header h1 {{ font-size:20px; margin:0; font-weight:600; }}
  header .sub {{ color:var(--muted); font-size:12px; }}
  header .credit {{ color:var(--accent); font-size:11px; font-weight:600;
                   letter-spacing:.3px; opacity:.85; margin-top:1px; }}
  .badge {{ background:var(--panel2); padding:3px 10px; border-radius:12px;
           font-size:12px; color:var(--muted); }}
  .toggle-wrap {{ margin-left:auto; display:flex; gap:8px; align-items:center; }}
  .toggle {{ display:flex; background:var(--panel); border-radius:8px;
            padding:3px; border:1px solid var(--line); }}
  .toggle button {{ background:none; border:none; color:var(--muted);
            padding:6px 16px; cursor:pointer; border-radius:6px;
            font-size:13px; font-weight:500; transition:all .15s; }}
  .toggle button.active {{ background:var(--accent); color:#fff; }}
  nav {{ display:flex; gap:4px; padding:0 24px; background:var(--panel);
        border-bottom:1px solid var(--line); overflow-x:auto; }}
  nav button {{ background:none; border:none; color:var(--muted);
        padding:12px 18px; cursor:pointer; font-size:13px;
        border-bottom:2px solid transparent; white-space:nowrap; }}
  nav button.active {{ color:var(--ink); border-bottom-color:var(--accent); }}
  main {{ padding:24px; max-width:1400px; margin:0 auto; }}
  .grid {{ display:grid; gap:16px; }}
  .cards {{ grid-template-columns:repeat(auto-fit,minmax(160px,1fr)); margin-bottom:20px; }}
  .card {{ background:var(--panel); border:1px solid var(--line);
          border-radius:10px; padding:16px; }}
  .card .v {{ font-size:26px; font-weight:700; color:var(--accent); }}
  .card .l {{ font-size:12px; color:var(--muted); margin-top:4px; }}
  .panel {{ background:var(--panel); border:1px solid var(--line);
           border-radius:10px; padding:18px; margin-bottom:16px; }}
  .panel h3 {{ margin:0 0 14px; font-size:15px; font-weight:600; }}
  .row {{ display:flex; gap:16px; flex-wrap:wrap; align-items:flex-start; }}
  .ctrl {{ display:flex; gap:8px; align-items:center; margin-bottom:12px;
          flex-wrap:wrap; }}
  select, input[type=text] {{ background:var(--panel2); color:var(--ink);
          border:1px solid var(--line); border-radius:6px; padding:6px 10px;
          font-size:13px; }}
  label.lbl {{ color:var(--muted); font-size:12px; }}
  .flip-btn {{ background:var(--panel2); color:var(--ink); border:1px solid var(--line);
          border-radius:6px; padding:6px 12px; font-size:13px; cursor:pointer;
          transition:all .15s; }}
  .flip-btn:hover {{ background:var(--accent); color:#fff; border-color:var(--accent); }}
  .flip-btn.active {{ background:var(--accent); color:#fff; border-color:var(--accent); }}
  canvas {{ max-width:100%; }}
  table {{ width:100%; border-collapse:collapse; font-size:12px; }}
  th, td {{ text-align:left; padding:7px 10px; border-bottom:1px solid var(--line); }}
  th {{ color:var(--muted); cursor:pointer; user-select:none; position:sticky;
       top:0; background:var(--panel); }}
  th:hover {{ color:var(--ink); }}
  tr:hover td {{ background:var(--panel2); }}
  .tbl-wrap {{ max-height:520px; overflow:auto; border:1px solid var(--line);
              border-radius:8px; }}
  .pill {{ padding:2px 8px; border-radius:10px; font-size:11px; font-weight:600; }}
  .pill.up {{ background:rgba(231,76,60,.18); color:#ff8a7a; }}
  .pill.down {{ background:rgba(52,152,219,.18); color:#7fc1f0; }}
  .pill.ns {{ background:rgba(74,85,104,.25); color:var(--muted); }}
  .hide {{ display:none !important; }}
  .tab-page {{ display:none; }}
  .tab-page.active {{ display:block; }}
  .legend {{ display:flex; gap:14px; font-size:12px; color:var(--muted);
            margin-top:8px; flex-wrap:wrap; }}
  .legend i {{ width:10px; height:10px; border-radius:2px; display:inline-block;
              margin-right:5px; vertical-align:middle; }}
  .note {{ color:var(--muted); font-size:12px; line-height:1.5; }}
  .seq {{ font-family:'Consolas',monospace; font-size:11px;
         letter-spacing:.5px; word-break:break-all; }}
  td.desc {{ max-width:240px; font-size:11px; color:var(--muted);
         overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }}
  td.desc:hover {{ white-space:normal; }}
</style>
</head>
<body>
<header>
  <svg viewBox="0 0 120 120" width="42" height="42" style="flex-shrink:0;filter:drop-shadow(0 2px 4px rgba(0,0,0,.5));" xmlns="http://www.w3.org/2000/svg">
    <defs>
      <radialGradient id="ph-surf" cx="38%" cy="35%" r="68%">
        <stop offset="0%" stop-color="#cdbda3"/><stop offset="45%" stop-color="#917c66"/><stop offset="100%" stop-color="#403229"/>
      </radialGradient>
      <radialGradient id="ph-term" cx="70%" cy="62%" r="58%">
        <stop offset="0%" stop-color="#160d06" stop-opacity="0.65"/><stop offset="100%" stop-color="#160d06" stop-opacity="0"/>
      </radialGradient>
      <radialGradient id="ph-spec" cx="32%" cy="28%" r="32%">
        <stop offset="0%" stop-color="#ece1cd" stop-opacity="0.5"/><stop offset="100%" stop-color="#ece1cd" stop-opacity="0"/>
      </radialGradient>
      <radialGradient id="ph-stick" cx="42%" cy="38%" r="62%">
        <stop offset="0%" stop-color="#2a1f15"/><stop offset="60%" stop-color="#3d2f22"/><stop offset="100%" stop-color="#5c4a35"/>
      </radialGradient>
      <clipPath id="ph-clip">
        <path d="M 60 14 C 74 13,90 17,98 28 C 106 40,105 55,100 67 C 94 80,82 90,68 96 C 54 101,38 100,27 92 C 15 83,11 68,13 54 C 15 40,23 27,34 20 C 43 14,51 14,60 14 Z"/>
      </clipPath>
    </defs>
    <path d="M 60 14 C 74 13,90 17,98 28 C 106 40,105 55,100 67 C 94 80,82 90,68 96 C 54 101,38 100,27 92 C 15 83,11 68,13 54 C 15 40,23 27,34 20 C 43 14,51 14,60 14 Z" fill="url(#ph-surf)"/>
    <path d="M 60 14 C 74 13,90 17,98 28 C 106 40,105 55,100 67 C 94 80,82 90,68 96 C 54 101,38 100,27 92 C 15 83,11 68,13 54 C 15 40,23 27,34 20 C 43 14,51 14,60 14 Z" fill="url(#ph-term)"/>
    <path d="M 60 14 C 74 13,90 17,98 28 C 106 40,105 55,100 67 C 94 80,82 90,68 96 C 54 101,38 100,27 92 C 15 83,11 68,13 54 C 15 40,23 27,34 20 C 43 14,51 14,60 14 Z" fill="url(#ph-spec)"/>
    <g clip-path="url(#ph-clip)">
      <ellipse cx="42" cy="50" rx="17" ry="14" fill="url(#ph-stick)" opacity="0.9"/>
      <ellipse cx="42" cy="50" rx="17" ry="14" fill="none" stroke="#241a10" stroke-width="1.5" opacity="0.6"/>
      <ellipse cx="40" cy="48" rx="15.5" ry="12.5" fill="none" stroke="#c2b292" stroke-width="0.5" opacity="0.35"/>
      <circle cx="42" cy="50" r="2.3" fill="#1c130b" opacity="0.5"/>
      <g stroke="#2e2318" stroke-width="0.5" opacity="0.30" fill="none">
        <path d="M 58 48 Q 78 52, 95 60"/><path d="M 57 54 Q 74 62, 90 72"/><path d="M 56 44 Q 72 42, 90 42"/>
      </g>
      <ellipse cx="80" cy="78" rx="7" ry="5.5" fill="none" stroke="#2e2318" stroke-width="1" opacity="0.5"/>
      <ellipse cx="80" cy="78" rx="7" ry="5.5" fill="#5c4a35" opacity="0.2"/>
      <circle cx="74" cy="38" r="2.5" fill="#4a3a28" opacity="0.3"/>
      <circle cx="34" cy="78" r="2" fill="#4a3a28" opacity="0.22"/>
    </g>
  </svg>
  <div style="display:flex; flex-direction:column; gap:1px;">
    <h1 style="margin:0;">Phobos</h1>
    <span class="sub" title="Peptide-level Heuristics for Bottom-up Omics Suite">Peptide-level Heuristics for Bottom-up Omics Suite · Peaks</span>
    <span class="credit">Proteogen — Université de Caen</span>
  </div>
  <span class="badge" id="badge-mode">Peptide level</span>
  <div class="toggle-wrap">
    <span class="lbl">View:</span>
    <div class="toggle" id="level-toggle">
      <button data-level="peptide" class="active">Peptide</button>
      <button data-level="protein">Protein</button>
    </div>
  </div>
</header>

<nav id="nav">
  <button data-tab="overview" class="active">Overview</button>
  <button data-tab="volcano">Volcano</button>
  <button data-tab="table">DE Table</button>
  <button data-tab="intensity">Intensity</button>
  <button data-tab="pca">PCA / UMAP</button>
  <button data-tab="qc">QC</button>
  <button data-tab="facets">Peptide facets</button>
</nav>

<main>
  <!-- OVERVIEW -->
  <section class="tab-page active" data-page="overview">
    <div class="grid cards" id="cards"></div>
    <div class="panel">
      <h3>Differential peptides per contrast</h3>
      <canvas id="chart-overview" height="110"></canvas>
      <div class="legend">
        <span><i style="background:#E74C3C"></i>Up-regulated</span>
        <span><i style="background:#3498DB"></i>Down-regulated</span>
      </div>
    </div>
  </section>

  <!-- VOLCANO -->
  <section class="tab-page" data-page="volcano">
    <div class="panel">
      <h3>Volcano plot</h3>
      <div class="ctrl">
        <label class="lbl">Contrast</label>
        <select id="volc-contrast"></select>
        <label class="lbl">Colour by</label>
        <select id="volc-colorby">
          <option value="status">Regulation</option>
          <option value="charge">Charge state</option>
          <option value="mod">Modification</option>
        </select>
        <label class="lbl"><input type="checkbox" id="volc-robust"/> Robust only</label>
        <button id="volc-flip" class="flip-btn" title="Swap numerator/denominator (A vs B → B vs A)">⇄ Invert</button>
      </div>
      <canvas id="chart-volcano" height="150"></canvas>
      <div class="legend" id="volc-legend"></div>
      <p class="note" id="volc-note"></p>
    </div>
  </section>

  <!-- DE TABLE -->
  <section class="tab-page" data-page="table">
    <div class="panel">
      <h3>Differential expression table</h3>
      <div class="ctrl">
        <label class="lbl">Contrast</label>
        <select id="tbl-contrast"></select>
        <label class="lbl">Filter</label>
        <select id="tbl-filter">
          <option value="all">All</option>
          <option value="sig">Significant only</option>
          <option value="up">Up only</option>
          <option value="down">Down only</option>
        </select>
        <input type="text" id="tbl-search" placeholder="search sequence / accession…"/>
        <span class="badge" id="tbl-count"></span>
      </div>
      <div class="tbl-wrap">
        <table id="de-table"><thead></thead><tbody></tbody></table>
      </div>
    </div>
  </section>

  <!-- INTENSITY -->
  <section class="tab-page" data-page="intensity">
    <div class="panel">
      <h3>Intensity across samples</h3>
      <div class="ctrl">
        <input type="text" id="int-search"
               placeholder="type a sequence or accession then Enter…"/>
        <span class="badge" id="int-label">—</span>
      </div>
      <canvas id="chart-intensity" height="120"></canvas>
      <p class="note">Per-sample log2 intensity, grouped by condition.
         Peptide view shows the selected feature; protein view shows the
         median rollup.</p>
    </div>
  </section>

  <!-- PCA / UMAP -->
  <section class="tab-page" data-page="pca">
    <div class="row">
      <div class="panel" style="flex:1; min-width:340px;">
        <h3>Sample PCA</h3>
        <canvas id="chart-pca" height="150"></canvas>
      </div>
      <div class="panel" style="flex:1; min-width:340px;">
        <h3>UMAP projection</h3>
        <canvas id="chart-umap" height="150"></canvas>
      </div>
    </div>
    <p class="note">PCA computed live from the current view's intensity matrix
       (peptide or protein). UMAP coordinates are precomputed (peptide level).</p>
  </section>

  <!-- QC -->
  <section class="tab-page" data-page="qc">
    <div class="grid cards" id="qc-cards"></div>
    <div class="row">
      <div class="panel" style="flex:1; min-width:340px;">
        <h3>Peptides detected per sample</h3>
        <canvas id="qc-detected" height="150"></canvas>
      </div>
      <div class="panel" style="flex:1; min-width:340px;">
        <h3>Missing values per sample</h3>
        <canvas id="qc-missing" height="150"></canvas>
        <p class="note">Real gaps <b>before</b> imputation (not detected in the
           sample). High, condition-correlated missingness is the MNAR signature
           that justifies QRILC.</p>
      </div>
    </div>
    <div class="panel">
      <h3>Intensity distribution per sample (log2)</h3>
      <canvas id="qc-boxplot" height="120"></canvas>
      <p class="note">Box = median ± IQR, whiskers = 5th–95th percentile.
         Consistent medians across samples indicate good normalisation.</p>
    </div>
    <div class="panel">
      <h3>RLE — Relative Log Expression</h3>
      <canvas id="qc-rle" height="120"></canvas>
      <p class="note">Per-peptide deviation from its median across samples.
         Boxes centred on 0 with similar spread = no residual technical bias.
         A shifted or wider box flags a problematic run.</p>
    </div>
  </section>

  <!-- FACETS -->
  <section class="tab-page" data-page="facets">
    <div class="row">
      <div class="panel" style="flex:1; min-width:300px;">
        <h3>Charge state distribution</h3>
        <canvas id="chart-charge" height="160"></canvas>
      </div>
      <div class="panel" style="flex:1; min-width:300px;">
        <h3>Modified vs unmodified</h3>
        <canvas id="chart-mod" height="160"></canvas>
      </div>
    </div>
    <div class="panel">
      <h3>Sequence length distribution</h3>
      <canvas id="chart-length" height="120"></canvas>
      <p class="note">Peptide-specific views — not available in protein mode.</p>
    </div>
  </section>
</main>

<script>
const DATA = {payload_json};
{_DASHBOARD_JS}
</script>
</body>
</html>"""


# ------------------------------------------------------------------------------
# Dashboard JavaScript (kept separate for readability; injected above)
# ------------------------------------------------------------------------------
_DASHBOARD_JS = r"""
// ===== State =====
let LEVEL = "peptide";
let charts = {};
const T = DATA.thresholds;
const COND_COLORS = {};
(function(){
  const palette = ["#5B8FF9","#F6BD16","#5AD8A6","#E8684A","#9270CA",
                   "#FF9D4D","#269A99","#FF99C3","#6DC8EC","#FF6B3B"];
  (DATA.meta.conditions||[]).forEach((c,i)=>COND_COLORS[c]=palette[i%palette.length]);
})();

function condOf(label){
  const s = (DATA.samples||[]).find(x=>x.label===label);
  return s ? s.condition : "?";
}

// ===== Tabs & toggle =====
document.querySelectorAll('#nav button').forEach(b=>{
  b.onclick=()=>{
    document.querySelectorAll('#nav button').forEach(x=>x.classList.remove('active'));
    b.classList.add('active');
    const tab=b.dataset.tab;
    document.querySelectorAll('.tab-page').forEach(p=>
      p.classList.toggle('active', p.dataset.page===tab));
    renderTab(tab);
  };
});
document.querySelectorAll('#level-toggle button').forEach(b=>{
  b.onclick=()=>{
    document.querySelectorAll('#level-toggle button').forEach(x=>x.classList.remove('active'));
    b.classList.add('active');
    LEVEL=b.dataset.level;
    document.getElementById('badge-mode').textContent=
      LEVEL==='peptide'?'Peptide level':'Protein level';
    // peptide-only tabs disabled in protein mode
    const facetBtn=document.querySelector('#nav button[data-tab="facets"]');
    facetBtn.style.opacity = LEVEL==='protein'?0.4:1;
    const active=document.querySelector('#nav button.active').dataset.tab;
    renderTab(active);
  };
});

// ===== Helpers =====
function pcol(c){ return T.use_padj ? c+"_p.adj" : c+"_p.val"; }
function statusOf(row, c){
  const d=row[c+"_diff"], p=row[pcol(c)];
  if(d==null||p==null) return "ns";
  if(d> T.lfc_min && p<T.p_thresh) return "up";
  if(d<-T.lfc_min && p<T.p_thresh) return "down";
  return "ns";
}
function destroy(id){ if(charts[id]){ charts[id].destroy(); delete charts[id]; } }

// ===== Overview =====
function renderOverview(){
  const cards=document.getElementById('cards');
  const m=DATA.meta;
  const items = LEVEL==='peptide' ? [
    [m.n_peptides_total||0,"Peptides"],
    [m.n_proteins_total||0,"Proteins"],
    [m.n_samples||0,"Samples"],
    [m.n_contrasts||0,"Contrasts"],
    [m.n_sig_anova||0,"ANOVA sig."],
  ] : [
    [m.n_proteins_total||(DATA.proteins||[]).length,"Proteins"],
    [m.n_samples||0,"Samples"],
    [m.n_contrasts||0,"Contrasts"],
  ];
  cards.innerHTML=items.map(([v,l])=>
    `<div class="card"><div class="v">${v}</div><div class="l">${l}</div></div>`).join('');

  const ups=[], downs=[];
  DATA.contrasts.forEach(c=>{
    let u=0,d=0;
    DATA.peptides.forEach(r=>{ const s=statusOf(r,c); if(s==='up')u++; else if(s==='down')d++; });
    ups.push(u); downs.push(d);
  });
  destroy('overview');
  charts.overview=new Chart(document.getElementById('chart-overview'),{
    type:'bar',
    data:{labels:DATA.contrasts.map(c=>c.replace(/_vs_/,' vs ')),
      datasets:[
        {label:'Up',data:ups,backgroundColor:'#E74C3C'},
        {label:'Down',data:downs,backgroundColor:'#3498DB'}]},
    options:{responsive:true,scales:{x:{stacked:true,ticks:{color:'#8b97a6'}},
      y:{stacked:true,ticks:{color:'#8b97a6'}}},
      plugins:{legend:{labels:{color:'#e6edf3'}}}}
  });
}

// ===== Volcano =====
function renderVolcano(){
  const panel = document.querySelector('[data-page="volcano"] .panel');
  const overlay = document.getElementById('volc-protein-msg');
  const innerEls = panel ? panel.querySelectorAll('.ctrl, canvas, .legend, #volc-note') : [];
  if(LEVEL==='protein'){
    innerEls.forEach(e=>e.classList.add('hide'));
    if(!document.getElementById('volc-protein-msg')){
      const div=document.createElement('div');
      div.id='volc-protein-msg'; div.className='note';
      div.style.lineHeight='1.7';
      div.innerHTML='<b>Differential statistics are computed at the PEPTIDE '+
        'level only.</b><br>Phobos fits the limma model on individual peptide '+
        'features. The protein view aggregates intensities (median rollup) for '+
        '<i>exploratory</i> visualisation — it carries no per-protein p-value or '+
        'fold-change, so a protein-level volcano would be statistically '+
        'misleading.<br><br>Switch to <b>Peptide level</b> (top-right) to see the '+
        'volcano, or use the Intensity / PCA tabs for the protein rollup.';
      panel.appendChild(div);
    } else { overlay.classList.remove('hide'); }
    return;
  }
  if(overlay) overlay.classList.add('hide');
  innerEls.forEach(e=>e.classList.remove('hide'));
  const sel=document.getElementById('volc-contrast');
  const flipBtn=document.getElementById('volc-flip');
  if(!sel.options.length){
    DATA.contrasts.forEach(c=>sel.add(new Option(c.replace(/_vs_/,' vs '),c)));
    sel.onchange=renderVolcano;
    document.getElementById('volc-colorby').onchange=renderVolcano;
    document.getElementById('volc-robust').onchange=renderVolcano;
    flipBtn.onclick=()=>{ window._volcFlip=!window._volcFlip;
                          flipBtn.classList.toggle('active', window._volcFlip);
                          renderVolcano(); };
  }
  const c=sel.value||DATA.contrasts[0];
  const colorBy=document.getElementById('volc-colorby').value;
  const robustOnly=document.getElementById('volc-robust').checked;
  const flipped=!!window._volcFlip;
  const sign=flipped?-1:1;
  // Robustness threshold = 80% of iterations (if available)
  const robMax=DATA.meta.n_iter||100;
  const robCut=Math.max(1, Math.round(robMax*0.8));
  const pts={};
  const push=(k,col,x,y,label)=>{ (pts[k]=pts[k]||{c:col,d:[]}).d.push({x,y,label}); };

  DATA.peptides.forEach(r=>{
    let d=r[c+"_diff"], p=r[pcol(c)];
    if(d==null||p==null) return;
    if(robustOnly){
      const rob=r["Robustness_Score_"+c];
      if(rob==null || rob<robCut) return;   // garder seulement les robustes
    }
    d=d*sign;                                // inversion du fold-change
    const y=-Math.log10(Math.max(p,1e-10));
    let key,col;
    if(colorBy==='status'){
      // statut recalculé sur le FC inversé (up/down échangés)
      let s='ns';
      if(d> T.lfc_min && p<T.p_thresh) s='up';
      else if(d<-T.lfc_min && p<T.p_thresh) s='down';
      key=s; col={up:'#E74C3C',down:'#3498DB',ns:'#4a5568'}[s];
    } else if(colorBy==='charge'){
      key='z='+(r.Charge!=null?r.Charge:'?');
      const pal=['#5B8FF9','#F6BD16','#5AD8A6','#E8684A','#9270CA'];
      col=pal[(String(r.Charge).charCodeAt(0))%pal.length];
    } else {
      const mod=(r.Modifications||'').toString().trim();
      key= mod.length? 'Modified':'Unmodified';
      col= mod.length? '#E8684A':'#5AD8A6';
    }
    push(key,col,d,y,r.Sequence||r.peptide_id);
  });

  destroy('volcano');
  charts.volcano=new Chart(document.getElementById('chart-volcano'),{
    type:'scatter',
    data:{datasets:Object.entries(pts).map(([k,v])=>({
      label:k,data:v.d,backgroundColor:v.c,pointRadius:3,pointHoverRadius:5}))},
    options:{responsive:true,
      plugins:{legend:{labels:{color:'#e6edf3'}},
        tooltip:{callbacks:{label:(ctx)=>{
          const p=ctx.raw; return `${p.label}  (log2FC=${p.x.toFixed(2)}, -log10p=${p.y.toFixed(2)})`;}}}},
      scales:{x:{title:{display:true,text:'log2 Fold Change',color:'#8b97a6'},
          ticks:{color:'#8b97a6'},grid:{color:'#2d3742'}},
        y:{title:{display:true,text:'-log10(p)',color:'#8b97a6'},
          ticks:{color:'#8b97a6'},grid:{color:'#2d3742'}}}}
  });
  // Sens du contraste affiché (inversé si flip actif)
  const parts=c.split("_vs_");
  const shown = flipped ? `${parts[1]} vs ${parts[0]}` : `${parts[0]} vs ${parts[1]}`;
  document.getElementById('volc-note').textContent=
    `Contrast: ${shown}${flipped?'  (inverted)':''}  •  `+
    `Thresholds: ${T.use_padj?'p.adj':'p.value'} < ${T.p_thresh}, `+
    `|log2FC| > ${T.lfc_min.toFixed(2)} (ratio ${T.ratio_min}). `+
    `Positive log2FC = enriched in ${flipped?parts[1]:parts[0]}.`;
}

// ===== DE Table =====
function renderTable(){
  const panel = document.querySelector('[data-page="table"] .panel');
  const innerEls = panel ? panel.querySelectorAll('.ctrl, .tbl-wrap') : [];
  let overlay = document.getElementById('tbl-protein-msg');
  if(LEVEL==='protein'){
    innerEls.forEach(e=>e.classList.add('hide'));
    if(!overlay){
      overlay=document.createElement('div');
      overlay.id='tbl-protein-msg'; overlay.className='note';
      overlay.style.lineHeight='1.7';
      overlay.innerHTML='<b>The DE table lists PEPTIDE-level statistics.</b><br>'+
        'No per-protein differential test is performed (intensities are only '+
        'median-aggregated for visualisation). Switch to <b>Peptide level</b> '+
        'to browse fold-changes, p-values, Pi-scores and robustness per peptide.';
      panel.appendChild(overlay);
    } else { overlay.classList.remove('hide'); }
    return;
  }
  if(overlay) overlay.classList.add('hide');
  innerEls.forEach(e=>e.classList.remove('hide'));
  const sel=document.getElementById('tbl-contrast');
  if(!sel.options.length){
    DATA.contrasts.forEach(c=>sel.add(new Option(c.replace(/_vs_/,' vs '),c)));
    sel.onchange=renderTable;
    document.getElementById('tbl-filter').onchange=renderTable;
    document.getElementById('tbl-search').oninput=renderTable;
  }
  const c=sel.value||DATA.contrasts[0];
  const filt=document.getElementById('tbl-filter').value;
  const q=document.getElementById('tbl-search').value.toLowerCase();

  let rows=DATA.peptides.map(r=>({
    seq:r.Sequence||r.peptide_id, acc:r.Accession, gene:r.Gene,
    desc:r.Description, z:r.Charge,
    mod:r.Modifications, diff:r[c+"_diff"], p:r[pcol(c)],
    pi:r["Pi_Score_"+c], rob:r["Robustness_Score_"+c],
    status:statusOf(r,c)
  })).filter(r=>r.diff!=null);

  if(filt==='sig') rows=rows.filter(r=>r.status!=='ns');
  else if(filt==='up') rows=rows.filter(r=>r.status==='up');
  else if(filt==='down') rows=rows.filter(r=>r.status==='down');
  if(q) rows=rows.filter(r=>(r.seq||'').toLowerCase().includes(q)||
                            (r.acc||'').toString().toLowerCase().includes(q)||
                            (r.gene||'').toString().toLowerCase().includes(q)||
                            (r.desc||'').toString().toLowerCase().includes(q));
  rows.sort((a,b)=>(b.pi||0)-(a.pi||0));

  document.getElementById('tbl-count').textContent=`${rows.length} peptides`;
  const thead=document.querySelector('#de-table thead');
  const tbody=document.querySelector('#de-table tbody');
  thead.innerHTML=`<tr><th>Sequence</th><th>Accession</th><th>Gene</th>
    <th>Description</th><th>z</th>
    <th>Mod</th><th>log2FC</th><th>${T.use_padj?'p.adj':'p.val'}</th>
    <th>Pi</th><th>Robust</th><th>Status</th></tr>`;
  tbody.innerHTML=rows.slice(0,500).map(r=>{
    const pill=`<span class="pill ${r.status}">${r.status.toUpperCase()}</span>`;
    return `<tr><td class="seq">${r.seq||''}</td><td>${r.acc||''}</td>
      <td>${r.gene||''}</td><td class="desc">${r.desc||''}</td>
      <td>${r.z!=null?r.z:''}</td><td>${r.mod||''}</td>
      <td>${r.diff!=null?r.diff.toFixed(3):''}</td>
      <td>${r.p!=null?r.p.toExponential(2):''}</td>
      <td>${r.pi!=null?r.pi.toFixed(2):''}</td>
      <td>${r.rob!=null?r.rob:''}</td><td>${pill}</td></tr>`;
  }).join('');
}

// ===== Intensity =====
function renderIntensity(){
  const inp=document.getElementById('int-search');
  inp.onkeydown=(e)=>{ if(e.key==='Enter') drawIntensity(inp.value.trim()); };
  // default: first peptide / protein
  if(!charts._int_init){
    charts._int_init=true;
    const seed = LEVEL==='peptide'
      ? (DATA.peptides[0]?.Sequence||'')
      : (DATA.proteins[0]?.Accession||'');
    drawIntensity(seed);
  }
}
function drawIntensity(query){
  let vals=null, label='—', cols=null;
  if(LEVEL==='peptide'){
    cols=DATA.meta.sample_cols||[];
    const hit=DATA.peptides.find(r=>
      (r.Sequence||'').toLowerCase()===query.toLowerCase() ||
      (r.peptide_id||'').toLowerCase().includes(query.toLowerCase()));
    if(hit){ vals=DATA.peptide_intensities[hit.peptide_id]; label=hit.Sequence||hit.peptide_id; }
  } else {
    cols=DATA.meta.protein_sample_cols||[];
    const hit=DATA.proteins.find(p=>
      (p.Accession||'').toLowerCase()===query.toLowerCase() ||
      (p.Accession||'').toLowerCase().includes(query.toLowerCase()));
    if(hit){ vals=hit.intensities; label=hit.Accession; }
  }
  document.getElementById('int-label').textContent=label;
  destroy('intensity');
  if(!vals||!cols) return;
  charts.intensity=new Chart(document.getElementById('chart-intensity'),{
    type:'bar',
    data:{labels:cols,datasets:[{label:label,data:vals,
      backgroundColor:cols.map(c=>COND_COLORS[condOf(c)]||'#5B8FF9')}]},
    options:{responsive:true,
      plugins:{legend:{display:false}},
      scales:{x:{ticks:{color:'#8b97a6',maxRotation:90}},
        y:{title:{display:true,text:'log2 intensity',color:'#8b97a6'},
          ticks:{color:'#8b97a6'}}}}
  });
}

// ===== PCA (live, 2-PC power-iteration on covariance) =====
function computePCA(matrix){
  // matrix: rows=features, cols=samples. Returns Nx2 sample coords.
  const nF=matrix.length, nS=matrix[0].length;
  // center per feature, build sample covariance (S x S)
  const cen=matrix.map(r=>{const m=r.reduce((a,b)=>a+b,0)/nS; return r.map(v=>v-m);});
  const cov=Array.from({length:nS},()=>new Array(nS).fill(0));
  for(let f=0;f<nF;f++) for(let i=0;i<nS;i++){const vi=cen[f][i];
    for(let j=0;j<nS;j++) cov[i][j]+=vi*cen[f][j];}
  for(let i=0;i<nS;i++) for(let j=0;j<nS;j++) cov[i][j]/=Math.max(nF-1,1);
  const pc=(A,exclude)=>{
    let v=new Array(nS).fill(0).map(()=>Math.random());
    for(let it=0;it<60;it++){
      let w=new Array(nS).fill(0);
      for(let i=0;i<nS;i++) for(let j=0;j<nS;j++) w[i]+=A[i][j]*v[j];
      if(exclude){const d=w.reduce((a,b,i)=>a+b*exclude[i],0);
        for(let i=0;i<nS;i++) w[i]-=d*exclude[i];}
      const n=Math.hypot(...w)||1; v=w.map(x=>x/n);
    }
    return v;
  };
  const v1=pc(cov,null), v2=pc(cov,v1);
  const proj=(v)=>{const out=new Array(nS).fill(0);
    for(let i=0;i<nS;i++){let s=0; for(let j=0;j<nS;j++) s+=cov[i][j]*v[j]; out[i]=s;}
    return out;};
  // sample coordinates = projection onto v1,v2 (use cov·v as score proxy)
  const c1=proj(v1), c2=proj(v2);
  return c1.map((_,i)=>[c1[i],c2[i]]);
}
function renderPCA(){
  // Plugin inline : dessine le label de chaque point en permanence
  const labelPlugin={
    id:'pointLabels',
    afterDatasetsDraw(chart){
      const {ctx}=chart;
      ctx.save();
      ctx.font='10px Segoe UI, sans-serif';
      ctx.fillStyle='#e6edf3';
      ctx.textAlign='left'; ctx.textBaseline='middle';
      chart.data.datasets.forEach((ds,di)=>{
        const meta=chart.getDatasetMeta(di);
        if(meta.hidden) return;
        meta.data.forEach((pt,i)=>{
          const lbl=ds.data[i] && ds.data[i].label;
          if(lbl) ctx.fillText(lbl, pt.x+8, pt.y);
        });
      });
      ctx.restore();
    }
  };
  let cols, recs;
  if(LEVEL==='peptide'){
    cols=DATA.meta.sample_cols||[];
    recs=Object.values(DATA.peptide_intensities||{});
  } else {
    cols=DATA.meta.protein_sample_cols||[];
    recs=(DATA.proteins||[]).map(p=>p.intensities);
  }
  // keep complete rows only
  const mat=recs.filter(r=>r&&r.every(v=>v!=null));
  destroy('pca');
  if(mat.length>2 && cols.length>2){
    const coords=computePCA(mat);
    const groups={};
    cols.forEach((c,i)=>{const cond=condOf(c);
      (groups[cond]=groups[cond]||[]).push({x:coords[i][0],y:coords[i][1],label:c});});
    charts.pca=new Chart(document.getElementById('chart-pca'),{
      type:'scatter',
      data:{datasets:Object.entries(groups).map(([k,v])=>({
        label:k,data:v,backgroundColor:COND_COLORS[k]||'#5B8FF9',
        pointRadius:7,pointHoverRadius:9}))},
      plugins:[labelPlugin],
      options:{responsive:true,
        layout:{padding:{right:40}},
        plugins:{legend:{labels:{color:'#e6edf3'}},
          tooltip:{callbacks:{label:(c)=>c.raw.label}}},
        scales:{x:{title:{display:true,text:'PC1',color:'#8b97a6'},
            ticks:{color:'#8b97a6'},grid:{color:'#2d3742'}},
          y:{title:{display:true,text:'PC2',color:'#8b97a6'},
            ticks:{color:'#8b97a6'},grid:{color:'#2d3742'}}}}
    });
  }
  // UMAP (precomputed)
  destroy('umap');
  if((DATA.umap||[]).length){
    const groups={};
    DATA.umap.forEach(u=>{(groups[u.condition]=groups[u.condition]||[]).push(
      {x:u.UMAP1,y:u.UMAP2,label:u.label});});
    charts.umap=new Chart(document.getElementById('chart-umap'),{
      type:'scatter',
      plugins:[labelPlugin],
      data:{datasets:Object.entries(groups).map(([k,v])=>({
        label:k,data:v,backgroundColor:COND_COLORS[k]||'#5B8FF9',
        pointRadius:7,pointHoverRadius:9}))},
      options:{responsive:true,
        layout:{padding:{right:40}},
        plugins:{legend:{labels:{color:'#e6edf3'}},
          tooltip:{callbacks:{label:(c)=>c.raw.label}}},
        scales:{x:{title:{display:true,text:'UMAP1',color:'#8b97a6'},
            ticks:{color:'#8b97a6'},grid:{color:'#2d3742'}},
          y:{title:{display:true,text:'UMAP2',color:'#8b97a6'},
            ticks:{color:'#8b97a6'},grid:{color:'#2d3742'}}}}
    });
  }
}

// ===== Facets (peptide only) =====
function renderFacets(){
  // charge
  const chargeCounts={};
  DATA.peptides.forEach(r=>{const z=r.Charge!=null?('z='+r.Charge):'z=?';
    chargeCounts[z]=(chargeCounts[z]||0)+1;});
  destroy('charge');
  charts.charge=new Chart(document.getElementById('chart-charge'),{
    type:'bar',
    data:{labels:Object.keys(chargeCounts),
      datasets:[{data:Object.values(chargeCounts),backgroundColor:'#5B8FF9'}]},
    options:{plugins:{legend:{display:false}},responsive:true,
      scales:{x:{ticks:{color:'#8b97a6'}},y:{ticks:{color:'#8b97a6'}}}}
  });
  // modifications
  let mod=0,nomod=0;
  DATA.peptides.forEach(r=>{const m=(r.Modifications||'').toString().trim();
    if(m.length) mod++; else nomod++;});
  destroy('mod');
  charts.mod=new Chart(document.getElementById('chart-mod'),{
    type:'doughnut',
    data:{labels:['Modified','Unmodified'],
      datasets:[{data:[mod,nomod],backgroundColor:['#E8684A','#5AD8A6']}]},
    options:{responsive:true,plugins:{legend:{labels:{color:'#e6edf3'}}}}
  });
  // sequence length
  const lenCounts={};
  DATA.peptides.forEach(r=>{const L=(r.Sequence||'').length;
    if(L>0) lenCounts[L]=(lenCounts[L]||0)+1;});
  const lens=Object.keys(lenCounts).map(Number).sort((a,b)=>a-b);
  destroy('length');
  charts.length=new Chart(document.getElementById('chart-length'),{
    type:'bar',
    data:{labels:lens,datasets:[{data:lens.map(l=>lenCounts[l]),
      backgroundColor:'#9270CA'}]},
    options:{plugins:{legend:{display:false}},responsive:true,
      scales:{x:{title:{display:true,text:'Peptide length (aa)',color:'#8b97a6'},
          ticks:{color:'#8b97a6'}},y:{ticks:{color:'#8b97a6'}}}}
  });
}

// ===== QC (native, computed client-side from peptide intensities) =====
function _quantile(sorted, q){
  if(!sorted.length) return null;
  const pos=(sorted.length-1)*q, base=Math.floor(pos), rest=pos-base;
  return sorted[base+1]!==undefined ? sorted[base]+rest*(sorted[base+1]-sorted[base])
                                    : sorted[base];
}
function renderQC(){
  const cols=DATA.meta.sample_cols||[];
  const recs=Object.values(DATA.peptide_intensities||{});
  if(!cols.length || !recs.length){
    document.querySelector('[data-page="qc"]').querySelectorAll('canvas')
      .forEach(c=>{const ctx=c.getContext('2d'); ctx.clearRect(0,0,c.width,c.height);});
    return;
  }
  const nS=cols.length;
  const colColors=cols.map(c=>COND_COLORS[condOf(c)]||'#5B8FF9');

  // Per-sample detected/missing : utiliser la matrice PRÉ-imputation (vrais trous)
  // si disponible, sinon retomber sur la matrice imputée (0 trou).
  const pre=DATA.qc_preimpute;
  let detected, missing, missFromPre=false;
  if(pre && pre.sample_cols && pre.sample_cols.length){
    detected=cols.map(c=>pre.detected[c]!=null?pre.detected[c]:0);
    missing =cols.map(c=>pre.missing[c]!=null?pre.missing[c]:0);
    missFromPre=true;
  } else {
    detected=new Array(nS).fill(0);
    missing =new Array(nS).fill(0);
    recs.forEach(row=>{ for(let i=0;i<nS;i++){
      if(row[i]==null) missing[i]++; else detected[i]++; }});
  }
  // Listes de valeurs (matrice imputée) pour les distributions/boxplots
  const valsBy=Array.from({length:nS},()=>[]);
  recs.forEach(row=>{ for(let i=0;i<nS;i++){ if(row[i]!=null) valsBy[i].push(row[i]); }});

  // Cards
  const totalPep = (pre && pre.n_peptides) ? pre.n_peptides : recs.length;
  const avgDet=Math.round(detected.reduce((a,b)=>a+b,0)/nS);
  const totMiss=missing.reduce((a,b)=>a+b,0);
  const pctMiss=(100*totMiss/(nS*totalPep)).toFixed(1);
  document.getElementById('qc-cards').innerHTML=[
    [totalPep,'Peptides (filtered)'],
    [nS,'Samples'],
    [avgDet,'Avg detected / sample'],
    [pctMiss+'%', missFromPre?'Missing (pre-imputation)':'Missing values'],
  ].map(([v,l])=>`<div class="card"><div class="v">${v}</div><div class="l">${l}</div></div>`).join('');

  // 1) Detected per sample
  destroy('qcdet');
  charts.qcdet=new Chart(document.getElementById('qc-detected'),{
    type:'bar',
    data:{labels:cols,datasets:[{data:detected,backgroundColor:colColors}]},
    options:{responsive:true,plugins:{legend:{display:false}},
      scales:{x:{ticks:{color:'#8b97a6',maxRotation:90,minRotation:45}},
        y:{ticks:{color:'#8b97a6'},title:{display:true,text:'# peptides',color:'#8b97a6'}}}}
  });

  // 2) Missing per sample
  destroy('qcmiss');
  charts.qcmiss=new Chart(document.getElementById('qc-missing'),{
    type:'bar',
    data:{labels:cols,datasets:[{data:missing,backgroundColor:'#E8684A'}]},
    options:{responsive:true,plugins:{legend:{display:false}},
      scales:{x:{ticks:{color:'#8b97a6',maxRotation:90,minRotation:45}},
        y:{ticks:{color:'#8b97a6'},title:{display:true,text:'# missing',color:'#8b97a6'}}}}
  });

  // Box stats helper -> floating bars [q1,q3] + median/whisker overlay plugin
  const stats=valsBy.map(v=>{
    const s=v.slice().sort((a,b)=>a-b);
    return {q1:_quantile(s,.25),med:_quantile(s,.5),q3:_quantile(s,.75),
            lo:_quantile(s,.05),hi:_quantile(s,.95)};
  });
  const boxPlugin=(getStats)=>({
    id:'box'+Math.random(),
    afterDatasetsDraw(chart){
      const {ctx,scales:{x,y}}=chart, st=getStats();
      ctx.save(); ctx.strokeStyle='#e6edf3'; ctx.lineWidth=1.2;
      st.forEach((s,i)=>{
        if(s.med==null) return;
        const cx=x.getPixelForValue(i), w=Math.min(x.width/st.length*0.5,22);
        // median line
        ctx.beginPath(); ctx.moveTo(cx-w,y.getPixelForValue(s.med));
        ctx.lineTo(cx+w,y.getPixelForValue(s.med)); ctx.stroke();
        // whiskers
        ctx.beginPath();
        ctx.moveTo(cx,y.getPixelForValue(s.q3)); ctx.lineTo(cx,y.getPixelForValue(s.hi));
        ctx.moveTo(cx,y.getPixelForValue(s.q1)); ctx.lineTo(cx,y.getPixelForValue(s.lo));
        ctx.moveTo(cx-w*0.6,y.getPixelForValue(s.hi)); ctx.lineTo(cx+w*0.6,y.getPixelForValue(s.hi));
        ctx.moveTo(cx-w*0.6,y.getPixelForValue(s.lo)); ctx.lineTo(cx+w*0.6,y.getPixelForValue(s.lo));
        ctx.stroke();
      });
      ctx.restore();
    }
  });

  // 3) Intensity boxplot (floating bar q1..q3 + plugin overlay)
  destroy('qcbox');
  charts.qcbox=new Chart(document.getElementById('qc-boxplot'),{
    type:'bar',
    data:{labels:cols,datasets:[{
      data:stats.map(s=>[s.q1,s.q3]),
      backgroundColor:colColors.map(c=>c+'aa'),borderColor:colColors,borderWidth:1}]},
    plugins:[boxPlugin(()=>stats)],
    options:{responsive:true,plugins:{legend:{display:false},
      tooltip:{callbacks:{label:(ctx)=>{const s=stats[ctx.dataIndex];
        return `med=${s.med.toFixed(2)}  IQR=[${s.q1.toFixed(2)}, ${s.q3.toFixed(2)}]`;}}}},
      scales:{x:{ticks:{color:'#8b97a6',maxRotation:90,minRotation:45}},
        y:{ticks:{color:'#8b97a6'},title:{display:true,text:'log2 intensity',color:'#8b97a6'}}}}
  });

  // 4) RLE: per-peptide deviation from its across-sample median
  const rleBy=Array.from({length:nS},()=>[]);
  recs.forEach(row=>{
    const present=row.filter(v=>v!=null);
    if(present.length<2) return;
    const s=present.slice().sort((a,b)=>a-b);
    const med=_quantile(s,.5);
    for(let i=0;i<nS;i++) if(row[i]!=null) rleBy[i].push(row[i]-med);
  });
  const rleStats=rleBy.map(v=>{
    const s=v.slice().sort((a,b)=>a-b);
    return {q1:_quantile(s,.25),med:_quantile(s,.5),q3:_quantile(s,.75),
            lo:_quantile(s,.05),hi:_quantile(s,.95)};
  });
  destroy('qcrle');
  charts.qcrle=new Chart(document.getElementById('qc-rle'),{
    type:'bar',
    data:{labels:cols,datasets:[{
      data:rleStats.map(s=>[s.q1,s.q3]),
      backgroundColor:colColors.map(c=>c+'aa'),borderColor:colColors,borderWidth:1}]},
    plugins:[boxPlugin(()=>rleStats)],
    options:{responsive:true,plugins:{legend:{display:false},
      tooltip:{callbacks:{label:(ctx)=>{const s=rleStats[ctx.dataIndex];
        return `median RLE=${s.med.toFixed(3)}`;}}}},
      scales:{x:{ticks:{color:'#8b97a6',maxRotation:90,minRotation:45}},
        y:{ticks:{color:'#8b97a6'},title:{display:true,text:'RLE (log2)',color:'#8b97a6'},
           suggestedMin:-2,suggestedMax:2}}}
  });
}

// ===== Router =====
function renderTab(tab){
  if(tab==='overview') renderOverview();
  else if(tab==='volcano') renderVolcano();
  else if(tab==='table') renderTable();
  else if(tab==='intensity'){ charts._int_init=false; renderIntensity(); }
  else if(tab==='pca') renderPCA();
  else if(tab==='qc') renderQC();
  else if(tab==='facets'){
    if(LEVEL==='protein'){
      document.querySelector('[data-page="facets"]').innerHTML=
        '<div class="panel"><p class="note">Peptide facets are not available '+
        'in protein view. Switch to Peptide level.</p></div>';
    } else renderFacets();
  }
}

// init
renderOverview();
"""


# ------------------------------------------------------------------------------
# Public entry point
# ------------------------------------------------------------------------------

def build_dashboard(xlsx_path: str, design: pd.DataFrame,
                    out_html: str, params: dict = None) -> "str | None":
    """
    Generate the self-contained Phobos peptide dashboard.

    Parameters
    ----------
    xlsx_path : path to Phobos_Peptide_Results.xlsx
    design    : experimental design DataFrame (label, condition, ...)
    out_html  : output .html path
    params    : pipeline params dict (thresholds). Falls back to defaults.

    Returns the output path, or None on failure.
    """
    if params is None:
        params = {}
    if not os.path.exists(xlsx_path):
        print(f"  [WARN] Excel not found ({xlsx_path}) — dashboard skipped.")
        return None

    try:
        payload = _extract_payload(xlsx_path, design, params)
        payload_json = json.dumps(payload, ensure_ascii=False,
                                  separators=(",", ":"))
        base_dir = os.path.dirname(os.path.abspath(xlsx_path))
        chart_tag = _chartjs_tag(base_dir if base_dir else ".")
        # also look next to this script for the local chart.js
        if _CHARTJS_CDN in chart_tag:
            here = os.path.dirname(os.path.abspath(__file__))
            chart_tag = _chartjs_tag(here)
        html = _html_template(payload_json, chart_tag)
        with open(out_html, "w", encoding="utf-8") as f:
            f.write(html)
        size_kb = os.path.getsize(out_html) / 1024
        print(f"  [OK] Dashboard written: {out_html} ({size_kb:.0f} KB)")
        return out_html
    except Exception as e:
        import traceback
        print(f"  [WARN] Dashboard generation failed: {type(e).__name__}: {e}")
        traceback.print_exc()
        return None


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Build the Phobos peptide dashboard")
    ap.add_argument("xlsx", help="Phobos Excel results file")
    ap.add_argument("design", help="ExperimentalDesign.csv (sep=';')")
    ap.add_argument("-o", "--out", default="phobos_dashboard.html")
    args = ap.parse_args()
    design = pd.read_csv(args.design, sep=";")
    build_dashboard(args.xlsx, design, args.out, params={})
