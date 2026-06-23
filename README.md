# Phobos — Peptide-level Heuristics for Bottom-up Omics Suite

<p align="center">
  <img src="phobos_icon.svg" alt="Phobos" width="200"/>
</p>

**Companion to Deimos for DDA data.**  
Peptide-level differential analysis (sequence + charge + modifications)
from the PEAKS `protein-peptides.csv` export.

---

## Phobos vs Deimos

| Aspect                  | Deimos (DIA/DIA-NN)            | Phobos (DDA/PEAKS)                   |
|-------------------------|-------------------------------|--------------------------------------|
| Input                   | `report.pg_matrix.tsv`        | `protein-peptides.csv`               |
| Unit                    | Protein (LFQ)                 | Precursor peptide (seq+chg+mod)      |
| Intensity               | DIA-NN normalised LFQ         | Raw Area → median normalisation      |
| Filtering               | n_rep − thr                   | -10lgP score + n_rep − thr           |
| Normalisation           | None (DIA-NN cross-run)       | Median (per-run shift)               |
| Statistics              | limma + DEqMS (optional)      | limma eBayes (DEqMS N/A at pep. level) |
| Protein aggregation     | —                             | Median rollup, post-hoc (viz only)   |
| Dashboard               | Yes                           | Yes (peptide/protein toggle)         |
| WGCNA / GO              | Yes                           | No (v1)                              |

---

## Required files

| File                       | Required | Description                              |
|----------------------------|----------|------------------------------------------|
| `protein-peptides.csv`     | ✅       | PEAKS combined export (.csv or .xlsx)    |
| `ExperimentalDesign.csv`   | ✅       | `label;condition;replicate` (sep `;`)    |

**Column matching:** each design `label` must point to an `Area <label>`
column in the PEAKS export.  
Example: label `S01` ↔ column `Area S01`.

A helper script pre-fills the design from the Area columns:

```bash
python make_design_template.py protein-peptides.csv
# -> ExperimentalDesign_template.csv (review condition/replicate before use)
```

---

## Installation

Same dependencies as Deimos:

```bash
pip install pandas numpy scipy scikit-learn umap-learn matplotlib seaborn \
  adjustText openpyxl pillow PyComplexHeatmap pyyaml
```

Phobos shares `limma_ebayes.py` and `config.py` with Deimos — **keep all
files in the same folder**.

For a fully offline dashboard, place `chart.umd.min.js` next to the scripts:
<https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js>

---

## Usage

```bash
# Interactive mode
python phobos.py

# YAML config (no prompts)
python phobos.py --config phobos_config_example.yaml

# Interactive + save config for future runs
python phobos.py --save-config
```

---

## Repository structure

```
phobos/
├── phobos.py                   # Main orchestrator
├── build_dashboard_phobos.py   # Interactive HTML dashboard generator
├── make_design_template.py     # Design template helper (from Area columns)
├── phobos_config_example.yaml  # YAML config template
├── limma_ebayes.py             # Statistical engine (shared with Deimos)
├── config.py                   # YAML config management (shared with Deimos)
└── ExperimentalDesign.csv      # Experimental design (to adapt)
```

---

## Excel output — `PeptideAnalysis_Results.xlsx`

| Sheet                    | Content                                          |
|--------------------------|--------------------------------------------------|
| `Methods`                | Documented parameters and methods                |
| `raw_peptides`           | Raw PEAKS data (all columns)                     |
| `Log2_Impute`            | Normalised + imputed log2 matrix (QRILC)         |
| `QC`                     | Detection frequency, boxplots, RLE, missing val. |
| `PCA_UMAP`               | PCA ± ellipses, UMAP, Pearson correlation        |
| `UMAP`                   | Per-sample UMAP coordinates                      |
| `Scatter_Plots`          | Per-contrast scatter plots                       |
| `Differential_Expression`| logFC, p.val, p.adj, Pi-score, robustness        |
| `Volcano_Plots`          | Facet volcano + individual volcanos              |
| `UpSet_Intersections`    | DEP peptides shared across contrasts             |
| `Zscore_Heatmap`         | Z-score of ANOVA-significant peptides            |
| `ANOVA_Results`          | F-stat, p.value, p.adj, per-peptide cluster      |
| `ANOVA_Clusters`         | Cluster heatmap + expression violins             |
| `Protein_Aggregation`    | Peptide → protein median rollup (post-hoc)       |
| `PTM_<mod>_DE` / `_ANOVA`| Per-PTM differential + ANOVA (if PTM analysis on)|

---

## Protein descriptions from a FASTA

PEAKS exports often omit a protein `Description` column. If a UniProt-style
FASTA file is present in the input folder, Phobos **auto-detects** it (the
filename does not matter) and recovers descriptions by joining on the
Accession — no configuration needed.

- **Join key:** the first UniProt ID, before the first `|`, of the first
  protein in the Accession group (e.g. `Q503D3|..:A0A..|..` → `Q503D3`).
- **Multi-protein groups:** only the first protein's description is kept.
- Gene names absent from the export are filled from the FASTA too.
- Unmatched accessions keep an empty description (no error).

Standalone use:

```bash
python fasta_descriptions.py peptides.csv [Danio.fasta] -o peptides_with_desc.csv
```

(Pure Python — no Biopython required.)

---

## PTM-targeted sub-analyses

In addition to the global peptide analysis, Phobos can re-run the **full
pipeline** (volcano, ANOVA, heatmap, UpSet, scatter) on **PTM-filtered
subsets** — one independent analysis per selected modification.

**Selection** is interactive at launch (a menu lists each PTM with its peptide
count), or fixed via `ptm_keys` in the YAML config.

**Available PTMs:** Phosphorylation (+79.97), Oxidation (+15.99), Methylation
mono (+14.02) / di (+28.03) / tri (+42.05) — also combined as `Methyl_all` —
and Acetylation (+42.01).

**Combined detection:** a peptide carries a PTM if matched by **either** the
PTM text column (e.g. `Oxidation (M)`) **or** an inline delta mass in the
sequence (e.g. `M(+15.99)`), within a tight tolerance. Acetylation (+42.011)
and trimethylation (+42.047) differ by only 0.036 Da: the tolerance is tight
and the text annotation takes precedence to disambiguate.

Each analysed PTM adds two Excel sheets (`PTM_<mod>_DE`, `PTM_<mod>_ANOVA`).
A subset with fewer than 12 carrying peptides is skipped.

---

## Interactive dashboard — `phobos_dashboard.html`

Self-contained HTML (Chart.js inlined, works offline). Six tabs: Overview,
Volcano, DE Table, Intensity, PCA/UMAP, Peptide facets (charge, modifications,
length).

A **Peptide ⇄ Protein** toggle switches the exploratory views (Intensity, PCA)
to the protein median rollup. Differential statistics (volcano, DE table) are
computed at the **peptide level only** — in protein view these tabs display an
explicit notice rather than misleading per-protein statistics.

---

## Methodological notes

**Why no DEqMS here?**  
DEqMS moderates residual variance as a function of the peptide count per
protein. At the peptide level (Phobos's unit), this model is redundant: each
feature is already a peptide. limma eBayes alone is the reference for
label-free peptide data.

**Median normalisation vs DIA-NN LFQ:**  
PEAKS exports raw Areas without cross-run normalisation. Median normalisation
corrects global per-run offsets (equivalent to `normalizeMedianValues()` in
DEP/R). It is robust to outliers and assumes no specific distribution. The RLE
plot lets you check its quality.

**QRILC imputation:**  
DDA produces predominantly MNAR missing values (peptides below the
detection/fragmentation threshold). QRILC (Lazar 2016) draws from the
truncated left tail of the log2 distribution, guaranteeing imputed values
below the estimated detection limit. The MNAR/MAR diagnostic confirms model
adequacy.

---

## Name

**Phobos** — moon of Mars, companion of Deimos.  
Acronym: **P**eptide-level **H**euristics for **B**ottom-up **O**mics **S**uite.  
Peptide-level DDA pipeline (*bottom-up*), twin of the Deimos LFQ DIA suite.

---

## Trademarks

PEAKS® is a trademark of Bioinformatics Solutions Inc. Phobos is an
independent pipeline that reads PEAKS export files; it is not affiliated
with, endorsed by, or sponsored by Bioinformatics Solutions Inc.
