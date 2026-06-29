# Changelog

All notable changes to Phobos are documented here.
Format based on [Keep a Changelog](https://keepachangelog.com/).

## [1.0.0] — 2025

### Added
- Initial release of **Phobos** — Peptide-level Heuristics for Bottom-up
  Omics Suite, companion to Deimos for DDA data.
- PEAKS `protein-peptides.csv` / `.xlsx` parser (precursor-level features:
  sequence + charge + modifications, inline `(+x.xx)` modifications preserved).
- `-10lgP` score filtering and `Area <label>` ↔ design column matching
  (exact-suffix priority, fuzzy fallback).
- Median normalisation of raw Areas + QRILC imputation (MNAR, DDA-suited)
  with MNAR/MAR diagnostic.
- Per-peptide imputation tracking columns (imputed flag + num_imputed
  count), to cross-check fragile hits against the robustness score.
- limma eBayes differential analysis at the peptide level (all pairwise
  contrasts), Pi-score and robustness score.
- QC suite: detection frequency, boxplots, RLE plot, missing-value heatmap,
  score and charge-state distributions.
- Exploratory analyses: PCA (± ellipses), UMAP, Pearson correlation, ANOVA
  with Z-score heatmaps and K-means clustering, UpSet intersections.
- Post-hoc protein aggregation (median rollup, visualisation only).
- **PTM-targeted sub-analyses**: the full pipeline (volcano/ANOVA/heatmap/
  UpSet) re-run on PTM-filtered peptide subsets. Interactive or YAML-driven
  selection (phospho, oxidation, methylation mono/di/tri, acetylation).
  Combined text + inline-mass detection with acetyl/trimethyl disambiguation.
- Multi-sheet Excel report (`PeptideAnalysis_Results.xlsx`).
- Self-contained interactive HTML dashboard with peptide/protein toggle,
  a volcano contrast-invert button, permanent sample labels on PCA/UMAP, and an
  interactive QC tab (detection, real pre-imputation missingness, intensity
  boxplots, RLE). A Log2_PreImpute sheet exposes the real missing values.
  peptide facets (charge, modifications, length), and Proteogen branding.
- `make_design_template.py` helper to pre-fill the experimental design from
  the Area columns.
- Phobos icon and mark (SVG), companion to the Deimos visual identity.
- **FASTA-based validation columns**: `amidation_Gflank` validates C-terminal
  amidation by checking for a +1 glycine in the precursor (PAM mechanism), and
  `signal_peptide` flags a signal peptide on the parent protein (von Heijne-style
  heuristic, yes/no). Both pure Python, computed when a FASTA is present.
- **Protein description recovery from a FASTA**: auto-detects a UniProt-style
  FASTA in the input folder and backfills the Description column (and missing
  gene names) by joining on the Accession. Pure Python, no Biopython.

### Notes
- Shares `limma_ebayes.py` and `config.py` with the Deimos pipeline
  (copy them alongside `phobos.py`).
- DEqMS is intentionally not used (redundant at the peptide level).
