# ==============================================================================
# limma_ebayes.py — Implémentation Python de l'eBayes de limma (Smyth 2004)
# Pipeline protéomique — Proteogen
# ==============================================================================
# Reproduit :
#   lmFit()         → régression OLS par protéine sur design matrix
#   makeContrasts() → combinaisons linéaires de coefficients
#   contrasts.fit() → recompute sur les contrastes
#   eBayes()        → modération empirique (emprunt de variance entre features)
#   topTable()      → statistiques finales avec p-values ajustées
# ==============================================================================

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma
from itertools import combinations


# ------------------------------------------------------------------------------
# 1. lmFit — régression OLS sur la design matrix
# ------------------------------------------------------------------------------
def lm_fit(expr_matrix: np.ndarray, design: np.ndarray):
    """
    Régression OLS protéine par protéine.

    Parameters
    ----------
    expr_matrix : np.ndarray (n_proteins × n_samples)
    design      : np.ndarray (n_samples × n_coefficients)

    Returns
    -------
    dict avec : coefficients, stdev_unscaled, sigma, df_residual
    """
    n_proteins, n_samples = expr_matrix.shape
    n_coef = design.shape[1]
    df_residual = n_samples - n_coef

    # Pseudo-inverse de la design matrix (stable même si colonnes corrélées)
    XtX_inv = np.linalg.pinv(design.T @ design)
    hat = XtX_inv @ design.T  # (n_coef × n_samples)

    coefficients = expr_matrix @ hat.T  # (n_proteins × n_coef)
    fitted = coefficients @ design.T    # (n_proteins × n_samples)
    residuals = expr_matrix - fitted    # (n_proteins × n_samples)

    # sigma² par protéine
    rss = np.sum(residuals ** 2, axis=1)
    sigma2 = rss / df_residual
    sigma = np.sqrt(np.maximum(sigma2, 0))

    # stdev_unscaled : écart-type non-scalé des coefficients (sans sigma)
    stdev_unscaled = np.sqrt(np.diag(XtX_inv))  # (n_coef,)
    stdev_unscaled_mat = np.tile(stdev_unscaled, (n_proteins, 1))  # broadcast

    return {
        "coefficients": coefficients,
        "stdev_unscaled": stdev_unscaled_mat,
        "sigma": sigma,
        "df_residual": df_residual,
        "design": design,
    }


# ------------------------------------------------------------------------------
# 2. contrasts_fit — projette le fit sur les contrastes
# ------------------------------------------------------------------------------
def contrasts_fit(fit: dict, contrast_matrix: np.ndarray):
    """
    Recalcule coefficients et stdev_unscaled pour chaque contraste.

    Parameters
    ----------
    fit             : dict retourné par lm_fit
    contrast_matrix : np.ndarray (n_coef × n_contrasts)

    Returns
    -------
    dict enrichi avec les statistiques par contraste
    """
    coef = fit["coefficients"]       # (n_proteins × n_coef)
    su   = fit["stdev_unscaled"]     # (n_proteins × n_coef)
    C    = contrast_matrix           # (n_coef × n_contrasts)

    # Nouveaux coefficients (logFC par contraste)
    new_coef = coef @ C              # (n_proteins × n_contrasts)

    # Variance non-scalée par contraste : diag(C^T * XtX_inv * C)
    XtX_inv = np.linalg.pinv(fit["design"].T @ fit["design"])
    stdev_contrasts = np.sqrt(np.diag(C.T @ XtX_inv @ C))  # (n_contrasts,)
    new_su = np.tile(stdev_contrasts, (coef.shape[0], 1))   # broadcast

    fit_c = fit.copy()
    fit_c["coefficients"] = new_coef
    fit_c["stdev_unscaled"] = new_su
    return fit_c


# ------------------------------------------------------------------------------
# 3. eBayes — modération empirique de Bayes (Smyth 2004)
# ------------------------------------------------------------------------------
def ebayes(fit: dict, fdr_global: bool = False):
    """
    Modération empirique de Bayes des variances résiduelles.
    Estime les hyperparamètres (d0, s0²) par fitting de la distribution
    a priori sur les sigma² observés.

    fdr_global : si True, la correction BH est appliquée sur l'ENSEMBLE des
        p-values (tous contrastes confondus, une seule famille de tests) au
        lieu d'une correction indépendante par contraste. Plus conservateur :
        contrôle le FDR au niveau de toute l'étude.

    Returns
    -------
    dict enrichi : t_stat, p_value, p_adj (BH), logFC, s2_post, df_total
    """
    sigma   = fit["sigma"]          # (n_proteins,)
    su      = fit["stdev_unscaled"] # (n_proteins × n_contrasts)
    df_res  = fit["df_residual"]
    coef    = fit["coefficients"]   # (n_proteins × n_contrasts)

    n_proteins, n_contrasts = coef.shape

    # --- Estimation des hyperparamètres d0, s0² (fitFDist exact de limma) ---
    # Référence : Smyth (2004), limma::fitFDist.
    # Transformation : e = log(s²) - digamma(df_res/2) + log(df_res/2)
    #   E[e]   = log(s0²) - digamma(d0/2) + log(d0/2)
    #   Var[e] = trigamma(df_res/2) + trigamma(d0/2)
    # => trigamma(d0/2) = Var[e] - trigamma(df_res/2)
    # => d0 = 2 * trigamma^{-1}(excess)
    # => s0² = exp( E[e] + digamma(d0/2) - log(d0/2) )
    s2 = sigma ** 2
    valid = s2 > 0
    z = np.log(s2[valid])

    e = z - digamma(df_res / 2.0) + np.log(df_res / 2.0)
    emean = np.mean(e)
    evar  = np.var(e, ddof=1) - polygamma(1, df_res / 2.0)

    if evar > 0:
        d0 = 2.0 * _inverse_trigamma(evar)
        d0 = float(min(d0, 1e6))
        s0_sq = np.exp(emean + digamma(d0 / 2.0) - np.log(d0 / 2.0))
    else:
        # Variances homogènes : df.prior = Inf (modération complète vers s0²)
        d0 = np.inf
        s0_sq = np.exp(emean)
    s0_sq = max(float(s0_sq), 1e-12)

    # --- Variance postérieure modérée ---
    if np.isinf(d0):
        s2_post = np.full_like(s2, s0_sq)
        df_total = float(df_res + 1e6)   # ddl quasi-infinis pour le test t
    else:
        df_total = df_res + d0
        s2_post = (df_res * s2 + d0 * s0_sq) / df_total  # shrinkage vers s0²
    s2_post = np.maximum(s2_post, 1e-12)

    # --- t-statistiques modérées et p-values ---
    t_stat  = np.zeros((n_proteins, n_contrasts))
    p_value = np.zeros((n_proteins, n_contrasts))
    p_adj   = np.zeros((n_proteins, n_contrasts))

    for j in range(n_contrasts):
        se_j = su[:, j] * np.sqrt(s2_post)
        se_j = np.maximum(se_j, 1e-10)
        t_j = coef[:, j] / se_j
        t_stat[:, j] = t_j

        # p-value bilatérale (distribution t avec df_total)
        p_j = 2 * stats.t.sf(np.abs(t_j), df=df_total)
        p_value[:, j] = p_j

        if not fdr_global:
            # Correction BH par contraste (une famille de tests par comparaison)
            p_adj[:, j] = _bh_correction(p_j)

    if fdr_global:
        # Correction BH GLOBALE : toutes les p-values empilées en une seule
        # famille de tests, puis re-dispatchées dans leur colonne d'origine.
        flat = p_value.flatten()
        valid = ~np.isnan(flat)
        adj_flat = np.full_like(flat, np.nan)
        adj_flat[valid] = _bh_correction(flat[valid])
        p_adj = adj_flat.reshape(p_value.shape)

    fit_eb = fit.copy()
    fit_eb.update({
        "t_stat":   t_stat,
        "p_value":  p_value,
        "p_adj":    p_adj,
        "s2_post":  s2_post,
        "df_total": df_total,
        "d0":       d0,
        "s0_sq":    s0_sq,
    })
    return fit_eb


# ------------------------------------------------------------------------------
# 3bis. spectra_count_ebayes — modération DEqMS (Zhu et al. 2020)
# ------------------------------------------------------------------------------
def spectra_count_ebayes(fit: dict, pep_count: np.ndarray,
                         loess_span: float = 0.75, fdr_global: bool = False):
    """
    Réplique DEqMS::spectraCounteBayes.

    Au lieu de modérer la variance vers une constante a priori (eBayes classique),
    DEqMS modélise la variance résiduelle EN FONCTION du nombre de peptides
    quantifiés par protéine. La régression se fait sur log2(count), suivant
    l'IBMT de Sartor et al. (2006) adapté à la protéomique.

    Parameters
    ----------
    fit       : dict retourné par lm_fit + contrasts_fit (avant ou après ebayes)
    pep_count : np.ndarray (n_proteins,) — nombre MINIMUM de peptides quantifiés
                par protéine across échantillons (depuis le pr_matrix)
    loess_span: fraction pour le lissage loess de la relation variance~count

    Returns
    -------
    dict enrichi : t_stat, p_value, p_adj, s2_post (DEqMS), df_prior, df_total
    """
    from scipy.interpolate import UnivariateSpline
    from scipy.special import digamma, polygamma

    sigma   = fit["sigma"]            # écart-type résiduel par protéine
    su      = fit["stdev_unscaled"]   # (n_proteins × n_contrasts)
    df_res  = fit["df_residual"]
    coef    = fit["coefficients"]     # (n_proteins × n_contrasts)
    n_proteins, n_contrasts = coef.shape

    s2 = sigma ** 2

    # --- Filtrer les protéines valides (variance > 0, count >= 1) ---
    pep_count = np.asarray(pep_count, dtype=float)
    valid = (s2 > 0) & np.isfinite(s2) & (pep_count >= 1) & np.isfinite(pep_count)

    log_count = np.log2(pep_count)
    log_var   = np.log(s2)

    # --- Régression de la variance log(s²) sur log2(count) ---
    # DEqMS ajuste une loess ; on l'approxime par un spline lissé monotone-friendly.
    x_fit = log_count[valid]
    y_fit = log_var[valid]
    order = np.argsort(x_fit)
    x_sorted, y_sorted = x_fit[order], y_fit[order]

    # Moyenne des y par valeur unique de x (stabilise le spline sur comptages discrets)
    x_unique = np.unique(x_sorted)
    y_means = np.array([y_sorted[x_sorted == xu].mean() for xu in x_unique])

    if len(x_unique) >= 4:
        # Spline lissé (équivalent loess)
        n_smooth = max(3, int(len(x_unique) * (1 - loess_span)) + 3)
        try:
            spl = UnivariateSpline(x_unique, y_means, k=min(3, len(x_unique) - 1),
                                   s=len(x_unique))
            pred_log_var = spl(log_count)
        except Exception:
            # Fallback : régression linéaire log(var) ~ log2(count)
            b = np.polyfit(x_fit, y_fit, 1)
            pred_log_var = np.polyval(b, log_count)
    else:
        # Trop peu de niveaux de comptage → régression linéaire
        b = np.polyfit(x_fit, y_fit, 1)
        pred_log_var = np.polyval(b, log_count)

    # Variance a priori prédite par le comptage
    s2_prior = np.exp(pred_log_var)
    s2_prior = np.where(np.isfinite(s2_prior) & (s2_prior > 0), s2_prior, np.nanmedian(s2[valid]))

    # --- Estimation du df.prior (degrés de liberté a priori) ---
    # DEqMS : basé sur la variance des résidus log(s²) - log(s²_prior) (méthode IBMT)
    residuals = log_var[valid] - pred_log_var[valid]
    var_resid = np.var(residuals, ddof=1)
    trig_df = polygamma(1, df_res / 2.0)
    excess = var_resid - trig_df
    if excess <= 1e-4:
        # Résidus quasi nuls : la variance est presque entièrement expliquée par le
        # comptage. On borne df_prior à une valeur réaliste (DEqMS ne dépasse
        # quasiment jamais quelques dizaines) pour éviter une sur-confiance.
        df_prior = float(df_res * 4)   # plafond conservateur
    else:
        df_prior = 2.0 * _inverse_trigamma(excess)
        df_prior = float(np.clip(df_prior, 0.1, df_res * 4))

    # --- Variance postérieure DEqMS ---
    df_total = df_res + df_prior
    s2_post = (df_res * s2 + df_prior * s2_prior) / df_total
    s2_post = np.maximum(s2_post, 1e-10)

    # --- t-statistiques et p-values modérées ---
    t_stat  = np.zeros((n_proteins, n_contrasts))
    p_value = np.full((n_proteins, n_contrasts), np.nan)
    p_adj   = np.full((n_proteins, n_contrasts), np.nan)

    for j in range(n_contrasts):
        se_j = su[:, j] * np.sqrt(s2_post)
        se_j = np.maximum(se_j, 1e-10)
        t_j  = coef[:, j] / se_j
        t_stat[:, j] = t_j
        p_j = 2 * stats.t.sf(np.abs(t_j), df=df_total)
        p_value[:, j] = p_j
        if not fdr_global:
            # BH par contraste sur les protéines valides uniquement
            p_adj_col = np.full(n_proteins, np.nan)
            p_adj_col[valid] = _bh_correction(p_j[valid])
            p_adj[:, j] = p_adj_col

    if fdr_global:
        # BH global sur toutes les p-values valides (protéines valides × contrastes)
        flat = p_value.copy()
        flat[~valid, :] = np.nan
        flat = flat.flatten()
        ok = ~np.isnan(flat)
        adj_flat = np.full_like(flat, np.nan)
        adj_flat[ok] = _bh_correction(flat[ok])
        p_adj = adj_flat.reshape(p_value.shape)

    # Invalider les protéines sans comptage exploitable
    t_stat[~valid, :]  = np.nan
    p_value[~valid, :] = np.nan

    fit_dq = fit.copy()
    fit_dq.update({
        "t_stat":   t_stat,
        "p_value":  p_value,
        "p_adj":    p_adj,
        "s2_post":  s2_post,
        "s2_prior": s2_prior,
        "df_prior": df_prior,
        "df_total": df_total,
        "pep_count": pep_count,
        "count_valid": valid,
    })
    return fit_dq


# ------------------------------------------------------------------------------
# 4. top_table — extraction des résultats pour un contraste donné
# ------------------------------------------------------------------------------
def top_table(fit_eb: dict, coef_idx: int, protein_names=None) -> pd.DataFrame:
    """
    Retourne un DataFrame trié pour le contraste coef_idx.

    Parameters
    ----------
    fit_eb       : dict retourné par ebayes()
    coef_idx     : index du contraste (0-based)
    protein_names: list/array de noms (optionnel)

    Returns
    -------
    DataFrame : logFC, t, P.Value, adj.P.Val
    """
    n = fit_eb["coefficients"].shape[0]
    names = protein_names if protein_names is not None else np.arange(n)

    df = pd.DataFrame({
        "name":      names,
        "logFC":     fit_eb["coefficients"][:, coef_idx],
        "t":         fit_eb["t_stat"][:, coef_idx],
        "P.Value":   fit_eb["p_value"][:, coef_idx],
        "adj.P.Val": fit_eb["p_adj"][:, coef_idx],
    })
    return df


# ------------------------------------------------------------------------------
# 5. make_contrasts — génère la matrice de contrastes (toutes paires)
# ------------------------------------------------------------------------------
def make_all_contrasts(group_names: list) -> tuple[np.ndarray, list[str]]:
    """
    Génère toutes les paires de groupes sous forme de matrice de contrastes.

    Returns
    -------
    contrast_matrix : np.ndarray (n_groups × n_contrasts)
    contrast_names  : list de str, ex: ["GroupA_vs_GroupB", ...]
    """
    n = len(group_names)
    pairs = list(combinations(range(n), 2))
    contrast_matrix = np.zeros((n, len(pairs)))
    contrast_names = []

    for k, (i, j) in enumerate(pairs):
        contrast_matrix[i, k] = 1
        contrast_matrix[j, k] = -1
        g1 = group_names[i].replace(".", "_")
        g2 = group_names[j].replace(".", "_")
        contrast_names.append(f"{g1}_vs_{g2}")

    return contrast_matrix, contrast_names


# ------------------------------------------------------------------------------
# 6. make_design_matrix — matrice design ~0 + condition
# ------------------------------------------------------------------------------
def make_design_matrix(conditions: list) -> tuple[np.ndarray, list[str]]:
    """
    Crée la design matrix indicator ~0 + condition (one-hot encoding).

    Returns
    -------
    design      : np.ndarray (n_samples × n_groups)
    group_names : list des noms de groupes (colonnes)
    """
    import re
    # Nettoyage similaire à make.names() de R
    clean = [re.sub(r"[^A-Za-z0-9_]", ".", c) for c in conditions]
    clean = ["X" + c if c[0].isdigit() else c for c in clean]

    group_names = sorted(set(clean), key=lambda x: clean.index(x))
    # Ordre stable = ordre d'apparition
    seen = []
    for c in clean:
        if c not in seen:
            seen.append(c)
    group_names = seen

    n_samples = len(clean)
    n_groups = len(group_names)
    design = np.zeros((n_samples, n_groups))
    for i, c in enumerate(clean):
        j = group_names.index(c)
        design[i, j] = 1.0

    return design, group_names


# ------------------------------------------------------------------------------
# Helpers internes
# ------------------------------------------------------------------------------
def _inverse_trigamma(x: float, max_iter: int = 50) -> float:
    """
    Inverse de la fonction trigamma : retourne y tel que trigamma(y) = x.
    Algorithme de Newton sur 1/y (Smyth, limma::trigammaInverse).
    """
    x = float(x)
    if np.isnan(x):
        return np.nan
    if x > 1e7:
        return 1.0 / np.sqrt(x)
    if x < 1e-6:
        return 1.0 / x
    # Initialisation (approximation de départ de Smyth)
    y = 0.5 + 1.0 / x
    for _ in range(max_iter):
        tri = float(polygamma(1, y))          # trigamma(y)
        tetra = float(polygamma(2, y))        # tetragamma(y) = trigamma'(y)
        # Newton sur g(y)=trigamma(y)-x, mais reparamétré pour stabilité :
        # delta = tri * (1 - tri/x) / tetra
        if tetra == 0:
            break
        delta = tri * (1.0 - tri / x) / tetra
        y += delta
        if -delta / y < 1e-8:
            break
    return max(y, 1e-6)


def _bh_correction(p_values: np.ndarray) -> np.ndarray:
    """Correction Benjamini-Hochberg (FDR)."""
    n = len(p_values)
    order = np.argsort(p_values)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)

    adj = p_values * n / ranked
    # Monotonie décroissante depuis la fin
    adj_sorted = adj[order]
    for i in range(n - 2, -1, -1):
        adj_sorted[i] = min(adj_sorted[i], adj_sorted[i + 1])
    adj[order] = adj_sorted
    return np.clip(adj, 0, 1)
