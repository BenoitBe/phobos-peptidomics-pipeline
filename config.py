# ==============================================================================
# config.py — Configuration management for Phobos (standalone)
# ==============================================================================
# Dedicated to Phobos: no WGCNA / GO / DEqMS prompts (unlike the Deimos config).
#
# Three modes, mirroring Deimos UX without the extra questions:
#   python phobos.py --config myproject.yaml   -> zero prompts (YAML-driven)
#   python phobos.py --save-config [name.yaml] -> ask params, then save them
#   python phobos.py                           -> reuse last_config or ask
#
# After every run, last_config_phobos.yaml is written into out_dir so the next
# run can offer to reuse it.
#
# resolve_config() keeps the SAME signature Phobos calls it with:
#   resolve_config(ask_params_fn, ask_go_params_fn=None, go_available=False,
#                  dash_available=False, pr_path="")
# The Deimos-specific arguments are accepted for compatibility but IGNORED.
# ==============================================================================

import os
import argparse


def _load_yaml(path: str) -> dict:
    import yaml
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config {path} is not a key/value mapping.")
    return data


def _save_yaml(params: dict, path: str):
    import yaml
    skip = {"volcano_lfc_min"}   # derived, recomputed at load
    clean = {k: v for k, v in params.items() if k not in skip}
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(clean, f, allow_unicode=True, sort_keys=False)


def _ensure_paths(params: dict, interactive: bool = True):
    """Prompt for input/design/output paths only if missing/non-existent.
    If interactive=False (e.g. --config provided everything), never prompt."""
    cur = params.get("peaks_csv") or params.get("tsv_path") or "protein-peptides.csv"
    if not os.path.exists(cur) and interactive:
        raw = input(f"\n  PEAKS peptide file [{cur}] -> ").strip()
        if raw:
            cur = raw
    params["peaks_csv"] = cur

    cur_d = params.get("design_path", "ExperimentalDesign.csv")
    if not os.path.exists(cur_d) and interactive:
        raw = input(f"  Experimental design CSV [{cur_d}] -> ").strip()
        if raw:
            cur_d = raw
    params["design_path"] = cur_d

    cur_o = params.get("out_dir", "phobos_output")
    # Only prompt for out_dir when it was NOT already specified in config/CLI
    if "out_dir" not in params and interactive:
        raw = input(f"  Output directory [{cur_o}] -> ").strip()
        if raw:
            cur_o = raw
    params["out_dir"] = cur_o


def _finalize(params: dict) -> dict:
    import numpy as np
    ratio = params.get("volcano_ratio_min")
    if ratio:
        params["volcano_lfc_min"] = float(np.log2(ratio))
    if "peaks_csv" not in params and "tsv_path" in params:
        params["peaks_csv"] = params["tsv_path"]
    return params


def resolve_config(ask_params_fn,
                   ask_go_params_fn=None,   # ignored (Deimos compat)
                   go_available=False,      # ignored
                   dash_available=False,    # ignored
                   pr_path=""):             # ignored
    """
    Resolve the Phobos configuration from CLI flags, a YAML file, a reusable
    last_config, or interactive prompts. Deimos-style extra args are unused.
    """
    parser = argparse.ArgumentParser(description="Phobos configuration",
                                     add_help=False)
    parser.add_argument("--config", "-c", metavar="FILE.yaml", default=None)
    parser.add_argument("--save-config", metavar="NAME.yaml", nargs="?",
                        const="__AUTO__", default=None)
    parser.add_argument("--peaks", "--csv", dest="peaks", default=None)
    parser.add_argument("--design", default=None)
    parser.add_argument("--out-dir", dest="out_dir", default=None)
    args, _unknown = parser.parse_known_args()

    params = None

    # 1) Explicit --config FILE.yaml -> zero prompts
    if args.config:
        try:
            params = _load_yaml(args.config)
            print(f"  [CONFIG] Loaded from {args.config}")
        except Exception as e:
            print(f"  [WARN] Could not load {args.config} ({e}). Interactive mode.")
            params = None

    # 2) Offer to reuse a previous last_config_phobos.yaml
    if params is None:
        guess_out = args.out_dir or "phobos_output"
        last = os.path.join(guess_out, "last_config_phobos.yaml")
        if os.path.exists(last):
            ans = input(f"\n  Found a previous config ({last}). Reuse it? "
                        f"[Y/n] -> ").strip().lower()
            if ans in ("", "y", "yes", "o", "oui"):
                try:
                    params = _load_yaml(last)
                    print("  [CONFIG] Previous config reloaded.")
                except Exception as e:
                    print(f"  [WARN] Reload failed ({e}). Asking parameters.")
                    params = None

    # 3) Interactive parameters
    if params is None:
        params = ask_params_fn()

    # CLI path overrides
    if args.peaks:
        params["peaks_csv"] = args.peaks
    if args.design:
        params["design_path"] = args.design
    if args.out_dir:
        params["out_dir"] = args.out_dir

    _ensure_paths(params, interactive=(args.config is None))
    params = _finalize(params)

    # Persist last_config in out_dir
    try:
        _save_yaml(params, os.path.join(params["out_dir"],
                                        "last_config_phobos.yaml"))
    except Exception:
        pass

    # Explicit --save-config NAME.yaml
    if args.save_config:
        try:
            import pandas as pd
            name = (args.save_config if args.save_config != "__AUTO__"
                    else f"phobos_config_"
                         f"{pd.Timestamp.now().strftime('%Y%m%d_%H%M')}.yaml")
            _save_yaml(params, name)
            print(f"  [CONFIG] Saved to {name}")
        except Exception as e:
            print(f"  [WARN] Could not save config: {e}")

    return params
