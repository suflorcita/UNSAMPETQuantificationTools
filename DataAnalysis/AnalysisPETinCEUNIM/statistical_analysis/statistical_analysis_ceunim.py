import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from statsmodels.formula.api import ols
from patsy.builtins import Q  # para nombres de ROIs con '-' u otros caracteres raros

# =======================
# CONFIG
# =======================
alpha = 0.05

base_dir = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ResultadosQuantification/group_analysis"
total_norm_path = os.path.join(base_dir, "PET_total_LC_Control_with_info.csv")

plots_dir = os.path.join(base_dir, "ANCOVA_total_norm_boxplots")
os.makedirs(plots_dir, exist_ok=True)

# =======================
# Helpers
# =======================
hemi_re = re.compile(r"-(L|R)$")

def split_hemi(roi_name: str):
    """
    Split ROI name into (base, hemi) where hemi is 'L' or 'R'.
    If no '-L'/'-R' suffix, returns (roi_name, None).
    """
    m = hemi_re.search(roi_name)
    if not m:
        return roi_name, None
    hemi = m.group(1)
    base = roi_name[:-2]  # remove '-L' or '-R'
    return base, hemi


def safe_filename(name: str) -> str:
    """Make ROI name safe for filenames."""
    name = name.replace("/", "_").replace("\\", "_")
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    return name[:200]


def ancova_all_rois(df: pd.DataFrame):
    """
    Run ANCOVA per ROI:
        ROI ~ Grupo (LC vs Control) + edad + sexo
    Returns:
        results_all: DataFrame with ROI, p, beta, means, ns
        df_clean: cleaned df used for plotting
    """
    required = {"ID", "Grupo", "edad", "sexo"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Faltan columnas: {missing}")

    df = df.copy()

    # Clean covariates / factors
    df["edad"] = pd.to_numeric(df["edad"], errors="coerce")
    df["sexo"] = df["sexo"].astype(str).str.strip()
    df["Grupo"] = df["Grupo"].astype(str).str.strip()

    # Keep only LC/Control and valid covariates
    df = df[df["Grupo"].isin(["LC", "Control"])].dropna(subset=["edad", "sexo", "Grupo"])

    # Reference = Control
    df["Grupo"] = pd.Categorical(df["Grupo"], categories=["Control", "LC"], ordered=True)

    non_roi = {"ID", "Grupo", "edad", "sexo"}
    roi_cols = [c for c in df.columns if c not in non_roi]

    term_group = 'C(Grupo, Treatment(reference="Control"))[T.LC]'

    rows = []
    for roi in roi_cols:
        sub = df.dropna(subset=[roi]).copy()

        # basic guards
        if sub.shape[0] < 10:
            continue
        if np.isclose(sub[roi].std(ddof=1), 0):
            continue

        # Raw group means + n
        gstats = sub.groupby("Grupo")[roi].agg(["mean", "count"])
        mean_ctrl = float(gstats.loc["Control", "mean"]) if "Control" in gstats.index else np.nan
        mean_lc = float(gstats.loc["LC", "mean"]) if "LC" in gstats.index else np.nan
        n_ctrl = int(gstats.loc["Control", "count"]) if "Control" in gstats.index else 0
        n_lc = int(gstats.loc["LC", "count"]) if "LC" in gstats.index else 0

        # ANCOVA
        formula = f'Q("{roi}") ~ C(Grupo, Treatment(reference="Control")) + edad + C(sexo)'
        try:
            model = ols(formula, data=sub).fit()
            p = float(model.pvalues.get(term_group, np.nan))
            beta = float(model.params.get(term_group, np.nan))

            rows.append({
                "ROI": roi,
                "p": p,
                "beta_LC_minus_Control": beta,
                "mean_Control": mean_ctrl,
                "mean_LC": mean_lc,
                "n_Control": n_ctrl,
                "n_LC": n_lc,
            })
        except Exception:
            continue

    results_all = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    return results_all, df


def plot_lr_pair(df_clean: pd.DataFrame, sig_df: pd.DataFrame, base: str, roi_L: str, roi_R: str, outpath: str):
    # más achatado
    fig, axes = plt.subplots(1, 2, figsize=(6,6), sharey=True)
    fig.suptitle(base, fontsize=12, fontweight="bold", y=1.02)

    # colores estilo tu ejemplo (Control azul, LC rojo)
    colors = {"Control": "#88BDE6", "LC": "#F08A84"}  # pastel-like

    def _panel(ax, roi, hemi):
        if roi is None or roi not in df_clean.columns:
            ax.axis("off")
            ax.set_title(f"{hemi}: missing")
            return

        data = df_clean[["Grupo", roi]].dropna()
        vals_control = data.loc[data["Grupo"] == "Control", roi].values
        vals_lc = data.loc[data["Grupo"] == "LC", roi].values

        # --- boxplot (matplotlib puro)
        bp = ax.boxplot(
            [vals_control, vals_lc],
            labels=["CONTROL", "LC"],
            widths=0.55,
            showfliers=False,
            patch_artist=True
        )

        # colorear cajas + líneas
        for patch, grp in zip(bp["boxes"], ["Control", "LC"]):
            patch.set_facecolor(colors[grp])
            patch.set_edgecolor("#444444")
            patch.set_alpha(0.85)
            patch.set_linewidth(0.9)

        for k in ["whiskers", "caps", "medians"]:
            for item in bp[k]:
                item.set_color("#444444")
                item.set_linewidth(0.9)

        # --- puntos (jitter)
        rng = np.random.default_rng(0)
        x1 = 1 + rng.normal(0, 0.06, size=len(vals_control))
        x2 = 2 + rng.normal(0, 0.06, size=len(vals_lc))
        ax.scatter(x1, vals_control, s=8, alpha=0.35, linewidths=0)
        ax.scatter(x2, vals_lc, s=8, alpha=0.35, linewidths=0)

        # estética
        ax.set_title(f"{base}_{hemi}", fontsize=10, fontweight="bold")
        ax.grid(True, axis="y", alpha=0.18)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # --- texto abajo (media + p + beta)
        row = sig_df.loc[sig_df["ROI"] == roi]
        p = float(row["p"].values[0]) if len(row) else np.nan
        beta = float(row["beta_LC_minus_Control"].values[0]) if len(row) else np.nan

        mC = np.mean(vals_control) if len(vals_control) else np.nan
        mL = np.mean(vals_lc) if len(vals_lc) else np.nan

        txt = (
            f"mean(Control)={mC:.3g} | mean(LC)={mL:.3g}\n"
            f"p={p:.3g} | beta(LC-Control)={beta:.3g}"
        )
        ax.text(0.5, -0.32, txt, transform=ax.transAxes, ha="center", va="top", fontsize=8.8)

    _panel(axes[0], roi_L, "L")
    _panel(axes[1], roi_R, "R")

    plt.tight_layout()
    plt.savefig(outpath, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_single_roi(df_clean: pd.DataFrame, sig_df: pd.DataFrame, roi: str, outpath: str):
    fig, ax = plt.subplots(1, 1, figsize=(4.4, 6))
    ax.set_title(roi, fontsize=11, fontweight="bold")

    colors = {"Control": "#88BDE6", "LC": "#F08A84"}

    data = df_clean[["Grupo", roi]].dropna()
    vals_control = data.loc[data["Grupo"] == "Control", roi].values
    vals_lc = data.loc[data["Grupo"] == "LC", roi].values

    bp = ax.boxplot(
        [vals_control, vals_lc],
        labels=["CONTROL", "LC"],
        widths=0.55,
        showfliers=False,
        patch_artist=True
    )

    for patch, grp in zip(bp["boxes"], ["Control", "LC"]):
        patch.set_facecolor(colors[grp])
        patch.set_edgecolor("#444444")
        patch.set_alpha(0.85)
        patch.set_linewidth(0.9)

    for k in ["whiskers", "caps", "medians"]:
        for item in bp[k]:
            item.set_color("#444444")
            item.set_linewidth(0.9)

    # puntos jitter
    rng = np.random.default_rng(0)
    x1 = 1 + rng.normal(0, 0.06, size=len(vals_control))
    x2 = 2 + rng.normal(0, 0.06, size=len(vals_lc))
    ax.scatter(x1, vals_control, s=8, alpha=0.35, linewidths=0)
    ax.scatter(x2, vals_lc, s=8, alpha=0.35, linewidths=0)

    ax.grid(True, axis="y", alpha=0.18)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    row = sig_df.loc[sig_df["ROI"] == roi]
    p = float(row["p"].values[0]) if len(row) else np.nan
    beta = float(row["beta_LC_minus_Control"].values[0]) if len(row) else np.nan

    mC = np.mean(vals_control) if len(vals_control) else np.nan
    mL = np.mean(vals_lc) if len(vals_lc) else np.nan

    txt = (
        f"mean(Control)={mC:.3g} | mean(LC)={mL:.3g}\n"
        f"p={p:.3g} | beta(LC-Control)={beta:.3g}"
    )
    ax.text(0.5, -0.30, txt, transform=ax.transAxes, ha="center", va="top", fontsize=8.8)

    plt.tight_layout()
    plt.savefig(outpath, dpi=220, bbox_inches="tight")
    plt.close(fig)



# =======================
# RUN (solo total_norm)
# =======================
if not os.path.exists(total_norm_path):
    raise FileNotFoundError(f"No existe: {total_norm_path}")

df_raw = pd.read_csv(total_norm_path)

results_all, df_clean = ancova_all_rois(df_raw)

# significativas (sin múltiples comparaciones)
sig = results_all[(results_all["p"].notna()) & (results_all["p"] < alpha)].copy()

# =======================
# Conteo por prefijo (TL/FL/PL/CG/...)
# =======================

# 1) extraer el prefijo antes del primer "-"
sig["prefix"] = sig["ROI"].str.extract(r"^([A-Za-z]+)-", expand=False)

print("\n" + "="*80)
print("📌 Conteo de ROIs SIGNIFICATIVAS por prefijo")
print("="*80)
print(sig["prefix"].value_counts(dropna=False))

# 2) (opcional) mapear prefijo a lóbulo/área “humana”
prefix_map = {
    "FL": "Frontal",
    "TL": "Temporal",
    "PL": "Parietal",
    "OL": "Occipital",
    "CG": "Cingulate",
}

sig["lobe"] = sig["prefix"].map(prefix_map).fillna("Other/Unknown")

print("\n" + "="*80)
print("📌 Conteo de ROIs SIGNIFICATIVAS por lóbulo (mapeado)")
print("="*80)
print(sig["lobe"].value_counts())


print("\n" + "=" * 80)
print(f"TOTAL_NORM | α={alpha} | ROIs significativas: {len(sig)}")
print("=" * 80)

if sig.empty:
    print("No hay ROIs significativas ajustando por edad y sexo.")
    raise SystemExit



# detectar hemisferio por sufijo -L / -R
sig["hemi"] = sig["ROI"].str.extract(r"-(L|R)$", expand=False)

# base = ROI sin el sufijo -L/-R (si no tiene, queda igual)
sig["base_roi"] = sig["ROI"].str.replace(r"-(L|R)$", "", regex=True)

# agrupar por base y ver qué hemisferios tiene
hemi_by_base = sig.groupby("base_roi")["hemi"].apply(lambda x: set(x.dropna())).to_dict()

both = sorted([b for b, hs in hemi_by_base.items() if hs == {"L","R"}])
only_L = sorted([b for b, hs in hemi_by_base.items() if hs == {"L"}])
only_R = sorted([b for b, hs in hemi_by_base.items() if hs == {"R"}])

# ROIs sin hemisferio (no terminan en -L/-R)
no_hemi = sorted(sig.loc[sig["hemi"].isna(), "ROI"].unique().tolist())

print("\n" + "="*80)
print(f"📌 Bases con ambos hemisferios significativos (L y R): {len(both)}")
print("="*80)
for b in both:
    print("-", b)

print("\n" + "="*80)
print(f"📌 Bases significativas solo IZQUIERDA (solo L): {len(only_L)}")
print("="*80)
for b in only_L:
    print("-", b + "-L")

print("\n" + "="*80)
print(f"📌 Bases significativas solo DERECHA (solo R): {len(only_R)}")
print("="*80)
for b in only_R:
    print("-", b + "-R")

if no_hemi:
    print("\n" + "="*80)
    print(f"📌 ROIs significativas sin hemisferio detectado: {len(no_hemi)}")
    print("="*80)
    for r in no_hemi:
        print("-", r)



# -----------------------
# Dividir hiper vs hipo (dirección ajustada por covariables)
# -----------------------
hyper = sig[sig["beta_LC_minus_Control"] > 0].copy()
hypo = sig[sig["beta_LC_minus_Control"] < 0].copy()

print("\n=== HIPERPERFUSIÓN (LC > Control) ===")
if hyper.empty:
    print("(ninguna)")
else:
    for _, r in hyper.sort_values("p").iterrows():
        print(f"- {r['ROI']} | p={r['p']:.4g} | beta={r['beta_LC_minus_Control']:.4g}")

print("\n=== HIPOPERFUSIÓN (LC < Control) ===")
if hypo.empty:
    print("(ninguna)")
else:
    for _, r in hypo.sort_values("p").iterrows():
        print(f"- {r['ROI']} | p={r['p']:.4g} | beta={r['beta_LC_minus_Control']:.4g}")

# -----------------------
# Agrupar por base ROI para generar plots:
# - Si L y R son significativas: un plot LR
# - Si solo una (L o R): plot single
# - Si no tiene hemi: plot single
# -----------------------
pairs = {}
for roi in sig["ROI"]:
    base, hemi = split_hemi(roi)
    if hemi is None:
        pairs.setdefault(base, {})["S"] = roi  # single (no hemi)
    else:
        pairs.setdefault(base, {})[hemi] = roi

n_lr = 0
n_single = 0

for base, d in pairs.items():
    has_L = "L" in d
    has_R = "R" in d
    has_S = "S" in d

    if has_L and has_R:
        n_lr += 1
        outpath = os.path.join(plots_dir, f"{safe_filename(base)}_LR_boxplot.png")
        plot_lr_pair(df_clean, sig, base, d["L"], d["R"], outpath)
    else:
        n_single += 1
        roi_single = d.get("L") or d.get("R") or d.get("S")
        outpath = os.path.join(plots_dir, f"{safe_filename(roi_single)}_boxplot.png")
        plot_single_roi(df_clean, sig, roi_single, outpath)

print("\n" + "=" * 80)
print(f"✅ Boxplots guardados en: {plots_dir}")
print(f"   - Pares L/R: {n_lr}")
print(f"   - Singles:  {n_single}")
print("=" * 80)

# Guardar CSV con significativas (opcional)
sig_out = os.path.join(plots_dir, "ANCOVA_total_norm_significativas.csv")
sig.to_csv(sig_out, index=False)
print(f"✅ CSV guardado: {sig_out}")
