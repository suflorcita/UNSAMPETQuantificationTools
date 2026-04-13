"""
=====================================================
ANCOVA: PET comparison (FLENI vs COVID)
Controlling for Age and Sex
=====================================================

This script performs an ANCOVA (Analysis of Covariance)
for each brain region in the unified PET dataset,
comparing FLENI and COVID groups while controlling
for age and sex.

Author: [Your Name]
Date: [YYYY-MM-DD]
=====================================================
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols

# =====================================================
# 1. Define paths
# =====================================================
base_dir = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025"
# pet_csv = os.path.join(base_dir, "group_analysis_unificado", "PET_Fleni_COVID_resampled_unificado_mean_value_final.csv")
pet_csv = os.path.join(base_dir, "group_analysis_unificado", "PET_Fleni_COVID_resampled_unificado_cerebellum.csv")
info_csv = os.path.join(base_dir, "Datos-UNSAM-FLENI-DatosParticipantes.csv")
out_csv = os.path.join(base_dir, "group_analysis_unificado", "ANCOVA_Fleni_vs_COVID_results.csv")

# =====================================================
# 2. Load datasets
# =====================================================
print("📂 Loading datasets...")



# PET atlas data (regional values)
df_pet = pd.read_csv(pet_csv)

# =====================================================
# Limpieza de IDs en df_pet
# =====================================================
df_pet["ID"] = df_pet["ID"].astype(str).str.strip()

# 1️⃣ Quitar prefijo 'sub-' de todos los IDs
df_pet["ID"] = df_pet["ID"].str.replace("^sub-", "", regex=True)

# 2️⃣ Quitar sufijo '_old' solo si existe
df_pet["ID"] = df_pet["ID"].str.replace("_old$", "", regex=True)

# =====================================================

# Clean region column names: replace "-" with "_"
df_pet.columns = df_pet.columns.str.replace("-", "_", regex=False)

# Eliminar columnas con nombre "_"
df_pet = df_pet.loc[:, df_pet.columns != "_"]


#df_pet.to_csv(os.path.join(base_dir, "group_analysis_unificado", "PET_Fleni_COVID_resampled_unificado_cerebellum_final.csv"))

# Read CSV, skip empty header rows automatically
df_info = pd.read_csv(info_csv, skip_blank_lines=True)

# Drop completely empty rows and columns
df_info = df_info.dropna(how="all").dropna(axis=1, how="all")

# Find the row that contains the real headers
header_row = df_info[df_info.iloc[:, 1] == "ID"].index
if len(header_row) > 0:
    header_idx = header_row[0]
    df_info.columns = df_info.iloc[header_idx]       # set that row as header
    df_info = df_info.iloc[header_idx + 1:]          # drop everything above
else:
    print("⚠️ Header row with 'ID' not found — check file manually.")

# Clean column names
df_info.columns = df_info.columns.str.strip()

# Keep relevant columns and rename to English
df_info = df_info.rename(columns={
    "ID": "ID",
    "Sexo": "Sex",
    "Edad": "Age"
})[["ID", "Sex", "Age"]]

# Convert types
df_info["Age"] = pd.to_numeric(df_info["Age"], errors="coerce")
df_info["Sex"] = df_info["Sex"].astype(str).str.strip().str.upper()

# Standardize sex labels to 'M' and 'F'
df_info["Sex"] = (
    df_info["Sex"]
    .str.strip()
    .str.lower()
    .replace({
        "m": "M", "f": "F",
        "masculino": "M", "femenino": "F",
        "male": "M", "female": "F"
    })
)




# =====================================================
# 3. Merge datasets
# =====================================================
df = df_pet.merge(df_info, on="ID", how="left")

# Remove rows with missing group or covariates
df = df.dropna(subset=["Grupo", "Age", "Sex"])

# Convert categorical variables
df["Grupo"] = df["Grupo"].astype("category")
df["Sex"] = df["Sex"].astype("category")




print(f"✅ Total subjects: {len(df)}")
print(df["Grupo"].value_counts())

# =====================================================
# 4. Identify numeric (regional) columns
# =====================================================
region_cols = df.select_dtypes(include=["number"]).columns.tolist()
region_cols = [col for col in region_cols if col not in ["Age"]]  # exclude covariates

# =====================================================
# 5. ANCOVA for each region
# =====================================================
print("\n⚙️  Running ANCOVA for each region (Group ~ Age + Sex)...")

results = []

region = "FL_subgenual_frontal_cortex_R"
print(df[region].describe())
print(df[region].isna().sum())
print(df[region].unique())

for region in region_cols:
    try:
        # Build ANCOVA model
        model = ols(f"{region} ~ Grupo + Age + C(Sex)", data=df).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        print(model)
        # Extract F and p for Group effect
        F_value = anova_table.loc["Grupo", "F"]
        p_value = anova_table.loc["Grupo", "PR(>F)"]

        # Mean values per group
        means = df.groupby("Grupo")[region].mean()
        mean_fleni = means.get("FLENI", float("nan"))
        mean_covid = means.get("COVID", float("nan"))

        # Append results
        results.append({
            "Region": region,
            "Mean_FLENI": mean_fleni,
            "Mean_COVID": mean_covid,
            "F_value": F_value,
            "p_value": p_value
        })
    except Exception as e:
        print(f"⚠️ Skipped region {region}: {e}")

res_df = pd.DataFrame(results).sort_values("p_value")
res_df.to_csv(out_csv, index=False)

# =====================================================
# 7. Print significant results and summary
# =====================================================
alpha = 0.05
significant = res_df[res_df["p_value"] < alpha]
n_total = len(res_df)
n_signif = len(significant)
percent_signif = (n_signif / n_total) * 100 if n_total > 0 else 0

print(f"\n🔢 Total regions analyzed: {n_total}")
print(f"🔥 Significant regions (p < {alpha}): {n_signif} ({percent_signif:.1f}%)\n")

if not significant.empty:
    print(significant.to_string(index=False, formatters={
        "Mean_FLENI": "{:.3f}".format,
        "Mean_COVID": "{:.3f}".format,
        "Diff": "{:.3f}".format,
        "F_value": "{:.3f}".format,
        "p_value": "{:.4f}".format
    }))
else:
    print("❌ No significant regions found.")

print(f"\n✅ Full ANCOVA results saved at: {out_csv}")

alpha = 0.05
sig_regions = res_df[res_df["p_value"] < alpha]["Region"].tolist()

if len(sig_regions) == 0:
    print("❌ No hay regiones significativas para graficar.")
else:
    plot_dir = os.path.join(base_dir, "group_analysis_unificado", "boxplots_cerebellum")
    os.makedirs(plot_dir, exist_ok=True)
    print(f"📦 Generando boxplots en: {plot_dir}\n")

    for region in sig_regions:
        plt.figure(figsize=(3, 7))
        sns.boxplot(
            data=df,
            x="Grupo",
            y=region,
            hue="Grupo",
            palette={"FLENI": "#2E86C1", "COVID": "#E74C3C"},
            dodge=False
        )
        sns.stripplot(
            data=df,
            x="Grupo",
            y=region,
            color="black",
            alpha=0.5,
            jitter=0.15
        )

        plt.title(f"{region}\n(p = {res_df.loc[res_df['Region'] == region, 'p_value'].values[0]:.4f})")
        plt.xlabel("Grupo")
        plt.ylabel("Mean uptake")
        plt.tight_layout()

        out_path = os.path.join(plot_dir, f"{region}_boxplot.png")
        plt.savefig(out_path, dpi=300)
        plt.close()

    print(f"✅ {len(sig_regions)} boxplots generados.")