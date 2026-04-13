#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PET_analysis_expanded.py
Análisis estadístico y visualización de datos PET (COVID, Control, ADNI)
Fecha: 2025-10-14 (Expansión)
"""

import os
import pandas as pd
from scipy.stats import ttest_ind, chi2_contingency
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_boxplot_global(df, region, title_suffix, output_dir, file_name_suffix, order=["CN", "Control", "COVID"]):
    """
    Genera y guarda un boxplot con strip plot para una región específica, mostrando 3 grupos.

    Args:
        df (pd.DataFrame): DataFrame que contiene los datos de PET y la columna 'Grupo'.
        region (str): Nombre de la columna de la región a plotear.
        title_suffix (str): Texto adicional para el título del plot (ej. resultados estadísticos).
        output_dir (str): Directorio donde se guardará el archivo.
        file_name_suffix (str): Sufijo para el nombre del archivo (ej. "TTest_Sig").
        order (list, optional): Orden de los grupos en el eje X.
    """
    plt.figure(figsize=(4, 6))

    # Colores consistentes para CN, Control, COVID
    palette = {"CN": "#7cb8ff", "Control": "#ff9999", "COVID": "#ffc85a"}

    # Boxplot
    sns.boxplot(
        x="Grupo",
        y=region,
        data=df,
        palette=palette,
        order=order,
        width=0.6,
        showmeans=True,
        meanprops={"marker": "o", "markerfacecolor": "black", "markeredgecolor": "black"}
    )

    # Strip plot (puntos individuales)
    sns.stripplot(
        x="Grupo",
        y=region,
        data=df,
        color="black",
        order=order,
        size=4,
        jitter=True,
        alpha=0.6
    )

    # Título global
    plt.title(f"{region}\n{title_suffix}", fontsize=10)
    plt.ylabel("Uptake normalizado")
    plt.xlabel("")
    plt.grid(axis="y", alpha=0.3)

    plt.tight_layout()

    # Nombre de archivo dinámico
    out_path = os.path.join(output_dir, f"{region}_Boxplot_{file_name_suffix}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()

    return out_path


# =======================
# 1. Paths
# =======================
# Utiliza los mismos paths de tu script original
base_pet = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/PET"
plots_dir = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/Plots/PET_7_OCT"
os.makedirs(plots_dir, exist_ok=True)

# Cargar datasets limpios
# Usaremos df_total para la consistencia con el código existente
df_total = pd.read_csv(os.path.join(base_pet, "group_analysis", "PET_totalnorm_ALL_with_info.csv"))
# df_cereb = pd.read_csv(os.path.join(base_pet, "group_analysis", "PET_cerebellum_norm_ALL_with_info.csv")) # No se usa en la expansión

# ===============================================
# 2-5. (Edad y Sexo t-test/Chi2 y Histogramas) - Se mantienen igual
# ===============================================

# 2. Test de edad (t-test)
covid_edad = df_total.loc[df_total["Grupo"] == "COVID", "edad"].dropna()
control_edad = df_total.loc[df_total["Grupo"] == "Control", "edad"].dropna()
t_stat, p_ttest = ttest_ind(covid_edad, control_edad, equal_var=False)
print("\n🧠 Test t para edad (Control vs COVID)")
print(f"Media edad Control: {control_edad.mean():.2f} ± {control_edad.std():.2f}")
print(f"Media edad COVID:   {covid_edad.mean():.2f} ± {covid_edad.std():.2f}")
print(f"T = {t_stat:.3f}, p = {p_ttest:.4f}")
if p_ttest < 0.05:
    print("➡️ Diferencia significativa en edad entre grupos.")
else:
    print("➡️ No hay diferencia significativa en edad.")

# 3. Chi-cuadrado sexo (Control vs COVID)
contingency = pd.crosstab(df_total["Grupo"], df_total["sexo"])
chi2, p_chi, dof, expected = chi2_contingency(contingency)
print("\n🚻 Test Chi-cuadrado para sexo (Control vs COVID)")
print(contingency)
print(f"Chi² = {chi2:.3f}, p = {p_chi:.4f}")
if p_chi < 0.05:
    print("➡️ Diferencia significativa en distribución de sexo entre grupos.")
    for sexo in contingency.columns:
        grupo_con_mas = contingency[sexo].idxmax()
        print(f"   Hay más {sexo.lower()}s en el grupo {grupo_con_mas}.")
else:
    print("➡️ No hay diferencia significativa en la distribución de sexo.")

# 4. Histogramas de edad
plt.figure(figsize=(7, 5))
plt.hist(control_edad, bins=10, color="skyblue", edgecolor="black", alpha=0.7, label="Control")
plt.hist(covid_edad, bins=10, color="lightcoral", edgecolor="black", alpha=0.6, label="COVID")
plt.title("Distribución de edad - Grupos Control y COVID", fontsize=13)
plt.xlabel("Edad (años)")
plt.ylabel("Frecuencia")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
out_path_hist = os.path.join(plots_dir, "Histograma_edad_Control_vs_COVID.png")
plt.savefig(out_path_hist, dpi=300)
plt.close()
print(f"\n📊 Histograma comparativo de edad guardado en: {out_path_hist}")

# 5. Boxplots de Cerebelo - COPIA DEL ORIGINAL
# Esta sección se ejecutará de todas formas para mantener la consistencia con el script original,
# pero luego la sección 6. ANCOVA General será la más importante.

df_sub_orig = df_total[df_total["Grupo"].isin(["CN", "Control", "COVID"])].copy()
regions_orig = ["cerebellum-R", "cerebellum-L"]
results_orig = []
boxplot_dir_orig = os.path.join(plots_dir, "Cerebellum_CN_Control_COVID")
os.makedirs(boxplot_dir_orig, exist_ok=True)
print("\n📊 Conteo de sujetos por grupo (CN, Control, COVID) - Sección Cerebelo:")
print(df_sub_orig["Grupo"].value_counts())

for region in regions_orig:
    # ANCOVA
    formula = f'Q("{region}") ~ C(Grupo) + edad + C(sexo)'
    model = ols(formula, data=df_sub_orig).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    p_value = anova_table.loc["C(Grupo)", "PR(>F)"]
    f_value = anova_table.loc["C(Grupo)", "F"]
    results_orig.append(
        {"Region": region, "F": f_value, "p": p_value, "Significativo": "Sí" if p_value < 0.05 else "No"})
    print(f"\n🔍 ANCOVA para {region} (CN, Control, COVID): p={p_value:.4f}")

    if p_value < 0.05:
        # POST HOC (Tukey HSD)
        formula_adj = f'Q("{region}") ~ edad + C(sexo)'
        model_adj = ols(formula_adj, data=df_sub_orig).fit()
        df_sub_orig["residuales"] = model_adj.resid
        tukey = pairwise_tukeyhsd(endog=df_sub_orig["residuales"], groups=df_sub_orig["Grupo"], alpha=0.05)
        print("📈 Test Post Hoc (Tukey HSD) sobre residuales ajustados:")
        print(tukey.summary())
        tukey_df = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])
        tukey_path = os.path.join(boxplot_dir_orig, f"Tukey_{region}_CN_Control_COVID.csv")
        tukey_df.to_csv(tukey_path, index=False)
        print(f"💾 Resultados Tukey guardados en: {tukey_path}")

    # BOXPlot - Usando la función temporalmente (luego se define bien)
    # plt.figure(figsize=(4, 6))
    # sns.boxplot(x="Grupo", y=region, data=df_sub_orig, palette=["#7cb8ff", "#ff9999", "#ffc85a"], width=0.6, showmeans=True, meanprops={"marker": "o", "markerfacecolor": "black", "markeredgecolor": "black"})
    # sns.stripplot(x="Grupo", y=region, data=df_sub_orig, color="black", size=4, jitter=True, alpha=0.6)
    # plt.title(f"{region}  (p = {p_value:.4f})", fontsize=12)
    # plt.ylabel("Uptake normalizado")
    # plt.xlabel("")
    # plt.grid(axis="y", alpha=0.3)
    # plt.tight_layout()
    # out_path_box = os.path.join(boxplot_dir_orig, f"{region}_Boxplot_CN_Control_COVID.png")
    # plt.savefig(out_path_box, dpi=300)
    # plt.close()
    # print(f"💾 Boxplot guardado en: {out_path_box}")

# Guardar tabla general ANCOVA original
out_path_orig_ancova = os.path.join(boxplot_dir_orig, "ANCOVA_Cerebellum_CN_Control_COVID.csv")
pd.DataFrame(results_orig).to_csv(out_path_orig_ancova, index=False)
print(f"\n💾 Resultados ANCOVA Cerebelo (original) guardados en: {out_path_orig_ancova}")

# =====================================================================
# 1B. FILTRADO DE DATOS: CN de menos de 70 años
# =====================================================================

print("\n⚙️ INICIO FILTRADO DE DATOS...")

# Crear una máscara de filtrado
# CN menores de 70 AÑOS
filtro_cn_edad = (df_total['Grupo'] == 'CN') & (df_total['edad'] >= 70)

# 1. Contar cuántos sujetos CN serán eliminados
n_cn_eliminados = filtro_cn_edad.sum()

# 2. Aplicar el filtro: mantener todos los sujetos que NO cumplen el criterio de eliminación
# (es decir, mantener Control, COVID, y CN < 70 años)
df_total_filtrado = df_total[~filtro_cn_edad].copy()

# Actualizar el DataFrame principal que usarán los análisis siguientes
df_total = df_total_filtrado

# 3. Reporte de la operación
print(f"   - Se encontraron y eliminaron {n_cn_eliminados} sujetos CN con 70 o más años.")
print(f"   - El número total de sujetos en el nuevo dataset es: {len(df_total)}")
print("   - Grupos restantes y sus tamaños:")
print(df_total['Grupo'].value_counts().to_string())
print("✅ FILTRADO COMPLETADO. Los análisis posteriores usarán este dataset.")

# Opcional: Volver a calcular y mostrar las estadísticas de edad para el nuevo grupo CN
df_cn_nuevo = df_total[df_total['Grupo'] == 'CN']
if not df_cn_nuevo.empty:
    print("\n   - Nueva Media de Edad para CN:")
    print(f"     {df_cn_nuevo['edad'].mean():.2f} ± {df_cn_nuevo['edad'].std():.2f} años (N={len(df_cn_nuevo)})")
else:
    print("\n   - ¡Advertencia! No quedan sujetos CN después del filtrado.")

# =====================================================================
# 2. ESTADÍSTICAS DESCRIPTIVAS FINALES (Conteo y Edad)
# =====================================================================
print("\n" + "=" * 80)
print("  ESTADÍSTICAS DESCRIPTIVAS DESPUÉS DEL FILTRADO (CN < 70 años)")
print("=" * 80)

# Filtrar los datos solo para los tres grupos de interés
df_edad_stats = df_total[df_total["Grupo"].isin(["CN", "Control", "COVID"])].copy()

# Calcular media, desviación estándar y conteo
edad_stats = df_edad_stats.groupby("Grupo")["edad"].agg(['mean', 'std', 'count']).reset_index()

print("\n--- Conteo de Sujetos y Distribución de Edad (Media ± DE) ---")
print("--------------------------------------------------")
print("  Grupo | Media (años) | Desv. Est. | N")
print("--------------------------------------------------")
for index, row in edad_stats.iterrows():
    grupo = row['Grupo'].ljust(7)
    mean_std = f"{row['mean']:.2f} ± {row['std']:.2f}"
    count = int(row['count'])
    print(f" {grupo}| {mean_std.ljust(15)}| {count}")
print("--------------------------------------------------")

# Opcional: Ejecutar un ANOVA simple para ver si hay diferencias de edad entre los 3 grupos
try:
    from scipy.stats import f_oneway

    edad_cn = df_edad_stats[df_edad_stats['Grupo'] == 'CN']['edad'].dropna()
    edad_control = df_edad_stats[df_edad_stats['Grupo'] == 'Control']['edad'].dropna()
    edad_covid = df_edad_stats[df_edad_stats['Grupo'] == 'COVID']['edad'].dropna()

    f_stat, p_anova_edad = f_oneway(edad_cn, edad_control, edad_covid)

    print(f"\nANOVA de Edad (3 Grupos): F={f_stat:.2f}, p={p_anova_edad:.4f}")
    if p_anova_edad < 0.05:
        print("  ⚠️ Advertencia: La edad aún difiere significativamente entre los grupos.")
    else:
        print("  ✅ La edad es homogénea entre los 3 grupos (no significativa).")

except Exception as e:
    print(f"\n⚠️ Error al calcular ANOVA de Edad: {e}")

# =================================================================================
# 3. ANCOVA y filtro: No sig entre controles, sí sig con COVID
# =================================================================================

from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

df_sub = df_total[df_total["Grupo"].isin(["CN", "Control", "COVID"])].copy()

exclusion_cols = ["ID", "Grupo", "edad", "sexo"]
all_regions = [col for col in df_sub.columns if col not in exclusion_cols]

results_ancova = []

output_dir = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/Plots/PET_Filtro_COVID"
os.makedirs(output_dir, exist_ok=True)

print("\n" + "=" * 80)
print("   ANCOVA (CN, Control, COVID) - No Sig entre controles / Sí Sig con COVID")
print("=" * 80)

for region in all_regions:
    try:
        # ANCOVA general
        formula = f'Q("{region}") ~ C(Grupo) + edad + C(sexo)'
        model = ols(formula, data=df_sub).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        p_ancova = anova_table.loc["C(Grupo)", "PR(>F)"]
    except Exception:
        continue

    if p_ancova >= 0.05:
        continue  # si la ANCOVA no es significativa, pasamos

    # Post Hoc (Tukey sobre residuales ajustados)
    formula_adj = f'Q("{region}") ~ edad + C(sexo)'
    model_adj = ols(formula_adj, data=df_sub).fit()
    df_sub["residuales"] = model_adj.resid

    tukey = pairwise_tukeyhsd(endog=df_sub["residuales"], groups=df_sub["Grupo"], alpha=0.05)
    tukey_df = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])

    # Extraer p-values
    def p_val(g1, g2):
        subset = tukey_df[
            (tukey_df["group1"].isin([g1, g2])) &
            (tukey_df["group2"].isin([g1, g2])) &
            (tukey_df["group1"] != tukey_df["group2"])
        ]
        return subset["p-adj"].iloc[0] if not subset.empty else 1.0

    p_cn_control = p_val("CN", "Control")
    p_covid_cn = p_val("COVID", "CN")
    p_covid_control = p_val("COVID", "Control")

    # Filtrar: no sig entre controles, sí sig con COVID
    if (p_cn_control >= 0.05) and ((p_covid_cn < 0.05) or (p_covid_control < 0.05)):
        # Medias por grupo
        medias = (
            df_sub.groupby("Grupo")[region]
            .mean()
            .reindex(["CN", "Control", "COVID"])
            .to_dict()
        )

        # Dirección del efecto
        rels = []
        if p_covid_cn < 0.05:
            if medias["COVID"] > medias["CN"]:
                rels.append("COVID > CN")
            else:
                rels.append("COVID < CN")
        if p_covid_control < 0.05:
            if medias["COVID"] > medias["Control"]:
                rels.append("COVID > Control")
            else:
                rels.append("COVID < Control")

        relation = " y ".join(rels) if rels else "—"

        results_ancova.append({
            "Region": region,
            "p_ANCOVA": p_ancova,
            "p_CN_vs_Control": p_cn_control,
            "p_COVID_vs_CN": p_covid_cn,
            "p_COVID_vs_Control": p_covid_control,
            "Media_CN": medias["CN"],
            "Media_Control": medias["Control"],
            "Media_COVID": medias["COVID"],
            "Relación_Significativa": relation
        })

# Exportar resultados
df_results = pd.DataFrame(results_ancova)
out_path = os.path.join(output_dir, "Regiones_NoSig_CNvsCtrl_PeroSig_COVID_conMedias.csv")
df_results.to_csv(out_path, index=False)

print(f"\n💾 Resultados guardados en: {out_path}")
if not df_results.empty:
    print(f"✅ {len(df_results)} regiones cumplen el criterio.")
    print(df_results[["Region", "Relación_Significativa"]].to_string(index=False))
else:
    print("❌ No se encontraron regiones que cumplan el criterio.")

# =================================================================================
# 4. ANCOVA y filtro: Significativo entre controles y también con COVID
# =================================================================================

results_ancova_sig_all = []

print("\n" + "=" * 80)
print("   ANCOVA (CN, Control, COVID) - Sig entre controles Y Sig con COVID")
print("=" * 80)

for region in all_regions:
    try:
        formula = f'Q("{region}") ~ C(Grupo) + edad + C(sexo)'
        model = ols(formula, data=df_sub).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        p_ancova = anova_table.loc["C(Grupo)", "PR(>F)"]
    except Exception:
        continue

    if p_ancova >= 0.05:
        continue

    # Post Hoc (Tukey)
    formula_adj = f'Q("{region}") ~ edad + C(sexo)'
    model_adj = ols(formula_adj, data=df_sub).fit()
    df_sub["residuales"] = model_adj.resid

    tukey = pairwise_tukeyhsd(endog=df_sub["residuales"], groups=df_sub["Grupo"], alpha=0.05)
    tukey_df = pd.DataFrame(data=tukey.summary().data[1:], columns=tukey.summary().data[0])

    def p_val(g1, g2):
        subset = tukey_df[
            (tukey_df["group1"].isin([g1, g2])) &
            (tukey_df["group2"].isin([g1, g2])) &
            (tukey_df["group1"] != tukey_df["group2"])
        ]
        return subset["p-adj"].iloc[0] if not subset.empty else 1.0

    p_cn_control = p_val("CN", "Control")
    p_covid_cn = p_val("COVID", "CN")
    p_covid_control = p_val("COVID", "Control")

    # Filtro: significativos entre controles y con COVID
    if (p_cn_control < 0.05) and ((p_covid_cn < 0.05) or (p_covid_control < 0.05)):
        medias = (
            df_sub.groupby("Grupo")[region]
            .mean()
            .reindex(["CN", "Control", "COVID"])
            .to_dict()
        )

        rels = []
        if p_covid_cn < 0.05:
            if medias["COVID"] > medias["CN"]:
                rels.append("COVID > CN")
            else:
                rels.append("COVID < CN")
        if p_covid_control < 0.05:
            if medias["COVID"] > medias["Control"]:
                rels.append("COVID > Control")
            else:
                rels.append("COVID < Control")
        if p_cn_control < 0.05:
            if medias["Control"] > medias["CN"]:
                rels.append("Control > CN")
            else:
                rels.append("Control < CN")

        relation = " y ".join(rels) if rels else "—"

        results_ancova_sig_all.append({
            "Region": region,
            "p_ANCOVA": p_ancova,
            "p_CN_vs_Control": p_cn_control,
            "p_COVID_vs_CN": p_covid_cn,
            "p_COVID_vs_Control": p_covid_control,
            "Media_CN": medias["CN"],
            "Media_Control": medias["Control"],
            "Media_COVID": medias["COVID"],
            "Relación_Significativa": relation
        })

# Exportar resultados
df_results_sig_all = pd.DataFrame(results_ancova_sig_all)
out_path_all = os.path.join(output_dir, "Regiones_Sig_CNvsCtrl_Y_Sig_COVID_conMedias.csv")
df_results_sig_all.to_csv(out_path_all, float_format="%.6e", index=False)

print(f"\n💾 Resultados guardados en: {out_path_all}")
if not df_results_sig_all.empty:
    print(f"✅ {len(df_results_sig_all)} regiones cumplen el criterio.")
    print(df_results_sig_all[["Region", "Relación_Significativa"]].to_string(index=False))
else:
    print("✅ No se encontraron regiones con significancia entre controles y COVID.")

import matplotlib.pyplot as plt
import seaborn as sns

# Crear carpeta para guardar plots del cerebelo
plots_cereb_dir = os.path.join(plots_dir, "Cerebellum_Plots_with_IDs")
os.makedirs(plots_cereb_dir, exist_ok=True)

for region in ["cerebellum-R", "cerebellum-L"]:
    print(f"\n==================== {region} ====================")
    df_tmp = df_total[df_total["Grupo"].isin(["CN", "Control", "COVID"])].copy()

    # ANCOVA
    model = ols(f'Q("{region}") ~ C(Grupo) + edad + C(sexo)', data=df_tmp).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    print("\n🔹 ANCOVA results:")
    print(anova_table)

    # Post Hoc (Tukey sobre residuales ajustados)
    formula_adj = f'Q("{region}") ~ edad + C(sexo)'
    model_adj = ols(formula_adj, data=df_tmp).fit()
    df_tmp["residuales_tmp"] = model_adj.resid
    tukey = pairwise_tukeyhsd(endog=df_tmp["residuales_tmp"], groups=df_tmp["Grupo"], alpha=0.05)
    print("\n🔹 Post Hoc (Tukey HSD) sobre residuales ajustados:")
    print(tukey.summary())

    # === BOX PLOT CON IDs ===
    plt.figure(figsize=(5, 6))
    palette = {"CN": "#7cb8ff", "Control": "#ff9999", "COVID": "#ffc85a"}
    order = ["CN", "Control", "COVID"]

    sns.boxplot(
        x="Grupo",
        y=region,
        data=df_tmp,
        order=order,
        palette=palette,
        showmeans=True,
        meanprops={"marker": "o", "markerfacecolor": "black", "markeredgecolor": "black"},
        width=0.6
    )

    # Puntos individuales con anotaciones
    sns.stripplot(
        x="Grupo",
        y=region,
        data=df_tmp,
        order=order,
        color="black",
        size=5,
        jitter=True,
        alpha=0.6
    )

    # Añadir IDs sobre los puntos
    for i, row in df_tmp.iterrows():
        x_pos = order.index(row["Grupo"])
        plt.text(
            x_pos + 0.05,
            row[region],
            str(row["ID"]),
            fontsize=6,
            alpha=0.8,
            rotation=30
        )

    plt.title(f"{region}\nANCOVA p = {anova_table.loc['C(Grupo)', 'PR(>F)']:.6g}", fontsize=10)
    plt.ylabel("Uptake normalizado")
    plt.xlabel("")
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    # Guardar figura
    out_path_plot = os.path.join(plots_cereb_dir, f"{region}_Boxplot_with_IDs.png")
    plt.savefig(out_path_plot, dpi=300)
    plt.close()

    print(f"💾 Boxplot con IDs guardado en: {out_path_plot}")
