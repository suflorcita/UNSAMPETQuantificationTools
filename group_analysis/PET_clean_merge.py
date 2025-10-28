import os
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.formula.api import ols

# =======================
# 1. Paths
# =======================
base_pet = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/PET"
plots_dir = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/Plots/PET_7_OCT"
adni_data = "/media/sol/Expansion/ADNI/FDG/ADNI2_FDG_T1_6_05_2023.csv"
ods_controls = os.path.join(base_pet, "Controles_CEUNIM", "lista_estudios_pet_v01_anon_procesados.ods")
csv_covid = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/Freesurfer/DataAnalysis/ResumenTotal_2_06.csv"


os.makedirs(plots_dir, exist_ok=True)

# COVID
covid_total = os.path.join(base_pet, "ProcesadosCOVID_PET_09_25", "group_analysis", "PET_mean_values_by_subject_Hammers.csv")
covid_cereb = os.path.join(base_pet, "ProcesadosCOVID_PET_09_25", "group_analysis", "PET_cerebellum_norm_by_subject_Hammers.csv")

# Controles
cont_total = os.path.join(base_pet, "ProcesadosControlesCOVID_PET_09_25", "group_analysis", "PET_mean_values_by_subject_Hammers.csv")
cont_cereb = os.path.join(base_pet, "ProcesadosControlesCOVID_PET_09_25", "group_analysis", "PET_cerebellum_norm_by_subject_Hammers.csv")

# ADNI
adni_total = os.path.join(base_pet, "ProcesadasADNI2025", "group_analysis", "PET_mean_values_by_subject_Hammers.csv")
adni_cereb = os.path.join(base_pet, "ProcesadasADNI2025", "group_analysis", "PET_cerebellum_norm_by_subject_Hammers.csv")

# =======================
# 2. Leer y agregar Grupo (COVID Y Control)
# =======================
df_covid_total = pd.read_csv(covid_total, index_col=0); df_covid_total["Grupo"] = "COVID"
df_covid_cereb = pd.read_csv(covid_cereb, index_col=0); df_covid_cereb["Grupo"] = "COVID"

df_cont_total = pd.read_csv(cont_total, index_col=0); df_cont_total["Grupo"] = "Control"
df_cont_cereb = pd.read_csv(cont_cereb, index_col=0); df_cont_cereb["Grupo"] = "Control"

# b. ADNI
df_adni = pd.read_csv(adni_data)
df_adni_total = pd.read_csv(adni_total, index_col=0)
df_adni_cereb = pd.read_csv(adni_cereb, index_col=0)

# Creamos un diccionario {Subject: Group} para mapear
group_dict = dict(zip(df_adni["Subject"], df_adni["Group"]))

# Agregamos la columna Grupo a los dataframes de ADNI
df_adni_total["Grupo"] = df_adni_total.index.map(group_dict)
df_adni_cereb["Grupo"] = df_adni_cereb.index.map(group_dict)

# =======================
# 3. Concatenar
# =======================
df_total = pd.concat([df_covid_total, df_cont_total, df_adni_total], axis=0).reset_index().rename(columns={"index": "ID"})
df_cereb = pd.concat([df_covid_cereb, df_cont_cereb, df_adni_cereb], axis=0).reset_index().rename(columns={"index": "ID"})

# =======================
# 4. Info demográfica
# =======================
df_controls = pd.read_excel(ods_controls, engine="odf")
df_controls.columns = df_controls.columns.str.strip().str.lower()
df_controls = df_controls.rename(columns={"anonimizado": "ID"})[["ID", "sexo", "edad"]]

df_covid = pd.read_csv(csv_covid)
df_covid = df_covid.rename(columns={"Genero": "sexo", "Edad": "edad"})[["ID", "sexo", "edad"]]

df_adni = df_adni.rename(columns={"Sex": "sexo", "Age": "edad", "Subject": "ID"})[["ID", "sexo", "edad"]]

df_info = pd.concat([df_controls, df_covid, df_adni], axis=0)

# =======================
# 5. Limpiar IDs
# =======================
for df in [df_total, df_cereb]:
    df["ID"] = df["ID"].replace({"CP0107old": "CP0107"})
    df.drop(df[df["ID"] == "CP0106"].index, inplace=True)

# =======================
# 6. Merge con info y guardar
# =======================
df_total = df_total.merge(df_info, on="ID", how="left")
df_cereb = df_cereb.merge(df_info, on="ID", how="left")

# Limpiar nombres de columnas: sacar espacios raros, \xa0, etc.
df_total.columns = df_total.columns.str.replace("\xa0", "", regex=True).str.strip()
df_cereb.columns = df_cereb.columns.str.replace("\xa0", "", regex=True).str.strip()

# Eliminar filas sin grupo asignado
df_total = df_total.dropna(subset=["Grupo"])
df_cereb = df_cereb.dropna(subset=["Grupo"])

# Mostrar conteo por grupo
print("Distribución de grupos (df_total):")
print(df_total["Grupo"].value_counts())
print("\nDistribución de grupos (df_cereb):")
print(df_cereb["Grupo"].value_counts())


outdir = os.path.join(base_pet, "group_analysis")
os.makedirs(outdir, exist_ok=True)

df_total.to_csv(os.path.join(outdir, "PET_totalnorm_ALL_with_info.csv"), index=False)
df_cereb.to_csv(os.path.join(outdir, "PET_cerebellum_norm_ALL_with_info.csv"), index=False)




# =======================
# 7. Regiones para análisis
# =======================
cols_analysis = [c for c in df_total.columns if c not in ['Grupo', 'sexo', 'edad', "ID", "-"]]


import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, chi2_contingency

# =======================
# 1. T-test de edad entre Control y COVID
# =======================
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


# =======================
# 2. Chi-cuadrado sexo (Control vs COVID)
# =======================
contingency = pd.crosstab(df_total["Grupo"], df_total["sexo"])
chi2, p_chi, dof, expected = chi2_contingency(contingency)

print("\n🚻 Test Chi-cuadrado para sexo (Control vs COVID)")
print(contingency)
print(f"Chi² = {chi2:.3f}, p = {p_chi:.4f}")

if p_chi < 0.05:
    print("➡️ Diferencia significativa en distribución de sexo entre grupos.")
    # Determinar en cuál hay más mujeres u hombres
    for sexo in contingency.columns:
        grupo_con_mas = contingency[sexo].idxmax()
        print(f"   Hay más {sexo.lower()}s en el grupo {grupo_con_mas}.")
else:
    print("➡️ No hay diferencia significativa en la distribución de sexo.")


# =======================
# 3. Histograma de edad para Control
# =======================
plt.figure(figsize=(6,4))
plt.hist(control_edad, bins=10, color="skyblue", edgecolor="black")
plt.title("Distribución de edad - Grupo Control")
plt.xlabel("Edad")
plt.ylabel("Frecuencia")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "Histograma_edad_Control.png"), dpi=300)
plt.show()
print("\n📊 Histograma de edad del grupo Control guardado.")

# =======================
# 3. Histograma comparativo de edad (Control vs COVID)
# =======================
plt.figure(figsize=(7,5))

# Histograma de ambos grupos
plt.hist(control_edad, bins=10, color="skyblue", edgecolor="black", alpha=0.7, label="Control")
plt.hist(covid_edad, bins=10, color="lightcoral", edgecolor="black", alpha=0.6, label="COVID")

# Títulos y ejes
plt.title("Distribución de edad - Grupos Control y COVID", fontsize=13)
plt.xlabel("Edad (años)")
plt.ylabel("Frecuencia")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

# Guardar y mostrar
out_path = os.path.join(plots_dir, "Histograma_edad_Control_vs_COVID.png")
plt.savefig(out_path, dpi=300)
plt.show()

print(f"\n📊 Histograma comparativo de edad guardado en: {out_path}")


# =======================
# 4. Media global de edad (todos los grupos)
# =======================
print("\n📈 Edad promedio por grupo:")
media_grupo = df_total.groupby("Grupo")["edad"].mean().round(2)
print(media_grupo)

print(f"\n🧮 Media general (todos los grupos, incluyendo ADNI): {df_total['edad'].mean():.2f} años")



# # =======================
# # 8. ANCOVA + Boxplots (4 grupos)
# # =======================
# for df, label, cols in [
#     (df_total, "meanbrain", cols_analysis),
#     (df_cereb, "cerebellum", cols_analysis)
# ]:
#     results = []
#     skipped = []
#
#     print(f"\n🔎 Analizando dataset: {label} ({len(cols)} regiones)")
#
#     for region in cols:
#         try:
#             # Modelo ANCOVA con 4 grupos
#             formula = f'Q("{region}") ~ C(Grupo) + edad + C(sexo)'
#             model = ols(formula, data=df).fit()
#             anova_table = sm.stats.anova_lm(model, typ=2)
#
#             # Extraer resultados del factor Grupo
#             p_value = anova_table.loc["C(Grupo)", "PR(>F)"]
#             f_value = anova_table.loc["C(Grupo)", "F"]
#
#             # Medias por grupo
#             means = df.groupby("Grupo")[region].mean().to_dict()
#
#             # Guardar en tabla de resultados
#             results.append({
#                 "Region": region,
#                 "F": f_value,
#                 "p": p_value,
#                 "Significativo": "Sí" if p_value < 0.05 else "No",
#                 **{f"Mean_{g}": means.get(g, float("nan")) for g in df["Grupo"].unique()}
#             })
#
#             # Si es significativo → imprimir y graficar
#             if p_value < 0.05:
#                 print(f"✅ {label} | {region} SIGNIFICATIVO: "
#                       f"F={f_value:.3f}, p={p_value:.4f}")
#
#             # Boxplot vertical para todos los grupos
#             plt.figure(figsize=(6, 8))
#             df.boxplot(column=region, by="Grupo", grid=False)
#             plt.title(f"{region} (p={p_value:.3f})")
#             plt.suptitle("")
#             plt.ylabel("Uptake normalizado")
#             plt.xticks(rotation=45)
#             plt.tight_layout()
#             out_path = os.path.join(plots_dir, f"{region}_PET_boxplot_{label}.png")
#             plt.savefig(out_path, dpi=300)
#             plt.close()
#             print(f"   📊 Boxplot guardado en: {out_path}")
#
#         except Exception as e:
#             print(f"⚠️ {label} | Región {region} salteada por error: {e}")
#             skipped.append({"Region": region, "Error": str(e)})
#             continue
#
#     # Guardar tabla de resultados
#     results_df = pd.DataFrame(results)
#     res_path = os.path.join(plots_dir, f"ANCOVA_Regional_PET_{label}_results.csv")
#     results_df.to_csv(res_path, index=False)
#     print(f"💾 Resultados ANCOVA guardados en: {res_path}")
#
#     # Guardar log de regiones salteadas
#     if skipped:
#         skipped_df = pd.DataFrame(skipped)
#         log_path = os.path.join(plots_dir, f"ANCOVA_skipped_regions_{label}.csv")
#         skipped_df.to_csv(log_path, index=False)
#         print(f"⚠️ Log de regiones salteadas guardado en: {log_path}")
