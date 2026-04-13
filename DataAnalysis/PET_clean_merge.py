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

