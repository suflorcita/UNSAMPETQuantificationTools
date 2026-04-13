import os
import pandas as pd

# Ruta principal
base_dir = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ProcesadosCOVID_PET_10_25"

# Listas para guardar resultados
subjects = []
uptake_values = []
mean_list = []
cereb_list = []

for subject in os.listdir(base_dir):

    # if subject.startswith("CP") and subject != "CP0160":
    if not subject.startswith("CP"):
        pass

    # subj_path = os.path.join(base_dir, subject,"Smoothed_quantification", "CSV", "subject_normalization_values_hammers.csv")
    subj_path = os.path.join(base_dir, subject, "CSV", "subject_normalization_values_Hammers.csv")

    if os.path.exists(subj_path):
            # Leer CSV
            df = pd.read_csv(subj_path, sep="\t|,", engine="python")
            df.columns = df.columns.str.strip()  # limpiar espacios

            # Guardar como fila
            subjects.append(subject)
            mean_list.append(df.set_index("Structure")["Normalization to total brain mean value"])
            cereb_list.append(df.set_index("Structure")["Normalization to cerebellum uptake values"])
            uptake_values.append(df.set_index("Structure")["Regional uptake mean values"])

# Concatenar
df_mean = pd.DataFrame(mean_list, index=subjects)
df_cereb = pd.DataFrame(cereb_list, index=subjects)
df_uptake_values = pd.DataFrame(uptake_values, index=subjects)

# Quitar el prefijo "sub-" del índice
#df_mean.index = df_mean.index.str.replace("^sub-", "", regex=True)
#df_cereb.index = df_cereb.index.str.replace("^sub-", "", regex=True)
#df_uptake_values.index = df_uptake_values.index.str.replace("^sub-", "", regex=True)

# Renombrar el índice
df_mean.index.name = "Subject"
df_cereb.index.name = "Subject"
df_uptake_values.index.name = "Subject"

# Quitar columna "-"
df_mean.drop(columns=['-'], inplace=True)
df_cereb.drop(columns=['-'], inplace=True)
df_uptake_values.drop(columns=['-'], inplace=True)


# # Guardar resultados
# outdir = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/Analisis_España"
outdir = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ResultadosQuantification"

os.makedirs(outdir, exist_ok=True)

# df_mean.to_csv(os.path.join(outdir, "PET_MRI_mean_values_by_subject_Hammers_ANTs_smoothed.csv"))
# df_cereb.to_csv(os.path.join(outdir, "PET_MRI_cerebellum_norm_by_subject_Hammers_mri_ANTs_smoothed.csv"))
# df_uptake_values.to_csv(os.path.join(outdir, "PET_uptake_values_by_subject_Hammers_mri_ANTs_smoothedq.csv"))

df_mean.to_csv(os.path.join(outdir, "PET_quantification_normalized_total_brain_LC.csv"))
df_cereb.to_csv(os.path.join(outdir, "PET_quantification_normalized_cerebellum_LC.csv"))
df_uptake_values.to_csv(os.path.join(outdir, "PET_quantification_uptake_LC.csv"))



print("✅ CSV generados en", outdir)

