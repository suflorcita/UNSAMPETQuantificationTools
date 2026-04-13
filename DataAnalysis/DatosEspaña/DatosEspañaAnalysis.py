import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress

# Ruta base
base_path = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/Analisis_España"

# Leer los CSV en DataFrames
df_uptake_ceunim = pd.read_csv(f"{base_path}/PET_uptake_values_by_subject_Hammers.csv")
df_neurocloud = pd.read_csv(f"{base_path}/CEUNIM_neurocloud_rois_analysis_data.csv", sep=";")

# CEUNIM
regiones_pares = {
    "middle_frontal_gyrus": ["FL-middle-frontal-gyrus-L", "FL-middle-frontal-gyrus-R"],
    "anterior_temporal_lobe_lateral": ["TL-anterior-temporal-lobe-lateral-part-L", "TL-anterior-temporal-lobe-lateral-part-R"],
    "anterior_temporal_lobe_medial": ["TL-anterior-temporal-lobe-medial-part-L", "TL-anterior-temporal-lobe-medial-part-R"],
    "subgenual_anterior_cingulate": ["FL-subgenual-frontal-cortex-L", "FL-subgenual-frontal-cortex-R"],
    "anterior_cingulate": ["CG-anterior-cingulate-gyrus-L", "CG-anterior-cingulate-gyrus-R"],
    "superior_frontal_gyrus": ["FL-superior-frontal-gyrus-L", "FL-superior-frontal-gyrus-R"],
    "middle_inferior_temporal_gyrus": ["TL-middle-and-inferior-temporal-gyrus-L", "TL-middle-and-inferior-temporal-gyrus-R"],
    "precentral_gyrus": ["FL-precentral-gyrus-L", "FL-precentral-gyrus-R"],
    "inferior_frontal_gyrus": ["FL-inferior-frontal-gyrus-L", "FL-inferior-frontal-gyrus-R"],
    "cuneus": ["OL-cuneus-L", "OL-cuneus-R"],
    "parahippocampal_gyrus": ["TL-parahippocampal-and-ambient-gyrus-L", "TL-parahippocampal-and-ambient-gyrus-R"],
    "fusiform_gyrus": ["TL-fusiform-gyrus-L", "TL-fusiform-gyrus-R"],
    "posterior_temporal_lobe": ["TL-posterior-temporal-lobe-L", "TL-posterior-temporal-lobe-R"],
    "lateral_occipital": ["OL-lateral-remainder-occipital-lobe-L", "OL-lateral-remainder-occipital-lobe-R"],
    "posterior_orbital_gyrus": ["FL-posterior-orbital-gyrus-L", "FL-posterior-orbital-gyrus-R"]
}

# -------------------------------------------------------
# Para la ínsula, combinar TODAS sus subdivisiones
# -------------------------------------------------------
insula_cols = [
    "insula-anterior-short-gyrus-L", "insula-anterior-short-gyrus-R",
    "insula-middle-short-gyrus-L", "insula-middle-short-gyrus-R",
    "insula-posterior-short-gyrus-L", "insula-posterior-short-gyrus-R",
    "insula-anterior-inferior-cortex-L", "insula-anterior-inferior-cortex-R",
    "insula-anterior-long-gyrus-L", "insula-anterior-long-gyrus-R",
    "insula-posterior-long-gyrus-L", "insula-posterior-long-gyrus-R"
]

# -------------------------------------------------------
# Calcular el promedio por región
# -------------------------------------------------------
df_mean_ceunim = pd.DataFrame(index=df_uptake_ceunim.index)

for region, cols in regiones_pares.items():
    df_mean_ceunim[region] = df_uptake_ceunim[cols].mean(axis=1)

# Promedio total de todas las subregiones de la ínsula
df_mean_ceunim["insula"] = df_uptake_ceunim[insula_cols].mean(axis=1)


# -------------------------------------------------------
# Diccionario de renombrado según nomenclatura CEUNIM
# -------------------------------------------------------
rename_dict = {
    "middle_frontal_gyrus": "Middle frontal gyrus",
    "anterior_temporal_lobe_lateral": "Anterior temporal lobe lateral part",
    "anterior_temporal_lobe_medial": "Anterior temporal lobe medial part",
    "subgenual_anterior_cingulate": "Subgenual anterior cingulate gyrus",
    "anterior_cingulate": "Gyrus cinguli anterior part",
    "superior_frontal_gyrus": "Superior frontal gyrus",
    "middle_inferior_temporal_gyrus": "Middle and inferior temporal gyrus",
    "precentral_gyrus": "Precentral gyrus",
    "inferior_frontal_gyrus": "Inferior frontal gyrus Brocca's area",
    "cuneus": "Cuneus",
    "parahippocampal_gyrus": "Parahippocampal gyrus",
    "fusiform_gyrus": "Fusiform gyrus",
    "posterior_temporal_lobe": "Posterior temporal lobe",
    "lateral_occipital": "Lateral remainder of occipital lobe",
    "posterior_orbital_gyrus": "Posterior orbital gyrus",
    "insula": "Insula"
}

# -------------------------------------------------------
# Aplicar el renombrado
# -------------------------------------------------------
df_mean_ceunim = df_mean_ceunim.rename(columns=rename_dict)
print(df_mean_ceunim)

# ESPAÑA


# --- 1. Calcular el promedio entre hemisferios ---
df_neurocloud["mean_value"] = df_neurocloud[["MEAN_VALUE_LH", "MEAN_VALUE_RH"]].mean(axis=1)

# --- 2. Agrupar por sujeto y región ---
# Si cada región aparece más de una vez por hemisferio, se promedian esas repeticiones
df_mean_neurocloud = (
    df_neurocloud.groupby(["PATIENT_ID", "ROI_NAME"])["mean_value"]
    .mean()
    .reset_index()
)

# --- 3. Pasar a formato ancho (una columna por región) ---
df_mean_neurocloud = df_mean_neurocloud.pivot(index="PATIENT_ID", columns="ROI_NAME", values="mean_value")

# --- 4. (Opcional) ordenar columnas alfabéticamente y resetear índice ---
df_mean_neurocloud = df_mean_neurocloud.reset_index().sort_index(axis=1)


# Graficos


# --- 1️⃣ Agregar columna de sujeto a df_mean_ceunim ---
# Asumiendo que el orden de df_uptake_ceunim coincide con los sujetos (una fila por sujeto)
# Si tenés una columna con el ID, reemplazá esta línea con: df_mean_ceunim["PATIENT_ID"] = df_uptake_ceunim["PATIENT_ID"]
df_mean_ceunim["Subject"] = df_uptake_ceunim.index

# Mover PATIENT_ID al inicio
df_mean_ceunim = df_mean_ceunim[["Subject"] + [c for c in df_mean_ceunim.columns if c != "Subject"]]

# --- 2️⃣ Verificar sujetos en ambos datasets ---
subjects_ceunim = set(df_mean_ceunim["Subject"])
subjects_neurocloud = set(df_mean_neurocloud["PATIENT_ID"])

only_in_ceunim = subjects_ceunim - subjects_neurocloud
only_in_neurocloud = subjects_neurocloud - subjects_ceunim

print(f"🧩 Sujetos en CEUNIM: {len(subjects_ceunim)}")
print(f"🧩 Sujetos en Neurocloud: {len(subjects_neurocloud)}")
print(f"🧠 Sujetos comunes: {len(subjects_ceunim & subjects_neurocloud)}")

if only_in_ceunim:
    print(f"⚠️ Sujetos solo en CEUNIM: {sorted(only_in_ceunim)}")
if only_in_neurocloud:
    print(f"⚠️ Sujetos solo en Neurocloud: {sorted(only_in_neurocloud)}")

# --- 3️⃣ Mantener solo los sujetos comunes y alinear ---
common_subjects = list(subjects_ceunim & subjects_neurocloud)

df_mean_ceunim = df_mean_ceunim[df_mean_ceunim["Subject"].isin(common_subjects)].set_index("Subject")
df_mean_neurocloud = df_mean_neurocloud[df_mean_neurocloud["PATIENT_ID"].isin(common_subjects)].set_index("PATIENT_ID")

# Asegurar mismo orden
df_mean_ceunim = df_mean_ceunim.loc[sorted(common_subjects)]
df_mean_neurocloud = df_mean_neurocloud.loc[sorted(common_subjects)]

# --- 4️⃣ Graficar comparaciones ---
output_dir = f"{base_path}/comparacion_ceunim_neurocloud"
os.makedirs(output_dir, exist_ok=True)

common_regions = set(df_mean_ceunim.columns).intersection(set(df_mean_neurocloud.columns))
print(f"\n📊 Regiones en común: {len(common_regions)} -> {sorted(common_regions)}")

for region in sorted(common_regions):
    x = df_mean_ceunim[region]
    y = df_mean_neurocloud[region]

    slope, intercept, r_value, p_value, std_err = linregress(x, y)

    plt.figure(figsize=(5, 5))
    plt.scatter(x, y, color='royalblue', alpha=0.7, label='Subjects')
    plt.plot(x, slope * x + intercept, color='red', label=f"y={slope:.2f}x+{intercept:.2f}")
    plt.title(f"{region}\nR²={r_value ** 2:.3f}  p={p_value:.3e}")
    plt.xlabel("CEUNIM uptake (Hammers)")
    plt.ylabel("Neurocloud mean value")
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/{region.replace(' ', '_')}.png", dpi=300)
    plt.close()

print(f"\n✅ Gráficos guardados en: {output_dir}")
