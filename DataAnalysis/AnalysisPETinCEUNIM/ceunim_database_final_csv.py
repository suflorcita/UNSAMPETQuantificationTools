import os
import pandas as pd

# =======================
# 1. Paths
# =======================
base_dir = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ResultadosQuantification"
outdir = os.path.join(base_dir, "group_analysis")
os.makedirs(outdir, exist_ok=True)

# CSVs PET (según imagen)
files = {
    "uptake": {
        "LC": "PET_quantification_uptake_LC.csv",
        "Control": "PET_quantification_uptake_controles.csv"
    },
    "total": {
        "LC": "PET_quantification_normalized_total_brain_LC.csv",
        "Control": "PET_quantification_normalized_total_brain_controles.csv"
    },
    "cereb": {
        "LC": "PET_quantification_normalized_cerebellum_LC.csv",
        "Control": "PET_quantification_normalized_cerebellum_controles.csv"
    }
}

# =======================
# 2. Info demográfica
# =======================
base_pet = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID"

ods_controls = os.path.join(
    base_pet, "PET", "Controles_CEUNIM", "lista_estudios_pet_v01_anon_procesados.ods"
)
csv_covid = os.path.join(
    base_pet, "Freesurfer", "DataAnalysis", "ResumenTotal_2_06.csv"
)

# Controles
df_controls = pd.read_excel(ods_controls, engine="odf")
df_controls.columns = df_controls.columns.str.strip().str.lower()
df_controls = df_controls.rename(columns={"anonimizado": "ID"})[["ID", "sexo", "edad"]]
df_controls["Grupo"] = "Control"

# Long COVID
df_lc = pd.read_csv(csv_covid)
df_lc = df_lc.rename(columns={"Genero": "sexo", "Edad": "edad"})[["ID", "sexo", "edad"]]
df_lc["Grupo"] = "LC"

# Info unificada
df_info = pd.concat([df_controls, df_lc], axis=0)

# =======================
# 3. Loop por tipo de métrica
# =======================
for metric, paths in files.items():

    dfs = []

    for group, fname in paths.items():
        path = os.path.join(base_dir, fname)
        # df = pd.read_csv(path, index_col=0)

        df = pd.read_csv(path)

        # normalizar el nombre del ID
        if "Subject" in df.columns:
            df = df.rename(columns={"Subject": "ID"})
        elif "subject" in df.columns:
            df = df.rename(columns={"subject": "ID"})
        else:
            # por si algún CSV viene con el ID como índice
            df = df.reset_index().rename(columns={"index": "ID"})

        #
        # # Si ya existe una columna "subject"/"Subject", usarla como ID
        # for col in ["subject", "Subject", "SUBJECT"]:
        #     if col in df.columns:
        #         df = df.rename(columns={col: "ID"})
        #         break
        # else:
        #     # si no existe, el ID está en el índice
        #     df = df.reset_index().rename(columns={"index": "ID"})

        df["Grupo"] = group
        dfs.append(df)

    df_all = pd.concat(dfs, axis=0, ignore_index=True)
    # Limpieza IDs conocida
    df_all["ID"] = df_all["ID"].replace({"CP0107old": "CP0107"})
    # df_all = df_all[df_all["ID"] != "CP0106"]

    # Merge con edad y sexo
    df_all = df_all.merge(df_info, on=["ID", "Grupo"], how="left")

    # Limpieza de columnas
    df_all.columns = (
        df_all.columns
        .str.replace("\xa0", "", regex=True)
        .str.strip()
    )

    # Eliminar sujetos sin info
    df_all = df_all.dropna(subset=["edad", "sexo"])

    # Guardar
    outname = f"PET_{metric}_LC_Control_with_info.csv"
    df_all.to_csv(os.path.join(outdir, outname), index=False)

    print(f"✅ Guardado: {outname}")
    print(df_all["Grupo"].value_counts(), "\n")
