import os
import pandas as pd

# =======================
# 1. Paths
# =======================
base_dir = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025"
info_csv = os.path.join(base_dir, "Datos-UNSAM-FLENI-DatosParticipantes.csv")

flenipath = os.path.join(base_dir, "Fleni", "group_analysis", "PET_summary_GM.csv")
covidpath = os.path.join(base_dir, "ProcesadosCOVID_PET_10_25_smoothed", "group_analysis", "PET_summary_GM.csv")

outdir = os.path.join(base_dir, "group_analysis_unificado")
os.makedirs(outdir, exist_ok=True)

# =======================
# 2. Read PET summaries
# =======================
df_fleni = pd.read_csv(flenipath, index_col=0)
df_covid = pd.read_csv(covidpath, index_col=0)

# Assign group
df_fleni["Grupo"] = "FLENI"
df_covid["Grupo"] = "COVID"

# =======================
# 3. Combine datasets
# =======================
df_comb = pd.concat([df_fleni, df_covid], axis=0).reset_index().rename(columns={"index": "ID"})

# Clean columns
df_comb.columns = df_comb.columns.str.replace("\xa0", "", regex=True).str.strip()
df_comb = df_comb.rename(columns={"Subject": "ID"})

df_comb = df_comb.dropna(how="all")

# Remove prefix/suffix from IDs to match info CSV
df_comb["ID"] = (
    df_comb["ID"].astype(str)
    .str.strip()
    .str.replace("^sub-", "", regex=True)
    .str.replace("_old$", "", regex=True)
)

# =======================
# 4. Read participant info (Edad + Sexo)
# =======================
df_info = pd.read_csv(info_csv, skip_blank_lines=True)
df_info = df_info.dropna(how="all").dropna(axis=1, how="all")

# Find header row (where second column == 'ID')
header_row = df_info[df_info.iloc[:, 1] == "ID"].index
if len(header_row) > 0:
    header_idx = header_row[0]
    df_info.columns = df_info.iloc[header_idx]
    df_info = df_info.iloc[header_idx + 1:]
else:
    print("⚠️ Header row with 'ID' not found — check manually.")

df_info.columns = df_info.columns.str.strip()

# Keep relevant columns
df_info = df_info.rename(columns={
    "ID": "ID",
    "Sexo": "Sex",
    "Edad": "Age"
})[["ID", "Sex", "Age"]]

# Clean types
df_info["Age"] = pd.to_numeric(df_info["Age"], errors="coerce")
df_info["Sex"] = df_info["Sex"].astype(str).str.strip().str.lower()
df_info["Sex"] = df_info["Sex"].replace({
    "m": "M", "f": "F",
    "masculino": "M", "femenino": "F",
    "male": "M", "female": "F"
})

# =======================
# 5. Merge PET + participant info
# =======================
df_final = df_comb.merge(df_info, on="ID", how="left")

# Drop subjects without age info
df_final = df_final.dropna(subset=["Age"])

print("📊 Group distribution:")
print(df_final["Grupo"].value_counts())
print(f"\nTotal subjects with age info: {len(df_final)}")

# =======================
# 6. Save merged dataset
# =======================
outfile = os.path.join(outdir, "PET_Fleni_COVID_GM_with_age.csv")
df_final.to_csv(outfile, index=False)
print(f"\n✅ Combined file with age saved at: {outfile}")
