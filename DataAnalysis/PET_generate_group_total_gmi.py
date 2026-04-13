import os
import pandas as pd

# === Base directory ===
base_dir = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ProcesadosCOVID_PET_10_25_smoothed"

# === Lists for results ===
rows = []

for subject in os.listdir(base_dir):
    subj_path = os.path.join(base_dir, subject, "CSV", "subject_normalization_values_GM.csv")
    if os.path.exists(subj_path):
        df = pd.read_csv(subj_path, sep="\t|,", engine="python")
        df.columns = df.columns.str.strip()


        # Extract values (e.g. TotalGM or global mean)
        mean_uptake = df.loc[df["Structure"] == "Total_GM-", "Regional uptake mean values"].values
        norm_total = df.loc[df["Structure"] == "Total_GM-", "Normalization to total brain mean value"].values

        if len(mean_uptake) and len(norm_total):
            rows.append({
                "Subject": subject,
                "Regional uptake mean value GM": mean_uptake[0],
                "Normalization to total brain mean value GM": norm_total[0]
            })

# === Create DataFrame ===
df_result = pd.DataFrame(rows)

# === Save output ===
outdir = os.path.join(base_dir, "group_analysis")
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, "PET_summary_GM.csv")
df_result.to_csv(outfile, index=False)

print(f"✅ File saved to {outfile}")
