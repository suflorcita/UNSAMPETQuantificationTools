import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import linregress

# ==============================
# 3️⃣ Función auxiliar para graficar (solo regresiones por grupo)
# ==============================
def scatter_with_regression(df, y_col, y_label, filename):
    # Crear figura
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=df,
        x="Age",
        y=y_col,
        hue="Grupo",
        sizes=(40, 200),
        alpha=0.75,
        palette=palette,
    )

    # Líneas por grupo (sin la global)
    for group, color in palette.items():
        sns.regplot(
            data=df[df["Grupo"] == group],
            x="Age",
            y=y_col,
            scatter=False,
            color=color,
            line_kws={
                "lw": 1.8,
                "alpha": 0.9,
                "ls": "-",
                "label": f"Tendencia {group}"
            }
        )

    # Personalización
    plt.title(f"Relación entre {y_label.lower()} y edad por grupo", fontsize=14)
    plt.xlabel("Edad (años)", fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Guardar
    outpath = f"{base_dir}/{filename}"
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.show()


# ==============================
# 1️⃣ Cargar datos
# ==============================
base_dir = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/group_analysis_unificado"
csv_path = f"{base_dir}/PET_Fleni_COVID_GM_with_age.csv"

df = pd.read_csv(csv_path)

# ==============================
# 2️⃣ Limpieza básica
# ==============================
df["Age"] = pd.to_numeric(df["Age"], errors="coerce")
df = df.dropna(subset=["Age", "Regional uptake mean value GM", "Normalization to total brain mean value GM"])

# --- Paleta de colores ---
palette = {"FLENI": "#2E86C1", "COVID": "#E74C3C"}

# ==============================
# ==============================
# 4️⃣ Generar ambos gráficos
# ==============================

# A. Uptake regional sin normalizar
scatter_with_regression(
    df,
    y_col="Regional uptake mean value GM",
    y_label="Valor medio de captación en el total de la Materia Gris",
    filename="scatter_edad_vs_uptake_no_normalizado.png"
)

# B. Uptake normalizado al valor medio cerebral
scatter_with_regression(
    df,
    y_col="Normalization to total brain mean value GM",
    y_label="Valor medio de captación en el total de la Materia Gris normalizado a Captacion total",
    filename="scatter_edad_vs_uptake_normalizado.png"
)

# ==============================
# 6️⃣ Gráfico adicional con IDs sobre los puntos
# ==============================

plt.figure(figsize=(10, 8))
sns.scatterplot(
    data=df,
    x="Age",
    y="Regional uptake mean value",
    hue="Grupo",
    size="Normalization to total brain mean value",
    sizes=(40, 200),
    alpha=0.75,
    palette=palette
)

# --- Línea global ---
sns.regplot(
    data=df,
    x="Age",
    y="Regional uptake mean value",
    scatter=False,
    color="red",
    line_kws={"lw": 2, "label": "Regresión global"}
)

# --- Líneas por grupo ---
for group, color in palette.items():
    sns.regplot(
        data=df[df["Grupo"] == group],
        x="Age",
        y="Regional uptake mean value",
        scatter=False,
        color=color,
        line_kws={"lw": 1.5, "alpha": 0.7, "ls": "--", "label": f"Tendencia {group}"}
    )

# --- Superponer IDs ---
for i, row in df.iterrows():
    plt.text(
        row["Age"] + 0.2,  # pequeño desplazamiento horizontal
        row["Regional uptake mean value"],
        str(row["ID"]),
        fontsize=8,
        alpha=0.8
    )

# --- Personalización ---
plt.title("Captación regional vs edad con IDs de sujetos", fontsize=14)
plt.xlabel("Edad (años)", fontsize=12)
plt.ylabel("Valor medio de captación regional (SUVr)", fontsize=12)
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.grid(True, alpha=0.3)
plt.tight_layout()

# --- Guardar ---
outpath = f"{base_dir}/scatter_edad_vs_uptake_con_IDs.png"
plt.savefig(outpath, dpi=300, bbox_inches="tight")
plt.show()

print(f"✅ Gráfico con IDs guardado en: {outpath}")
