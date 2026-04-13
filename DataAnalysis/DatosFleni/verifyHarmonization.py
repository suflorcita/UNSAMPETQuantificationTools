import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

# =========================================================================
# ⚙️ CONFIGURACIÓN
# =========================================================================

# Directorio donde se encuentran los archivos armonizados
CONSOLIDATED_CSV_DIR = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/Quantification_Output/DatosEspaña"

# Sufijo de los archivos armonizados
HARMONIZED_CSV_SUFFIX = "_NEUROHARMONIZED.csv"

# Columna que identifica los grupos/sitios (la columna 'Group' que creamos)
GROUP_COLUMN = 'SITE'

# Columna que identifica a los sujetos
ID_COLUMN = 'Subject'

# Directorio donde se guardarán las imágenes de los boxplots
OUTPUT_DIR = os.path.join(CONSOLIDATED_CSV_DIR, 'Boxplots_Sitios_Armonizados')


# =========================================================================

def generar_boxplots_por_region(df, filename):
    """
    Genera boxplots para todas las ROIs en el DataFrame, comparando la
    distribución entre los sitios de adquisición.
    """

    # Identificar las columnas de las ROIs
    data_cols = [col for col in df.columns if col not in [ID_COLUMN, GROUP_COLUMN]]

    if not data_cols:
        print("❌ ADVERTENCIA: No se encontraron columnas de datos (ROIs) en el archivo.")
        return

    # Creamos un DataFrame largo (tidy format) para Seaborn para todas las ROIs
    df_long = pd.melt(df,
                      id_vars=[ID_COLUMN, GROUP_COLUMN],
                      value_vars=data_cols,
                      var_name='ROI',
                      value_name='Valor Armonizado')

    # Nombre base para la carpeta de resultados (ej: PET_cerebellum_norm)
    base_name = filename.replace(HARMONIZED_CSV_SUFFIX, '')
    output_subdir = os.path.join(OUTPUT_DIR, base_name)

    if not os.path.exists(output_subdir):
        os.makedirs(output_subdir)

    print(f"🔎 Generando {len(data_cols)} boxplots en la carpeta: {output_subdir}")

    for roi in data_cols:
        plt.figure(figsize=(6, 5))

        # Filtrar datos para la ROI actual
        df_roi = df_long[df_long['ROI'] == roi]

        # Generar el Boxplot comparando los grupos (sitios)
        sns.boxplot(x=GROUP_COLUMN, y='Valor Armonizado', data=df_roi, palette="Set1")
        # Añadir los puntos individuales para mejor visualización
        sns.stripplot(x=GROUP_COLUMN, y='Valor Armonizado', data=df_roi, color=".3", size=4, jitter=True)

        plt.title(f'Armonizado: Distribución por Sitio para {roi}', fontsize=12)
        plt.xlabel('Sitio de Adquisición', fontsize=10)
        plt.ylabel('Valor Armonizado', fontsize=10)

        # Guardar la figura
        file_name = os.path.join(output_subdir, f'boxplot_{roi}.png')
        plt.savefig(file_name, bbox_inches='tight')
        plt.close()

    print(f"✅ Boxplots completados para {filename} y guardados en {output_subdir}.\n")


# --- 🏃 FUNCIÓN PRINCIPAL DE EJECUCIÓN ---
def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print(f"Buscando archivos armonizados en: {CONSOLIDATED_CSV_DIR}")

    # Encontrar todos los archivos armonizados en el directorio
    harmonized_files = [f for f in os.listdir(CONSOLIDATED_CSV_DIR) if f.endswith(HARMONIZED_CSV_SUFFIX)]

    if not harmonized_files:
        print(f"❌ ERROR: No se encontraron archivos con el sufijo '{HARMONIZED_CSV_SUFFIX}'.")
        print("Asegúrate de ejecutar primero el script de armonización.")
        return

    print(f"Archivos armonizados encontrados: {len(harmonized_files)}")
    print("-" * 50)

    for filename in harmonized_files:
        file_path = os.path.join(CONSOLIDATED_CSV_DIR, filename)

        try:
            df_harmonized = pd.read_csv(file_path)
            print(f"Procesando: {filename}. Dimensiones: {df_harmonized.shape}")
            generar_boxplots_por_region(df_harmonized, filename)

        except Exception as e:
            print(f"❌ ERROR al procesar {filename}: {e}")


if __name__ == '__main__':
    main()