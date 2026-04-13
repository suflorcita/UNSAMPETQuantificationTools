import os
import re
import pandas as pd
import numpy as np
from neuroHarmonize.harmonizationLearn import harmonizationLearn
import sys

# =====================================================
# 1. ⚙️ CONFIGURACIÓN
# =====================================================

CONSOLIDATED_CSV_DIR = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/Quantification_Output/DatosEspaña"

CSV_FILES_TO_PROCESS = [
    "PET_cerebellum_norm_by_subject_Hammers_mri_fleni.csv",
    "PET_mean_values_by_subject_Hammers_fleni.csv",
    "PET_uptake_values_by_subject_Hammers_mri_fleni.csv",
]

SUBJECT_ID_COLUMN = "Subject"
FLENI_REGEX = re.compile(r'^(NC\d+)')
# Usamos un sufijo diferente para distinguir el método de armonización
HARMONIZED_CSV_SUFFIX = "_NEUROHARMONIZED.csv"


# =====================================================
# 2. 🧠 FUNCIONES DE CLASIFICACIÓN Y ARMONIZACIÓN
# =====================================================

def classify_subject_id(subject_id):
    """Clasifica el ID del sujeto como FLENI o CEUNIM."""
    if isinstance(subject_id, str):
        if FLENI_REGEX.match(subject_id):
            return "FLENI"
        elif subject_id.startswith('CP') or subject_id.startswith('s'):
            return "CEUNIM"
    return "UNKNOWN"


def apply_harmonization_neuroHarmonize(df, id_column, batch_column='Group'):
    """
    Aplica el algoritmo ComBat usando neuroHarmonize (función harmonizationLearn)
    y reconstruye el DataFrame usando pd.concat para seguridad.
    """

    # 1. PREPARACIÓN DE DATOS

    # Columnas de identificación (Subject y Group)
    identifier_cols = [id_column, batch_column]

    # Columnas de datos (ROIs)
    data_cols = [col for col in df.columns if col not in identifier_cols]

    if not data_cols:
        raise ValueError("No se encontraron columnas de datos numéricos (ROIs) para armonizar.")

    # Crear la Matriz de datos (Sujetos x ROIs) - Formato esperado
    data_array = df[data_cols].values.astype('float64')

    # Matriz de covariables (Batch/Sitio)
    covariates_df = df[[batch_column]].copy()
    # Solución al error 'SITE'
    covariates_df.rename(columns={batch_column: 'SITE'}, inplace=True)

    # 2. ARMONIZACIÓN

    print(f"🔎 Iniciando neuroHarmonize. Data: {data_array.shape[0]} Sujetos x {data_array.shape[1]} ROIs.")

    # ⚠️ USAMOS EL NOMBRE CORREGIDO: harmonizationLearn
    _, data_harmonized_array = harmonizationLearn(
        data_array,
        covariates_df
    )

    print("🔎 Armonización completada.")

    # 3. RECONSTRUCCIÓN (Soluciona el error 'length 92; 2 is required')
    # 3a. Convertir el array armonizado de vuelta a DataFrame de solo datos
    # Mantenemos el índice original del DataFrame (df.index) para asegurar la alineación
    df_harmonized_data = pd.DataFrame(data_harmonized_array, columns=data_cols, index=df.index)


    # 3b. Aislar las columnas de identificación del DataFrame original (Subject, Group)
    df_identifiers = df[identifier_cols]

    # 3c. Concatenar de forma segura (alineando por el índice)
    df_harmonized = pd.concat([df_identifiers, df_harmonized_data], axis=1)

    return df_harmonized


# =====================================================
# 3. 📌 FUNCIÓN PRINCIPAL PARA PROCESAR CADA CSV
# =====================================================

def process_single_csv(file_path, id_column, suffix):
    filename = os.path.basename(file_path)
    print(f"\n--- 🚀 Procesando: {filename} ---")

    # --------------------------
    # 1. CARGA Y LIMPIEZA
    # --------------------------
    df = pd.read_csv(file_path, header=0)

    cols_to_drop = []

    # Limpieza: Unnamed, Guion ('-'), y 100% NaN
    unnamed_cols = [col for col in df.columns if 'Unnamed:' in str(col)]
    cols_to_drop.extend(unnamed_cols)

    if '-' in df.columns:
        cols_to_drop.append('-')
        print("   ⚠️ Columna '-' eliminada.")

    nan_cols = df.columns[df.isnull().all()].tolist()
    cols_to_drop.extend(nan_cols)

    if cols_to_drop:
        df = df.drop(columns=list(set(cols_to_drop)))
        print(f"   ⚠️ Columnas eliminadas en limpieza: {list(set(cols_to_drop))}")

    if id_column not in df.columns:
        print(f"   ❌ ERROR: La columna '{id_column}' no existe en el CSV.")
        return False

    # --------------------------
    # 2. CLASIFICAR SUJETOS (BATCH) Y CONVERTIR A NUMÉRICO
    # --------------------------
    df['SITE'] = df[id_column].apply(classify_subject_id)
    print("   ✅ Columna 'Group' añadida.")

    # Convertir columnas de ROI a float, ignorando errores (se convierten a NaN)
    data_cols = [col for col in df.columns if col not in [id_column, 'SITE']]
    df[data_cols] = df[data_cols].apply(pd.to_numeric, errors='coerce')

    # --------------------------
    # 3. FILTRAR UNKNOWN Y FILAS CON NaN (Importante para neuroHarmonize)
    # --------------------------

    # Eliminar sujetos 'UNKNOWN'
    df_filtered = df[df['SITE'] != 'UNKNOWN'].copy()
    num_unknown_dropped = len(df) - len(df_filtered)

    # Eliminar filas (sujetos) que contengan cualquier NaN en los datos de ROI (ComBat no maneja NaNs)
    df_cleaned = df_filtered.dropna(subset=data_cols).copy()
    num_nan_dropped = len(df_filtered) - len(df_cleaned)

    if num_unknown_dropped > 0:
        print(f"   ⚠️ Se eliminaron {num_unknown_dropped} sujetos 'UNKNOWN'.")
    if num_nan_dropped > 0:
        print(f"   ⚠️ Se eliminaron {num_nan_dropped} sujetos con valores NaN en los datos de ROI.")

    if len(df_cleaned) == 0:
        print("   ❌ No quedan sujetos válidos para la armonización.")
        return False

    print(f"   ✔️ Sujetos válidos finales para armonización: {len(df_cleaned)}")

    # --------------------------
    # 4. EJECUTAR NEUROHARMONIZE
    # --------------------------
    df_harmonized = apply_harmonization_neuroHarmonize(df_cleaned, id_column, batch_column='SITE')

    # --------------------------
    # 5. GUARDAR
    # --------------------------
    out_path = os.path.join(os.path.dirname(file_path),
                            os.path.splitext(filename)[0] + suffix)

    df_harmonized.to_csv(out_path, index=False)
    print(f"   ✔️ Guardado archivo armonizado en:\n     {out_path}")

    return True


# =====================================================
# 4. 🚀 EJECUCIÓN
# =====================================================

if __name__ == "__main__":
    try:
        total_processed = 0

        for filename in CSV_FILES_TO_PROCESS:
            file_path = os.path.join(CONSOLIDATED_CSV_DIR, filename)

            if not os.path.exists(file_path):
                print(f"❌ ADVERTENCIA: Archivo no encontrado: {filename}. Saltando.")
                continue

            if process_single_csv(file_path, SUBJECT_ID_COLUMN, HARMONIZED_CSV_SUFFIX):
                total_processed += 1

        print("\n--- 🏁 PROCESO COMPLETO FINALIZADO ---")
        print(f"✔️ Archivos armonizados exitosamente: {total_processed} de {len(CSV_FILES_TO_PROCESS)}")

    except ImportError as e:
        if 'neuroHarmonize' in str(e):
            print("\n❌ ERROR CRÍTICO: Falta la librería neuroHarmonize.")
            print("   Instala con: pip install neuroHarmonize")
        else:
            print(f"\n❌ ERROR CRÍTICO de importación: {e}")
        sys.exit(1)

    except Exception as e:
        print(f"\n❌ ERROR FATAL INESPERADO: {e}")
        sys.exit(1)