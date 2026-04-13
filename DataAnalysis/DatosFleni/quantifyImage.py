import os
import sys
import SimpleITK as sitk
import pandas as pd
import numpy as np
import re  # Importar el módulo de expresiones regulares

# =====================================================
# 1. 📂 CONFIGURACIÓN DE RUTAS DE MÓDULOS (Importación)
# =====================================================

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, '../..'))

if project_root not in sys.path:
    sys.path.append(project_root)
    print(f"✅ Añadida la raíz del proyecto para la importación: {project_root}")

try:
    from quantification_pipeline import run_quantification
    import PETquantification as quant
    import normalization as norm

    print("✅ Módulos de cuantificación importados correctamente.")
except ImportError as e:
    print(f"❌ ERROR de Importación: No se pudo importar la función o sus módulos: {e}")
    sys.exit(1)

# =====================================================
# 2. ⚙️ CONFIGURACIÓN DE EJECUCIÓN Y RUTAS DE DATOS
# =====================================================

# Directorio base para los archivos de imagen PET (ruta absoluta - sin cambios)
BASE_DIR = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025"

# --- Definición de Patrones y Rutas ---
# Definimos los patrones de búsqueda y dónde se encuentran.
# Usamos re.compile para compilar los patrones una sola vez.

SUBJECT_PATTERNS = [
    # 1. Sujetos NC (asumimos Fleni)
    {"folder_base": os.path.join(BASE_DIR, "Fleni"),
     "pattern_regex": re.compile(r'(NC\d+)')},  # Busca 'NC' seguido de dígitos

    # 2. Sujetos CP (asumimos COVID - ejemplo: CP0024)
    {"folder_base": os.path.join(BASE_DIR, "ProcesadosCOVID_PET_10_25"),
     "pattern_regex": re.compile(r'(CP\d+)')},  # Busca 'CP' seguido de dígitos

    # 3. Sujetos S (Control/Otros - ejemplo: s93895370F-5462-00001-000001)
    {"folder_base": os.path.join(BASE_DIR, "ProcesadosControlesCOVID_PET_10_25"),
     # ⭐️ PATRÓN CORREGIDO: Captura cualquier cadena que empiece con 's' y contenga dígitos,
     # letras mayúsculas, guiones y números (que es el nombre completo de la carpeta).
     "pattern_regex": re.compile(r'(s[a-zA-Z0-9-]+)')},
    # NOTA: Si los sujetos 'S' están en una carpeta diferente a 'ProcesadosCOVID...',
    #       cambia la variable 'folder_base' en este diccionario.
]

# --- Rutas de Atlas (Usando project_root/data/) ---
data_dir = os.path.join(project_root, 'data')

path_MNI_152_T1_gm_mask = os.path.join(data_dir, 'atlas', 'mni_icbm152_gm_mask.nii.gz')
path_Hammers = os.path.join(data_dir, 'atlas', "Hammers_mith-n30r95-MaxProbMap-gm-MNI152-SPM12.nii.gz")
path_DKT = os.path.join(data_dir, 'atlas', "Desikan_Killiany_MNI_SPM12.nii.gz")
path_CerebrA = os.path.join(data_dir, 'atlas', "CerebrA.nii")

# --- Archivos CSV de Etiquetas ---
LABELS_CSV_PATHS = {
    "Hammers": os.path.join(data_dir, 'labels', "labels_Hammers.csv"),
    "DKT": os.path.join(data_dir, 'labels', "labels_DKT.csv"),
    "CerebrA": os.path.join(data_dir, 'labels', "labels_CerebrA.csv"),
    "GM": os.path.join(data_dir, 'labels', "labels_gm_mask.csv"),
}

# Carpeta de salida final (se crea dentro de la carpeta actual 'datasets/')
OUTPUT_ROOT = os.path.join(BASE_DIR, "Quantification_Output")
os.makedirs(OUTPUT_ROOT, exist_ok=True)

# Configuramos el atlas a usar
ATLAS_NAME = "Hammers"

# --- Selección de Atlas ---
atlas_to_use = ATLAS_NAME
path_atlas = path_Hammers
path_labels_atlas_csv = LABELS_CSV_PATHS["Hammers"]

if atlas_to_use == "DKT":
    path_atlas = path_DKT
    path_labels_atlas_csv = LABELS_CSV_PATHS["DKT"]
elif atlas_to_use == "CerebrA":
    path_atlas = path_CerebrA
    path_labels_atlas_csv = LABELS_CSV_PATHS["CerebrA"]
elif atlas_to_use == "GM":
    path_atlas = path_MNI_152_T1_gm_mask
    path_labels_atlas_csv = LABELS_CSV_PATHS["GM"]

## 📁 Cargar Datos del Atlas
try:
    df_labels = pd.read_csv(path_labels_atlas_csv)
    print(f"✅ Atlas CSV cargado: {os.path.basename(path_labels_atlas_csv)}")
except FileNotFoundError:
    print(f"❌ ERROR: No se encontró el archivo de etiquetas del atlas en: {path_labels_atlas_csv}")
    sys.exit(1)


# =====================================================
# 4. 🚀 BUCLE PRINCIPAL DE EJECUCIÓN (Modificado)
# =====================================================

def get_subject_id(path, regex_pattern):
    """Extrae el ID del sujeto usando un patrón de regex."""
    match = regex_pattern.search(path)
    return match.group(1) if match else None


if __name__ == "__main__":

    print(f"\n--- INICIANDO CUANTIFICACIÓN con ATLAS: {ATLAS_NAME} ---")

    processed_count = 0

    for pattern_config in SUBJECT_PATTERNS:
        folder_base = pattern_config["folder_base"]
        regex = pattern_config["pattern_regex"]

        print(f"\n--- Buscando sujetos en: {os.path.basename(folder_base)} ---")

        # Iterar sobre las carpetas de sujetos dentro de la carpeta base
        try:
            for subject_folder in os.listdir(folder_base):
                # Intentar encontrar el ID del sujeto basado en el nombre de la carpeta
                subject_id = get_subject_id(subject_folder, regex)

                if subject_id:
                    # Construir la ruta completa del archivo PET (ESTO PUEDE NECESITAR AJUSTE)
                    # Asumo que el archivo PET se encuentra dentro de la carpeta del sujeto.

                    # EJEMPLO 1 (Para sujetos NC o CP):
                    # pet_path = os.path.join(folder_base, subject_folder, "PET_image_to_template_ANT.nii.gz")

                    # EJEMPLO 2 (Si la estructura es más compleja como el ejemplo anterior):
                    if subject_id.startswith('NC'):
                        pet_path = os.path.join(folder_base, subject_folder,
                                                "PET_image_to_template_ANT.nii.gz")  # Suponemos que el archivo está aquí
                    elif subject_id.startswith('CP'):
                        pet_path = os.path.join(folder_base, subject_folder,
                                                "PET_smoothed_image_to_template_ANTWarped.nii.gz")
                    elif subject_id.startswith('s'):
                        # Se necesita saber el nombre exacto del archivo PET para sujetos S
                        pet_path = os.path.join(folder_base, subject_folder,
                                                "PET_smoothed_image_to_template_ANTWarped.nii.gz")
                    else:
                        continue  # No coincide con los patrones definidos

                    if not os.path.exists(pet_path):
                        print(f"❌ Saltando {subject_id}: Archivo PET no encontrado en {pet_path}")
                        continue

                    print(f"\n=================================================")
                    print(f"  PROCESSING: {subject_id}")
                    print(f"=================================================")

                    # 1. Usar solo el subject_id como nombre de carpeta de salida
                    path_output_subject = os.path.join(OUTPUT_ROOT, subject_id)

                    # 2. Ejecutar la función de cuantificación
                    run_quantification(
                        path_PET_image=pet_path,
                        path_atlas=path_atlas,
                        df_labels_atlas=df_labels,
                        path_output=path_output_subject,
                        atlas_name=ATLAS_NAME,
                    )
                    processed_count += 1

        except FileNotFoundError:
            print(f"❌ ERROR: La carpeta base {folder_base} no fue encontrada. Revisa la ruta.")
            continue
        except Exception as e:
            print(f"❌ ERROR GENERAL procesando la carpeta {folder_base}: {e}")
            continue

    print("\n--- PROCESO FINALIZADO ---")
    print(f"Total de sujetos procesados: {processed_count}")
    print(f"Resultados guardados en: {os.path.abspath(OUTPUT_ROOT)}")