import os
import sys
import SimpleITK as sitk
import processPET4D as reg
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

try:
    from quantification_pipeline import run_quantification
    import PETquantification as quant
    import normalization as norm
    import registration_tool as regtools

    print("✅ Módulos de cuantificación importados correctamente.")
except ImportError as e:
    print(f"❌ ERROR de Importación: No se pudo importar la función o sus módulos: {e}")
    sys.exit(1)

path_my_images = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ProcesadosCOVID_PET_MRI_09_2025"
path_MNI_152_T1 = os.path.join(project_root, 'data', 'atlas', 'MNI152_T1_1mm.nii.gz')
path_MNI_152_T1_brain = os.path.join(project_root, 'data', 'atlas', 'MNI152_T1_1mm_brain.nii.gz')

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

# Configuramos el atlas a usar
ATLAS_NAME = "Hammers"
path_atlas = path_Hammers
path_labels_atlas_csv = LABELS_CSV_PATHS["Hammers"]

try:
    df_labels = pd.read_csv(path_labels_atlas_csv)
    print(f"✅ Atlas CSV cargado: {os.path.basename(path_labels_atlas_csv)}")
except FileNotFoundError:
    print(f"❌ ERROR: No se encontró el archivo de etiquetas del atlas en: {path_labels_atlas_csv}")
    sys.exit(1)


for subj in os.listdir(path_my_images):
    if not subj.startswith("CP"):
        continue



    subj_dir = os.path.join(path_my_images, subj)

    if not os.path.isdir(subj_dir):
        continue



    path_pet_image = os.path.join(subj_dir, f"Nifti/{subj}_rec-acstat-PET.nii.gz")
    path_mri_image = os.path.join(subj_dir, "Nifti", f"{subj}_t1mprage.nii.gz")


    if not os.path.exists(path_pet_image) or not os.path.exists(path_mri_image):
        continue

    output_smoothed_dir = os.path.join(subj_dir, "Smoothed_quantification")
    os.makedirs(output_smoothed_dir, exist_ok=True)

    smoothed_pet_path = regtools.gaussian_smooth(
        input_image=path_pet_image,
        output_dir=output_smoothed_dir,
        prefix=f"{subj}_PET_smoothed_5mm",
        fwhm=5.0
    )

    path_tx_pet_to_mri = os.path.join(subj_dir, f"Registration/tx_PET_2_RMN.tfm")
    path_tx_mri_to_mni152_flirt = os.path.join(subj_dir, f"FLIRT/T1_Norm_MNI_152_FLIRT.mat")
    path_tx_mri_to_mni152_ANT_1 = os.path.join(subj_dir, f"ANTs/T1_Norm_MNI_152_ANT1Warp.nii.gz")
    path_tx_mri_to_mni152_ANT_2 = os.path.join(subj_dir, f"ANTs/T1_Norm_MNI_152_ANT0GenericAffine.mat")



    # ================================
    # (1) PET smoothed → MRI
    # ================================
    img_pet_smoothed = sitk.ReadImage(smoothed_pet_path)
    img_mri_original = sitk.ReadImage(path_mri_image)

    pet_smoothed_in_mri = reg.read_and_apply_tx(
        tx_file=path_tx_pet_to_mri,
        image=img_pet_smoothed,
        ref_image=img_mri_original
    )

    out_smoothed_pet_mri = os.path.join(output_smoothed_dir, f"{subj}_PET_smoothed_5mm_in_MRI.nii.gz")
    sitk.WriteImage(pet_smoothed_in_mri, out_smoothed_pet_mri)

    print(f"✓ PET smoothed → MRI: {out_smoothed_pet_mri}")

    # ================================
    # (2) PET MRI → PET MNI 152
    # ================================

    # Add ANTs to PATH inside the script
    os.environ["ANTSPATH"] = "/usr/local/ANTs/bin"
    os.environ["PATH"] = f"{os.environ['ANTSPATH']}:" + os.environ["PATH"]

    if os.path.exists(path_tx_mri_to_mni152_ANT_1) and os.path.exists(path_tx_mri_to_mni152_ANT_2):
        print("✓ Applying ANTs transform")

        out_smoothed_pet_mni_ants = regtools.apply_ants(
            input_image=out_smoothed_pet_mri,  # PET en espacio MRI
            ref_image=path_MNI_152_T1_brain,  # referencia MNI152 T1 brain
            tx1=path_tx_mri_to_mni152_ANT_1,  # affine
            tx2=path_tx_mri_to_mni152_ANT_2,  # warp
            output_dir=output_smoothed_dir,
            prefix=f"{subj}_PET_smoothed_5mm_in_MNI152_ANTs"
        )

        print(f"✓ PET → MRI → MNI152 (ANTs): {out_smoothed_pet_mni_ants}")

        try:
            run_quantification(
                path_PET_image=out_smoothed_pet_mni_ants,
                path_atlas=path_atlas,
                df_labels_atlas=df_labels,
                path_output=output_smoothed_dir,
                atlas_name=ATLAS_NAME,
            )
            print(f"Quantification OK")
        except Exception as e:
            print(f"❌ ERROR GENERAL Cuantificando: {e}")
        continue