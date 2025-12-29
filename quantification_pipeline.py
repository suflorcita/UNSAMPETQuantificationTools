#!/usr/bin/env python3
import os
import SimpleITK as sitk
import pandas as pd
import normalization as norm
import PETquantification as quant


def run_quantification(
    path_PET_image,
    path_atlas,
    df_labels_atlas,
    path_output,
    atlas_name="Hammers",
    path_MNI152_PET=None,
):
    """
    Ejecuta la cuantificación de FDG-PET para un sujeto dado.
    Incluye normalización por valor medio y por cerebelo (si aplica),
    comparación con el atlas MNI152 y generación de CSVs e imágenes sintéticas.
    """

    os.makedirs(path_output, exist_ok=True)
    path_CSV = os.path.join(path_output, "CSV")
    os.makedirs(path_CSV, exist_ok=True)
    path_norm = os.path.join(path_output, "Normalization")
    os.makedirs(path_norm, exist_ok=True)
    path_synth = os.path.join(path_output, "Synthetic_Image")
    os.makedirs(path_synth, exist_ok=True)

    image_PET = sitk.ReadImage(path_PET_image)
    image_atlas = sitk.ReadImage(path_atlas)

    # --- Cuantificación base ---
    df_intensity, image_intensity = quant.PET_FDG_quantification(
        image_PET, image_atlas, df_labels_atlas, atlas=atlas_name
    )

    df_intensity.to_csv(os.path.join(path_CSV, f"subject_uptake_{atlas_name}.csv"))
    values_mean = df_intensity["mean_PET"]
    n_labels = df_intensity["n_label"]
    name_structures = (
        df_intensity["structure"]
        if atlas_name == "DKT"
        else df_intensity["structure"] + "-" + df_intensity["hemisphere"]
    )

    # --- Normalización al valor medio ---
    norm_img_mean = norm.intensity_normalization(image_PET)
    sitk.WriteImage(norm_img_mean, os.path.join(path_norm, "intensity_normalization_mean.nii.gz"))
    df_mean, _ = quant.PET_FDG_quantification(norm_img_mean, image_atlas, df_labels_atlas, atlas=atlas_name)
    values_mean_norm = df_mean["mean_PET"]

    # --- Normalización al cerebelo (si aplica) ---
    if atlas_name != "GM":
        if atlas_name == "DKT":
            cerebellum_name = "Cerebellum-Cortex"
        elif atlas_name == "CerebrA":
            cerebellum_name = "Cerebellum-Gray-Matter"
        else:
            cerebellum_name = "cerebellum"

        if atlas_name == "DKT":
            cereb_val = float(
                df_intensity.loc[df_intensity["structure"] == cerebellum_name, "mean_PET"].iloc[0]
            )
        else:
            L = df_intensity.query("structure == @cerebellum_name and hemisphere == 'L'")["mean_PET"].iloc[0]
            R = df_intensity.query("structure == @cerebellum_name and hemisphere == 'R'")["mean_PET"].iloc[0]
            cereb_val = (L + R) / 2

        norm_img_cereb = norm.intensity_normalization(image_PET, mode="scalar", scalar=cereb_val)
        sitk.WriteImage(norm_img_cereb, os.path.join(path_norm, f"intensity_normalization_cerebellum_{atlas_name}.nii.gz"))
        df_cereb, _ = quant.PET_FDG_quantification(norm_img_cereb, image_atlas, df_labels_atlas, atlas=atlas_name)
        values_cereb_norm = df_cereb["mean_PET"]
    else:
        values_cereb_norm = None

    # --- Crear DataFrame final de normalización ---
    df_subject = pd.DataFrame({
        "Structure": name_structures,
        "n_label": n_labels,
        "Regional uptake mean values": values_mean,
        "Normalization to total brain mean value": values_mean_norm,
    })
    if values_cereb_norm is not None:
        df_subject["Normalization to cerebellum uptake values"] = values_cereb_norm

    df_subject.to_csv(os.path.join(path_CSV, f"subject_normalization_values_{atlas_name}.csv"))

    # --- Guardar imágenes sintéticas ---
    sitk.WriteImage(image_intensity, os.path.join(path_synth, "synthetic_image_subject.nii.gz"))

    print(f"✅ Quantificación completada para {atlas_name} en {path_output}")
    return df_subject
