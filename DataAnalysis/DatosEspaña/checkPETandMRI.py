import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import SimpleITK as sitk

# Ruta base
path_my_images = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ProcesadosCOVID_PET_MRI_09_2025"
path_images_martin = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/PET_MRI_dataset_Martin"

output_dir = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ChequearMartin"

# --- FUNCIONES ---
def get_views(image):
    arr = sitk.GetArrayFromImage(image)
    if arr.ndim == 4:
        arr = arr[0]
    z, y, x = arr.shape[0] // 2, arr.shape[1] // 2, arr.shape[2] // 2
    axial = np.flipud(arr[z, :, :])
    coronal = np.flipud(np.fliplr(arr[:, y, :]))
    sagittal = np.flipud(arr[:, :, x])
    return [axial, coronal, sagittal]

def auto_contrast(slice_array):
    data = slice_array[slice_array > 0]
    return np.percentile(data, [1, 99]) if data.size else (0, 1)

def plot_compare(subj, img1_path, img2_path, modality):
    try:
        img1, img2 = sitk.ReadImage(img1_path), sitk.ReadImage(img2_path)
        slices1, slices2 = get_views(img1), get_views(img2)
        cmap = "gray" if modality == "T1" else "turbo"

        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        fig.patch.set_facecolor("black")
        fig.suptitle(f"{subj} – {modality}", color="white", fontsize=14)

        # --- Arriba: Sol ---
        for i, slc in enumerate(slices1):
            vmin, vmax = auto_contrast(slc)
            axes[0, i].imshow(slc, cmap=cmap, vmin=vmin, vmax=vmax)
            axes[0, i].set_title(f"{['Axial','Coronal','Sagittal'][i]} – Sol", color="white", fontsize=10)
            axes[0, i].axis("off")

        # --- Abajo: Martín ---
        for i, slc in enumerate(slices2):
            vmin, vmax = auto_contrast(slc)
            axes[1, i].imshow(slc, cmap=cmap, vmin=vmin, vmax=vmax)
            axes[1, i].set_title(f"{['Axial','Coronal','Sagittal'][i]} – Martín", color="white", fontsize=10)
            axes[1, i].axis("off")

        plt.subplots_adjust(wspace=0.02, hspace=0.08)
        out_path = os.path.join(output_dir, f"{subj}_{modality}_Compare.png")
        plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close()
        print(f"✅ Saved comparison for {subj} ({modality})")

    except Exception as e:
        print(f"⚠️ Error with {subj} ({modality}): {e}")

# --- MAIN ---
my_files = {f.split("_")[0] for f in os.listdir(path_my_images) if f.startswith("CP")}
martin_files = {f.split("_")[0] for f in os.listdir(path_images_martin) if f.startswith("CP")}

subjects = sorted(my_files & martin_files)
print(f"🔍 Found {len(subjects)} subjects in both folders.")

# --- NUEVO: Chequeo de diferencias ---
missing_in_martin = sorted(my_files - martin_files)
missing_in_me = sorted(martin_files - my_files)

if missing_in_martin:
    print(f"\n Sujetos que están en MIS imágenes pero no en las de Martín:\n{', '.join(missing_in_martin)}")
if missing_in_me:
    print(f"\n Sujetos que están en las de Martín pero no en MIS imágenes:\n{', '.join(missing_in_me)}")

if not missing_in_martin and not missing_in_me:
    print("\n✅ Todos los sujetos coinciden entre ambas carpetas.\n")

for subj in os.listdir(path_my_images):
    if not subj.startswith("CP"):
        continue

    subj_dir = os.path.join(path_my_images, subj)
    if not os.path.isdir(subj_dir):
        continue

    # --- Tus archivos (dentro de subcarpetas) ---
    t1_mine = os.path.join(subj_dir, f"Nifti/{subj}_t1mprage.nii.gz")
    pet_mine = os.path.join(subj_dir, f"Nifti/{subj}_rec-acstat-PET.nii.gz")

    # --- Archivos de Martín (carpeta plana) ---
    t1_martin = os.path.join(path_images_martin, f"{subj}_t1mprage.nii.gz")
    pet_martin = os.path.join(path_images_martin, f"{subj}_FDG_pet.nii.gz")

    # --- Comparación T1 ---
    if os.path.exists(t1_mine) and os.path.exists(t1_martin):
        plot_compare(subj, t1_mine, t1_martin, "T1")
    else:
        if not os.path.exists(t1_mine):
            print(f"Falta T1 mío para {subj}: {t1_mine}")
        if not os.path.exists(t1_martin):
            print(f"Falta T1 de Martín para {subj}: {t1_martin}")

    # --- Comparación PET ---
    if os.path.exists(pet_mine) and os.path.exists(pet_martin):
        plot_compare(subj, pet_mine, pet_martin, "PET")
    else:
        if not os.path.exists(pet_mine):
            print(f"Falta PET mío para {subj}: {pet_mine}")
        if not os.path.exists(pet_martin):
            print(f"Falta PET de Martín para {subj}: {pet_martin}")

