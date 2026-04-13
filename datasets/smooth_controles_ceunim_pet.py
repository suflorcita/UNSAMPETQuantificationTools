#!/usr/bin/env python3
# ============================================================
# Gaussian Smoothing and QC for Long COVID PET dataset
# ============================================================

import os, sys, ants
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Allow imports from parent directory (project root)
# ------------------------------------------------------------
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))
from registration_tool import gaussian_smooth


# ============================================================
# CONFIGURATION
# ============================================================
input_root = "/mnt/d87cc26d-5470-443c-81c1-e09b68ee4730/Sol/COVID/PET/BIDS_PET"
output_root = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/LongCOVID_Smoothed"
qc_root = os.path.join(output_root, "QC")

os.makedirs(output_root, exist_ok=True)
os.makedirs(qc_root, exist_ok=True)

cmap = "turbo"  # Clinical PET color map
dpi = 200





# ============================================================
# FUNCTIONS
# ============================================================


def resample_to_2x2x3p27(input_path, output_dir, prefix="resampled_2x2x3p27"):
    """
    Resample a NIfTI image to 2x2x3.27 mm voxel spacing using ANTsPy.
    """
    output_path = os.path.join(output_dir, f"{prefix}.nii.gz")

    # Load image
    img = ants.image_read(input_path)

    # Target spacing
    new_spacing = (2.0, 2.0, 3.27)

    # Resample image
    resampled = ants.resample_image(img, resample_params=new_spacing, use_voxels=False, interp_type=1)
    # interp_type=1 → linear interpolation

    # Save result
    ants.image_write(resampled, output_path)

    return output_path

def get_views(image):
    """Return axial, coronal and sagittal mid slices."""
    arr = sitk.GetArrayFromImage(image)
    z, y, x = arr.shape[0] // 2, arr.shape[1] // 2, arr.shape[2] // 2
    axial = np.flipud(arr[z, :, :])
    coronal = np.flipud(arr[:, y, :])
    sagittal = np.flipud(arr[:, :, x])
    return [axial, coronal, sagittal]


def generate_qc_plot(image_path, subj, fwhm_label, qc_dir, low_margin=1.2):
    """
    Generate 3-view QC plot for a PET image.
    Low-intensity voxels and background appear black (extended margin).
    Adaptive contrast with +10% upper range extension.
    """
    img = sitk.ReadImage(image_path)
    data = sitk.GetArrayFromImage(img).astype(np.float32)
    data[data < 0] = 0

    # Máscara cerebral
    mask = data > np.percentile(data, 20)
    brain_vals = data[mask]

    # Escala de intensidades
    vmin = np.percentile(brain_vals, 1)
    vmax = np.percentile(brain_vals, 99) * 1.10

    # 🔧 Margen adicional para el fondo negro
    cutoff = vmin * low_margin

    # Cortes centrales
    z, y, x = data.shape[0] // 2, data.shape[1] // 2, data.shape[2] // 2
    slices = [np.flipud(data[z, :, :]),
              np.flipud(data[:, y, :]),
              np.flipud(data[:, :, x])]

    fig, axes = plt.subplots(1, 3, figsize=(9, 4))
    fig.patch.set_facecolor("black")

    for ax, slc, title in zip(axes, slices, ["Axial", "Coronal", "Sagittal"]):
        # Enmascarar valores bajos → fondo negro con margen
        masked = np.ma.masked_where(slc <= cutoff, slc)
        im = ax.imshow(masked, cmap="turbo", vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_title(title, color="white", fontsize=9)
        ax.axis("off")
        ax.set_facecolor("black")

    # Barra de color
    cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(colors="white", labelsize=8)
    cbar.set_label("PET Intensity", color="white", fontsize=9)

    fig.suptitle(f"{subj} – Smoothed PET ({fwhm_label})",
                 color="white", fontsize=12, weight="bold")
    plt.subplots_adjust(left=0.05, right=0.9, top=0.88, bottom=0.05, wspace=0.05)

    out_png = os.path.join(qc_dir, f"{subj}_Smoothed_{fwhm_label}_QC.png")
    plt.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="black")
    plt.close()
    print(f"✅ Saved QC plot → {out_png}")

def generate_qc_plot_physical(image_path, subj, fwhm_label, qc_dir, low_margin=1.2):
    """
    Generate 3-view QC plot preserving voxel aspect ratio (e.g. for resampled PET 2x2x3.27 mm).
    """
    img = sitk.ReadImage(image_path)
    data = sitk.GetArrayFromImage(img).astype(np.float32)
    spacing = img.GetSpacing()  # (x, y, z) voxel size
    data[data < 0] = 0

    # Máscara cerebral
    mask = data > np.percentile(data, 20)
    brain_vals = data[mask]

    # Escala de intensidades
    vmin = np.percentile(brain_vals, 1)
    vmax = np.percentile(brain_vals, 99) * 1.10
    cutoff = vmin * low_margin

    # Cortes centrales
    z, y, x = data.shape[0] // 2, data.shape[1] // 2, data.shape[2] // 2
    slices = [np.flipud(data[z, :, :]), np.flipud(data[:, y, :]), np.flipud(data[:, :, x])]
    titles = ["Axial", "Coronal", "Sagittal"]

    # Aspect ratios basados en spacing (x, y, z)
    aspect_ratios = [
        spacing[1] / spacing[0],  # Axial
        spacing[2] / spacing[0],  # Coronal
        spacing[2] / spacing[1]   # Sagittal
    ]

    fig, axes = plt.subplots(1, 3, figsize=(9, 4))
    fig.patch.set_facecolor("black")

    for ax, slc, title, aspect in zip(axes, slices, titles, aspect_ratios):
        masked = np.ma.masked_where(slc <= cutoff, slc)
        im = ax.imshow(masked, cmap="turbo", vmin=vmin, vmax=vmax, interpolation="nearest", aspect=aspect)
        ax.set_title(title, color="white", fontsize=9)
        ax.axis("off")
        ax.set_facecolor("black")

    # Barra de color
    cbar_ax = fig.add_axes([0.92, 0.25, 0.02, 0.5])
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(colors="white", labelsize=8)
    cbar.set_label("PET Intensity", color="white", fontsize=9)

    fig.suptitle(f"{subj} – Smoothed PET ({fwhm_label})",
                 color="white", fontsize=12, weight="bold")
    plt.subplots_adjust(left=0.05, right=0.9, top=0.88, bottom=0.05, wspace=0.05)

    out_png = os.path.join(qc_dir, f"{subj}_{fwhm_label}_QC_physical.png")
    plt.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="black")
    plt.close()
    print(f"✅ Saved QC plot with physical proportions → {out_png}")





# ============================================================
# MAIN LOOP
# ============================================================
for subj in sorted(os.listdir(input_root)):
    subj_dir = os.path.join(input_root, subj, "pet")  # 👈 bajar al nivel de 'pet'
    if not os.path.isdir(subj_dir):
        continue

    # Buscar el PET reconstruido con corrección de atenuación (acstat)
    pet_files = [f for f in os.listdir(subj_dir) if f.endswith("rec-acstat_pet.nii.gz")]
    if not pet_files:
        print(f"⚠️ No PET found for {subj}")
        continue

    pet_file = os.path.join(subj_dir, pet_files[0])
    subj_out = os.path.join(output_root, subj)
    os.makedirs(subj_out, exist_ok=True)

    print(f"\n🧠 Processing {subj}")

    resampled_pet = resample_to_2x2x3p27(pet_file, subj_out, prefix="Raw_PET_resampled_2x2x3p27")

    # generate_qc_plot_physical(resampled_pet, subj, "Raw_resampled", qc_root, low_margin=100)

    # # Smooth with FWHM = 5 mm
    # smoothed5 = gaussian_smooth(pet_file, subj_out, prefix="smoothed_FWHM5", fwhm=5.0)
    # generate_qc_plot(smoothed5, subj, "FWHM5", qc_root, low_margin=100)

    # Smooth with FWHM = 5 mm
    smoothed5resampled = gaussian_smooth(resampled_pet, subj_out, prefix="smoothed_FWHM5_resampled", fwhm=5.0)
    generate_qc_plot_physical(smoothed5resampled, subj, "smoothed_FWHM5resampled", qc_root, low_margin=100)

    # # Smooth with FWHM = 6 mm
    # smoothed6 = gaussian_smooth(pet_file, subj_out, prefix="smoothed_FWHM6", fwhm=6.0)
    # generate_qc_plot(smoothed6, subj, "FWHM6", qc_root, low_margin=100)
    #
    # smoothed6resampled = gaussian_smooth(resampled_pet, subj_out, prefix="smoothed_FWHM6_resampled", fwhm=6.0)
    # generate_qc_plot_physical(smoothed6resampled, subj, "FWHM6resampled", qc_root, low_margin=100)
    #
    # # Smooth with FWHM = 8 mm
    # smoothed8 = gaussian_smooth(pet_file, subj_out, prefix="smoothed_FWHM8", fwhm=8.0)
    # generate_qc_plot(smoothed8, subj, "FWHM8", qc_root, low_margin=100)
    #
    # # Smooth with FWHM = 8 mm
    # smoothed8resampled = gaussian_smooth(resampled_pet, subj_out, prefix="smoothed_FWHM8_resampled", fwhm=8.0)
    # generate_qc_plot_physical(smoothed8resampled, subj, "FWHM8resampled", qc_root, low_margin=100)


print("\n✅ All Long COVID PET subjects processed successfully.")
