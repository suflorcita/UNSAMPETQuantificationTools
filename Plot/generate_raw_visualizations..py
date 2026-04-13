import os
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt

# =====================================================
# CONFIGURATION
# =====================================================
datasets = {
    "FLENI": "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/Fleni",
    "COVID": "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ProcesadosCOVID_PET_10_25",
}

raw_filename = "raw_pet_image.nii.gz"
colormap = "turbo"   # o "turbo" para un estilo más moderno
dpi = 300


# =====================================================
# FUNCTIONS
# =====================================================


def get_views(image):
    """Return axial, coronal and sagittal mid-slices from a 3D or 4D PET image."""
    arr = sitk.GetArrayFromImage(image)
    if arr.ndim == 4:
        arr = arr[0]
    z, y, x = arr.shape[0] // 2, arr.shape[1] // 2, arr.shape[2] // 2
    axial = np.flipud(arr[z, :, :])
    coronal = np.flipud(np.fliplr(arr[:, y, :]))
    sagittal = np.flipud(arr[:, :, x])
    return [axial, coronal, sagittal]


def auto_contrast(slice_array):
    """Return contrast limits based on 1st–99th percentile."""
    data = slice_array[slice_array > 0]
    if data.size == 0:
        return 0, 1
    return np.percentile(data, [1, 99])


def generate_raw_qc(subj, subj_dir, output_folder):
    """Generate QC figure for a single subject."""
    file_path = os.path.join(subj_dir, raw_filename)
    if not os.path.exists(file_path):
        print(f"❌ Missing raw PET for {subj}")
        return

    try:
        img = sitk.ReadImage(file_path)
        spacing = img.GetSpacing()
        slices = get_views(img)

        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        fig.patch.set_facecolor("black")

        for i, (view, slc) in enumerate(zip(["Axial", "Coronal", "Sagittal"], slices)):
            vmin, vmax = auto_contrast(slc)

            # Aspect ratio from voxel spacing
            if view == "Axial":
                aspect = spacing[1] / spacing[0]
            elif view == "Coronal":
                aspect = spacing[2] / spacing[0]
            else:
                aspect = spacing[2] / spacing[1]

            # Clean background
            slc_clean = np.where(slc < 0.01 * np.max(slc), 0, slc)
            masked = np.ma.masked_where(slc_clean == 0, slc_clean)

            im = axes[i].imshow(masked, cmap=colormap, vmin=vmin, vmax=vmax, aspect=aspect)
            axes[i].set_facecolor("black")
            axes[i].set_title(view, color="white", fontsize=10)
            axes[i].axis("off")

        # Colorbar
        cbar_ax = fig.add_axes([0.93, 0.25, 0.02, 0.5])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.ax.tick_params(colors="white")
        cbar.set_label("PET Intensity", color="white")

        fig.suptitle(f"{subj} – Raw PET (QC views)", color="white", fontsize=12)
        plt.subplots_adjust(left=0.03, right=0.9, top=0.88, bottom=0.05, wspace=0.05)

        os.makedirs(output_folder, exist_ok=True)
        out_path = os.path.join(output_folder, f"{subj}_RawPET_QC.png")
        plt.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close()
        print(f"✅ Saved QC image for {subj}: {out_path}")

    except Exception as e:
        print(f"⚠️ Error loading {subj}: {e}")


# =====================================================
# MAIN LOOP FOR BOTH DATASETS
# =====================================================
for dataset_name, input_folder in datasets.items():
    print(f"\n================ {dataset_name} =================")
    output_plot_folder = os.path.join(input_folder, "QC_raw_" + dataset_name)
    os.makedirs(output_plot_folder, exist_ok=True)

    for subj in sorted(os.listdir(input_folder)):
        subj_dir = os.path.join(input_folder, subj)
        if not os.path.isdir(subj_dir) or subj.startswith("QC"):
            continue
        print(f"📁 Processing {dataset_name} subject: {subj}")
        generate_raw_qc(subj, subj_dir, output_plot_folder)

print("\n✅ All QC plots completed for both datasets.")
