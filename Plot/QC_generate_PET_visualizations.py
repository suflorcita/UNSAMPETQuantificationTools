import os
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt

# =====================================================
# CONFIG
# =====================================================
input_folder = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ProcesadosControlesCOVID_PET_10_25"
output_plot_folder = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ProcesadosControlesCOVID_PET_10_25/QC"
os.makedirs(output_plot_folder, exist_ok=True)

# --- Base PETs ---
base_patterns = {
    "PET_image_to_template_flirt.nii.gz": "FLIRT (Raw PET)",
    "PET_image_to_template_ANT.nii.gz": "ANTs (Raw PET)",
    "PET_smoothed_image_to_template_flirt.nii.gz": "FLIRT (Smoothed)",
    "PET_smoothed_image_to_template_ANTWarped.nii.gz": "ANTs (Smoothed)",
}

# --- Normalized PETs ---
norm_patterns = {
    os.path.join("Normalization", "intensity_normalization_mean.nii.gz"): "Mean Normalized",

    # Cualquiera de estos nombres es válido
    os.path.join("Normalization", "intensity_normalization_cerebellum_Hammers.nii.gz"): "Cerebellum Normalized",
    os.path.join("Normalization", "intensity_normalization_cerebellum.nii.gz"): "Cerebellum Normalized",
}

# --- Unique + Summed PET ---
unique_patterns = {
    "unique_fdg_pet_image.nii.gz": "Unique FDG PET (Frame 1)",
    "Summed_PET_image.nii.gz": "Summed PET",
}

# Visualization parameters
cmap = "plasma"
vmin_default, vmax_default = 0, 10


# =====================================================
# FUNCTIONS
# =====================================================
def get_views(image, frame_index=0):
    """
    Return correctly oriented axial, coronal and sagittal slices.
    If 4D, show only one frame (default: first).
    """
    arr = sitk.GetArrayFromImage(image)

    # Handle 4D PET
    if arr.ndim == 4:
        n_frames = arr.shape[0]
        if frame_index < 0 or frame_index >= n_frames:
            frame_index = 0
        print(f"🧩 Using frame {frame_index+1}/{n_frames} (shape={arr.shape})")
        arr = arr[frame_index, :, :, :]

    # Middle slices
    z, y, x = arr.shape[0] // 2, arr.shape[1] // 2, arr.shape[2] // 2

    # Orientation corrections
    axial = np.flipud(arr[z, :, :])
    coronal = np.flipud(np.fliplr(arr[:, y, :]))
    sagittal = np.flipud(arr[:, :, x])

    return [axial, coronal, sagittal]


def auto_contrast(slice_array):
    """Compute percentile-based contrast limits."""
    data = slice_array[slice_array > 0]
    if data.size > 0:
        return np.percentile(data, [1, 99])
    else:
        return vmin_default, vmax_default


def generate_grid(subj, subj_dir, file_patterns, grid_name):
    """Generate a grid of PET views."""
    fig, axes = plt.subplots(len(file_patterns), 3, figsize=(9, 3 * len(file_patterns)))
    fig.patch.set_facecolor("black")

    for row, (fname, label) in enumerate(file_patterns.items()):
        file_path = os.path.join(subj_dir, fname)
        if not os.path.exists(file_path):
            for col in range(3):
                axes[row, col].axis("off")
                axes[row, col].set_title("Missing", color="gray")
            print(f"❌ Missing {fname}")
            continue

        try:
            img = sitk.ReadImage(file_path)
            # Show first frame only for 4D unique_fdg_pet_image
            if "unique_fdg_pet_image" in fname:
                slices = get_views(img, frame_index=0)
            else:
                slices = get_views(img)

            for col, slc in enumerate(slices):
                #vmin, vmax = auto_contrast(slc)
                vmin, vmax = 0, 6.5
                im = axes[row, col].imshow(
                    np.ma.masked_where(slc == 0, slc),
                    cmap=cmap, vmin=vmin, vmax=vmax
                )
                axes[row, col].axis("off")

                view_name = ["Axial", "Coronal", "Sagittal"][col]
                subtitle = f"{label} – {view_name}"
                axes[row, col].set_title(subtitle, color="white", fontsize=8)

        except Exception as e:
            print(f"⚠️ Error loading {fname}: {e}")
            for col in range(3):
                axes[row, col].axis("off")
                axes[row, col].set_title("Error", color="red")

    # Colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(colors="white")
    cbar.set_label("PET Intensity", color="white")

    # Title and save
    fig.suptitle(f"{subj} – {grid_name} (views)", color="white", fontsize=12)
    plt.subplots_adjust(left=0.05, right=0.9, top=0.93, bottom=0.05, wspace=0.05, hspace=0.2)

    output_path = os.path.join(output_plot_folder, f"{subj}_{grid_name}.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"✅ Saved: {output_path}")


# =====================================================
# MAIN LOOP
# =====================================================
for subj in sorted(os.listdir(input_folder)):
    print(subj)
    # if subj.startswith("CP0"):
    # subj_dir = os.path.join(input_folder, subj, "Smoothed_quantification")
    subj_dir = os.path.join(input_folder, subj)

    # if not os.path.isdir(subj_dir):
    #     continue

    print(f"\n📁 Processing subject: {subj}")

    # Grid 1: Base PETs (FLIRT + ANTs)
    # generate_grid(subj, subj_dir, base_patterns, "BasePET_Grid")

    # Grid 2: Normalization maps
    if os.path.exists(os.path.join(subj_dir, "Normalization")):
        generate_grid(subj, subj_dir, norm_patterns, "Normalization_Grid")

    # Grid 3: Unique and Summed PETs
    #generate_grid(subj, subj_dir, unique_patterns, "Unique_Summed_Grid")
