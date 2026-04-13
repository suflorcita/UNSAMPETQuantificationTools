import os
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt

# =====================================================
# CONFIG (Using necessary parts from the original script)
# =====================================================
output_plot_folder = "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ProcesadosCOVID_PET_MRI_09_2025/QC"
os.makedirs(output_plot_folder, exist_ok=True)

cmap = "plasma"
vmin_default, vmax_default = 0, 10  # Default from original script
# Use the Vmin, Vmax found in the original generate_grid
VMIN_CUSTOM, VMAX_CUSTOM = 0, 6.5

# --- Custom Image Paths ---
custom_images = {
    # The label will be used as the row title
    "CP0024 - PET Smoothed ANTWarped": "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ProcesadosCOVID_PET_10_25/CP0024/PET_smoothed_image_to_template_ANTWarped.nii.gz",
    "Fleni/NC5 - PET Raw ANT": "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/Fleni/NC5/PET_smoothed_image_to_template_ANT1Warp.nii.gz",
    "s93895370F... - PET Smoothed ANTWarp": "/media/sol/Expansion/Sol/Procesamiento-Imagenes-PET-10-2025/ProcesadosControlesCOVID_PET_10_25/s93895370F-5462-00001-000001/PET_smoothed_image_to_template_ANT1Warp.nii.gz",
}


# =====================================================
# FUNCTIONS (Adapted from original script)
# =====================================================
def auto_contrast(slice_array):
    """Compute percentile-based contrast limits."""
    data = slice_array[slice_array > 0] # Solo considera píxeles PET (mayor a 0)
    if data.size > 0:
        # ⭐️ Calcula el 1er y 99o percentil. Esto elimina outliers (valores extremos).
        return np.percentile(data, [1, 99])
    else:
        return vmin_default, vmax_default


def get_views(image, frame_index=0):
    """
    Return correctly oriented axial, coronal and sagittal slices.
    If 4D, show only one frame (default: first).
    """
    arr = sitk.GetArrayFromImage(image)

    # Handle 4D PET
    if arr.ndim == 4:
        n_frames = arr.shape[3]
        if frame_index < 0 or frame_index >= n_frames:
            frame_index = 0
        print(f"🧩 Using frame {frame_index + 1}/{n_frames} (shape={arr.shape})")
        arr = arr[frame_index, :, :, :]

    # Handle 3D PET (already done if 4D was processed)
    if arr.ndim == 3:
        # Middle slices
        z, y, x = arr.shape[0] // 2, arr.shape[1] // 2, arr.shape[2] // 2

        # Orientation corrections
        axial = np.flipud(arr[z, :, :])
        coronal = np.flipud(np.fliplr(arr[:, y, :]))
        sagittal = np.flipud(arr[:, :, x])

        return [axial, coronal, sagittal]

    # Handle unexpected dimensions
    print(f"⚠️ Image has unexpected dimensions: {arr.ndim}")
    return [np.zeros((1, 1)) for _ in range(3)]  # Return empty views


def generate_custom_grid(image_dict, grid_name="Custom_Grid"):
    """
    Generate a 3x3 grid of PET views for the specified image files.
    """

    n_rows = len(image_dict)
    fig, axes = plt.subplots(n_rows, 3, figsize=(9, 3 * n_rows))
    fig.patch.set_facecolor("black")

    # Ensure axes is 2D even for a single image case, though here it's 3
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    im = None  # To store the last image for the colorbar

    for row, (label, file_path) in enumerate(image_dict.items()):
        if not os.path.exists(file_path):
            for col in range(3):
                axes[row, col].axis("off")
                axes[row, col].set_title("Missing", color="gray")
            print(f"❌ Missing {file_path}")
            continue

        try:
            print(f"✅ Loading: {label}")
            img = sitk.ReadImage(file_path)
            # Use frame 0 if 4D, otherwise 3D processing is default
            slices = get_views(img, frame_index=0)

            for col, slc in enumerate(slices):

                # Use fixed vmin/vmax for consistent comparison
                vmin, vmax = auto_contrast(slc)
                im = axes[row, col].imshow(
                    # Mask zeros for black background
                    np.ma.masked_where(slc == 0, slc),
                    cmap=cmap, vmin=vmin, vmax=vmax
                )
                axes[row, col].axis("off")

                view_name = ["Axial", "Coronal", "Sagittal"][col]
                # Set row label on the first column title
                if col == 0:
                    title_text = f"{label}\n({view_name})"
                else:
                    title_text = view_name

                axes[row, col].set_title(title_text, color="white", fontsize=8)

        except Exception as e:
            print(f"⚠️ Error loading {file_path}: {e}")
            for col in range(3):
                axes[row, col].axis("off")
                axes[row, col].set_title("Error", color="red")

    # Colorbar (Must run only if at least one image was successfully loaded)
    if im:
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.ax.tick_params(colors="white")
        cbar.set_label("PET Intensity", color="white")
        cbar.ax.yaxis.set_tick_params(labelsize=8)  # Smaller text for white color

    # Title and save
    fig.suptitle(f"Custom PET Image Comparison Grid", color="white", fontsize=12)
    plt.subplots_adjust(left=0.05, right=0.9, top=0.93, bottom=0.05, wspace=0.05, hspace=0.2)

    output_path = os.path.join(output_plot_folder, f"{grid_name}.png")
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print(f"\n✅ Saved Custom Grid: {output_path}")


# =====================================================
# EXECUTION
# =====================================================
# Call the new function with your specific image paths
generate_custom_grid(custom_images, "3x3_Custom_PET_Grid")