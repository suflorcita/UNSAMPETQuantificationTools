import SimpleITK as sitk
import numpy as np
import pandas as pd
import os
import PETquantification as quant
import matplotlib.pyplot as plt

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

path_Hammers = os.path.join(script_dir, 'data', 'atlas', "Hammers_mith-n30r95-MaxProbMap-gm-MNI152-SPM12.nii.gz")
path_MNI_152_T1_brain = os.path.join(script_dir, 'data', 'atlas', 'MNI152_T1_1mm_brain.nii.gz')

labels_Hammers_csv_path =  os.path.join(script_dir, "Labels/labels_Hammers.csv")
df_labels_Hammers = pd.read_csv(labels_Hammers_csv_path)

# Read T1_brain image.
MNI_152_T1_brain_sitk = sitk.ReadImage(path_MNI_152_T1_brain)

# Read Hammers image.
Hammers_image = sitk.ReadImage(path_Hammers)
def resample_sitk(input_image, reference_image):
    """
    Resamples a SimpleITK image to match the spacing,
    size, and orientation of a reference image.

    Parameters
    ----------
    input_image  : SimpleITK.Image
        The input image to be resampled.
    reference_image : SimpleITK.Image
        The reference image whose spacing, size,
        and orientation will be used for resampling.

    Returns
    ---------
    resampled_image: simpleITK image object

    """
    # Create an identity transform
    identity_transform = sitk.Transform()

    # Define the interpolator to use
    interpolator = sitk.sitkLinear

    # Define the default pixel value for points outside the input image
    defaultPixelValue = 0

    # Define the output pixel type
    outputPixelType = sitk.sitkFloat32

    # Resample the input image
    resampled_image = sitk.Resample(input_image, reference_image, identity_transform, interpolator, defaultPixelValue,
                                    outputPixelType)
    return resampled_image


def show_normalized_image(image1, image2, slice):

    # Select the slice indices along the third dimension
    slice_index1 = slice
    slice_index2 = slice

    # Extract the 2D slices from the 3D array
    slice_data1 = image1[:, :, slice_index1]
    slice_data2 = image2[:, :, slice_index2]

    # Create a figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Show the first slice in the first subplot
    axes[1].imshow(slice_data1, cmap='gray')
    axes[1].axis('off')

    # Show the second slice in the second subplot
    axes[0].imshow(slice_data2, cmap='gray')
    axes[0].axis('off')

    # Adjust the spacing between subplots
    plt.tight_layout()

    # Show the figure
    plt.show()


def show_pet_image(mri_image, pet_image, slice):

    # Select the slice indices along the third dimension
    slice_index1 = slice
    slice_index2 = slice

    # Extract the 2D slices from the 3D array
    slice_data1 = mri_image[:, :, slice_index1]
    slice_data2 = pet_image[:, :, slice_index2]

    # Rotate the image 180 degrees
    rotated_slice_data1 = np.rot90(slice_data1, 2)
    rotated_slice_data2 = np.rot90(slice_data2, 2)

    # Create a figure with two subplots
    fig, ax = plt.subplots(figsize=(10, 5))

    # Show the first slice in the first subplot
    ax.imshow(rotated_slice_data1, cmap='gray')
    ax.axis('off')

    # Show the second slice in the second subplot
    ax.imshow(rotated_slice_data2, cmap='hot', alpha=0.5)
    ax.axis('off')

    # Adjust the spacing between subplots
    plt.tight_layout()

    # Show the figure
    plt.show()


def intensity_normalization(image, mode="avg total value", mask_only_brain = MNI_152_T1_brain_sitk,
                            atlas=None, atlas_image=None, df_labels=None, scalar=None):
    ''' image: Sitk object image '''
    original_image = image

    if mode  == "avg total value":

        # Resample mask only brain to input image
        if mask_only_brain.GetSize() != image.GetSize():
            mask_only_brain = resample_sitk(mask_only_brain, image)

        mask_image_array = sitk.GetArrayFromImage(mask_only_brain) # Convert mask to array
        image_only_brain = sitk.GetArrayFromImage(image)

        # for i in range(0, mask_image_array.shape[0], 10):
        #     show_pet_image(mask_image_array, image_only_brain, i)

        # Create Mask only brain
        image_only_brain[mask_image_array == 0] = 0

        # Compute avg value
        avg_total_value = np.sum(image_only_brain) / image_only_brain.size
        scalar = avg_total_value

    elif mode == "cerebellum":
        # Quantify using Hammers atlas
        df_atlas_intensity, image_atlas_intensity = quant.PET_FDG_quantification(image, atlas_image,
                                                                                 df_labels,
                                                                                 atlas=atlas)
        # Normalize using Cerebellum values
        if atlas == "Hammers":
            cerebellum_name = "cerebellum"
        else:
            cerebellum_name = 'Cerebellum-Cortex'

        # Extract cerebellum values
        cerebellum_R = float((df_atlas_intensity.loc[(df_atlas_intensity['structure'] == cerebellum_name)
                                                     & (df_atlas_intensity['hemisphere'] == 'R')]['mean_PET']).iloc[0])
        cerebellum_L = float((df_atlas_intensity.loc[(df_atlas_intensity['structure'] ==  cerebellum_name)
                                                     & (df_atlas_intensity['hemisphere'] == 'L')]['mean_PET']).iloc[0])
        cerebellum = (cerebellum_R + cerebellum_L) / 2

        print(cerebellum)


        scalar = cerebellum

    image_array = sitk.GetArrayFromImage(image)
    normalized_image = image_array / scalar

    # if mode == "avg total value":
    #     # show image
    #     for i in range(0, image_array.shape[0], 10):
    #       show_normalized_image(image_array, normalized_image, i)

    normalized_image = sitk.GetImageFromArray(normalized_image)
    normalized_image.CopyInformation(image)

    return normalized_image


def suv_image(image, weight, total_dose_bq):
    '''image: sitk image
    suv: csv'''
    image_array = sitk.GetArrayFromImage(image)

    image_suv = (image_array * weight)/total_dose_bq

    image_suv = sitk.GetImageFromArray(image_suv)
    image_suv.CopyInformation(image)
    return image_suv



