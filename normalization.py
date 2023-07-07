import SimpleITK as sitk
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def show_image(image1, image2, slice): 

    # Select the slice indices along the third dimension
    slice_index1 = slice
    slice_index2 = slice

    # Extract the 2D slices from the 3D array
    slice_data1 = image1[:, :, slice_index1]
    slice_data2 = image2[:, :, slice_index2]

    # Create a figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Show the first slice in the first subplot
    axes[0].imshow(slice_data1, cmap='gray')
    axes[0].axis('off')

    # Show the second slice in the second subplot
    axes[1].imshow(slice_data2, cmap='gray')
    axes[1].axis('off')

    # Adjust the spacing between subplots
    plt.tight_layout()

    # Show the figure
    plt.show()    


def intensity_normalization(image, mode="avg total value"): 
    image_array = sitk.GetArrayFromImage(image)

    if mode =="avg total value":
        avg_total_value = np.sum(image_array) / image_array.size
    #elif mode == "cerebellum":
    
    normalized_image = image_array / avg_total_value

    # show image 
    # for i in range(0, image_array.shape[0], 50): 
    #     show_image(image_array, normalized_image, i)

    normalized_image = sitk.GetImageFromArray(normalized_image)
    normalized_image.CopyInformation(image)
    return normalized_image




if __name__ == '__main__':
    path_brain_image  = "/home/solcat/PET/002_S_4270/ANTs/PET_Norm_MNI_152_ANT.nii.gz"
    brain_image = sitk.ReadImage(path_brain_image)
    norm_image = intensity_normalization(brain_image)


