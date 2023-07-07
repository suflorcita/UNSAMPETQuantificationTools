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


def intensity_normalization(image, mode=("avg total value", 0)):
    image_array = sitk.GetArrayFromImage(image)

    if mode[0] == "avg total value":
        avg_total_value = np.sum(image_array) / image_array.size
        scalar = avg_total_value

    elif mode[0] == "cerebellum":
        scalar = mode[1]

    normalized_image = image_array / scalar

    # show image 
    # for i in range(0, image_array.shape[0], 50): 
    #     show_image(image_array, normalized_image, i)

    normalized_image = sitk.GetImageFromArray(normalized_image)
    normalized_image.CopyInformation(image)
    return normalized_image


