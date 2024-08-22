import os.path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
import SimpleITK as sitk
import sys
from scipy import ndimage as ndi
from nilearn import plotting, datasets, surface


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

def plot_hypometabolism_maps(segmented_image, brain_image, output_path, coord=(-3, -40, 15)):
    # Check if spacing, size, or orientation differ
    if (segmented_image.GetSpacing() != brain_image.GetSpacing() or
            segmented_image.GetSize() != brain_image.GetSize() or
            segmented_image.GetDirection() != brain_image.GetDirection()):
                brain_image = resample_sitk(brain_image, segmented_image)

    ax = plt.subplot()
    plotting.plot_stat_map(segmented_image, bg_img=brain_image,
                           annotate=False,
                           cmap="bwr", black_bg=False, cut_coords=coord, display_mode="ortho", axes=ax,
                           threshold=10, vmax=50)

    ax.remove()

    plt.savefig(output_path)
    return


def plot_bar_charts(segmented_image, brain_image, output_path):

    return
