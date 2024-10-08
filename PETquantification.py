import SimpleITK as sitk
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

path_MNI_152_T1_brain = os.path.join(script_dir, 'data', 'atlas', 'MNI152_T1_1mm_brain.nii.gz')
label_Hammers_right_cerebellum = 17
label_Hammers_left_cerebellum = 18
aseg_right_cerebellum = 47
aseg_left_cerebellum = 8

# Read T1_brain image.
MNI_152_T1_brain_sitk = sitk.ReadImage(path_MNI_152_T1_brain)


def new_mask_without_ventricles(image_only_brain, aseg_segmentation):
    '''

    Parameters
    ----------
    image_only_brain = sitk image mask only brain (with ventricles)
    aseg_segmentation = sitk image freesurfer subcortical segmentation

    Returns
    -------

    '''

    left_lateral_vent_label = 4
    left_inf_lateral_vent_label = 5
    right_lateral_vent_label = 43
    right_inf_lateral_vent_label = 44
    third_vent_label = 43
    fourth_vent_label = 44

    new_brain_without_ventricles_array = sitk.GetArrayFromImage(image_only_brain)
    aseg_segmentation_array = sitk.GetArrayFromImage(aseg_segmentation)

    new_brain_without_ventricles_array[aseg_segmentation_array == left_lateral_vent_label] = 0
    new_brain_without_ventricles_array[aseg_segmentation_array == left_inf_lateral_vent_label] = 0
    new_brain_without_ventricles_array[aseg_segmentation_array == right_inf_lateral_vent_label] = 0
    new_brain_without_ventricles_array[aseg_segmentation_array == right_lateral_vent_label] = 0
    new_brain_without_ventricles_array[aseg_segmentation_array == fourth_vent_label] = 0
    new_brain_without_ventricles_array[aseg_segmentation_array == third_vent_label] = 0

    new_brain_without_ventricles = sitk.GetImageFromArray(new_brain_without_ventricles_array)
    new_brain_without_ventricles.CopyInformation(image_only_brain)

    return new_brain_without_ventricles


def intensity_cerebellum_normalization(df_intensity, atlas="Hammers"):
    """
    Add a column to the dataframe
    with the normalized intensity value per area, using the cerebellum as the reference.


    Parameters
    ----------
    df_intensity: pandas.DataFrame

    Returns
    ---------
    df_intensity: pandas.DataFrame

    """

    # cerebellum values
    if atlas == "Hammers":
        cerebellum_R = float((df_intensity.loc[(df_intensity['structure'] == 'cerebellum')
                                               & (df_intensity['hemisphere'] == 'R')]['mean_PET']).iloc[0])
        cerebellum_L = float((df_intensity.loc[(df_intensity['structure'] == 'cerebellum')
                                               & (df_intensity['hemisphere'] == 'L')]['mean_PET']).iloc[0])
        cerebellum = (cerebellum_R + cerebellum_L) / 2
    elif atlas == "DKT":
        cerebellum_R = float((df_intensity.loc[(df_intensity['structure'] == 'Cerebellum-Cortex')
                                               & (df_intensity['hemisphere'] == 'R')]['mean_PET']).iloc[0])
        cerebellum_L = float((df_intensity.loc[(df_intensity['structure'] == 'Cerebellum-Cortex')
                                               & (df_intensity['hemisphere'] == 'L')]['mean_PET']).iloc[0])
        cerebellum = (cerebellum_R + cerebellum_L) / 2
    # normalize with cerebellum
    normalization_values = []

    for index, row in df_intensity.iterrows():
        normalization_values.append(row["mean_PET"] / cerebellum)

    df_intensity["normalization_cerebellum"] = normalization_values

    return df_intensity


def intensity_mean_total(df_intensity, brain_image, path_mask_only_brain):
    """
    Add a column to the dataframe
    with the normalized intensity value per area, using the mean as the reference.


    Parameters
    ----------
    df_intensity: pandas.DataFrame
    brain_image: Array of a brain image

    Returns
    ---------
    df_intensity: pandas.DataFrame

    """
    # Read MNI152 brain to mask
    mask_only_brain = sitk.ReadImage(path_mask_only_brain)
    mask_only_brain_array = sitk.GetArrayFromImage(mask_only_brain)

    brain_image[mask_only_brain_array == 0] = 0

    # avg_total_value
    mean_total_value = np.sum(brain_image) / brain_image.size
    scalar = mean_total_value

    # normalize with mean total value
    normalization_values = []

    for index, row in df_intensity.iterrows():
        normalization_values.append(row["mean_PET"] / scalar)

    df_intensity["normalization_mean"] = normalization_values

    return df_intensity


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


def PET_FDG_quantification(brain_image, brain_segmentation,
                           name_labels, atlas="Hammers", second_segmentation=False, normalization=None,
                           normalization_mean=False, mask_only_brain=MNI_152_T1_brain_sitk):
    """
    Generates a dataframe with signal intensity per area and
    a synthetic image with the intensity value per area.

    Parameters
    ----------
    brain_image : sitk object
        A path of brain image of a nifti image
    brain_segmentation : sitk object
        A path of brain segmentation of a nifti image
    name_labels: pandas.Dataframe
        A dataframe with the name of labels [format: n_label -- structure]
        Index:
            RangeIndex
        Columns:
            Name: n_label, dtype: int64
            Name: structure, dtype: str
            Name:
    atlas : str, optional
        Atlas used for the segmentation (default is "Hammers").
    path_second_segmentation : str, optional
        Path to a FreeSurfer aseg segmentation (default is None).
    normalization : boolean, optional
        Specifies whether to perform normalization with cerebellum values.

    Returns
    ---------
    new_df: pandas.dataframe
        Columns:
            Name: , dytpe:
            Name: , dytpe:
    new_image: SimpleITK image object
        Signal intensty values in each voxel.

    """

    # Resample brain image to segmentation
    if brain_image.GetSize() != brain_segmentation.GetSize():
        brain_image_not_resample = brain_image
        brain_image = resample_sitk(brain_image, brain_segmentation)

    # Get numpy array:
    array_segmentation = sitk.GetArrayFromImage(brain_segmentation)
    array_brain = sitk.GetArrayFromImage(brain_image)

    # Create new image
    new_image = np.zeros_like(array_segmentation)

    # labels
    labels = np.unique(array_segmentation)

    # create df
    new_df = pd.DataFrame()

    # second image cerebellum
    if second_segmentation:
        brain_segmentation2 = resample_sitk(second_segmentation, brain_segmentation)
        array_segmentation2 = sitk.GetArrayFromImage(brain_segmentation2)

    for i, label in enumerate(labels):

        # match label with name of structure
        label_row = name_labels.loc[name_labels['n_label'] == label]  # find row with label in csv
        structure_name = str(label_row['structure'].values)[2:-2]  # match with the structure name

        if atlas == "Hammers":
            # split structure name
            structure_name = structure_name.split("-")

            if structure_name[-1] in ["L", "R"]:
                hemisphere = structure_name[-1]
                structure_name = "-".join(structure_name[:-1])
            else:
                hemisphere = ""
                structure_name = "-".join(structure_name)

        elif atlas == "DKT":
            # split structure name
            structure_name = structure_name.split("-")

            if structure_name[0] in ["Left", "Right"]:
                if structure_name[0] == "Left":
                    hemisphere = "L"
                else:
                    hemisphere = "R"
                structure_name = "-".join(structure_name[1:])

            elif len(structure_name) > 1:
                if structure_name[1] in ["lh", "rh"]:
                    if structure_name[1] == "lh":
                        hemisphere = "L"
                        structure_name = "-".join(structure_name[2:])
                    else:
                        hemisphere = "R"
                        structure_name = "-".join(structure_name[2:])
                else:
                    structure_name = "-".join(structure_name[1:])

            else:
                structure_name = "-".join(structure_name)
                hemisphere = ""

        mask_label = array_segmentation == label  # create a mask

        # right cerebellum: 17 is Hammers right cerebellum, 47 is aseg right cerebellum
        if second_segmentation and atlas == "Hammers" and label == label_Hammers_right_cerebellum:
            mask_label = array_segmentation2 == aseg_right_cerebellum

        # left cerebellum: 18 is Hammers left cerebellum , 8 is aseg left cerebellum
        if second_segmentation and atlas == "Hammers" and label == label_Hammers_left_cerebellum:
            mask_label = array_segmentation2 == aseg_left_cerebellum

        # signal intensity
        mean = np.mean(array_brain[mask_label])

        # copy signal intensty in new image
        new_image[mask_label] = mean

        # append new row to df
        row = pd.Series(
            {'n_label': int(label), 'mean_PET': mean, 'structure': structure_name, 'hemisphere': hemisphere})
        new_df = pd.concat([new_df, row.to_frame().T], ignore_index=True)

    # Normalization cerebellum
    if normalization:
        new_df = intensity_cerebellum_normalization(new_df, atlas)

    # Normalization mean total value
    if normalization_mean:
        if atlas == "Hammers":
            array_brain_not_resample = sitk.GetArrayFromImage(brain_image_not_resample)
            new_df = intensity_mean_total(new_df, array_brain_not_resample, mask_only_brain)
        else:
            new_df = intensity_mean_total(new_df, array_brain, mask_only_brain)

    new_image = sitk.GetImageFromArray(new_image)
    new_image.CopyInformation(brain_segmentation)
    return new_df, new_image


def image_change(df_subject, df_atlas, segmented_brain, mode="Mean"):
    """
    Generates a dataframe with the percentage change per area, comparing the image with an atlas,
    and a synthetic image with the corresponding percentage change.

    Parameters
    ----------
    df_subject : pandas.Dataframe
    df_atlas : pandas.Dataframe
    segmented_brain: simple
    cerebellum:



    new_df: pandas.dataframe
        Columns:
            Name: , dytpe:
            Name: , dytpe:
    new_image: SimpleITK image object
        Signal intensty values in each voxel.

    """

    # create df
    new_df = pd.DataFrame()

    # Get numpy array:
    array_segmented_brain = sitk.GetArrayFromImage(segmented_brain)

    # Create new image
    new_image = np.zeros_like(array_segmented_brain)

    # labels
    labels = df_atlas["n_label"]

    for label in labels:

        if mode == "cerebellum":
            name_column = "Normalization to cerebellum uptake values"
        elif mode == "mean":
            name_column = "Normalization to total brain mean value"

        intensity_subject = float((df_subject.loc[df_subject["n_label"] == label][name_column]).iloc[0])
        intensity_atlas = float((df_atlas.loc[df_atlas["n_label"] == label][name_column]).iloc[0])

        # structure
        structure = str(df_subject.loc[df_subject["n_label"] == label]["Structure"].values)[2:-2]

        # hemispheres
        # hemisphere = str(df_subject.loc[df_subject["n_label"] == label]["hemisphere"].values)[2:-2]

        change = ((intensity_subject - intensity_atlas) / intensity_atlas) * 100

        if label == 0:
            change = 0

        # create mask
        mask_label = array_segmented_brain == label

        new_image[mask_label] = change

        # append new row to df
        row = pd.Series({'n_label': int(label), 'change': change, 'structure': structure})

        # , 'hemisphere': hemisphere
        new_df = pd.concat([new_df, row.to_frame().T], ignore_index=True)

    new_image = sitk.GetImageFromArray(new_image)
    new_image.CopyInformation(segmented_brain)

    return new_image, new_df


def image_change_2(df_subject, df_atlas, segmented_brain, second_segmentation=None):
    """
    Generates a dataframe with the percentage change per area, comparing the image with an atlas,
    and a synthetic image with the corresponding percentage change.

    Parameters
    ----------
    df_subject : pandas.Dataframe
    df_atlas : pandas.Dataframe
    segmented_brain: simple
    cerebellum:



    new_df: pandas.dataframe
        Columns:
            Name: , dytpe:
            Name: , dytpe:
    new_image: SimpleITK image object
        Signal intensty values in each voxel.

    """

    # create df
    new_df = pd.DataFrame()

    # second image cerebellum
    if second_segmentation:
        brain_segmentation2 = resample_sitk(second_segmentation, segmented_brain)
        array_segmentation2 = sitk.GetArrayFromImage(brain_segmentation2)

    # Get numpy array:
    array_segmented_brain = sitk.GetArrayFromImage(segmented_brain)

    # Create new image
    new_image = np.zeros_like(array_segmented_brain)

    # labels
    labels = df_atlas["n_label"]

    for label in labels:

        intensity_subject = float((df_subject.loc[df_subject["n_label"] == label]["mean_PET"]).iloc[0])
        intensity_atlas = float((df_atlas.loc[df_atlas["n_label"] == label]["mean_PET"]).iloc[0])

        # structure
        structure = str(df_subject.loc[df_subject["n_label"] == label]["structure"].values)[2:-2]

        # hemispheres
        hemisphere = str(df_subject.loc[df_subject["n_label"] == label]["hemisphere"].values)[2:-2]

        change = ((intensity_subject - intensity_atlas) / intensity_atlas) * 100

        if label == 0:
            change = 0

        # append new row to df
        row = pd.Series({'n_label': int(label), 'cambio': change, 'structure': structure, 'hemisphere': hemisphere})
        new_df = pd.concat([new_df, row.to_frame().T], ignore_index=True)

        # change mask for left and right cerebellum
        # right
        if second_segmentation and label == label_Hammers_right_cerebellum:
            mask_right_cerebellum = array_segmentation2 == aseg_right_cerebellum

            # create new image
            new_image[mask_right_cerebellum] = change
            continue

        # left
        if second_segmentation and label == label_Hammers_left_cerebellum:
            # mask
            mask_left_cerebellum = array_segmentation2 == aseg_left_cerebellum

            # create new image
            new_image[mask_left_cerebellum] = change
            continue

        # create mask
        mask_label = array_segmented_brain == label

        new_image[mask_label] = change

    new_image = sitk.GetImageFromArray(new_image)
    new_image.CopyInformation(segmented_brain)

    return new_image, new_df
