import SimpleITK as sitk
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def intensity_cerebellum_normalization(df_intensity):
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
    cerebellum_R = float((df_intensity.loc[(df_intensity['structure'] == 'cerebellum')
                                        & (df_intensity['hemisphere'] == 'R')]['mean_PET']).iloc[0])
    cerebellum_L = float((df_intensity.loc[(df_intensity['structure'] == 'cerebellum')
                                        & (df_intensity['hemisphere'] == 'L')]['mean_PET']).iloc[0])
    cerebellum = (cerebellum_R + cerebellum_L) / 2



    # normalize with cerebellum
    normalization_values = []

    for index, row in df_intensity.iterrows():
        normalization_values.append(row["mean_PET"] / cerebellum)

    df_intensity["normalization"] = normalization_values

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


def PET_FDG_quantification(path_brain_image, path_brain_segmentation,
                          name_labels, atlas="Hammers", path_second_segmentation=False, normalization=None):
    """ 
    Generates a dataframe with signal intensity per area and 
    a synthetic image with the intensity value per area.
    
    Parameters
    ----------
    path_brain_image : str 
        A path of brain image of a nifti image 
    path_brain_segmentation : str
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
    brain_image = sitk.ReadImage(path_brain_image)
    brain_segmentation = sitk.ReadImage(path_brain_segmentation)

    # Resample brain image to segmentation
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
    if path_second_segmentation:
        brain_segmentation2 = sitk.ReadImage(path_second_segmentation)
        brain_segmentation2 = resample_sitk(brain_segmentation2, brain_segmentation)
        array_segmentation2 = sitk.GetArrayFromImage(brain_segmentation2)

    for i, label in enumerate(labels):

        # match label with name of structure
        label_row = name_labels.loc[name_labels['n_label'] == label]  # find row with label in csv
        structure_name = str(label_row['structure'].values)[2:-2]  # match with the structure name

        if atlas == "Hammers":
            # split between two hemisphere
            structure_name = structure_name.split("-")

            if structure_name[-1] in ["L", "R"]:
                hemisphere = structure_name[-1]
                structure_name = "-".join(structure_name[:-1])
            else:
                hemisphere = ""
                structure_name = "-".join(structure_name)

        mask_label = array_segmentation == label  # create a mask
        

        # right cerebellum: 17 is Hammers right cerebellum, 47 is aseg right cerebellum
        if path_second_segmentation and atlas == "Hammers" and label == 17:
            mask_label1 = mask_label
            mask_label = array_segmentation2 == 47


        # left cerebellum: 18 is Hammers left cerebellum , 8 is aseg left cerebellum
        if path_second_segmentation and atlas == "Hammers" and label == 18:         
            mask_label = array_segmentation2 == 8      

        # signal intensity
        mean = np.mean(array_brain[mask_label])

        # copy signal intensty in new image 
        new_image[mask_label] = mean

        # append new row to df
        row = pd.Series(
            {'n_label': int(label), 'mean_PET': mean, 'structure': structure_name, 'hemisphere': hemisphere})
        new_df = pd.concat([new_df, row.to_frame().T], ignore_index=True)


    # Normalization
    if normalization:
        new_df = intensity_cerebellum_normalization(new_df) 

    new_image = sitk.GetImageFromArray(new_image)
    new_image.CopyInformation(brain_segmentation)
    return new_df, new_image

def image_change(df_subject, df_atlas, path_segmented_brain, path_second_segmentation=None):
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
    if path_second_segmentation:
        brain_segmentation2 = sitk.ReadImage(path_second_segmentation)
        array_segmentation2 = sitk.GetArrayFromImage(path_second_segmentation)

    segmented_brain = sitk.ReadImage(path_segmented_brain)

    # Get numpy array:
    array_segmented_brain = sitk.GetArrayFromImage(segmented_brain)

    # Create new image
    new_image = np.zeros_like(array_segmented_brain)

    # labels
    labels = df_atlas["n_label"]

    for label in labels:

        intensity_subject = float((df_subject.loc[df_subject["n_label"] == label]["normalization"]).iloc[0])
        intensity_atlas = float((df_atlas.loc[df_atlas["n_label"] == label]["normalization"]).iloc[0])

        # structure
        structure = str(df_subject.loc[df_subject["n_label"] == label]["structure"].values)[2:-2]

        # hemispheres
        hemisphere = str(df_subject.loc[df_subject["n_label"] == label]["hemisphere"].values)[2:-2]

        change = ((intensity_subject - intensity_atlas) / intensity_atlas) * 100

        if label == 0:
            change = 0

        # change mask for left and right cerebellum
        # right
        if path_second_segmentation and label == 17:
            mask_right_cerebellum = array_segmentation2 == 47

            # create new image
            new_image[mask_right_cerebellum] = change
            continue

        # left
        if path_second_segmentation and label == 18:
            # mask
            mask_left_cerebellum = array_segmentation2 == 8

            # create new image
            new_image[mask_left_cerebellum] = change
            continue

        # create mask
        mask_label = array_segmented_brain == label

        new_image[mask_label] = change

        # append new row to df
        row = pd.Series({'n_label': int(label), 'cambio': change, 'structure': structure, 'hemisphere': hemisphere})
        new_df = pd.concat([new_df, row.to_frame().T], ignore_index=True)

    new_image = sitk.GetImageFromArray(new_image)
    new_image.CopyInformation(segmented_brain)

    return new_image, new_df
