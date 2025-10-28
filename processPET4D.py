import os

import SimpleITK as sitk
import matplotlib.pyplot as plt

import ImageRegistration as reg


def copy_information_4d_3d(image4d, image3d):
    spacing_4d = image4d.GetSpacing()
    origin_4d = image4d.GetOrigin()
    direction_4d = image4d.GetDirection()
    direction_3d = [direction_4d[i] for i in (0, 1, 2, 4, 5, 6, 8, 9, 10)]
    image3d.SetSpacing(spacing_4d[:])
    image3d.SetOrigin(origin_4d)
    image3d.SetDirection(direction_3d)

    return image3d


def image4D_to_3D(image_4D):
    images_3D = []
    # Array from image. Four dimensions: (x,y,z,t)
    image_4D_array = sitk.GetArrayFromImage(image_4D)

    for i in range(image_4D_array.shape[0]):
        image_3d = sitk.GetImageFromArray(image_4D_array[i, :, :, :])
        image_3d = sitk.Cast(image_3d, sitk.sitkFloat32)  # cast
        image_3d = copy_information_4d_3d(image_4D, image_3d)
        images_3D.append(image_3d)

    return images_3D


def register_and_sum_pet_images(PET_images):
    """
    Register all PET frames rigidly to the first one
    and return the aligned frames plus their summed 3D image.
    """

    # Use the first PET frame as reference
    reference_PET_image = PET_images[0]

    # Rigidly register all other frames to the reference
    for i in range(1, len(PET_images)):
        resultReg = reg.RigidImageRegistration(PET_images[i], reference_PET_image, printLog=True)
        PET_images[i] = resultReg['image']

    # Initialize an empty image with the same metadata as the reference
    sum_pet_3d_image = sitk.Image(reference_PET_image.GetSize(), reference_PET_image.GetPixelID())
    sum_pet_3d_image.CopyInformation(reference_PET_image)

    # Add all registered frames to obtain one single PET image
    for image in PET_images:
        sum_pet_3d_image = sitk.Add(sum_pet_3d_image, image)

    return sum_pet_3d_image, PET_images


def register_PET_MRI(PET_image, MRI_image):
    """
    Main function to register PET to MRI.
    - Loads PET (directory, 4D or 3D)
    - Registers PET frames to each other and sums them
    - Registers the summed PET to the MRI
    - Returns: individual frames, summed PET, registered PET, and transformation
    """

    # --- Load PET ---
    if os.path.isdir(PET_image):
        PET_images = []
        for image in os.listdir(PET_image):
            path_image = os.path.join(PET_image, image)
            img = sitk.ReadImage(path_image)
            img = sitk.Cast(img, sitk.sitkFloat32)
            PET_images.append(img)
    else:
        image = sitk.ReadImage(PET_image)
        if len(image.GetSize()) == 4:  # PET is 4D
            PET_images = image4D_to_3D(image)
        else:  # PET is already 3D
            PET_images = [image]

    # --- Load MRI ---
    MRI_image = sitk.ReadImage(MRI_image)
    MRI_image = sitk.Cast(MRI_image, sitk.sitkFloat32)

    # --- Register and sum PET frames ---
    sum_pet_3d_image, PET_images = register_and_sum_pet_images(PET_images)

    # --- Register summed PET to MRI ---
    result_registration = reg.RigidImageRegistration(sum_pet_3d_image, MRI_image, printLog=True)
    register_pet_t1 = result_registration["image"]
    txPET2MRI = result_registration["tx"]

    return PET_images, sum_pet_3d_image, register_pet_t1, txPET2MRI

def read_and_apply_tx(tx_file, image, ref_image = []):

    # Read the transform file
    with open(tx_file, 'r') as f:
        lines = f.readlines()

    # Initialize variables to store the parameters
    parameters = None
    fixed_parameters = None

    # Loop through the lines and find the line starting with "Parameters"
    for line in lines:
        if line.startswith("Parameters:"):
            parameters = list(map(float, line.strip().split(": ")[1].split()))
        elif line.startswith("FixedParameters:"):
            fixed_parameters = list(map(float, line.strip().split(": ")[1].split()))

    # Create the tx
    transform = sitk.Euler3DTransform()
    transform.SetParameters(parameters)
    transform.SetFixedParameters(fixed_parameters)

    output_image = reg.ApplyRegTransform(image, transform, refImage=ref_image, interpolator=sitk.sitkLinear)
    return output_image


if __name__ == '__main__':
    path_PET_image = "/home/sol/PET_MRI/CODIGOS/Images/CEUNIM/PET"
    path_MRI_image = "/home/sol/PET_MRI/CODIGOS/Images/CEUNIM/T1.nii.gz"

    sum_pet, register_image, tx_PET_2_RMN = register_PET_MRI(path_PET_image, path_MRI_image)


