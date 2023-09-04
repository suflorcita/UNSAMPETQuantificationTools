import os

import SimpleITK as sitk
import matplotlib.pyplot as plt
import ImageRegistration as reg


def copy_information_4d_3d(image4d, image3d):
    spacing_4d = image4d.GetSpacing()
    origin_4d = image4d.GetOrigin()

    image3d.SetSpacing(spacing_4d[:])
    image3d.SetOrigin(origin_4d)

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


def register_PET_MRI(PET_image, MRI_image):

    # If 4D Convert PET 4D image in multiple 3D image
    if os.path.isdir(PET_image):
        PET_images = []
        images = os.listdir(PET_image)

        for image in images:
            path_image = PET_image + "/" + image
            image_pet = sitk.ReadImage(path_image)
            image_pet = sitk.Cast(image_pet, sitk.sitkFloat32)  # cast
            PET_images.append(image_pet)
    else:
        image = sitk.ReadImage(PET_image)
        if len(image.GetSize()) == 4:
            PET_images = image4D_to_3D(image)
        else:
            PET_images = [image]


    # Read MRI image
    MRI_image = sitk.ReadImage(MRI_image)
    MRI_image = sitk.Cast(MRI_image, sitk.sitkFloat32)

    # # Register PET images between them and sum them:
    # Registration between frames
    reference_PET_image = PET_images[0]
    for i in range(1, len(PET_images)):
        resultReg = reg.RigidImageRegistration(PET_images[i], reference_PET_image, printLog=True)
        PET_images[i] = resultReg['image']

    # Compute the sum
    sum_pet_3d_image  = PET_images[0]
    for image in PET_images[1:]:
        sum_pet_3d_image = sitk.Add(sum_pet_3d_image , image)

    # Register to t1 using the sum
    result_registration = reg.RigidImageRegistration(sum_pet_3d_image, MRI_image, printLog=True)

    register_pet_t1 = result_registration["image"]
    txPET2MRI = result_registration["tx"]


    return PET_images, sum_pet_3d_image, register_pet_t1, txPET2MRI


if __name__ == '__main__':
    path_PET_image = "/home/sol/PET_MRI/CODIGOS/Images/CEUNIM/PET"
    path_MRI_image = "/home/sol/PET_MRI/CODIGOS/Images/CEUNIM/T1.nii.gz"

    sum_pet, register_image, tx_PET_2_RMN = register_PET_MRI(path_PET_image, path_MRI_image)


