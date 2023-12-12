import processPET4D as reg
import ImageRegistration as imreg
import SimpleITK as sitk
import subprocess
import pandas as pd
import PETquantification as quant
import normalization as norm
import sys
import os
import shutil, math

if __name__ == '__main__':

    # images
    if len(sys.argv) == 5:
        path_PET_image = sys.argv[1]
        path_MRI_image = sys.argv[2]
        subject = sys.argv[3]
        output_path = sys.argv[4]

    elif len(sys.argv) == 4:
        path_PET_image = sys.argv[1]
        path_MRI_image = None
        subject = sys.argv[2]
        output_path = sys.argv[3]
    else:
        path_PET_image = "./FDG-PET.nii.gz"
        path_MRI_image = "./T1_MPRAGE.nii.gz"
        subject = "anonymous"
        output_path = "./"

    # Path ATLAS
    path_MNI_152_T1 = "./data/atlas/MNI152_T1_1mm.nii.gz"
    path_MNI_152_T1_brain = "./data/atlas/MNI152_T1_1mm_brain.nii.gz"
    path_MNI_152_PET = "./data/atlas/MNI152_PET_1mm.nii"
    path_Hammers = "./data/atlas/Hammers_mith-n30r95-MaxProbMap-gm-MNI152-SPM12.nii.gz"
    path_DKT_MNI152 = "./data/atlas/"
    path_PET_template_CN = "ATLAS/AtlasPETFDG_FiltOutliers_CN.nii"


    # labels segmentation
    labels_FS_csv_path = "Labels/FS_labels.csv"
    labels_Hammers_csv_path = "Labels/labels_Hammers.csv"

    # output
    output_path = output_path + "/" + subject

    if not os.path.exists(output_path):
        os.mkdir(output_path)


    if not path_MRI_image:
        aseg = False

        # Copy Image to Processed dir
        PET_image_processed_dir = os.path.join(output_path, "Raw_PET_image.nii.gz")
        shutil.copy(path_PET_image, PET_image_processed_dir)

        PET_image = sitk.ReadImage(path_PET_image)

        if len(PET_image.GetSize()) == 4:

            # Creates FRAMES dir
            frames_dir = os.path.join(output_path, "FRAMES")

            if not os.path.exists(frames_dir):
                os.mkdir(frames_dir)

            PET_images = reg.image4D_to_3D(PET_image)



            for i, PET_image_frame in enumerate(PET_images):
                new_dir_PET_frame = os.path.join(frames_dir, f"FRAME_{i}")

                if not os.path.exists(new_dir_PET_frame):
                    os.mkdir(new_dir_PET_frame)

                # Apply Gaussian Filter to PET Image
                # Calculate the standard deviation for an FWHM of 6 mm
                fwhm = 6.0
                sigma = fwhm / (2 * math.sqrt(2 * math.log(2)))

                # Create Gaussian Filter
                gaussian_filter = sitk.SmoothingRecursiveGaussianImageFilter()

                # Set the filter parameters
                gaussian_filter.SetSigma(sigma)

                # Apply the Gaussian filter to the input image
                output_image = gaussian_filter.Execute(PET_image_frame)

                # Save the output image
                path_smoothed_image = os.path.join(new_dir_PET_frame, "Smoothed_PET_image.nii.gz")
                sitk.WriteImage(output_image, path_smoothed_image)
                PET_image_frame = path_smoothed_image

                # Normalize to Atlas PET

                # FLIRT
                path_PET_Template_Flirt = os.path.join(new_dir_PET_frame, "PETimage_template_flirt.nii.gz")
                path_FLIRT_tx = os.path.join(new_dir_PET_frame, "PETimage_template_flirt.mat")
                flirt_command = f"flirt -in {PET_image_frame} -ref {path_PET_template_CN} -out {path_PET_Template_Flirt} -omat {path_FLIRT_tx} -interp trilinear -dof 12"
                subprocess.run([flirt_command], shell=True)


                # ANT

                # Use 16 threads
                subprocess.run(["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=16"], shell=True)


                ANT_path = "/home/sol/PET_MRI/CODIGOS/antsRegistrationSyN.sh"
                path_ANT_image = os.path.join(new_dir_PET_frame, "PETimage_template_ANT")
                ANT_command = f"{ANT_path} -d 3 -f {path_PET_template_CN} -m {path_PET_Template_Flirt} -o {path_ANT_image}"
                subprocess.run([ANT_command], shell=True)


                print("ANTs:OK")

                path_PET_final = path_ANT_image + "Warped.nii.gz"

                print(path_PET_final)

        else:
            # Apply Gaussian Filter to PET Image
            # Calculate the standard deviation for an FWHM of 6 mm
            fwhm = 6.0
            sigma = fwhm / (2 * math.sqrt(2 * math.log(2)))

            # Create Gaussian Filter
            gaussian_filter = sitk.SmoothingRecursiveGaussianImageFilter()

            # Set the filter parameters
            gaussian_filter.SetSigma(sigma)

            # Apply the Gaussian filter to the input image
            output_image = gaussian_filter.Execute(PET_image)

            # Save the output image
            path_smoothed_image = os.path.join(output_path, "Smoothed_PET_image.nii.gz")
            sitk.WriteImage(output_image, path_smoothed_image)
            PET_image_frame = path_smoothed_image

            # Normalize to Atlas PET

            # FLIRT
            path_PET_Template_Flirt = os.path.join(output_path, "PETimage_template_flirt.nii.gz")
            path_FLIRT_tx = os.path.join(output_path, "PETimage_template_flirt.mat")
            flirt_command = f"flirt -in {PET_image_frame} -ref {path_PET_template_CN} -out {path_PET_Template_Flirt} -omat {path_FLIRT_tx} -interp trilinear -dof 12"
            subprocess.run([flirt_command], shell=True)

            # ANT

            # Use 16 threads
            subprocess.run(["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=16"], shell=True)

            ANT_path = "/home/sol/PET_MRI/CODIGOS/antsRegistrationSyN.sh"
            path_ANT_image = os.path.join(output_path, "PETimage_template_ANT")
            ANT_command = f"{ANT_path} -d 3 -f {path_PET_template_CN} -m {path_PET_Template_Flirt} -o {path_ANT_image}"
            subprocess.run([ANT_command], shell=True)

            print("ANTs:OK")

            path_PET_final = path_ANT_image + "Warped.nii.gz"

    else:
        freesurfer = True
        if freesurfer:
            subject_dir = output_path + "/Freesurfer"
            if not os.path.exists(subject_dir):
                os.mkdir(subject_dir)

            recon_all_command = f"recon-all -s {subject} -i {path_MRI_image} -sd {subject_dir} -all -parallel -openmp 16"

            subprocess.run([recon_all_command], shell=True)

        path_T1_FS_mgz = subject_dir + "/" + subject + "/mri/T1.mgz"

        if os.path.exists(path_T1_FS_mgz):
            print("FS: ok. T1 is FS T1.mgz")
            path_T1 = output_path + "/T1_FS.nii"
            mri_convert_t1_FS = f"mri_convert {path_T1_FS_mgz} {path_T1}"
            subprocess.run([mri_convert_t1_FS], shell=True)

            path_T1_brain_mgz = subject_dir + "/" + subject + "/mri/brainmask.mgz"
            path_T1_brain = output_path + "/T1_FS_onlybrain.nii"
            mri_convert_t1_only_brain = f"mri_convert {path_T1_brain_mgz} {path_T1_brain}"
            subprocess.run([mri_convert_t1_only_brain], shell=True)

            # Segmentation
            aseg = True
            path_segmentation_mgz = subject_dir + "/" + subject + "/mri/aseg.mgz"
            path_segmentation = output_path + "/aseg.nii"
            mri_convert_segmentation = f"mri_convert {path_segmentation_mgz} {path_segmentation}"
            subprocess.run([mri_convert_segmentation], shell=True)

            # segmentation2
            path_aparc_dkt_mgz = subject_dir + "/" + subject + "/mri/aparc.DKTatlas+aseg.mgz"
            path_aparc = output_path + "/aparc.DKTatlas+aseg.nii"
            mri_convert_parcellation = f"mri_convert {path_aparc_dkt_mgz} {path_aparc}"
            subprocess.run([mri_convert_parcellation], shell=True)


        else:
            print("FS: error. Doing BET.")
            # Extract skull with BET
            path_BET_images = output_path + "/BET/"
            if not os.path.exists(path_BET_images):
                os.mkdir(path_BET_images)

            path_BET_image = path_BET_images + "/BET_image.nii.gz"
            bet_command = f"bet {path_MRI_image} {path_BET_image} -f 0.3 -R "
            subprocess.run([bet_command], shell=True)

            print("BET: ok")

            path_T1 = path_MRI_image
            path_T1_brain = path_BET_image

            aseg = False

        # Registration: PET + RMN

        # Do PET + RMN registration when dir doesnt exist
        path_PET_registration = output_path + "/Registration/"
        registration_image_path = path_PET_registration + "/Registration_image_PET_T1.nii.gz"
        pet_3d_image_path = path_PET_registration + "/PET_sum_image.nii.gz"

        if not os.path.exists(path_PET_registration):
            os.mkdir(path_PET_registration)

        PET_frames, pet_3d_image, register_image, tx_PET_2_RMN = reg.register_PET_MRI(path_PET_image, path_T1)

        # Write Euler transform (PET to RMN)
        path_transform = os.path.join(path_PET_registration, 'tx_PET_2_RMN.tfm')
        sitk.WriteTransform(tx_PET_2_RMN, path_transform)

        # Write registered image (RMN + PET) and PET sum image
        sitk.WriteImage(register_image, registration_image_path)
        sitk.WriteImage(pet_3d_image, pet_3d_image_path)

        path_individual_frames = []

        # Write every individual frame
        for i, PET_image in enumerate(PET_frames, start=1):
            path_frame = path_PET_registration + f'FRAME_{i}/'

            if not os.path.exists(path_frame):
                os.mkdir(path_frame)

            path_individual_frames.append(path_frame)

            path_PET_image = path_frame + "PET_image_register.nii.gz"
            sitk.WriteImage(PET_image, path_PET_image)

        print("Registration: ok")

        # Normalization:

        # Flirt
        # Apply FLIRT to T1
        path_flirt_images = output_path + "/FLIRT/"
        if not os.path.exists(path_flirt_images):
            os.mkdir(path_flirt_images)

        path_FLIRT_image = path_flirt_images + "/T1_Norm_MNI_152_FLIRT.nii.gz"
        path_FLIRT_tx = path_flirt_images + "/T1_Norm_MNI_152_FLIRT.mat"

        flirt_command = f"flirt -in {path_T1} -ref {path_MNI_152_T1} -out {path_FLIRT_image} -omat {path_FLIRT_tx} -interp trilinear -dof 12"
        subprocess.run([flirt_command], shell=True)

        print("FLIRT: ok")

        # Apply FLIRT transform to T1 only brain
        path_FLIRT_t1 = path_flirt_images + "/T1_brain_MNI_152_FLIRT.nii.gz"
        flirt_apply_transform_t1_brain_command = f"flirt -in {path_T1_brain} -applyxfm -init {path_FLIRT_tx} " \
                                            f"-out {path_FLIRT_t1} " \
                                            f"-paddingsize 0.0 -interp trilinear -ref {path_MNI_152_T1_brain}"
        subprocess.run([flirt_apply_transform_t1_brain_command], shell=True)

        # Apply FLIRT transform to PET Image
        path_FLIRT_PET = path_flirt_images + "/PET_Norm_MNI_152_FLIRT.nii.gz"
        flirt_apply_transform_PET_command = f"flirt -in {registration_image_path} -applyxfm -init {path_FLIRT_tx} " \
                                            f"-out {path_FLIRT_PET} " \
                                            f"-paddingsize 0.0 -interp trilinear -ref {path_MNI_152_T1_brain}"
        subprocess.run([flirt_apply_transform_PET_command], shell=True)

        # If exists, apply FLIRT transform to segmentation image
        if aseg:
            path_FLIRT_segmentation = path_flirt_images + "/aseg_Norm_MNI_152_FLIRT.nii.gz"
            flirt_apply_transform_aseg_command = f"flirt -in {path_segmentation} -applyxfm -init {path_FLIRT_tx} " \
                                            f"-out {path_FLIRT_segmentation} " \
                                            f"-paddingsize 0.0 -interp nearestneighbour -ref {path_MNI_152_T1_brain}"
            subprocess.run([flirt_apply_transform_aseg_command], shell=True)

        # ANT

        # Install ANTs and set environment variables in antsRegistrationSyNQuick.sh script
        ANT_path = "/home/sol/PET_MRI/CODIGOS/antsRegistrationSyN.sh"
        path_ANT_images = output_path + "/ANTs/"
        if not os.path.exists(path_ANT_images):
            os.mkdir(path_ANT_images)

        path_ANT_image = path_ANT_images + "/T1_Norm_MNI_152_ANT"
        ANT_command = f"{ANT_path} -d 3 -f {path_MNI_152_T1_brain} -m {path_FLIRT_t1} -o {path_ANT_image}"
        subprocess.run([ANT_command], shell=True)

        print("ANT: ok")

        # ANT transform
        # path_ANT_apply_transform = "/usr/local/ANTs/bin/antsApplyTransforms"
        path_ANT_transform1 = path_ANT_image + "1Warp.nii.gz"
        path_ANT_transform2 = path_ANT_image + "0GenericAffine.mat"

        # Apply ANT transform to PET
        path_ANT_PET = path_ANT_images + "/PET_Norm_MNI_152_ANT.nii.gz"

        ANT_transform_PET_command = f"antsApplyTransforms -d 3 -i {path_FLIRT_PET} " \
                                    f"-r {path_MNI_152_T1_brain} -o {path_ANT_PET}  " \
                                    f"-t {path_ANT_transform1} -t {path_ANT_transform2}"

        subprocess.run([ANT_transform_PET_command], shell=True)

        # If exists, apply FLIRT transform to segmentation image
        if aseg:
            path_ANT_segmentation = path_ANT_images + "/Aseg_Norm_MNI_152_ANT.nii.gz"

            ANT_transform_segmentation_command = f"antsApplyTransforms -d 3 -i {path_FLIRT_segmentation} " \
                                        f"-r {path_MNI_152_T1_brain} " \
                                        f"-o {path_ANT_segmentation}  " \
                                        f"-t {path_ANT_transform1} -t {path_ANT_transform2} -n NearestNeighbor "

            subprocess.run([ANT_transform_segmentation_command], shell=True)


    # Dataframe of name of labels
    df_labels_FS = pd.read_csv(labels_FS_csv_path)
    df_labels_Hammers = pd.read_csv(labels_Hammers_csv_path)

    path_CSV_files = output_path + "/CSV/"
    if not os.path.exists(path_CSV_files):
        os.mkdir(path_CSV_files)

    # Quantification FDG-PET in subject to process
    path_PET_final = path_ANT_PET
    if aseg:
        path_aseg_segmentation = path_ANT_segmentation
    else:
        path_aseg_segmentation = False

    df_subject_intensity, image_subject_intensity = quant.PET_FDG_quantification(path_PET_final, path_Hammers,
                                                                                 df_labels_Hammers,
                                                                                 path_second_segmentation=path_aseg_segmentation,
                                                                                 atlas="Hammers",
                                                                                 normalization=True)
    df_subject_intensity.to_csv(path_CSV_files + "subject_intensity.csv")

    # Quantification FDG-PET in ATLAS MNI152
    df_MNI152_intensity, image_MNI152_intensity = quant.PET_FDG_quantification(path_MNI_152_PET,
                                                                               path_Hammers,
                                                                               df_labels_Hammers,
                                                                               atlas="Hammers", normalization=True)
    df_MNI152_intensity.to_csv(path_CSV_files + "MNI152_intensity.csv")

    path_images = output_path + "/Synthetic_Image"
    if not os.path.exists(path_images):
        os.mkdir(path_images)

    sitk.WriteImage(image_subject_intensity, path_images + "/synthetic_image_subject.nii.gz")
    sitk.WriteImage(image_MNI152_intensity, path_images + "/synthetic_image_MNI152.nii.gz")

    # Compare subject with Atlas MNI152 dataset

    image_diff, df_diff = quant.image_change(df_subject_intensity, df_MNI152_intensity, path_Hammers)
    df_diff.to_csv(path_CSV_files + "Intensity_image_changes.csv")
    sitk.WriteImage(image_diff, path_images + "/synthetic_image_changes.nii.gz")

    # Intensity Normalization

    # Cerebellum values
    cerebellum_R = float((df_subject_intensity.loc[(df_subject_intensity['structure'] == 'cerebellum')
                                                   & (df_subject_intensity['hemisphere'] == 'R')]['mean_PET']).iloc[0])
    cerebellum_L = float((df_subject_intensity.loc[(df_subject_intensity['structure'] == 'cerebellum')
                                                   & (df_subject_intensity['hemisphere'] == 'L')]['mean_PET']).iloc[0])
    cerebellum = (cerebellum_R + cerebellum_L) / 2

    brain_image = sitk.ReadImage(path_PET_final)

    norm_image = norm.intensity_normalization(brain_image, mode=("cerebellum", cerebellum))
    sitk.WriteImage(norm_image, output_path + "/intensity_normalization_image.nii.gz")

    # Bar charts
    path_chart = output_path + "/charts"

    if not os.path.exists(path_chart):
        os.mkdir(path_chart)
