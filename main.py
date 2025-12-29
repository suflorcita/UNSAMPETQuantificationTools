#!/usr/bin/env python3

import processPET4D as reg
import argparse
import SimpleITK as sitk
import subprocess
import pandas as pd
import PETquantification as quant
import normalization as norm
import registration_tool as regtools
# import plots
import os
import shutil, math

if __name__ == '__main__':
    # Create the parser
    parser = argparse.ArgumentParser(description="A PET Quantification Tools")

    # MRI and PET images
    parser.add_argument('-m', '--mri', help='MRI image')  # MRI image
    parser.add_argument('-p', '--pet', help='PET image')  # PET image
    parser.add_argument("-o", "--output", help='Name of output dir')  # Output path
    parser.add_argument("-s", "--subjet", help='Name of subject')  # Subject name
    parser.add_argument("-t", "--pet-template", help='PET template in case there is no MRI image')
    parser.add_argument("-f", "--freesurfer-dir", help='path of Freesurfer files')
    parser.add_argument("-a", "--atlas", help='Atlas for quantification. Default: Hammers')


    # Add the --no-flirt argument
    parser.add_argument(
        '--no-flirt',
        action='store_false',
        dest='flirt',
        help='Disable the FLIRT command. Use this flag to skip the FLIRT operation.'
    )

    # By default, args.flirt will be True unless --no-flirt is specified
    parser.set_defaults(flirt=True)

    # Add the --no-ant argument
    parser.add_argument(
        '--no-ant',
        action='store_false',
        dest='ants',
        help='Disable the ANT command. Use this flag to skip the ANT operation.'
    )

    # By default, args.ants will be True unless --no-flirt is specified
    parser.set_defaults(ants=True)

    # Add the --no-freesurfer argument
    parser.add_argument(
        '--no-freesurfer',
        action='store_false',
        dest='freesurfer',
        help='Disable the Freesurfer command. Use this flag to skip the Freesurfer operation.'
    )

    # By default, args.flirt will be True unless --no-flirt is specified
    parser.set_defaults(freesurfer=True)

    # Add the --no-smoothed argument
    parser.add_argument(
        '--no-smoothed',
        action='store_false',
        dest='smoothed',
        help='Disable the use of smoothed PET images. Use this flag to process raw PET instead.'
    )

    # By default, smoothed=True unless --no-smoothed is specified
    parser.set_defaults(smoothed=True)

    # Parse the arguments
    args = parser.parse_args()

    path_PET_image = args.pet
    path_MRI_image = args.mri
    path_PET_template = args.pet_template
    subject = args.subjet
    output_path = args.output
    freesurfer_dir = args.freesurfer_dir
    atlas = args.atlas

    freesurfer = args.freesurfer



    # Get the directory where the script is located
    # Check if the current file is a symbolic link
    if os.path.islink(__file__):
        # Resolve the symbolic link
        script_path = os.readlink(__file__)
        # Get the directory containing the symbolic link
        script_dir = os.path.dirname(os.path.abspath(script_path))
    else:
        # Use the directory of the original file
        script_dir = os.path.dirname(os.path.abspath(__file__))

    if path_PET_image == None:
        path_PET_image = os.path.join(script_dir, "./FDG-PET.nii.gz")

    if os.path.exists(os.path.join(os.path.join(script_dir, "./T1_MPRAGE.nii.gz"))):
        path_MRI_image = os.path.join(script_dir, "./T1_MPRAGE.nii.gz")

    if subject == None:
        subject = os.path.join(script_dir,  "anonymous")

    if output_path == None:
        output_path = os.path.join(script_dir,  "./anonymous")

    if freesurfer_dir == None:
        freesurfer_dir = ""

    if atlas == None:
        atlas = "Hammers"



    # Path ATLAS

    path_MNI_152_T1 = os.path.join(script_dir, 'data', 'atlas', 'MNI152_T1_1mm.nii.gz')
    path_MNI_152_T1_brain = os.path.join(script_dir, 'data', 'atlas', 'MNI152_T1_1mm_brain.nii.gz')
    path_MNI_152_T1_gm_mask = os.path.join(script_dir, 'data', 'atlas', 'mni_icbm152_gm_mask.nii.gz')

    path_MNI_152_PET = os.path.join(script_dir, 'data', 'atlas', "MNI152_PET_1mm.nii")
    path_Hammers = os.path.join(script_dir, 'data', 'atlas', "Hammers_mith-n30r95-MaxProbMap-gm-MNI152-SPM12.nii.gz")
    path_DKT = os.path.join(script_dir, 'data', 'atlas', "Desikan_Killiany_MNI_SPM12.nii.gz")
    path_CerebrA = os.path.join(script_dir, 'data', 'atlas', "CerebrA.nii")
    path_PET_template_CN =  os.path.join(script_dir, 'data', 'atlas',"ATLAS/AtlasPETFDG_FiltOutliers_CN.nii")

    if atlas == "DKT":
        path_atlas = path_DKT
    elif atlas == "CerebrA":
        path_atlas = path_CerebrA
    elif atlas == "GM":
        path_atlas = path_MNI_152_T1_gm_mask
    else:
        path_atlas = path_Hammers

    # labels segmentation
    path_labels_FS_csv =  os.path.join(script_dir, 'data', 'labels', "labels_FS.csv")
    path_labels_DKT_csv = os.path.join(script_dir,  'data', 'labels', "labels_DKT.csv")
    path_labels_Hammers_csv = os.path.join(script_dir, 'data', 'labels', "labels_Hammers.csv")
    path_labels_CerebrA_csv = os.path.join(script_dir, 'data', 'labels', "labels_CerebrA.csv")
    path_labels_gm_csv = os.path.join(script_dir, 'data', 'labels', "labels_gm_mask.csv")

    if path_atlas == path_DKT:
        path_labels = path_labels_DKT_csv
    elif path_atlas == path_CerebrA:
        path_labels = path_labels_CerebrA_csv
    elif path_atlas == path_MNI_152_T1_gm_mask:
        path_labels = path_labels_gm_csv
    else:
        path_labels = path_labels_Hammers_csv


    # Output
    output_path = os.path.join(output_path, subject)

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # Use all disponible threads
    os.environ["TF_NUM_INTRAOP_THREADS"] = str(os.cpu_count())
    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(os.cpu_count())
    path_ANT = "antsRegistrationSyN.sh"
    path_ANT_apply_transform = "antsApplyTransforms"

    process_frames = False


    if path_PET_template is None:
        # If not provided through the terminal, it will attempt to use the MNI152 by default
        path_PET_template = path_MNI_152_PET


    # If PET is a dir convert to an image 4D
    if os.path.isdir(path_PET_image):
        output_fdg_path = os.path.join(output_path, "unique_fdg_pet_image.nii.gz")

        # List all PET NIfTI files in the directory
        fdg_files = sorted([
            os.path.join(path_PET_image, f)
            for f in os.listdir(path_PET_image)
            if f.endswith('.nii') or f.endswith('.nii.gz')
        ])

        # If no PET files found
        if len(fdg_files) == 0:
            path_PET_image = None

        # If only one PET file found, use it directly
        elif len(fdg_files) == 1:
            path_PET_image = fdg_files[0]

        # If multiple PET files found, join them into a 4D image
        else:
            pet_images = [sitk.ReadImage(fdg_file) for fdg_file in fdg_files]
            combined_fdg_image = sitk.JoinSeries(pet_images)
            sitk.WriteImage(combined_fdg_image, output_fdg_path)
            path_PET_image = output_fdg_path
            process_frames = True

    path_PET_final = None
    path_normalized_MRI_image = None
    path_normalized_PET_image = None

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    if not path_MRI_image:
        aseg = False
        freesurfer = False

        PET_image = sitk.ReadImage(path_PET_image)

        # --- PET 4D case: split, register & sum ---
        if len(PET_image.GetSize()) == 4:
            # Split into frames
            PET_images = reg.image4D_to_3D(PET_image)

            # Register frames among themselves and get summed PET
            sum_pet_3d_image, PET_images = reg.register_and_sum_pet_images(PET_images)

            # Save frames
            frames_dir = os.path.join(output_path, "Frames")
            os.makedirs(frames_dir, exist_ok=True)
            for i, frame in enumerate(PET_images):
                frame_path = os.path.join(frames_dir, f"PET_frame_{i}.nii.gz")
                sitk.WriteImage(frame, frame_path)

            # Save summed PET
            sum_path = os.path.join(output_path, "Summed_PET_image.nii.gz")
            sitk.WriteImage(sum_pet_3d_image, sum_path)

            # Work with summed PET as main image
            path_PET_image = sum_path

        # --- Run pipeline on summed PET ---
        path_PET_final, path_normalized_PET_image, transform_files = regtools.pet_registration_pipeline(
            path_PET_image=path_PET_image,
            path_PET_template=path_PET_template,
            output_path=output_path,
            args=args,
            path_ANT="antsRegistrationSyN.sh",
            path_ANT_apply="antsApplyTransforms"
        )

        # --- Ensure transforms exist even if pipeline skipped them ---

        transform_files = regtools.ensure_transform_files(output_path, transform_files)

        # --- Normalize each frame ---
        if process_frames == True:

            # --- Apply the same transforms to each individual frame ---
            frames_dir = os.path.join(output_path, "Frames")
            frames_to_template_dir = os.path.join(output_path, "Frames_to_template")
            os.makedirs(frames_to_template_dir, exist_ok=True)

            # Get all PET frame files (registered to the summed PET previously)
            PET_frame_paths = [os.path.join(frames_dir, f) for f in os.listdir(frames_dir) if f.endswith(".nii.gz")]

            # Normalize each frame to the template
            PET_frames_normalized = regtools.normalize_pet_frames_to_template(
                PET_frames=PET_frame_paths,
                path_PET_template=path_PET_template,
                transform_files=transform_files,
                output_dir=frames_to_template_dir
            )


    else:
        #freesurfer_dir = os.path.join(output_path, "Freesurfer")
        if freesurfer:
            if not freesurfer_dir:
                freesurfer_dir = os.path.join(output_path, "Freesurfer")
                if not os.path.exists(freesurfer_dir):
                    os.mkdir(freesurfer_dir)

                recon_all_command = f"recon-all -s {subject} -i {path_MRI_image} -sd {freesurfer_dir} -all -parallel -openmp 16"

                subprocess.run([recon_all_command], shell=True)
                freesurfer_dir = os.path.join(freesurfer_dir, subject)

        path_T1_FS_mgz = os.path.join(freesurfer_dir, "mri", "T1.mgz")
        print(path_T1_FS_mgz)

        if os.path.exists(path_T1_FS_mgz):
            print("FS: ok. T1 is FS T1.mgz")
            path_T1 = output_path + "/T1_FS.nii"
            mri_convert_t1_FS = f"mri_convert {path_T1_FS_mgz} {path_T1}"
            subprocess.run([mri_convert_t1_FS], shell=True)

            path_T1_brain_mgz = os.path.join(freesurfer_dir, "mri", "brainmask.mgz")

            path_T1_brain = output_path + "/T1_FS_onlybrain.nii"
            mri_convert_t1_only_brain = f"mri_convert {path_T1_brain_mgz} {path_T1_brain}"
            subprocess.run([mri_convert_t1_only_brain], shell=True)

            # Segmentation
            aseg = True
            path_segmentation_mgz = os.path.join(freesurfer_dir, "mri","aseg.mgz")
            path_segmentation = output_path + "/aseg.nii"
            mri_convert_segmentation = f"mri_convert {path_segmentation_mgz} {path_segmentation}"
            subprocess.run([mri_convert_segmentation], shell=True)

            # Parcellation
            path_aparc_mgz = os.path.join(freesurfer_dir, "mri", "aparc+aseg.mgz")
            path_aparc = os.path.join(output_path, "aparc+aseg.nii")

            # MRI convert: DKT atlas
            mri_convert_parcellation = f"mri_convert {path_aparc_mgz} {path_aparc}"
            if not os.path.exists(path_aparc):
                subprocess.run([mri_convert_parcellation], shell=True)

            # MRI convert: Destrieux atlas
            path_aparc_destrieux_mgz =  os.path.join(freesurfer_dir, "mri","aparc.a2009s+aseg.mgz")
            path_aparc_destrieux = os.path.join(freesurfer_dir, "aparc.a2009s+aseg.nii")

            mri_convert_parcellation_destrieux = f"mri_convert {path_aparc_destrieux_mgz} {path_aparc_destrieux}"

            if not os.path.exists(path_aparc_destrieux):
                subprocess.run([mri_convert_parcellation_destrieux], shell=True)





        else:
            print("FS: error. Doing BET.")

            # Extract skull with BET
            path_BET_images = os.path.join(output_path,"BET")
            if not os.path.exists(path_BET_images):
                os.mkdir(path_BET_images)

            path_BET_image =  os.path.join(path_BET_images,"BET_image.nii.gz")
            bet_command = f"bet {path_MRI_image} {path_BET_image} -R -f 0.35 -g 0 -o -m "
            subprocess.run([bet_command], shell=True)

            print("BET: ok")

            path_T1 = path_MRI_image
            path_T1_brain = path_BET_image

            aseg = False

        # Registration: PET + RMN

        # PET + RMN registration
        path_PET_registration_dir = os.path.join(output_path, "Registration")
        path_PET_registered_to_MRI = os.path.join(path_PET_registration_dir, "Registration_image_PET_T1.nii.gz")
        pet_3d_image_path = os.path.join(path_PET_registration_dir, "PET_sum_image.nii.gz")

        if not os.path.exists(path_PET_registration_dir):
            os.mkdir(path_PET_registration_dir)

        PET_frames, PET_3d_image, PET_registered_to_MRI, tx_PET_2_MRI = reg.register_PET_MRI(path_PET_image, path_T1)
        # PET_frames: List of registered PET frames to the first image
        # PET_3D_image: Sum of all PET frames (if more than one) registered to the first one.
        # PET_registered_to_MRI: Sum of all PET frames registered to MRI
        # tx_PET_2_RMN: Transformation of PET registered to MRI.

        # Write Euler transform (PET to MRI)
        path_transform = os.path.join(path_PET_registration_dir, 'tx_PET_2_RMN.tfm')
        sitk.WriteTransform(tx_PET_2_MRI, path_transform)

        # Write registered image (PET to MRI) and the sum of all PET frames
        sitk.WriteImage(PET_registered_to_MRI, path_PET_registered_to_MRI)
        sitk.WriteImage(PET_3d_image, pet_3d_image_path)

        path_individual_frames = []

        if len(PET_frames) > 1:

            # Write each individual frame
            # The frames are registered to a reference image, which is the first frame
            for i, registered_frame in enumerate(PET_frames, start=1):
                path_frame = os.path.join(path_PET_registration_dir, f"FRAME_{i}")

                if not os.path.exists(path_frame):
                    os.mkdir(path_frame)

                path_individual_frames.append(path_frame)

                # Write the registered frame
                path_PET_image = os.path.join(path_frame, f"PET_frame_to_reference_{i}.nii.gz")
                sitk.WriteImage(registered_frame, path_PET_image)

                # Register the frame to the MRI image using the PET-to-MRI transform
                # then write the image.
                path_PET_frame_registered_to_MRI = os.path.join(path_frame, f"PET_frame_to_MRI_{i}.nii.gz")


                PET_frame_registered_to_MRI = reg.read_and_apply_tx(path_transform, registered_frame, ref_image=PET_registered_to_MRI)
                sitk.WriteImage(PET_frame_registered_to_MRI, path_PET_frame_registered_to_MRI)

        path_PET_final = path_PET_registered_to_MRI  # Set the reference PET as the PET image registered to MRI

        # Normalization:
        #
        # # Flirt
        # # Apply FLIRT to T1
        # if args.flirt:
        #     path_flirt_images = os.path.join(output_path, "FLIRT")
        #     if not os.path.exists(path_flirt_images):
        #         os.mkdir(path_flirt_images)
        #
        #     path_FLIRT_image = os.path.join(path_flirt_images, "T1_Norm_MNI_152_FLIRT.nii.gz")
        #     path_FLIRT_tx = os.path.join(path_flirt_images, "T1_Norm_MNI_152_FLIRT.mat")
        #
        #     flirt_command = f"flirt -in {path_T1} -ref {path_MNI_152_T1} -out {path_FLIRT_image} -omat {path_FLIRT_tx} -interp trilinear -dof 12"
        #     subprocess.run([flirt_command], shell=True)
        #
        #
        #     # Apply FLIRT transform to T1 only brain
        #     path_flirt_MRI = os.path.join(path_flirt_images, "T1_brain_MNI_152_FLIRT.nii.gz")
        #     flirt_apply_transform_t1_brain_command = f"flirt -in {path_T1_brain} -applyxfm -init {path_FLIRT_tx} " \
        #                                         f"-out {path_flirt_MRI} " \
        #                                         f"-paddingsize 0.0 -interp trilinear -ref {path_MNI_152_T1_brain}"
        #
        #
        #     path_normalized_MRI_image = path_flirt_MRI
        #     subprocess.run([flirt_apply_transform_t1_brain_command], shell=True)
        #     print(flirt_command)
        #
        #
        #     # Apply FLIRT transform to PET Image
        #     path_FLIRT_PET = os.path.join(path_flirt_images, "PET_Norm_MNI_152_FLIRT.nii.gz")
        #     flirt_apply_transform_PET_command = f"flirt -in {path_PET_registered_to_MRI} -applyxfm -init {path_FLIRT_tx} " \
        #                                         f"-out {path_FLIRT_PET} " \
        #                                         f"-paddingsize 0.0 -interp trilinear -ref {path_MNI_152_T1_brain}"
        #     subprocess.run([flirt_apply_transform_PET_command], shell=True)
        #
        #     path_PET_final = path_FLIRT_PET
        #     path_normalized_PET_image = path_FLIRT_PET
        #     # If exists, apply FLIRT transform to segmentation image
        #     if aseg:
        #         path_FLIRT_segmentation = path_flirt_images + "/aseg_Norm_MNI_152_FLIRT.nii.gz"
        #         flirt_apply_transform_aseg_command = f"flirt -in {path_segmentation} -applyxfm -init {path_FLIRT_tx} " \
        #                                         f"-out {path_FLIRT_segmentation} " \
        #                                         f"-paddingsize 0.0 -interp nearestneighbour -ref {path_MNI_152_T1_brain}"
        #         path_normalized_aseg_segmentation = path_FLIRT_segmentation
        #         subprocess.run([flirt_apply_transform_aseg_command], shell=True)
        # else:
        #     path_flirt_MRI = None
        #     path_FLIRT_PET = None
        #     path_FLIRT_segmentation = None

        # --- FLIRT removed ---
        path_FLIRT_MRI = None
        path_FLIRT_PET = None
        path_FLIRT_segmentation = None
        # ---------------------


        # ANT

        if args.ants:
            # T1
            if path_FLIRT_MRI == None:
                path_ANT_input_t1 = path_T1_brain
            else:
                path_ANT_input_t1 = path_FLIRT_MRI

            # PET
            if path_FLIRT_PET == None:
                path_ANT_input_PET = path_PET_registered_to_MRI
            else:
                path_ANT_input_PET = path_FLIRT_PET

            if aseg:
                if path_FLIRT_segmentation == None:
                    path_ANT_input_aseg = path_segmentation
                else:
                    path_ANT_input_aseg = path_FLIRT_segmentation


            # Install ANTs and set environment variables in antsRegistrationSyNQuick.sh script
            path_ANT_images = os.path.join(output_path, "ANTs")
            if not os.path.exists(path_ANT_images):
                os.mkdir(path_ANT_images)

            path_ANT_image = os.path.join(path_ANT_images, "T1_Norm_MNI_152_ANT")
            ANT_command = f"{path_ANT} -d 3 -f {path_MNI_152_T1_brain} -m {path_ANT_input_t1} -o {path_ANT_image}"
            print(ANT_command)
            subprocess.run([ANT_command], shell=True)

            path_normalized_MRI_image = path_ANT_image +  "Warped.nii.gz"

            print("ANT: ok")

            # ANT transform
            path_ANT_transform1 = path_ANT_image + "1Warp.nii.gz"
            path_ANT_transform2 = path_ANT_image + "0GenericAffine.mat"

            # Apply ANT transform to PET
            path_ANT_PET = os.path.join(path_ANT_images, "PET_Norm_MNI_152_ANT.nii.gz")

            ANT_transform_PET_command = f"{path_ANT_apply_transform} -d 3 -i {path_ANT_input_PET} " \
                                        f"-r {path_MNI_152_T1_brain} -o {path_ANT_PET}  " \
                                        f"-t {path_ANT_transform1} -t {path_ANT_transform2}"

            subprocess.run([ANT_transform_PET_command], shell=True)

            path_PET_final = path_ANT_PET
            path_normalized_PET_image = path_ANT_PET

            # If exists, apply ANT transform to segmentation image
            if aseg:
                path_ANT_segmentation = path_ANT_images + "/Aseg_Norm_MNI_152_ANT.nii.gz"

                ANT_transform_segmentation_command = f"{path_ANT_apply_transform} -d 3 -i {path_ANT_input_aseg} " \
                                            f"-r {path_MNI_152_T1_brain} " \
                                            f"-o {path_ANT_segmentation}  " \
                                            f"-t {path_ANT_transform1} -t {path_ANT_transform2} -n NearestNeighbor "
                path_normalized_aseg_segmentation = path_ANT_segmentation

                subprocess.run([ANT_transform_segmentation_command], shell=True)


    # Dataframe of name of labels
    df_labels_FS = pd.read_csv(path_labels_FS_csv)
    df_labels_atlas = pd.read_csv(path_labels)

    # --- Ensure normalized PET (and MRI if exists) ---
    normalized = regtools.ensure_normalized_image(
        output_path=output_path,
        path_MRI_image=path_MRI_image,
        freesurfer=freesurfer,
    )



    path_normalized_PET_image = normalized["PET"]
    path_normalized_MRI_image = normalized["MRI"]
    path_normalized_aseg_segmentation = normalized["aseg"]


    path_CSV_files = os.path.join(output_path, "CSV")

    if not os.path.exists(path_CSV_files):
        os.mkdir(path_CSV_files)

    path_normalization_images = os.path.join(output_path, "Normalization")

    if not os.path.exists(path_normalization_images):
        os.mkdir(path_normalization_images)


    # Quantification FDG-PET Using a Non FS atlas
    image_subject_PET = sitk.ReadImage(path_normalized_PET_image)

    if path_MRI_image:
        image_subject_MRI = sitk.ReadImage(path_normalized_MRI_image)


    image_atlas = sitk.ReadImage(path_atlas)
    image_aseg_segmentation = None

    if aseg:
        image_aseg_segmentation = sitk.ReadImage(path_normalized_aseg_segmentation)


    df_subject_intensity, image_subject_intensity = quant.PET_FDG_quantification(image_subject_PET, image_atlas,
                                                                                 df_labels_atlas,
                                                                                 second_segmentation=image_aseg_segmentation,
                                                                                 atlas=atlas)



    df_subject_intensity.to_csv(os.path.join(path_CSV_files,f"subject_uptake_{atlas}.csv"))

    values_subject_mean_activity_atlas = df_subject_intensity["mean_PET"]
    n_label_subject_mean_activity_atlas = df_subject_intensity["n_label"]

    if atlas == "DKT":
        name_structures = df_subject_intensity["structure"]
    else:
        name_structures = df_subject_intensity["structure"] + '-' + df_subject_intensity["hemisphere"]

    if not atlas == "GM":
        # Extract Cerebellum values
        if atlas == "DKT":
            cerebellum_name = "Cerebellum-Cortex"
        elif atlas == "CerebrA":
            cerebellum_name = "Cerebellum-Gray-Matter"
        else:
            cerebellum_name = "cerebellum"


        # Extract cerebellum values
        if atlas == "DKT":
            cerebellum_subject = float((df_subject_intensity.loc[(df_subject_intensity['structure'] == cerebellum_name)]
            ['mean_PET']).iloc[0])
        else:
            cerebellum_R_subject = float((df_subject_intensity.loc[(df_subject_intensity['structure'] == cerebellum_name)
                                                                   & (df_subject_intensity['hemisphere'] == 'R')][
                'mean_PET']).iloc[0])

            cerebellum_L_subject = float((df_subject_intensity.loc[(df_subject_intensity['structure'] == cerebellum_name)
                                                                   & (df_subject_intensity['hemisphere'] == 'L')][
                'mean_PET']).iloc[0])
            cerebellum_subject = (cerebellum_R_subject + cerebellum_L_subject) / 2


        norm_subject_image_cerebellum = norm.intensity_normalization(image_subject_PET, mode="scalar",
                                                                     scalar=cerebellum_subject)


        # Quantification Normalization Cerebellum
        df_subject_intensity_norm_cerebellum_atlas, image_subject_intensity_norm_cerebellum = quant.PET_FDG_quantification(
            norm_subject_image_cerebellum, image_atlas,
            df_labels_atlas,
            atlas=atlas)

        values_subject_normalization_cerebellum_atlas = df_subject_intensity_norm_cerebellum_atlas["mean_PET"]

        output_normalization_cerebellum = os.path.join(path_normalization_images,
                                                       f'intensity_normalization_cerebellum_{atlas}.nii.gz')

        sitk.WriteImage(norm_subject_image_cerebellum, output_normalization_cerebellum)

    ### Normalized to mean value
    norm_subject_image_mean = norm.intensity_normalization(image_subject_PET)
    output_normalization_mean = os.path.join(path_normalization_images,
                                             'intensity_normalization_mean.nii.gz')
    sitk.WriteImage(norm_subject_image_mean, output_normalization_mean)

    # Quantification Normalization mean Value
    df_subject_intensity_norm_mean_atlas, image_subject_intensity_norm_mean = quant.PET_FDG_quantification(
        norm_subject_image_mean, image_atlas,
        df_labels_atlas,
        atlas=atlas)

    values_subject_normalization_mean = df_subject_intensity_norm_mean_atlas["mean_PET"]

    # Subject Intensity Normalization CSV
    output_subject_csv_atlas = os.path.join(path_CSV_files, f"subject_normalization_values_{atlas}.csv")

    # --- Create the DataFrame depending on the atlas type ---
    if atlas == "GM":
        # If the atlas is GM, exclude normalization to cerebellum uptake values
        df_normalization_subject = pd.DataFrame({
            'Structure': name_structures,
            'n_label': n_label_subject_mean_activity_atlas,
            'Regional uptake mean values': values_subject_mean_activity_atlas,
            'Normalization to total brain mean value': values_subject_normalization_mean,
        })
    else:
        # For all other atlases, include both normalizations
        df_normalization_subject = pd.DataFrame({
            'Structure': name_structures,
            'n_label': n_label_subject_mean_activity_atlas,
            'Regional uptake mean values': values_subject_mean_activity_atlas,
            'Normalization to total brain mean value': values_subject_normalization_mean,
            'Normalization to cerebellum uptake values': values_subject_normalization_cerebellum_atlas,
        })

    df_normalization_subject.to_csv(output_subject_csv_atlas)

    # Quantification FDG-PET in ATLAS MNI152
    MNI_152_PET_image = sitk.ReadImage(path_MNI_152_PET)
    df_MNI152_intensity, image_MNI152_intensity = quant.PET_FDG_quantification(MNI_152_PET_image,
                                                                               image_atlas,
                                                                               df_labels_atlas,
                                                                               atlas=atlas)
    df_MNI152_intensity.to_csv(os.path.join(path_CSV_files, f"MNI152_intensity_{atlas}.csv"))

    values_MNI152_mean_activity_atlas = df_MNI152_intensity["mean_PET"]


    # MNI152 Normalization to cerebellum
    if not atlas == "GM":
        # Extract cerebellum values
        if atlas == "DKT":
            cerebellum_MNI152 = float((df_MNI152_intensity.loc[(df_MNI152_intensity['structure'] == cerebellum_name)
                                                      ]['mean_PET']).iloc[0])
        else:
            cerebellum_R_MNI152 = float((df_MNI152_intensity.loc[(df_MNI152_intensity['structure'] == cerebellum_name)
                                                           & (df_MNI152_intensity['hemisphere'] == 'R')]['mean_PET']).iloc[0])
            cerebellum_L_MNI152 = float((df_MNI152_intensity.loc[(df_MNI152_intensity['structure'] == cerebellum_name)
                                                           & (df_MNI152_intensity['hemisphere'] == 'L')]['mean_PET']).iloc[0])
            cerebellum_MNI152 = (cerebellum_R_MNI152 + cerebellum_L_MNI152) / 2

        norm_MNI152_image_cerebellum = norm.intensity_normalization(MNI_152_PET_image, mode="scalar",
                                                             scalar=cerebellum_MNI152)

        # Quantification Normalization Cerebellum in ATLAS MNI152
        df_MNI152_intensity_norm_cerebellum_atlas, image_MNI152_intensity_norm_cerebellum = quant.PET_FDG_quantification(
            norm_MNI152_image_cerebellum, image_atlas,
            df_labels_atlas,
            atlas=atlas)

        values_MNI152_normalization_cerebellum_atlas = df_MNI152_intensity_norm_cerebellum_atlas["mean_PET"]

    ### Normalized to mean value
    norm_MNI152_image_mean = norm.intensity_normalization(MNI_152_PET_image)

    # Quantification Normalization mean Value
    df_MNI152_intensity_norm_mean_atlas, image_MNI152_intensity_norm_mean = quant.PET_FDG_quantification(
        norm_MNI152_image_mean, image_atlas,
        df_labels_atlas,
        atlas=atlas)

    values_MNI152_normalization_mean = df_MNI152_intensity_norm_mean_atlas["mean_PET"]

    # MNI152 Intensity Normalization CSV
    output_MNI152_csv_name_atlas = os.path.join(path_CSV_files, f"MNI152_normalization_values_{atlas}.csv")

    # --- Create the DataFrame depending on the atlas type ---
    if atlas == "GM":
        # If the atlas is GM, skip normalization to cerebellum uptake values
        df_normalization_MNI152 = pd.DataFrame({
            'Structure': name_structures,
            'n_label': n_label_subject_mean_activity_atlas,
            'Regional uptake mean values': values_MNI152_mean_activity_atlas,
            'Normalization to total brain mean value': values_MNI152_normalization_mean,
        })
    else:
        # For all other atlases, include both normalizations
        df_normalization_MNI152 = pd.DataFrame({
            'Structure': name_structures,
            'n_label': n_label_subject_mean_activity_atlas,
            'Regional uptake mean values': values_MNI152_mean_activity_atlas,
            'Normalization to total brain mean value': values_MNI152_normalization_mean,
            'Normalization to cerebellum uptake values': values_MNI152_normalization_cerebellum_atlas,
        })

    # --- Save the resulting CSV file ---


    df_normalization_MNI152.to_csv(output_MNI152_csv_name_atlas)

    # Generate Synthetic image
    path_synthetic_images = os.path.join(output_path,"Synthetic_Image")

    if not os.path.exists(path_synthetic_images):
        os.mkdir(path_synthetic_images)

    sitk.WriteImage(image_subject_intensity, os.path.join(path_synthetic_images , "synthetic_image_subject.nii.gz"))
    sitk.WriteImage(image_MNI152_intensity, os.path.join(path_synthetic_images , "synthetic_image_MNI152.nii.gz"))


    # Compare subject with Atlas MNI152 dataset
    if not atlas == "GM":
        image_diff_cerebellum, df_diff_cerebellum = quant.image_change(df_normalization_subject,
                                                                       df_normalization_MNI152,
                                                                       image_atlas,
                                                                       mode="cerebellum")
        path_synthetic_image_norm_cerebellum = os.path.join(path_synthetic_images, f"synthetic_cerebellum_image_changes_{atlas}.nii.gz")
        sitk.WriteImage(image_diff_cerebellum, path_synthetic_image_norm_cerebellum)


    image_diff_mean, df_diff_mean = quant.image_change(df_normalization_subject,
                                                                   df_normalization_MNI152,
                                                                   image_atlas,
                                                                   mode="mean")
    # Images
    path_synthetic_image_norm_mean = os.path.join(path_synthetic_images, "synthetic_mean_image_changes.nii.gz")


    sitk.WriteImage(image_diff_mean, path_synthetic_image_norm_mean)

    # CSV
    # --- Create the DataFrame depending on the atlas type ---
    if atlas == "GM":
        # If the atlas is GM, skip cerebellum normalization
        df_changes = pd.DataFrame({
            'Structure': name_structures,
            'Change between images normalized to mean uptake values': df_diff_mean["change"],
        })
    else:
        # For all other atlases, include both normalizations
        df_changes = pd.DataFrame({
            'Structure': name_structures,
            'Change between images normalized to cerebellum': df_diff_cerebellum["change"],
            'Change between images normalized to mean uptake values': df_diff_mean["change"],
        })

    if freesurfer:
        output_synthetic_image_Destrieux = os.path.join(path_synthetic_images, "synthetic_image_subject_Destrieux.nii.gz")
        output_synthetic_image_DKT = os.path.join(path_synthetic_images, "synthetic_image_subject_DKT.nii.gz")

        # Read image
        aparc_dkt_image = sitk.ReadImage(path_aparc)
        aparc_destrieux_image = sitk.ReadImage(path_aparc_destrieux)


        # # New mask of only brain (without ventricles)
        path_T1_only_brain = os.path.join(output_path, "T1_FS_onlybrain.nii")
        T1_only_brain_with_ventricles = sitk.ReadImage(path_T1_only_brain)

        T1_only_brain_without_ventricles = quant.new_mask_without_ventricles(T1_only_brain_with_ventricles,
                                                                       aparc_dkt_image)

        output_brain_without_ventricles = os.path.join(output_path, "T1_FS_onlybrain_without_vent.nii")
        sitk.WriteImage(T1_only_brain_without_ventricles, output_brain_without_ventricles)

        output_csv_freesurfer_dkt = os.path.join(path_CSV_files, "normalization_values_freesurfer_dkt.csv")
        output_csv_freesurfer_destrieux = os.path.join(path_CSV_files, "normalization_values_freesurfer_destrieux.csv")


        registred_PET_image = sitk.ReadImage(path_PET_registered_to_MRI)


        # DKT
        df_subject_intensity_freesurfer_dkt, image_subject_intensity_freesurfer_dkt = quant.PET_FDG_quantification(
            registred_PET_image,
            aparc_dkt_image,
            df_labels_FS,
            atlas="FS")

        # Destrieux
        df_subject_intensity_freesurfer_destrieux, image_subject_intensity_freesurfer_destrieux = quant.PET_FDG_quantification(
            registred_PET_image,
            aparc_destrieux_image,
            df_labels_FS,
            atlas="FS")

        # Write synthetic images
        sitk.WriteImage(image_subject_intensity_freesurfer_dkt, output_synthetic_image_DKT)
        sitk.WriteImage(image_subject_intensity_freesurfer_destrieux, output_synthetic_image_Destrieux)

        # Name of the structures
        name_structures_dkt = df_subject_intensity_freesurfer_dkt["structure"] + '-' + \
                              df_subject_intensity_freesurfer_dkt["hemisphere"]
        name_structures_destrieux = df_subject_intensity_freesurfer_destrieux["structure"] + '-' + \
                                    df_subject_intensity_freesurfer_destrieux[
                                        "hemisphere"]

        # Mean Regional PET uptake values
        values_mean_activity_dkt = df_subject_intensity_freesurfer_dkt["mean_PET"]
        values_mean_activity_destrieux = df_subject_intensity_freesurfer_destrieux["mean_PET"]

        ### INTENSITY NORMALIZED IMAGES ###
        # Normalized to cerebellum
        # Extract Cerebellum values
        # Quantify using Hammers atlas

        cerebellum_freesurfer_name = 'Cerebellum-Cortex'

        # Extract cerebellum values (both DKT and aseg have the same value for cerebellum)
        cerebellum_R = float((df_subject_intensity_freesurfer_dkt.loc[
            (df_subject_intensity_freesurfer_dkt['structure'] == cerebellum_freesurfer_name)
            & (df_subject_intensity_freesurfer_dkt['hemisphere'] == 'R')]['mean_PET']).iloc[0])
        cerebellum_L = float((df_subject_intensity_freesurfer_dkt.loc[
            (df_subject_intensity_freesurfer_dkt['structure'] == cerebellum_freesurfer_name)
            & (df_subject_intensity_freesurfer_dkt['hemisphere'] == 'L')]['mean_PET']).iloc[0])
        cerebellum = (cerebellum_R + cerebellum_L) / 2

        # Normalize image
        norm_image_freesurfer_cerebellum_dkt = norm.intensity_normalization(registred_PET_image,
                                                                            mode="scalar",
                                                                            scalar=cerebellum)
        output_normalization_cerebellum_freesurfer = os.path.join(path_normalization_images,
                                                                  'intensity_normalization_freesurfer_cerebellum.nii.gz')

        sitk.WriteImage(norm_image_freesurfer_cerebellum_dkt, output_normalization_cerebellum_freesurfer)

        # Quantification Normalized to Cerebellum

        df_subject_intensity_freesurfer_cerebellum_dkt, image_subject_intensity_freesurfer_cerebellum_dkt = quant.PET_FDG_quantification(
            norm_image_freesurfer_cerebellum_dkt,
            aparc_dkt_image,
            df_labels_FS,
            atlas="FS")

        df_subject_intensity_freesurfer_cerebellum_destrieux, image_subject_intensity_freesurfer_cerebellum_destrieux = quant.PET_FDG_quantification(
            norm_image_freesurfer_cerebellum_dkt,
            aparc_destrieux_image,
            df_labels_FS,
            atlas="FS")

        # Mean Regional PET uptake values normalized to cerebellum
        values_normalization_cerebellum_dkt = df_subject_intensity_freesurfer_cerebellum_dkt["mean_PET"]
        values_normalization_cerebellum_destrieux = df_subject_intensity_freesurfer_cerebellum_destrieux["mean_PET"]


        ### Normalized to mean value
        norm_image_mean_freesurfer = norm.intensity_normalization(registred_PET_image,
                                                                  mask_only_brain=T1_only_brain_without_ventricles)
        output_normalization_mean_freesurfer = os.path.join(path_normalization_images,
                                                            'intensity_normalization_mean_freesurfer.nii.gz')
        sitk.WriteImage(norm_image_mean_freesurfer, output_normalization_mean_freesurfer)

        # Quantification Normalized to mean total uptake value
        df_subject_intensity_freesurfer_mean_dkt, image_subject_intensity_freesurfer_mean_dkt = quant.PET_FDG_quantification(
            norm_image_mean_freesurfer,
            aparc_dkt_image,
            df_labels_FS,
            atlas="FS")

        df_subject_intensity_freesurfer_mean_destrieux, image_subject_intensity_freesurfer_mean_destrieux = quant.PET_FDG_quantification(
            norm_image_mean_freesurfer,
            aparc_destrieux_image,
            df_labels_FS,
            atlas="FS")

        values_normalization_mean_dkt = df_subject_intensity_freesurfer_mean_dkt["mean_PET"]
        values_normalization_mean_destrieux = df_subject_intensity_freesurfer_mean_destrieux["mean_PET"]



        df_normalization_dkt = pd.DataFrame({'Structure': name_structures_dkt,
                                             'Regional uptake mean values': values_mean_activity_dkt,
                                             'Normalization to total brain mean value': values_normalization_mean_dkt,
                                             'Normalization to cerebellum uptake values': values_normalization_cerebellum_dkt,
                                             })

        df_normalization_destrieux = pd.DataFrame({'Structure': name_structures_destrieux,
                                                   'Regional uptake mean values': values_mean_activity_destrieux,
                                                   'Normalization to total brain mean value': values_normalization_mean_destrieux,
                                                   'Normalization to cerebellum uptake values': values_normalization_cerebellum_destrieux,
                                                   })

        df_normalization_dkt.to_csv(output_csv_freesurfer_dkt)
        df_normalization_destrieux.to_csv(output_csv_freesurfer_destrieux)



    # # Cerebellum values
    # cerebellum_R = float((df_subject_intensity.loc[(df_subject_intensity['structure'] == 'cerebellum')
    #                                                & (df_subject_intensity['hemisphere'] == 'R')]['mean_PET']).iloc[0])
    # cerebellum_L = float((df_subject_intensity.loc[(df_subject_intensity['structure'] == 'cerebellum')
    #                                                & (df_subject_intensity['hemisphere'] == 'L')]['mean_PET']).iloc[0])
    # cerebellum = (cerebellum_R + cerebellum_L) / 2
    #
    # brain_image = sitk.ReadImage(path_PET_final)
    #
    # norm_image = norm.intensity_normalization(brain_image, mode=("cerebellum", cerebellum))
    # sitk.WriteImage(norm_image, output_path + "/intensity_normalization_image.nii.gz")

    # Final Dir (Final files)
    path_final_dir = os.path.join(output_path, "Final")
    path_final_PET_normalized_image = os.path.join(path_final_dir, "spatial_normalized_pet_image.nii.gz")

    path_CSV_final_dir = os.path.join(path_final_dir, "CSV")
    path_Normalized_final_dir = os.path.join(path_final_dir, "Normalized")
    path_Registration_final_dir = os.path.join(path_final_dir, "Registration")

    if not os.path.exists(path_final_dir):
        os.mkdir(path_final_dir)

    # Copy PET and MRI if exists
    shutil.copy2(path_normalized_PET_image, path_final_PET_normalized_image)

    if path_MRI_image:
        path_final_MRI_normalized_image = os.path.join(path_final_dir, "spatial_normalized_mri_image.nii.gz")

        shutil.copy2(path_normalized_MRI_image, path_final_MRI_normalized_image)
        shutil.copytree(path_PET_registration_dir, path_Registration_final_dir, dirs_exist_ok=True)

    # Copy CSV and Normalized image
    shutil.copytree(path_CSV_files,path_CSV_final_dir, dirs_exist_ok=True)
    shutil.copytree(path_normalization_images, path_Normalized_final_dir, dirs_exist_ok=True)


    # Basic Data Analysis (Hipometabolism maps and Bars Charts)
    path_plots = os.path.join(output_path, "Plots")

    if not os.path.exists(path_plots):
        os.mkdir(path_plots)

    # Hypometabolism Maps
    # output_hypometabolism_map_cerebellum = os.path.join(path_plots, "Hypometabolism_maps_cerebellum.png")
    # output_hypometabolism_map_mean = os.path.join(path_plots, "Hypometabolism_maps_mean.png")
    #
    # plots.plot_hypometabolism_maps(path_synthetic_image_norm_cerebellum, path_normalized_MRI_image, output_path, output_hypometabolism_map_cerebellum, coord=(-3, -40, 15))
    # plots.plot_hypometabolism_maps(path_synthetic_image_norm_mean, path_normalized_MRI_image, output_path, output_hypometabolism_map_mean,
    #                                coord=(-3, -40, 15))
    #
    # # Bar charts
    # output_plot_cerebellum_top_ten = os.path.join(path_plots, f"{subject}_top_ten_cerebellum.png")
    # output_plot_mean_top_ten = os.path.join(path_plots, f"{subject}_top_ten_mean.png")
    #
    # plots.plot_top_ten_regions(df_changes, df_normalization_subject, df_normalization_MNI152,
    #                            output_plot_cerebellum_top_ten, mode="cerebellum")
    # plots.plot_top_ten_regions(df_changes, df_normalization_subject, df_normalization_MNI152,
    #                            output_plot_mean_top_ten)
