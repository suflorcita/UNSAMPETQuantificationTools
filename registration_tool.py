import os
import subprocess
import SimpleITK as sitk
import math


# -------------------
# Generic helpers
# -------------------

def gaussian_smooth(input_image, output_dir, prefix="Image", fwhm=6.0, name=None):
    """
    Apply Gaussian smoothing to an image and save.
    """
    sigma = fwhm / (2 * math.sqrt(2 * math.log(2)))
    gaussian_filter = sitk.SmoothingRecursiveGaussianImageFilter()
    gaussian_filter.SetSigma(sigma)

    img = sitk.ReadImage(input_image)

    # --- Extracción de la Información Espacial Original ---
    original_origin = img.GetOrigin()
    original_direction = img.GetDirection()

    smoothed = gaussian_filter.Execute(img)

    smoothed.SetOrigin(original_origin)
    smoothed.SetDirection(original_direction)


    base_name = name if name else prefix
    out_path = os.path.join(output_dir, f"{base_name}.nii.gz")

    smoothed_final = sitk.Cast(smoothed, sitk.sitkFloat32)

    sitk.WriteImage(smoothed, out_path)
    return out_path


def run_flirt(input_image, ref_image, output_dir, prefix="Image", name=None):
    """
    Run FLIRT between input_image and ref_image.
    """
    base_name = name if name else prefix
    path_registered = os.path.join(output_dir, f"{base_name}.nii.gz")
    path_mat = os.path.join(output_dir, f"{base_name}_TX.mat")

    flirt_cmd = (
        f"flirt -in {input_image} -ref {ref_image} "
        f"-out {path_registered} -omat {path_mat} "
        f"-interp trilinear -dof 12"
    )
    subprocess.run([flirt_cmd], shell=True, check=True)
    return path_registered, path_mat


def apply_flirt(input_image, ref_image, mat_file, output_dir, prefix="Image", name=None):
    """
    Apply a FLIRT transform to input_image.
    """
    base_name = name if name else prefix
    out_path = os.path.join(output_dir, f"{base_name}.nii.gz")

    apply_cmd = (
        f"flirt -in {input_image} -applyxfm -init {mat_file} "
        f"-out {out_path} -paddingsize 0.0 -interp trilinear -ref {ref_image}"
    )
    subprocess.run([apply_cmd], shell=True, check=True)
    return out_path


def run_ants(input_image, ref_image, output_dir, prefix="Image", name=None, path_ANT="antsRegistrationSyN.sh"):
    """
    Run ANTs registration (SyN).
    """
    base_name = name if name else prefix
    ants_prefix = os.path.join(output_dir, base_name)
    out_warped = ants_prefix + "Warped.nii.gz"

    ants_cmd = f"{path_ANT} -d 3 -f {ref_image} -m {input_image} -o {ants_prefix}"
    subprocess.run([ants_cmd], shell=True, check=True)

    return out_warped, ants_prefix + "1Warp.nii.gz", ants_prefix + "0GenericAffine.mat"


def apply_ants(input_image, ref_image, tx1, tx2, output_dir, prefix="Image", name=None, path_ANT_apply="antsApplyTransforms"):
    """
    Apply ANTs transforms to input_image.
    """
    base_name = name if name else prefix
    out_path = os.path.join(output_dir, f"{base_name}.nii.gz")

    apply_cmd = (
        f"{path_ANT_apply} -d 3 -i {input_image} -r {ref_image} "
        f"-o {out_path} -t {tx1} -t {tx2}"
    )
    subprocess.run([apply_cmd], shell=True, check=True)
    return out_path


import os
import glob


def ensure_transform_files(output_path, transform_files):
    """
    Ensure transform_files contains FLIRT or ANTs transforms.
    If missing, check the output_path for existing files.
    """
    # --- Check FLIRT ---
    if "flirt_mat" not in transform_files:
        flirt_mat = os.path.join(output_path, "PET_image_to_template_flirt.mat")
        if os.path.exists(flirt_mat):
            transform_files["flirt_mat"] = flirt_mat

    # --- Check ANTs ---
    if "ants_tx" not in transform_files:
        tx1 = os.path.join(output_path, "PET_image_to_template_ANT1Warp.nii.gz")
        tx2 = os.path.join(output_path, "PET_image_to_template_ANT0GenericAffine.mat")
        if os.path.exists(tx1) and os.path.exists(tx2):
            transform_files["ants_tx"] = (tx1, tx2)

    return transform_files


def ensure_normalized_image(output_path, path_MRI_image=None, freesurfer=False, quantify_using_flirt=False):
    """
    Ensure that a normalized PET (and optionally MRI, aseg) image exists.
    Prioritize ANTs, fallback to FLIRT. Raises if none exist.

    Args:
        output_path (str): directory where normalized images should exist.
        path_MRI_image (str|None): if MRI exists, check both PET and MRI normalization.
        freesurfer (bool): if aseg segmentation should be included.
        quantify_using_flirt (bool): force using FLIRT PET if available.

    Returns:
        dict: with keys { 'PET': path, 'MRI': path|None, 'aseg': path|None }
    """

    results = {"PET": None, "MRI": None, "aseg": None}

    if path_MRI_image:  # MRI case
        path_flirt = {
            "PET": os.path.join(output_path, "FLIRT", "PET_Norm_MNI_152_FLIRT.nii.gz"),
            "MRI": os.path.join(output_path, "FLIRT", "T1_brain_MNI_152_FLIRT.nii.gz"),
            "aseg": os.path.join(output_path, "FLIRT", "aseg_Norm_MNI_152_FLIRT.nii.gz"),
        }
        path_ants = {
            "PET": os.path.join(output_path, "ANTs", "PET_Norm_MNI_152_ANT.nii.gz"),
            "MRI": os.path.join(output_path, "ANTs", "T1_Norm_MNI_152_ANTWarped.nii.gz"),
            "aseg": os.path.join(output_path, "ANTs", "Aseg_Norm_MNI_152_ANT.nii.gz"),
        }

        if os.path.exists(path_ants["PET"]):
            results["PET"] = path_ants["PET"]
            results["MRI"] = path_ants["MRI"]
            if freesurfer:
                results["aseg"] = path_ants["aseg"]

        elif os.path.exists(path_flirt["PET"]):
            results["PET"] = path_flirt["PET"]
            results["MRI"] = path_flirt["MRI"]
            if freesurfer:
                results["aseg"] = path_flirt["aseg"]

    else:  # PET-only case
        path_flirt = os.path.join(output_path, "PET_image_to_template_flirt.nii.gz")
        path_ants = os.path.join(output_path, "PET_image_to_template_ANT.nii.gz")

        if quantify_using_flirt and os.path.exists(path_flirt):
            results["PET"] = path_flirt

        elif os.path.exists(path_ants):
            results["PET"] = path_ants

        elif os.path.exists(path_flirt):
            results["PET"] = path_flirt

    # --- Sanity check ---
    if not results["PET"]:
        raise FileNotFoundError("No Normalized PET Image found in output_path.")

    return results


# -------------------
# PET -> Template pipeline
# -------------------
def pet_registration_pipeline(path_PET_image, path_PET_template, output_path, args,
                              path_ANT="antsRegistrationSyN.sh", path_ANT_apply="antsApplyTransforms"):
    """
    Run PET -> Template registration using FLIRT and/or ANTs.
    Optionally skips Gaussian smoothing if --no-smoothed is used.
    Ensures transforms are applied only when relevant.
    """

    # 1. Selección de imagen base
    if getattr(args, "smoothed", True):
        path_input_for_registration = gaussian_smooth(
            input_image=path_PET_image,
            output_dir=output_path,
            prefix="Smoothed_PET_image"
        )
        input_label = "smoothed"
    else:
        path_input_for_registration = path_PET_image
        input_label = "image"

    path_PET_final = path_input_for_registration
    path_normalized_PET_image = None
    transform_files = {}

    # 2. FLIRT registration
    if args.flirt:
        path_PET_to_Template_Flirt, path_FLIRT_tx = run_flirt(
            input_image=path_input_for_registration,
            ref_image=path_PET_template,
            output_dir=output_path,
            prefix=f"PET_{input_label}_to_template_flirt"
        )

        # Solo aplicar la matriz si hubo smoothing
        if args.smoothed:
            path_PET_to_Template_Flirt_applied = apply_flirt(
                input_image=path_PET_image,
                ref_image=path_PET_template,
                mat_file=path_FLIRT_tx,
                output_dir=output_path,
                prefix=f"PET_image_to_template_flirt"
            )
            path_PET_final = path_PET_to_Template_Flirt_applied
            path_normalized_PET_image = path_PET_to_Template_Flirt_applied
        else:
            # No need to apply transform again if registration used raw directly
            path_PET_final = path_PET_to_Template_Flirt
            path_normalized_PET_image = path_PET_to_Template_Flirt

        transform_files["flirt_mat"] = path_FLIRT_tx

    # 3. ANTs registration
    if args.ants:
        ants_prefix = f"PET_{input_label}_to_template_ANT"

        warped, tx1, tx2 = run_ants(
            input_image=path_input_for_registration,
            ref_image=path_PET_template,
            output_dir=output_path,
            prefix=ants_prefix,
            path_ANT=path_ANT
        )

        # Si hubo smoothing, aplicar la transformada al raw
        if args.smoothed:
            path_PET_to_Template_ANT = apply_ants(
                input_image=path_PET_image,
                ref_image=path_PET_template,
                tx1=tx1,
                tx2=tx2,
                output_dir=output_path,
                prefix="PET_image_to_template_ANT",
                path_ANT_apply=path_ANT_apply
            )
            path_PET_final = path_PET_to_Template_ANT
        else:
            path_PET_final = warped

        path_normalized_PET_image = path_PET_final
        transform_files["ants_tx"] = (tx1, tx2)

    return path_PET_final, path_normalized_PET_image, transform_files



def normalize_pet_frames_to_template(PET_frames, path_PET_template, transform_files, output_dir):
    """
    Apply a PET→Template transform (FLIRT or ANTs) to each PET frame and save the results.

    Args:
        PET_frames (list[str]): List of paths to individual PET frames (e.g., PET_frame_0.nii.gz).
        path_PET_template (str): Path to the PET template image.
        transform_files (dict): Transformations returned by pet_registration_pipeline.
                                Example: {"flirt_mat": path} or {"ants_tx": (tx1, tx2)}.
        output_dir (str): Directory where normalized frames will be saved.

    Returns:
        normalized_frames (list[str]): List of paths to the normalized PET frames.
    """

    os.makedirs(output_dir, exist_ok=True)
    normalized_frames = []

    for i, frame_path in enumerate(PET_frames):
        # --- Case: FLIRT transform ---
        if "flirt_mat" in transform_files:
            out_path = apply_flirt(
                input_image=frame_path,
                ref_image=path_PET_template,
                mat_file=transform_files["flirt_mat"],
                output_dir=output_dir,
                prefix=f"PET_frame_{i}"
            )
            normalized_frames.append(out_path)

        # --- Case: ANTs transform ---
        if "ants_tx" in transform_files:
            tx1, tx2 = transform_files["ants_tx"]
            out_path = apply_ants(
                input_image=frame_path,
                ref_image=path_PET_template,
                tx1=tx1,
                tx2=tx2,
                output_dir=output_dir,
                prefix=f"PET_frame_{i}"
            )
            normalized_frames.append(out_path)

    return normalized_frames
