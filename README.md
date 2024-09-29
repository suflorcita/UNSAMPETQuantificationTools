# UNSAM PET Quantification Tools
Pipeline for cuantification of $[^{18}F] FDG$ PET Images 

Positron Emission Tomography (PET) and Magnetic Resonance Imaging (MRI) are important tools for the study of neurodegenerative diseases. Existing tools often perform individual steps of preprocessing, registration and segmentation separately. We propose an integrated approach that combines these steps into a unified framework for quantification of brain $[^{18}F] FDG$ PET images, using a Python script. 

## Input images 
The proposed pipeline receives dynamic or static PET images and a high-resolution anatomical MRI image .

## Output 
This pipeline generates a series of spatially and intensity-normalized images and a report containing quantitative metrics of regional uptake and hypometabolism maps.

## Systems requirements

The pipeline is developed for use on a Linux system. To use this tool, you will need to download and install the following software packages:

- FreeSurfer (https://surfer.nmr.mgh.harvard.edu/)
- FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
- ANTs (http://stnava.github.io/ANTs/)

Please follow the installation instructions provided by the developers for these tools.

## Python libraries requirements

This project relies on several Python libraries for its operation. Make sure you have these libraries installed in your Python environment:

- SimpleITK
- pandas
- subprocess
- os
- shutil 
- matplotlib
  
## Command Line Usage 

To run this project, clone the repository to your local machine and navigate to the project directory in the terminal. Then, use the following command to execute the main script:

Basic Command Structure

```bash
PETQuantification [-h] [-m MRI] [-p PET] [-o OUTPUT] [-s SUBJECT]
                  [-t PET_TEMPLATE] [-f FREESURFER_DIR] [--no-flirt]
                  [--no-ant] [--no-freesurfer]
```
### **Options:**

- `-h, --help`  
  Show the help message and exit.
  
- `-m MRI, --mri MRI`  
  Path to the MRI image. Use this flag if you have an MRI image.

- `-p PET, --pet PET`  
  Path to the PET image (required).

- `-o OUTPUT, --output OUTPUT`  
  Name of the output directory where processed data will be saved.

- `-s SUBJECT, --subject SUBJECT`  
  Name of the subject being processed.

- `-t PET_TEMPLATE, --pet-template PET_TEMPLATE`  
  PET template to be used when no MRI image is provided.

- `-f FREESURFER_DIR, --freesurfer-dir FREESURFER_DIR`  
  Path to the FreeSurfer directory (required for FreeSurfer-based processing).

- `--no-flirt`  
  Disable the FLIRT command. Use this flag to skip the FLIRT registration process.

- `--no-ant`  
  Disable the ANTs command. Use this flag to skip the ANTs registration process.

- `--no-freesurfer`  
  Disable the FreeSurfer command. Use this flag to skip FreeSurfer-based
