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
  matplotlib
  
## Command Line Usage 

To run this project, clone the repository to your local machine and navigate to the project directory in the terminal. Then, use the following command to execute the main script:

```bash
python3 main.py path_PET path_MRI subject_name output_folder
```

For example: 





