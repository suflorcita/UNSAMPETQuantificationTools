# UNSAM PET Quantification Tools
Pipeline for cuantification of $[^{18}F] FDG$ PET Images 

Positron Emission Tomography (PET) and Magnetic Resonance Imaging (MRI) are important tools for the study of neurodegenerative diseases. Existing tools often perform individual steps of preprocessing, registration and segmentation separately. We propose an integrated approach that combines these steps into a unified framework for quantification of brain $[^{18}F] FDG$ PET images, using a python script. 

## Input images 
The proposed pipeline receives dynamic or static PET images and a high-resolution anatomical MRI image .

## Output 
This pipeline generate a series of spatial and intensity normalized images and a repport of quantitative metrics of regional $[^{18}F] FDG$ uptake and hypometabolism maps. 

## Systems requirements

The pipeline is developed for use on a Linux system. In order to use this tool, you need to download FreeSurfer (https://surfer.nmr.mgh.harvard.edu/), FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) and ANTs (http://stnava.github.io/ANTs/). Install this tools according to the instructions given by the developers.

## Python libraries requirements






