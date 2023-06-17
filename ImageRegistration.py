#! python3
from __future__ import print_function

import SimpleITK as sitk
import numpy as np
import sys
import os

# Rigid registration using mutual information for multi modal. Mainly used to register a PET image to the MRI image of the same subject:
def RigidImageRegistration(movingImage, fixedImage, printLog = False, fixedMask = []):

    regMethod = sitk.ImageRegistrationMethod()

    regMethod.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50) # Mutual information for multi modality image registration:
    # Optimization:
    regMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=1.0,
                                            numberOfIterations=200,
                                            convergenceMinimumValue=1e-6)
    regMethod.SetMetricSamplingStrategy(regMethod.RANDOM)
    regMethod.SetOptimizerScalesFromPhysicalShift()
    regMethod.SetMetricSamplingPercentage(0.2)
    regMethod.SetShrinkFactorsPerLevel([4, 2, 1])
    regMethod.SetSmoothingSigmasPerLevel([2,1,0])

    # Initialization for Rigid registration:
    tx = sitk.CenteredTransformInitializer(fixedImage, movingImage,
                                           sitk.Euler3DTransform(),
                                           sitk.CenteredTransformInitializerFilter.GEOMETRY)

    regMethod.SetInitialTransform(tx)

    # Linear interpolator:
    regMethod.SetInterpolator(sitk.sitkLinear)

    # If mask, set it:
    if fixedMask:
        regMethod.SetMetricFixedMask(fixedMask)

    # Registration
    outTx = regMethod.Execute(fixedImage, movingImage)

    if printLog:
        print("-------")
        print(outTx)
        print(f"Optimizer stop condition: {regMethod.GetOptimizerStopConditionDescription()}")
        print(f" Iteration: {regMethod.GetOptimizerIteration()}")
        print(f" Metric value: {regMethod.GetMetricValue()}")

    # Apply transform:
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fixedImage)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(0)
    resampler.SetTransform(outTx)

    regImage = resampler.Execute(movingImage)

    regResults = {'image': regImage, 'tx': outTx}

    return regResults

# Affine registration using cross correlation. Mainly used as the first step to normalize MRI into MNI152
def AffineImageRegistration(movingImage, fixedImage, printLog = False, fixedMask = []):

    regMethod = sitk.ImageRegistrationMethod()

    regMethod.SetMetricAsMattesMutualInformation(50)
    # Optimization:
    regMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=1.0,
                                              numberOfIterations=200,
                                            convergenceMinimumValue=1e-6)
    regMethod.SetMetricSamplingStrategy(regMethod.RANDOM)
    regMethod.SetOptimizerScalesFromPhysicalShift()
    regMethod.SetMetricSamplingPercentage(0.2)
    regMethod.SetShrinkFactorsPerLevel([4, 2, 1])
    regMethod.SetSmoothingSigmasPerLevel([2,1,0])

    # Initialization for Rigid registration:
    tx = sitk.CenteredTransformInitializer(fixedImage, movingImage,
                                           sitk.AffineTransform(3),
                                           sitk.CenteredTransformInitializerFilter.GEOMETRY)

    regMethod.SetInitialTransform(tx)

    # Linear interpolator:
    regMethod.SetInterpolator(sitk.sitkLinear)

    # If mask, set it:
    if fixedMask:
        regMethod.SetMetricFixedMask(fixedMask)

    # Registration
    outTx = regMethod.Execute(fixedImage, movingImage)

    if printLog:
        print("-------")
        print(outTx)
        print(f"Optimizer stop condition: {regMethod.GetOptimizerStopConditionDescription()}")
        print(f" Iteration: {regMethod.GetOptimizerIteration()}")
        print(f" Metric value: {regMethod.GetMetricValue()}")

    # Apply transform:
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fixedImage)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(0)
    resampler.SetTransform(outTx)

    regImage = resampler.Execute(movingImage)

    regResults = {'image': regImage, 'tx': outTx}

    return regResults

# Nonlinear image registration using cross-correlation. Mainly to normalize MRI into NMI152.
def NonlinearImageRegistration(movingImage, fixedImage, printLog = False):

    regMethod = sitk.ImageRegistrationMethod()

    regMethod.SetMetricAsCorrelation()
    # Optimization:
    regMethod.SetOptimizerAsGradientDescentLineSearch(learningRate=1.0,
                                              numberOfIterations=200,
                                            convergenceMinimumValue=1e-6)
    regMethod.SetMetricSamplingStrategy(regMethod.RANDOM)
    regMethod.SetOptimizerScalesFromPhysicalShift()
    regMethod.SetMetricSamplingPercentage(0.1)
    regMethod.SetShrinkFactorsPerLevel([4, 2, 1])
    regMethod.SetSmoothingSigmasPerLevel([2,1,0])

    # Initialization for bspline image registration:
    transformDomainMeshSize = [12, 12, 6]
    tx = sitk.BSplineTransformInitializer(fixedImage, [12, 12, 6])

    regMethod.SetInitialTransform(tx,
                               inPlace=True)

    # Linear interpolator:
    regMethod.SetInterpolator(sitk.sitkLinear)
    # Registration
    outTx = regMethod.Execute(fixedImage, movingImage)

    if printLog:
        print("-------")
        print(outTx)
        print(f"Optimizer stop condition: {regMethod.GetOptimizerStopConditionDescription()}")
        print(f" Iteration: {regMethod.GetOptimizerIteration()}")
        print(f" Metric value: {regMethod.GetMetricValue()}")

    # Apply transform:
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fixedImage)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetDefaultPixelValue(0)
    resampler.SetTransform(outTx)

    regImage = resampler.Execute(movingImage)

    regResults = {'image': regImage, 'tx': outTx}

    return regResults

def ApplyRegTransform(inputImage, tx, refImage = [], interpolator = sitk.sitkLinear):
    # If not reference, stay on the input image space
    if not refImage:
        refImage = inputImage
    # Apply transform:
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(refImage)
    resampler.SetInterpolator(interpolator)
    resampler.SetDefaultPixelValue(0)
    resampler.SetTransform(tx)

    txImage = resampler.Execute(inputImage)
    return txImage