Code for object detection in astronomical images

Copyright(c) 2015, Caixia Zheng, Jesus Pulido, Paul Thorman, Bernd Hamann
All Rights Reserved.

-----------------------------------------------------------------------------------------------
This code is being made available for research purposes only, it shall not be used, rewritten,or adapted as the basis of a commercial software or hardware product. The authors make no representations about the suitability of this code for any purpose. It is provided "as is" without express or implied warranty.
------------------------------------------------------------------------------------------------

1、OVERVIEW

* Data access:
Download the data from "http://lsst.astro.washington.edu/data/deep/" for testing the global detetcion pipeline,
and download the data from "http://dls.physics.ucdavis.edu" for testing the local detecion pipline.

* Requirements: Matlab and the Image Processing Toolbox
                or
                GNU Octave with Image and FITS packages

* To execute the local method for astronomical object detection, run the script Main_Local_Detection.m.
* To execute the global method for astronomical object detection, run the script Main_Global_Detection.m.

BGsubtract.m: Function to estimate the background of an astronomical image.

imgnew.m: Function to clip the image intensity for computing the background.

global_detect.m: Function to detect all pixels set larger than a threshold in whole image as the objects.

select_seedpoints.m: Function to select the large objects in the global detection results as the seed points for the watershed partition.

local_detect.m: Function to divide whole image to sub-regions and detect objects in each region.

deblend_Multi.m: Function to deblending the touching or overlapping objects by using the multi-threshold method and the watershed method.

label_objects.m: Function to label each individual object by different value.

compt_statistics.m: Function to compute some statistics of detected objects.

NClean.m: Function to clean the artifacts (noise) around large objects for global detection results.

deblend_Wa.m: Function to deblend the touching or overlapping objects directly using the watershed method

* For noise level estimation, we use Masayuki Tanaka’s matlab codes “NoiseLevel.m”

* For ellipse fitting, we use the matlab codes “fitellipse.m” provided by Andrew Fitzgibbon, Maurizio Pilu and Bob Fisher.

* For computing the distances from the points to the ellipse, we use Hui Ma’s matlab codes “Residuals_ellipse.m”.

2、REFERENCE
This is an example code for the algorithm described in the paper:
 Caixia Zheng, Jesus Pulido, Paul Thorman, Bernd Hamann, “An Improved Method for Object Detection in Astronomical Images”.
