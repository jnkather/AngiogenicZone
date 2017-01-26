# AngiogenicZone
Distribution analysis of CD34-stained blood vessels in histopathological whole slide images

##What is this?

This is a general outline of the blood vessel analysis presented in this repository. It is described in more detail in the following paper:

"Identification of a characteristic vascular belt zone in human colorectal cancer"

For questions, please contact Dr. Jakob N. Kather of Heidelberg University, Germany (http://orcid.org/0000-0002-3730-5348). Please note that this is experimental software not designed for routine work. When you use it, you have to manually change a couple of variables in the code each time you run the code. Please also observe the disclaimer and the license. Please cite our paper if you re-use parts of the method.

##General usage

Preparation:

First, download external dependencies: 
- Yen.m and CumMeanVar.m from [1]. Place these files in the sub-folder "Yen"
- inpoly.m from [2]. Place this file in the sub-folder "poly_stuff".

[1] http://dismac.dii.unipg.it/paper/code.html from: Francesco Bianconi et al.: "A sequential machine vision procedure for assessing paper impurities" doi: 10.1016/j.compind.2013.12.001

[2] https://www.mathworks.com/matlabcentral/fileexchange/10391-fast-points-in-polygon-test/content/poly_stuff/inpoly.m 


The general usage is as follows:

Use "main_blockproc_2016.m" to analyze the distribution of blood vessels in histological whole slide images. These images (tiled TIFF, Scans of CD34 stained tissue slides) represent two batches.

- "01_CD34_full_cohort_TIFF_images/First Cohort" for the First cohort
- "01_CD34_full_cohort_TIFF_images/Validation Cohort" for the Validation cohort

Please note that the actual image files are not part of this repository. Please refer to the paper for a DOI to these images.

Before running, two folders have to be defined in the file "getConstants.m":

- the "cache folder" that is used to store all intermediate results (including segmentation masks etc.). Default for this folder is '../cache_multivariate/'
- the "results folder" that is used to store all final results (including heat maps for each image, thumbnails etc.). Default for this folder ist '../results_multivariate/'

At the end of the whole analysis, the final CSV table will be stored in the root folder.

##Manual settings

Because often it is desirable not to analyze all images at once, the image filenames have to be specified in the variable "imgArray" at the beginning of "main_blockproc_2016.m". The following actions are successively applied to each image:

1. Particle detection: Read the image using "blockproc()" and perform an automatic color deconvolution and thresholding. The output is a binary mask that is stored in the cache folder. Please note that this process can cause artifacts if the staining is too strong or too weak in some parts of the image. Generally, this is not a problem if the artifacts are not located within the tumor ROI.
2. Request ROI: The user is asked to manually draw a polygonal ROI in each image. The resulting coordinates are stored in the cache folder. All ROI coordinates for the images in our data set are already stored in the cache folder. For all cases with staining artifacts, the ROIs are drawn in such a way that the artifacts are excluded. Please note that the names of the ROIs have to be exactly "FullTumorValidated", "LumenValidated" and "AdjacentTissueValidated".
3. Analysis in ROI: All blood vessels within the ROI are identified and the hot spot analysis is performed (see our previous paper: http://dx.doi.org/10.18632/oncotarget.4383)
4. Create a thumbnail: A thumbnail is saved as an image file. This thumbnail shows a reduced version of the original whole slide image with all ROIs as an overlay. (saved in the cache folder)
5. Analysis of distances: The distance of all blood vessels within the ROI "FullTumorValidated" to the other ROIs is measured. Then, the distribution of distances is analyzed as described in our paper.

If you run this whole process, you will receive intermediate results in the cache folder and final heat maps in the results folder. The last step is to merge all these results into a large CSV table. This can be done by running the script "main_measureHotspotArea.m". This script reads all intermediate measurements from the cache folder and results folder and merges them in an output file. The name of the table can be specified within this file, by default it is "QuantificationHotspotTable.csv".

##Performance

For one image, this procedure takes between 30 and 60 minutes on a standard computer workstation. Often, it is desirable to perform one of these tasks only. Therefore, it can be specified which of these five steps should be performed. To change this, you have to assign 0 or 1 to the corresponding variables ("do") in "main_blockproc_2016.m").

##Random numbers

Please note that this program uses random numbers at different points. Generally, the observed distribution of blood vessels is compared to random distributions several times. For the main part of the program (detection of a vascular belt zone / angiogenic zone), this is done reproducibly by a) seeding the random number generator before each run and b) performing N=100 random experiments and using the average. However, these precautions might not be present in all other parts of the program so please be aware that repeated experimental runs can yield slightly different results. 
 
