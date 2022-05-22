# Generating cell counts from images

This directory holds the code used in (Susswein et al., 
2021)[https://www.biorxiv.org/content/10.1101/2022.02.18.481041v1] to process 
images of fluorescent cells and extract counts of each population. This 
pipeline is optimized for cells fluorescing red and green, particularly with 
the mCherry and egfp. This pipeline has been validated against these 
particular fluorescent markers and there are no guarantees it will perform as 
well in other conditions (see the validation figure or the code to produce it 
in the repository).

## Suggestions for use

If one wishes to use this image processing pipeline -- understanding that it 
is not validated outside of the experimental setting described in the 
manuscript -- they are welcome and encouraged to do so under the MIT license.

The functions to perform the actual image processing are located in the file 
`cell_counting_functions.py`. The actual workhorse functions are meant to be 
called from a wrapper function `total_cell_count()`, but may be accessed 
directly if one wishes. The `total_cell_count()` function expects the image to 
have red-labeled cells, green-labeled cells, or both. It takes as its first 
argument the path to an image as well as some optional keyword arguments and 
returns two integers: the number of observed red cells and the number of 
observed green cells. 

In this function, the optional keyword argument `show_plot` is passed to the 
workhorse `cell_count()` function, where, if true, the function plots the 
image with each identified cell labeled with a red circle. This keyword 
argument is by default true. It serves as a useful diagnostic that the 
algorithm is correctly identifying the cells in the image.

## Methods

This pipeline has two main steps: cleaning and denoising the image and 
counting the number of cells in the given RGB color channel.

### Image Denoising

The denoising algorithm operates under the assumption that each pixel can only 
be occupied by a red cell or a green cell, not both. It iterates through each 
pixel in the image, comparing the pixel intensity in the red and green 
channels, and zeroing out the channel with the lower pixel intensity. The 
algorithm then linearly scales the pixels to be in [0, 1] using the 
`scikit-image` `adjust_gamma()` function as a convenience function by passing 
the image to it with `gamma` set to 1. 


# Cell counting

The number of cells in a given color channel are counted by the `cell_count()` 
function, with the red and green color channels counted separately (the blue 
channel is ignored). The cells in the channel are counted as part of a three 
step process.

First, the channel exposure is adjusted with a sigmoid correction, the channel 
is thresholded with a local adaptive filter, and the remaining small objects 
in the image are removed. We use a local adaptive filter because of washout 
along the edges of the image leading to systematic changes in mean intensity 
that reduce the effectiveness of non-local methods. Then, we replace the 
remaining pixels in the channel above the threshold with their original 
intensity.

Second, we repeat the process from the first step, but with higher thresholds.  
This process reduces the remaining noise in the image, but is better able to 
discern cells from residual noise because of the first, more gentle,  cleaning 
step. We finish this step by increasing the contrast with a power law 
correction with gamma set to 5.

Third, and finally, we actually count the cells in the image using using the 
(`skimage.feature.blob_log()` 
function)[https://scikit-image.org/docs/stable/auto_examples/features_detection/plot_blob.html].  
This function computes the Laplacian of the Gaussian objects by successively 
computing the standard deviations of the identified objects until the radius 
of the object is obtained. This function returns the position (x and 
y coordinates) and the radius of the detected object. The `cell_count()` 
function takes this object and returns the number of objects detected -- the 
number of cells in the channel. This is the number that is returned for the 
number of red or green cells detected by the algorithm.
