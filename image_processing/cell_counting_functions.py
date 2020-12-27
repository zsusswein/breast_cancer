#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:45:54 2019

@author: zsusswein
"""

## This file houses the cell counting functions to import to other files


###########################
# Libraries

from skimage import filters, io, morphology, feature, exposure
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#############################

def cell_count(image, show_plot):
# this function takes an rgb channel of a cell image and returns the cell 
# count based on the light expressed from that channel and an image with each
# cell found circled in red. The function equalizes and thresholds the image 
# twice (to preserve the dimmer cells) before applying a power law filter to
# increase contrast between cells before counting.
    
                
    #########################
    ## FIRST EQUALIZE/THRESHOLD -- more gentle
    
    # adjust exposure using sigmoid correction
    rescale_low = exposure.adjust_sigmoid(image = image,
                                          cutoff=.2, gain=15, inv=False)

    # threshold the image with local adaptive filter (because washout on edges)
    thresh_init = filters.threshold_local(rescale_low, 201)
    thresh_new = rescale_low >= thresh_init
    
    # remove remaining crap from image
    cleaned_init = morphology.remove_small_objects(thresh_new, 20)

    # replace  pixels w/ a T in the thresholded image with the actual intensity 
    filtered_low = rescale_low

    for x in range(len(cleaned_init)):
        for y in range(len(cleaned_init[x])):
            if not cleaned_init[x,y]:
                filtered_low[x,y] = 0
                
    ########################
    ## SECOND EQUALIZE/THRESHOLD -- less gentle

    # use an adaptive method to locally contrast the image
    equalized_hist = exposure.equalize_adapthist(filtered_low, clip_limit=.02)
        
    # adjust exposure using sigmoid correction
    rescale_high = exposure.adjust_sigmoid(image = equalized_hist, 
                                           cutoff=.3, gain=10, inv=False)    
    
    # threshold the image with an local adaptive filter
    thresh = filters.threshold_local(rescale_high, 501)
    thresh_image = rescale_high >= thresh
        
    # clean away the junk
    cleaned = morphology.remove_small_objects(thresh_image, 150)
    
    # replace  pixels w/ a T in the thresholded image with the actual intensity 
    filtered_max = rescale_high
    
    for x in range(len(cleaned)):
        for y in range(len(cleaned[x])):
            if not cleaned[x,y]:
                filtered_max[x,y] = 0
                
    # power law correction to increase contrast
    gamma_corrected = exposure.adjust_gamma(filtered_max, gamma=5, gain=10)

    # find all the cells
    blobs = feature.blob_log(gamma_corrected, min_sigma=5, max_sigma=15, 
                             overlap=0.8, threshold=.5, num_sigma=30)
 
    # plot all the cells found as red circles on the original image
    if show_plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.imshow(image)
        for blob in blobs:
            y, x, r = blob
            c = plt.Circle((x, y), r, color='red', linewidth=1, fill=False)
            ax.add_patch(c)
        plt.tight_layout()
        plt.show()
    
    num_cells = len(blobs)

    ##########################
    # return the total count
    return num_cells

########################
    
def clean_channel(image):
# this function takes the full image as a numpy array the rgb channels and 
# returns (cleaned) single red and green channels. The channels check if the
# or green channel has a higher intensity for each pixel and zeroes out the one
# with lower intensity. This procedure removes the background red and green
# light from coculture images
    
    # channel image
    clean_red = image[:,:,0]
    clean_green = image[:,:,1]
    
    # iterate through each pixel and zero lower value
    for x in range(len(image)):
        for y in range(len(image[x])):
            if image[x,y,0] > image[x,y,1]: # red > green
                clean_green[x,y] = 0
            elif image[x,y,0] < image[x,y,1]: # green < red
                clean_red[x,y] = 0
            else:                     # ties have both zeroed out
                clean_red[x,y] = 0
                clean_green[x,y] = 0
    
    # threshold the red image with an ISODATA filter to remove noise
    clean_red = exposure.adjust_gamma(clean_red, gamma=1)
    clean_green = exposure.adjust_gamma(clean_green, gamma=1)
                
     ###################### 
     
    # return cleaned and channeled red and green intensities
    return clean_red, clean_green
  
#########################

def total_cell_count(image_path, specify_image = False, specified_image = False, show_plot = True):
# this function takes an image path and returns the counts of red and green
# cells in an image.

    # convert the path to the image
    if not specify_image:
        image = io.imread(image_path) # image in MxNx3 where the 3rd dim is RBG
    else:
        image = specified_image
    
    clean_red, clean_green = clean_channel(image) # clean coculture image
    
    num_red = cell_count(clean_red, show_plot)        # count the number of cells of that color
    num_green = cell_count(clean_green, show_plot)
    
    ##################
    return num_red, num_green

###############
