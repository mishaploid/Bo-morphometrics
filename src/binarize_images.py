#!/Users/DontPanic/anaconda3/envs/plantcv/bin/python

# Brassica oleracea leaf segmentation
# adapted mostly from PlantCV Multi Plant Workflow:
# https://plantcv.readthedocs.io/en/stable/multi-plant_tutorial/

import sys, traceback
import cv2
import os
import re
import numpy as np
import argparse
import string
from plantcv import plantcv as pcv

### Parse command-line arguments
def options():
    parser = argparse.ArgumentParser(description="Imaging processing with opencv")
    parser.add_argument("-i", "--image", help="Input image file.", required=True)
    parser.add_argument("-o", "--outdir", help="Output directory for image files.", required=True)
    parser.add_argument("-n", "--names", help="path to txt file with names of genotypes to split images into", required=False)
    parser.add_argument("-D", "--debug", help="Turn on debug, prints intermediate images.", action=None)
    parser.add_argument("-r", "--result", help="Specify a json file to store results", required=False)
    args = parser.parse_args()
    return args

### Main workflow
def main():
    # Get options
    args = options()

    # Read image
    img, path, filename = pcv.readimage(args.image)

    pcv.params.debug=args.debug #set debug mode

    # STEP 1: white balance (for comparison across images)
    # inputs:
    #   img = image object, RGB colorspace
    #   roi = region for white reference (position of ColorChecker Passport)
    img1 = pcv.white_balance(img, roi = (910, 3555, 30, 30))

    # STEP 2: Mask out color card and stake
    # inputs:
    #   img = grayscale image ('a' channel)
    #   p1 = (x,y) coordinates for top left corner of rectangle
    #   p2 = (x,y) coordinates for bottom right corner of rectangle
    #   color = color to make the mask (white here to match background)
    masked, binary, contours, hierarchy = pcv.rectangle_mask(img1, (0, 2000), (1300, 4000), color = "white")
    masked2, binary, contours, hierarchy = pcv.rectangle_mask(masked, (0, 3600), (4000, 4000), color = "white")

    # STEP 3: Convert from RGB colorspace to LAB colorspace
    # Keep green-magenta channel (a)
    # inputs:
    #   img = image object, RGB colorspace
    #   channel = color subchannel ('l' = lightness, 'a' = green-magenta, 'b' = blue-yellow)
    a = pcv.rgb2gray_lab(masked2, 'l')

    # STEP 4: Set a binary threshold on the saturation channel image
    # inputs:
    #   img = img object, grayscale
    #   threshold = treshold value (0-255) - need to adjust this
    #   max_value = value to apply above treshold (255 = white)
    #   object_type = light or dark
    img_binary = pcv.threshold.binary(a, 118, 255, object_type = "dark")

    # STEP 5: Apply Gaussian blur to binary image (reduced noise)
    # inputs:
    #   img = img object, binary
    #   ksize = tuple of kernel dimensions, e.g. (5,5)
    blur_image = pcv.median_blur(img_binary, 10)

    # STEP 6: Fill small objects (speckles)
    # inputs:
    #   img = img object, binary
    #   size = minimum object area size in pixels
    fill_image1 = pcv.fill(blur_image, 150000)

    # STEP 7: Invert image to fill gaps
    # inputs:
    #   img = img object, binary
    inv_image = pcv.invert(fill_image1)
    # rerun fill on inverted image
    inv_fill = pcv.fill(inv_image, 25000)
    # invert image again
    fill_image = pcv.invert(inv_fill)

    # STEP 8: Dilate to avoid losing detail
    # inputs:
    #   img = img object, binary
    #   ksize = kernel size
    #   i = iterations (number of consecutive filtering passes)
    dilated = pcv.dilate(fill_image, 2, 1)

    # STEP 9: Find objects (contours: black-white boundaries)
    # inputs:
    #   img = img object, RGB colorspace
    #   mask = binary image used for object detection
    id_objects, obj_hierarchy = pcv.find_objects(img1, dilated)

    # STEP 10: Define region of interest (ROI)
    # inputs:
    #   img = img object to overlay ROI
    #   x = x-coordinate of upper left corner for rectangle
    #   y = y-coordinate of upper left corner for rectangle
    #   h = height of rectangle
    #   w = width of rectangle
    roi_contour, roi_hierarchy = pcv.roi.rectangle(img=img1, x=20, y=10, h=3000, w=3000)

    # STEP 11: Keep objects that overlap with the ROI
    # inputs:
    #   img = img where selected objects will be displayed
    #   roi_type = options are 'cutto', 'partial' (objects are partially inside roi), or 'largest' (keep only the biggest boi)
    #   roi_countour = contour of roi, output from 'view and adjust roi' function (STEP 10)
    #   roi_hierarchy = contour of roi, output from 'view and adjust roi' function (STEP 10)
    #   object_contour = contours of objects, output from 'identifying objects' function (STEP 9)
    #   obj_hierarchy = hierarchy of objects, output from 'identifying objects' function (STEP 9)
    roi_objects, roi_obj_hierarchy, kept_mask, obj_area = pcv.roi_objects(img1, 'partial', roi_contour, roi_hierarchy, id_objects, obj_hierarchy)

    # STEP 12: Cluster multiple contours in an image based on user input of rows/columns
    # inputs:
    #   img = img object (RGB colorspace)
    #   roi_objects = object contours in an image that will be clustered (output from STEP 11)
    #   roi_obj_hierarchy = object hierarchy (also from STEP 11)
    #   nrow = number of rows for clustering (desired rows in image even if no leaf present in all)
    #   ncol = number of columns to cluster (desired columns in image even if no leaf present in all)
    clusters_i, contours, hierarchies = pcv.cluster_contours(img1, roi_objects, roi_obj_hierarchy, 3, 3)

    # STEP 13: select and split clustered contours to export into multiple images
    # also checks if number of inputted filenames matches number of clustered contours
    # if no filenames, objects are numbered in order
    # inputs:
    #   img = masked RGB image
    #   grouped_contour_indexes = indexes of clustered contours, output of 'cluster_contours' (STEP 12)
    #   contours = contours of cluster, output of 'cluster_contours' (STEP 12)
    #   hierarchy = object hierarchy (from STEP 12)
    #   outdir = directory to export output images
    #   file = name of input image to use as basename (uses filename from 'readimage')
    #   filenames = (optional) txt file with list of filenames ordered from top to bottom/left to right

    # Set global debug behavior to None (default), "print" (to file), or "plot" (Jupyter Notebooks or X11)
    pcv.params.debug = None
    out = args.outdir
    # names = args.namesout = "./"

    output_path, imgs, masks = pcv.cluster_contour_splitimg(img1, clusters_i, contours, hierarchies, out, file = filename, filenames = None)

# Call program
if __name__ == '__main__':
    main()
