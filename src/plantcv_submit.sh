#!/bin/bash

~/Box/BoleraceaLeafScans/plantcv-workflow.py \
-d ~/Box/BoleraceaLeafScans/testdir \
-a filename \
-p ~/Box/BoleraceaLeafScans/binarize_images.py \
-j ~/Box/BoleraceaLeafScans/testrun.json \
-i ~/Box/BoleraceaLeafScans/plantcv_test \
-f camera_id_treatment \
-t tif \
-T 1 

# -d  directory of images
# -a  adaptor to indicate structure to grab metadata from (e.g. filename)
# -p  workflow to run over directory of images
# -j  json database name
# -i  output directory for images
# -f  meta data format map (e.g. 'imgtype_camera_frame_zoom_id')
# -t  type extension for output (e.g. 'tif', 'jpg' (default))
# -T  number of threads to use
# -w  writeimg option (True will write images, default = False)
