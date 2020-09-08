################################################################################
## Load leaf outlines and accession info for B. oleracea leaf scans
## Authors: Sarah Turner-Hissong & Makenzie Mabry
## Date: 25 November 2019
################################################################################

# This script uses R/Momocs to load segmented leaf scans for a diversity panel
# of Brassica oleracea. The output is an RData 'Coo' object of outlines with 
# a common landmark and corresponding sample information. 

################################################################################
## setup
################################################################################

# load packages
# devtools::install_github("MomX/Momocs")
library(Momocs)
library(tidyverse)

################################################################################
## load data
## list image names and use Momocs::import_jpg() to read in outlines
## returns a 'Coo' object (collection of coordinates)
################################################################################

# list of images to exclude
# these are images where the leaf was torn or where a stake was segmented
exclude <- read.csv("~/Box/BoleraceaLeafScans/images_to_exclude.csv", 
                    header = FALSE)[,1] %>%
  as.character() 

# list segmented images and import outline coordinates 
img <- list.files("~/Box/BoleraceaLeafScans/segmented_images",
                  pattern = "_mask.jpg",
                  full.names = TRUE)[2151:2200] %>%
  .[!grepl(paste(exclude, collapse = "|"), .)] %>% # exclude problem images
  import_jpg() 

################################################################################
## import info for each image
## e.g. sample id, morphotype, leaf stage, etc. 
################################################################################

img_names <- data.frame(img = names(img))

leafInfo <- read.csv("~/Box/BoleraceaLeafScans/src/sample_info.csv", 
                     sep = ",", 
                     header = TRUE) %>%
  mutate_all(as.factor) %>% # set all variables to factors
  mutate(img = str_replace(img, ".jpg", "")) %>%
  left_join(img_names, ., by = "img") %>%
  mutate(morphotype = fct_recode(morphotype, "viridis" = "virdis"),
         morphotype = fct_recode(morphotype, "palmifolia" = "longata"))

leaves <- Out(img, fac = leafInfo) # add factor info to outline (Coo) object

leaves <- leaves %>%
  Momocs::filter(!img %in% exclude)

################################################################################
## plot leaf outlines to check for issues
################################################################################

leaves %>%
  coo_smooth(50) %>% 
  coo_center() %>% # center images by origin
  coo_alignxax() %>% # align coordinates along longest axis
  coo_rotate(theta = -80) %>% 
  panel(fac = "morphotype", palette = col_summer, names = "morphotype") 

# ggsave("~/Desktop/ghat/leaf_outlines.png", height = 20, width = 20)

################################################################################
# define landmarks and save outlines
################################################################################

leaves <- def_ldk(leaves, nb.ldk = 4)

save(leaves, file = "~/Box/BoleraceaLeafScans/data/processed/leaf_outlines2151_2200.RData")
