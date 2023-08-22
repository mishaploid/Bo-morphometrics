################################################################################
## Elliptical Fourier descriptors and PCA for B. oleracea leaves
## Authors: Sarah Turner-Hissong & Makenzie Mabry
## Date: 25 November 2019
## Updated: 27 May 2020
################################################################################

library(tidyverse)
library(Momocs)

# load data ---------------------------------------------------------
# scanned images
load("data/processed/aligned_scans.RData")
# images from Snogerup et al. 
load("data/processed/Snogerup_landmarked.RData")

# specify outlier images --------------------------------------------------
# images that are misaligned 

pca_outliers <- c("B_oleracea012_A_1_p1_mask", 
                  "B_oleracea023_B_2_p2_mask",
                  "B_oleracea034_B_leaf2_3_1_p1_mask", 
                  "B_oleracea117_B_2_p2_mask",
                  "B_oleracea157_B_0_p0_mask", 
                  "B_oleracea207_C_0_p0_mask",
                  "B_oleracea218_A_0_p0_mask", 
                  "B_oleracea232_D_0_p0_mask")

# align, center, and scale Snogerup images --------------------------------

Snogerup_aln <- Snogerup_outlines %>% 
  coo_align() %>% 
  coo_center() %>% 
  coo_scale() %>%
  coo_slide(ldk = 3) 

leaves_Snogerup <- leaves_aln %>% 
  Momocs::filter(!img %in% pca_outliers) %>% 
  Momocs::combine(., Snogerup_aln) %>% 
  coo_align() %>%
  coo_center() %>%
  coo_scale() %>%
  coo_slide(ldk = 3) 

# Elliptical Fourier Analysis ---------------------------------------------

# Quantitative calibration via harmonic power
# retain n harmonics to gather 99% of harmonic power
harmonic_power <- calibrate_harmonicpower_efourier(leaves_Snogerup, 
                                                   nb.h = 12)

ggsave("reports/harmonic_ranks.png")

# Calibration via deviations of reconstructed and original shapes 
# uses a range of harmonic numbers 
harmonic_deviations <- calibrate_deviations_efourier(leaves_Snogerup)

ggsave("reports/harmonic_deviations.png")

# Calibrate reconstructions
leaves_Snogerup %>%
  calibrate_reconstructions_efourier(id = 1, 
                                     range = 1:9)

ggsave("reports/calibrate_reconstructions.png")

# compute elliptical Fourier descriptors (EFDs)
# EFDs with symmetric & asymmetric variation 
leaves_f_asym <- efourier(leaves_Snogerup,
                          nb.h = 7,
                          norm = FALSE)

# EFDs with symmetric variation only 
leaves_f_sym <- efourier(leaves_Snogerup, 
                     nb.h = 7,
                     norm = FALSE) %>% 
  rm_asym() # remove asymmetry

# export EFD RData object -------------------------------------------------

save(leaves_f_sym, leaves_f_asym, leaves_aln, leaves_Snogerup,
     file = "data/processed/EFDs.RData")
