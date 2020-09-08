################################################################################
## Elliptical Fourier descriptors and PCA for B. oleracea leaves
## Authors: Sarah Turner-Hissong & Makenzie Mabry
## Date: 25 November 2019
## Updated: 27 May 2020
################################################################################

library(tidyverse)
library(Momocs)
library(viridis)

# load data ---------------------------------------------------------
load("~/Box/BoleraceaLeafScans/data/processed/aligned_scans.RData")
load("~/Box/BoleraceaLeafScans/data/processed/Snogerup_landmarked.RData")

# specify outlier images --------------------------------------------------

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

################################################################################
## Elliptical Fourier analysis
################################################################################

# Quantitative calibration via harmonic power
# retain n harmonics to gather 99% of harmonic power
harmonic_power <- calibrate_harmonicpower_efourier(leaves_Snogerup, 
                                                   nb.h = 12)

ggsave("~/Box/BoleraceaLeafScans/reports/harmonic_ranks.png")

# Calibration via deviations of reconstructed and original shapes 
# uses a range of harmonic numbers 
harmonic_deviations <- calibrate_deviations_efourier(leaves_Snogerup)

ggsave("~/Box/BoleraceaLeafScans/reports/harmonic_deviations.png")

# Calibrate reconstructions
leaves_Snogerup %>%
  calibrate_reconstructions_efourier(id = 1, 
                                     range = 1:9)

ggsave("~/Box/BoleraceaLeafScans/reports/calibrate_reconstructions.png")

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

################################################################################
## Principal component analysis
################################################################################

leaves.p <- PCA(leaves_f_sym) 

# export the principal components
leaves.p %>% 
  as_df() %>%
  select(img, id, rep, leafNum, morphotype, pheno, PC1, PC2, PC3, PC4) %>% 
  write_csv("~/Box/BoleraceaLeafScans/data/processed/leaf_pcs.csv")

leaves.p$fac <- leaves.p$fac %>% 
  mutate(source = ifelse(pheno %in% "Snogerup",
                         paste0(morphotype),
                         "Mabry")) 

leaves.s <- leaves.p %>% 
  Momocs::filter(pheno %in% "Snogerup") 

png("~/Box/BoleraceaLeafScans/reports/Snogerup_pca.png",
    height = 5, width = 5, units = "in", res = 300)
leaves.s %>% 
  plot_PCA(~source,
           chull = FALSE,
           points = FALSE,
           palette = pal_qual_Dark2) %>% 
  layer_points(cex = 2,
               transp = 0.2,
               pch = 17)
dev.off() 

png("~/Box/BoleraceaLeafScans/reports/pc2_pc3.png",
    height = 4,
    width = 10,
    units = "in",
    res = 300)

# PCA symmetric and asymmetric
par(mfrow = c(1,2))

plot_PCA(leaves.p,  
         ~morphotype,
         chull = FALSE,
         legend = FALSE, 
         axes = c(3,4),
         palette = pal_qual_Dark2,
         title = "symmetric and asymmetric variation") 

# PCA no asymmetric
leaves_f_sym %>% 
  PCA() %>% 
  plot_PCA(~morphotype,
           chull = FALSE, 
           legend = FALSE,
           axes = c(3,4),
           palette = pal_qual_Dark2,
           title = "symmetric variation only")

dev.off() 

# overplot measurements

leaves.p %>% 
  as_df() %>% 
  left_join(all, .) %>% 
  select(id:PC2) %>% 
  gather(key = measurement,
         value = value, 
         -c(id:rep, morphotype:pheno, PC1:PC2)) %>% 
  group_by(measurement) %>% 
  mutate(percentile = percent_rank(value)) %>% 
  ggplot(., aes(PC1, PC2, color = percentile)) +
  scale_color_viridis_c() +
  geom_point() +
  facet_wrap(~measurement,
             ncol = 4) +
  theme_bw() +
  theme(strip.text = element_text(size = 14))

ggsave("~/Box/BoleraceaLeafScans/reports/pca_traditional_metrics.png",
       height = 5, width = 10)





leaves.p %>% as_df() %>% 
  ggplot() +
  aes(x = PC1, y = PC2) + # col = leafNum) + 
  coord_equal() + 
  geom_point(size = 0.5) + 
  geom_density2d() + 
  theme_light() +
  facet_wrap(~morphotype)

# variance of PC scores by morphotype 
leaves.p %>%
  as_df() %>%
  select(img, id, rep, leafNum, morphotype, pheno, PC1, PC2, PC3) %>% 
  mutate(morphotype = fct_reorder(morphotype, PC1)) %>% 
  gather(key = pc,
         value = value,
         -img, -id, -rep, -leafNum, -morphotype, -pheno) %>% 
  ggplot(aes(x = value, y = morphotype, col = morphotype)) +
  geom_jitter(height = 0.25) +
  facet_grid(~pc,
             scales = "free_x") +
  theme_bw()

leaves.p %>%
  as_df() %>% 
  filter(morphotype %in% c("insularis", "villosa")) %>% 
  ggplot(aes(x = PC1, y = PC2, col = morphotype)) + 
  scale_color_manual(values = pal_qual_Dark2(22)[c(10,21)]) + 
  geom_point(col = "black", size = 1) + 
  geom_density2d(lwd = 1.5) + 
  # facet_wrap(~morphotype, 
  #            ncol = 6) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 18)) +
  xlim(-0.45, 0.45) + 
  ylim(-0.4, 0.4) + 
  ggtitle("B. insularis v. B. villosa")

ggsave("~/Box/BoleraceaLeafScans/reports/insularis_v_villosa.png",
       height = 3,
       width = 4)

################################################################################
## shape variation along PC axes 
################################################################################

gg <- PCcontrib(leaves.p,
                nax = 1:4,
                sd.r = c(-3, -1.5, 0, 1.5, 3))

# eigenleaf representations for PCs representing variance in EFDs 
sum_eig <- sum(leaves.p$eig)
percents <- leaves.p$eig/sum_eig*100
labs <- data.frame(pc = c("PC1:", "PC2:", "PC3:", "PC4:"),
                   percent = paste0(round(percents[1:4], 1), "%"))

x_facet_labs <- paste(labs$pc, labs$percent, sep = "\n")
names(x_facet_labs) <- c(1:4)
y_facet_labs <- c("-3 s.d.", "-1.5 s.d.", "Mean", "+1.5 s.d.", "3 s.d.")
names(y_facet_labs) <- c(-3, -1.5, 0, 1.5, 3)

col_pal <- rev(pal_div_BrBG(5))
col_pal[3] <- "dark gray"

gg$gg + 
  coord_flip() +
  scale_x_reverse() + 
  geom_polygon(aes(col = as.factor(shp1)),
               fill = "white",
               size = 1.5) + 
  scale_color_manual(values = col_pal) + 
  theme_void() +
  facet_grid(nax ~ shp,
             labeller = labeller(shp = y_facet_labs,
                                 nax = x_facet_labs),
             switch = "y") +
  theme(legend.position = "none",
        strip.text.y.left = element_text(angle = 0,
                                         hjust = 0),
        text = element_text(size = 24))

ggsave("~/Box/BoleraceaLeafScans/reports/PC_shape_variation.png")


# linear discriminant analysis of harmonics 
harmonic_deviations <- calibrate_deviations_efourier(leaves_Snogerup)


### SCRATCH 
PCcontrib(leaves.p, nax = 1:3)

leaves.p %>%
  LDA(~morphotype) %>%
  classification_metrics() 

leaves.l <- LDA(leaves.p, 'morphotype', retain = 0.99)

leaves.cv <- leaves.l$CV.tab

# plot confusion matrix 

counts <- leaves$fac %>%
  group_by(morphotype) %>% 
  count() 

leaves.cv %>% 
  as_tibble() %>% 
  left_join(., counts, by = c("actual" = "morphotype")) %>% 
  left_join(., prop_correct[,c("actual", "prop_correct")], by = "actual") %>% 
  mutate(prop = n.x/n.y,
         actual = fct_reorder(actual, -prop_correct),
         classified = fct_relevel(classified, levels(actual))) %>% 
  ggplot(aes(classified, actual, fill = prop)) + 
  geom_tile() +
  scale_fill_viridis() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("~/Box/BoleraceaLeafScans/reports/lda_confusion_matrix.png",
       height = 5,
       width = 6.5)

plot_LDA(leaves.l)

mean_species <- MSHAPES(leaves.f,
                        fac = "morphotype")$shp

coo_plot(mean_species$insularis,
         border = pal_qual_Dark2(22)[11],
         lwd = 3)
coo_plot(mean_species$villosa,
         border = pal_qual_Dark2(1),
         lwd = 3,
         plot.new = FALSE)


tps_iso(mean_species$insularis, 
        mean_species$villosa, 
        cont = FALSE, 
        amp = 1,
        grid = TRUE,
        shp.border=pal_qual_Dark2(22)[c(10,21)],
        palette = viridis_pal(option = "D", alpha = 0),
        shp.lwd = c(5,5))


# plot all leaves
palmifolia <- Momocs::filter(leaves_aln, 
                             morphotype %in% "palmifolia")

ms_palmifolia <- MSHAPES(leaves.f, "morphotype")


pile(palmifolia,
     transp = 0.9,
     col = pal_qual_Dark2(1))

coo_plot(ms_palmifolia$shp$palmifolia,
         border = col) 

for (i in 1:length(palmifolia)) {
  coo_plot(palmifolia[i], 
           plot.new = FALSE,
           border = col)
}
coo_plot(palmifolia[2], plot.new = FALSE)

