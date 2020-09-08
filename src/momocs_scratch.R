# devtools::install_github("MomX/Momocs")
library(Momocs)
library(tidyverse)
library(readxl)
library(RColorBrewer)


setwd("~/")

exclude <- read.csv("~/Box/BoleraceaLeafScans/images_to_exclude.csv", header = FALSE)[,1] %>%
  as.character() 

# list and read images into R
img <- list.files("~/Box/BoleraceaLeafScans/segmented_images",
                  pattern = "p0_mask.jpg",
                  full.names = TRUE) %>%
  .[!grepl(paste(exclude, collapse = "|"), .)] %>%
  import_jpg() 

img %>% 
  Out() %>%
  coo_center() %>%
  coo_scale() %>%
  coo_align() %>%
  panel(col="blue")

factors <- data.frame(image = names(img)) %>% # [c(-2383,-2413)]) %>%
  separate(image, into = c("genus", "id", "rep", "stage", "mask"), 
           remove = FALSE,
           sep = "_") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(id = gsub("oleracea", "", id)) %>%
  filter(image %in% names(img)) # [c(-2383,-2413)])

accession_info <- read_xlsx("~/Documents/B_oleracea/Bo_sample_ids.xlsx")  

leafInfo <- read.csv("~/Box/BoleraceaLeafScans/src/sample_info.csv", sep = ",", 
                     header = TRUE) %>%
  mutate_all(as.factor) %>%
  mutate(img = str_replace(img, ".jpg", ""))

factors <- left_join(factors, leafInfo, by = c("image" = "img"))

# factors <- left_join(factors, accession_info, by = c("id", "rep")) %>%
#   select(genus, id, rep, stage, Morphotype) %>%
#   mutate_if(is.character, as.factor) %>%
#   mutate(Morphotype = fct_recode(Morphotype, kohlrabi = "Purple Kohlrabi", kohlrabi = "Green Kohlrabi",
#                                  kohlrabi = "Kohlrabi", Cauliflower = "Caulflower",
#                                  Cabbage = "Red Cabbage",
#                                  "walking stick kale" = "Giant Jersey Kale"),
#          morph = fct_collapse(Morphotype,
#                              kale = c("Collard Greens", 
#                                       "walking stick kale",
#                                       "Lacinato kale",
#                                       "Perpetual Kale",
#                                       "Chinese White Kale"),
#                              cabbage = c("Cabbage",
#                                          "Marrow Cabbage",
#                                          "Savoy Cabbage",
#                                          "Tronchuda Kale",
#                                          "Ornamental cabbage"),
#                              cauliflower = c("Cauliflower",
#                                              "Romanesco")),
#         type = fct_collapse(morph,
#                             wild = c("Wild C", "Wild oleracea"),
#                             kale = c("kale"),
#                             tuber = c("kohlrabi"),
#                             heading = c("cabbage", "Broccoli",
#                                         "Brussels Sprouts", "cauliflower")))

# create collection of coordinates
# leaves <- Out(img[c(-2383,-2413)], fac = factors)

# leaves <- Out(img, fac = factors)
leaves <- Out(img, fac = factors) 

# c(1:841, 843:1114)
leaves %>%
  # coo_smooth(15) %>% # smooth to 10 points (probably need to adjust this)
  coo_center() %>% # center images by origin
  # coo_scale() %>%
  coo_align() %>% # align coordinates along longest axis - may not want to do this (some leaves are wider than they are long)
  coo_rotate(theta = -80) %>% 
  panel(fac = "morphotype", palette = col_summer) 

# save.image("Box/BoleraceaLeafScans/Processed Leaf Scans/src/momocs_test.R")
# load("Box/BoleraceaLeafScans/Processed Leaf Scans/src/momocs_test.RData")

# plot outlines
leaves %>% 
  coo_center() %>%
  coo_align() %>%
  pile()


# define landmarks
# leaves2 <- def_ldk(leaves, nb.ldk = 1)
leaves$ldk <- leaves2$ldk 

leaves %>%
  coo_align() %>%
  coo_center() %>%
  coo_slide(ldk = 1) %>%
  pile()

leaves3 <- leaves %>%
  coo_align() %>%
  coo_center() %>%
  coo_scale() %>%
  coo_slide(ldk = 1)

leaves4 <- Momocs::filter(leaves3, leafNum == "2")

# calibrate number of harmonics required
leaves4 %>%
  calibrate_harmonicpower_efourier(nb.h = 15) # 99% of harmonic power gathered by X harmonics

# need to think more about how to interpret harmonics 
calibrate_deviations_efourier(leaves4)

calibrate_reconstructions_efourier(leaves4) # neat

leaves.f <- efourier(leaves4, nb.h = 7, norm = FALSE) 

leaves.f 

# leaves.f <- rm_asym(leaves.f)

boxplot(leaves.f, drop = 1) 


# leaves4 <- Momocs::filter(leaves.f, !id == "112")

leaves.p <- PCA(leaves.f)
class(leaves.p)

plot_PCA(leaves.p, "morphotype", chull = FALSE)
plot_PCA(leaves.p, "leafNum")
plot_PCA(leaves.p, ~pheno, chull = FALSE) %>% 
  layer_stars()

leaves.p %>% as_df() %>% 
  # filter(Morphotype %in% c("Cabbage", "Brussels Sprouts")) %>%
  ggplot() +
  aes(x = PC1, y = PC2, col = morphotype) + 
  coord_equal() + 
  geom_point() + 
  geom_density2d() + 
  theme_light() +
  facet_wrap(~morphotype)

pc <- leaves.p %>% 
  as_df() %>%
  as_tibble()
  
pc %>%
  gather(., key = "pcs", value = "value", -c(1:12)) %>%
  filter(pcs %in% c("PC1", "PC2", "PC3", "PC4")) %>%
  ggplot(., aes(morphotype, value, col = morphotype)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  facet_grid(~pcs)


  # plot_ly(x = .$PC1, y = .$PC2, z = .$PC3, 
  #         type = "scatter3d", 
  #         mode = "markers", 
  #         color = ~stage)

KMEANS(leaves.p, centers = 6) # cool! 

ms <- MSHAPES(leaves.f, "morphotype")
ms <- ms$shp

coo_plot(ms$medullosa, border = "orange")
coo_plot(ms$sabellica, plot.new = FALSE, border = "blue")
coo_plot(ms$gemmifera, plot.new = FALSE, border = "green")
coo_plot(ms$gongylodes, plot.new = FALSE, border = "red")

leaves.f %>% mshapes("morphotype") %>% plot_mshapes()


circularity <- coo_circularity(leaves4) %>%
  map_df(~data_frame(circularity = .x), .id = ".id")

area <- coo_area(leaves4) %>% 
  map_df(~data_frame(area = .x), .id = ".id")

convex <- coo_convexity(leaves4) %>%
  map_df(~data_frame(area = .x), .id = ".id")

solidity <- coo_solidity(leaves4) %>%
  map_df(~data_frame(solidity = .x), .id = ".id")

w <- coo_width(leaves4) %>%
  map_df(~data_frame(width = .x), .id = ".id")

l <- coo_length(leaves4) %>%
  map_df(~data_frame(length = .x), .id = ".id")

all <- Reduce(function(x, y) merge(x, y, all=TRUE, by = ".id"), 
              list(pc, circularity, area, convex, solidity, w, l))

all <- all %>%
  mutate(ar = width/length)
         # morph = fct_collapse(Morphotype,
         #                    kale = c("Collard Greens", 
         #                             "walking stick kale",
         #                             "Lacinato kale",
         #                             "Perpetual Kale",
         #                             "Chinese White Kale"),
         #                    cabbage = c("Cabbage",
         #                                "Marrow Cabbage",
         #                                "Savoy Cabbage",
         #                                "Tronchuda Kale",
         #                                "Ornamental cabbage"),
         #                    cauliflower = c("Cauliflower",
         #                                    "Romanesco")),
         # type = fct_collapse(morph,
         #                     wild = c("Wild C", "Wild oleracea"),
         #                     kale = c("kale"),
         #                     tuber = c("kohlrabi"),
         #                     heading = c("cabbage", "Broccoli",
         #                                 "Brussels Sprouts", "cauliflower")))

fct_count(all$morphotype)


pca <- all %>%
  ggplot(., aes(PC1, PC2, colour = morphotype)) +
  theme_minimal() + 
  coord_fixed() +
  geom_density_2d() 
#   scale_colour_brewer(palette = "Dark2") 

p11 <- all %>%
  ggplot(., aes(ar, circularity/100)) +
  theme_minimal() + 
  coord_fixed() +
  geom_point(size = 0.75) +   
  # xlim(c(0, 1)) +
  # ylim(c(0, 1)) + 
  scale_colour_brewer(palette = "Dark2") +
  xlab("Aspect Ratio") +
  ylab("Circularity/5")

p12 <- all %>%
  ggplot(., aes(ar, circularity/100, colour = morphotype)) +
  theme_minimal() + 
  coord_fixed() +
  # xlim(c(0, 1)) +
  # ylim(c(0, 1)) + 
  geom_density_2d(size = 1.5) + 
  # scale_colour_brewer(palette = "Dark2") +
  xlab("Aspect Ratio") +
  ylab("Circularity/5")
  

p21 <- ggplot(all, aes(ar, solidity^8)) +
    geom_point(size = 0.75) + 
    theme_minimal() + 
    coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) + 
  xlab("Aspect Ratio") +
  ylab(expression(Solidity^8))

p22 <- ggplot(all, aes(ar, solidity^8, colour = morph)) +
  theme_minimal() + 
  coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) + 
  geom_density_2d(size = 1.5) + 
  scale_colour_brewer(palette = "Dark2") +
  xlab("Aspect Ratio") +
  ylab(expression(Solidity^8))

p31 <- ggplot(all, aes(solidity^8, circularity/5)) +
  geom_point(size = 0.75) + 
  theme_minimal() + 
  coord_fixed() +
  xlab(expression(Solidity^8)) +
  ylab("Circularity/5")
  
p32 <- ggplot(all, aes(solidity^8, circularity/5, colour = type)) +
  theme_minimal() + 
  coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1) + 
  geom_density_2d(size = 1.5) + 
  scale_colour_brewer(palette = "Dark2") +
  xlab(expression(Solidity^8)) +
  ylab("Circularity/5")

library(cowplot)

plot_grid(p11, p12, p21, p22, p31, p32, ncol = 2)

library(viridis)

p2.1 <- ggplot(all, aes(PC1, PC2, colour = ar)) +
  geom_point() + 
  theme_minimal() + 
  coord_fixed() +
  xlim(-0.35, 0.35) +
  ylim(-0.35, 0.35) + 
  scale_color_viridis(name = "Aspect Ratio") +
  xlab("PC1 (49.8%)") +
  ylab("PC2 (15.2%)")

p2.2 <- ggplot(all, aes(PC1, PC2, colour = circularity/100)) +
  geom_point() + 
  theme_minimal() + 
  coord_fixed() +
  xlim(-0.35, 0.35) +
  ylim(-0.35, 0.35) + 
  scale_color_viridis(name = "Circularity/5") +
  xlab("PC1 (49.8%)") +
  ylab("PC2 (15.2%)")

p2.3 <- ggplot(all, aes(PC1, PC2, colour = solidity^8)) +
  geom_point() + 
  theme_minimal() + 
  coord_fixed() +
  xlim(-0.35, 0.35) +
  ylim(-0.35, 0.35) + 
  scale_color_viridis(name = expression(Solidity^8)) +
  xlab("PC1 (49.8%)") +
  ylab("PC2 (15.2%)")


p2.4 <- ggplot(all, aes(PC1, PC2, colour = morphotype)) +
  geom_density_2d(size = 1.25)  +
  theme_minimal() + 
  # scale_colour_brewer(palette = "Dark2") +
  coord_fixed() +
  # xlim(-0.35, 0.35) +
  # ylim(-0.35, 0.35) + 
  xlab("PC1") +
  ylab("PC2") +
  facet_wrap(~morphotype)


plot_grid(p2.1, p2.2, p2.3, p2.4, ncol = 2)


## just kales 
kales <- which(!grepl("kale|Kale|Wild", leaves$fac$Morphotype))

kalesdf <- leaves %>%
  Momocs::slice(., kales)

kalesdf %>%
  calibrate_harmonicpower_efourier(nb.h = 15) # 99% of harmonic power gathered by X harmonics

kales.f <- efourier(kalesdf, nb.h = 13, norm = TRUE) 

kales.f 

# leaves.f <- rm_asym(leaves.f)

boxplot(kales.f, drop = 1) 

kales.p <- PCA(kales.f)

plot_PCA(kales.p, ~Morphotype)

kalespc <- kales.p %>% 
  as_df() %>%
  as_tibble()

all <- Reduce(function(x, y) merge(x, y, all=TRUE, by = ".id"), 
              list(pc, circularity, area, convex, solidity, w, l))

all <- all %>%
  mutate(ar = width/length,
         morph = fct_collapse(Morphotype,
                              kale = c("Collard Greens", 
                                       "walking stick kale",
                                       "Lacinato kale",
                                       "Perpetual Kale",
                                       "Chinese White Kale"),
                              cabbage = c("Cabbage",
                                          "Marrow Cabbage",
                                          "Savoy Cabbage",
                                          "Tronchuda Kale",
                                          "Ornamental cabbage"),
                              cauliflower = c("Cauliflower",
                                              "Romanesco")),
         type = fct_collapse(morph,
                             wild = c("Wild C", "Wild oleracea"),
                             kale = c("kale"),
                             tuber = c("kohlrabi"),
                             heading = c("cabbage", "Broccoli",
                                         "Brussels Sprouts", "cauliflower")))

all %>%
  ggplot(., aes(PC1, PC2, colour = morph)) +
  geom_density_2d(size = 1.25)  +
  theme_minimal() + 
  coord_fixed() +
  facet_wrap(~morph)

