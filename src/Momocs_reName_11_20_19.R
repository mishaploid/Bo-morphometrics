##11.19.19### Morophometerics script for renameing for momocs

# devtools::install_github("MomX/Momocs")
library(Momocs)
library(tidyverse)
library(readxl)
library(viridis)
library(RColorBrewer)


setwd("~/Box Sync/BoleraceaLeafScans/")

#1. read in list of files with problems that we do not want in the analysis, followed by reading in files

exclude <- read.csv("images_to_exclude.csv", header = FALSE)[,1] %>%
  as.character() 

img <- list.files("segmented_images",
                  pattern = "_mask.jpg",
                  full.names = TRUE) %>%
  .[!grepl(paste(exclude, collapse = "|"), .)] %>%
  import_jpg() 


#2. Read in file to assoicate file names with rep/samples/morphotype information

leafInfo <- read.csv("RenameFile.csv", sep = ",", header = TRUE, row.names = 1)

exclude2 <- read.csv("images_to_exclude2.csv", header = FALSE)[,1]

leafInfo2 <- leafInfo[ ! row.names(leafInfo) %in% exclude2, ]

leafInfo2[, "Morphotype"] <- as.factor(leafInfo2[, "Morphotype"])
leafInfo2[, "RepLetter"] <- as.factor(leafInfo2[, "RepLetter"])
leafInfo2[, "LeafNum"] <- as.factor(leafInfo2[, "LeafNum"])
leafInfo2[, "Phenotype"] <- as.factor(leafInfo2[, "Phenotype"])
leafInfo2[, "SampleNum"] <- as.factor(leafInfo2[, "SampleNum"])


##just checking that it is reading it correctly
img %>% 
  Out() %>%
  coo_center() %>%
  coo_scale() %>%
  coo_align() %>%
  panel(col="blue")

##example how to pull out specific files
row.names(leafInfo[leafInfo$Morphotype == "gemmifera", ])

leaves <- Out(img, fac = leafInfo2)


# c(1:841, 843:1114)
leaves %>%
  coo_center() %>%
  # coo_smooth(15) %>% # smooth to 10 points (probably need to adjust this)
  coo_center() %>% # center images by origin
  # coo_scale() %>%
  coo_align() %>% # align coordinates along longest axis - may not want to do this (some leaves are wider than they are long)
  panel(fac = "Morphotype", palette = col_sari, names=fac_dispatcher(leaves, 4)) 

# save.image("Box/BoleraceaLeafScans/Processed Leaf Scans/src/momocs_test.R")
# load("Box/BoleraceaLeafScans/Processed Leaf Scans/src/momocs_test.RData")

#panel(leaves, fac = "Morphoype")

# plot outlines
leaves %>% pile()

leaves2 <- leaves %>%
  # coo_smooth(15) %>% # smooth to 10 points (probably need to adjust this)
  coo_center() %>% # center images by origin
  # coo_scale() %>%
  coo_align()

# calibrate number of harmonics required
leaves2 %>%
  calibrate_harmonicpower_efourier(nb.h = 10) # 99% of harmonic power gathered by X harmonics

# need to think more about how to interpret harmonics 
calibrate_deviations_efourier(leaves2)

calibrate_reconstructions_efourier(leaves2) # neat

leaves.f <- efourier(leaves2, nb.h = 6, norm = FALSE) 

leaves.f 

# leaves.f <- rm_asym(leaves.f)

boxplot(leaves.f, drop = 1) 

leaves.p <- PCA(leaves.f)
class(leaves.p)


plot_PCA(leaves.p)
plot_PCA(leaves.p, "Morphotype", legend = TRUE)
plot_PCA(leaves.p, ~Morphotype, chull = FALSE) %>% 
  layer_stars()

leaves.p %>% as_df() %>% 
  # filter(Morphotype %in% c("Cabbage", "Brussels Sprouts")) %>%
  ggplot() +
  aes(x = PC1, y = PC2, col = Morphotype) + 
  coord_equal() + 
  geom_point() + 
  geom_density2d() + 
  theme_light() +
  facet_wrap(~Morphotype)

pc <- leaves.p %>% 
  as_df() %>%
  as_tibble()

# plot_ly(x = .$PC1, y = .$PC2, z = .$PC3, 
#         type = "scatter3d", 
#         mode = "markers", 
#         color = ~stage)

KMEANS(leaves.p, centers = 5) # cool! 

ms <- MSHAPES(leaves.f, "Morphotype")
ms <- ms$shp

coo_plot(ms$`Marrow Cabbage`, border = "orange")
coo_plot(ms$`Curly Kale`, plot.new = FALSE, border = "blue")
coo_plot(ms$`Brussels Sprouts`, plot.new = FALSE, border = "green")
coo_plot(ms$`kohlrabi`, plot.new = FALSE, border = "red")

leaves.f %>% mshapes("morph") %>% plot_mshapes()


circularity <- coo_circularity(leaves) %>%
  map_df(~data_frame(circularity = .x), .id = ".id")

area <- coo_area(leaves) %>% 
  map_df(~data_frame(area = .x), .id = ".id")

convex <- coo_convexity(leaves) %>%
  map_df(~data_frame(area = .x), .id = ".id")

solidity <- coo_solidity(leaves) %>%
  map_df(~data_frame(solidity = .x), .id = ".id")

w <- coo_width(leaves) %>%
  map_df(~data_frame(width = .x), .id = ".id")

l <- coo_length(leaves) %>%
  map_df(~data_frame(length = .x), .id = ".id")

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

fct_count(all$morph)


pcs <- all %>%
  ggplot(., aes(PC1, PC2, colour = type)) +
  theme_minimal() + 
  coord_fixed() +
  geom_density_2d() + 
  scale_colour_brewer(palette = "Dark2") 

p11 <- all %>%
  ggplot(., aes(ar, circularity/5)) +
  theme_minimal() + 
  coord_fixed() +
  geom_point(size = 0.75) +   
  xlim(c(0, 1)) +
  ylim(c(0, 1)) + 
  scale_colour_brewer(palette = "Dark2") +
  xlab("Aspect Ratio") +
  ylab("Circularity/5")

p12 <- all %>%
  ggplot(., aes(ar, circularity/5, colour = type)) +
  theme_minimal() + 
  coord_fixed() +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) + 
  geom_density_2d(size = 1.5) + 
  scale_colour_brewer(palette = "Dark2") +
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

p2.2 <- ggplot(all, aes(PC1, PC2, colour = circularity/5)) +
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


p2.4 <- ggplot(all, aes(PC1, PC2, colour = type)) +
  geom_density_2d(size = 1.25)  +
  theme_minimal() + 
  scale_colour_brewer(palette = "Dark2") +
  coord_fixed() +
  # xlim(-0.35, 0.35) +
  # ylim(-0.35, 0.35) + 
  xlab("PC1") +
  ylab("PC2") +
  facet_wrap(~type)


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


all %>%
  ggplot(., aes(PC1, PC2, colour = morph)) +
  geom_density_2d(size = 1.25)  +
  theme_minimal() + 
  coord_fixed() +
  facet_wrap(~morph)

