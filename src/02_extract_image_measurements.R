################################################################################
## Traditional shape descriptors for B. oleracea leaf scans
## Authors: Sarah Turner-Hissong & Makenzie Mabry
## Date: 25 November 2019
################################################################################

library(tidyverse)
library(Momocs)
library(GGally)

################################################################################
## load outlines with landmarks 
################################################################################
# list RData files 
files <- list.files("~/Box/BoleraceaLeafScans/data/processed",
                    pattern = "leaf_outlines",
                    full.names = TRUE)

# load RData files with unique names 
leaves_data <- map(files, function(x) {
  load(x)
  get(ls()[ls() != "filename"])
})

# combine RData files into single object of class 'Out' and 'Coo'
leaves <- Momocs::combine(leaves_data[1:length(leaves_data)])

class(leaves)

# drop unused factor levels
leaves$fac <- leaves$fac %>% 
  mutate_if(is.factor, fct_drop)

################################################################################
## align, center, and scale outlines
################################################################################

leaves_aln <- leaves %>%
  coo_align() %>%
  coo_center() %>%
  coo_scale() %>%
  coo_slide(ldk = 3) 

################################################################################
## plot outlines to check for issues 
################################################################################

leaves_aln %>% 
  Momocs::slice(1:1000) %>% 
  stack() 

png("~/Box/BoleraceaLeafScans/reports/panel_all_leaves.png",
    width = 10,
    height = 10,
    units = "in",
    res = 300)
leaves_aln %>% 
  coo_rotate(theta = pi/2) %>% 
  panel(., fac = "morphotype",
        palette = pal_qual_Dark2) 
dev.off() 

################################################################################
## export traditional shape descriptors 
################################################################################
  # area = area of the shape (px^2)
  # circularity = ratio of area to perimeter, 'compactness' or 'roundness'
  #               4*pi*(area/perimeter^2)
  # convex = area of a convex hull (smallest polygon that bounds the shape)
  # solidity = area/convex hull; distinguishes leaves with/without lobes, leaflets, etc.
  # width = width of the shape
  # length = length of the shape
  # ar = length-to-width ratio

# area
area <- coo_area(leaves) %>% 
  map_df(~data_frame(area = .x), .id = ".id")

# circularity
circularity <- coo_circularity(leaves) %>%
  map_df(~data_frame(circularity = .x), .id = ".id")

# convex hull area 
convex <- coo_convexity(leaves) %>%
  map_df(~data_frame(convex_area = .x), .id = ".id")

# solidity
solidity <- coo_solidity(leaves) %>%
  map_df(~data_frame(solidity = .x), .id = ".id")

# width
w <- coo_width(leaves) %>%
  map_df(~data_frame(width = .x), .id = ".id")

# length
l <- coo_length(leaves) %>%
  map_df(~data_frame(length = .x), .id = ".id")

# combine into a single data frame
all <- Reduce(function(x, y) 
  merge(x, y, all=TRUE, by = ".id"), 
  list(area, circularity, convex, solidity, w, l))

# aspect ratio 
all <- all %>%
  mutate(ar = width/length) %>%
  left_join(leaves$fac, ., by = c("img" = ".id")) %>%
  mutate_if(is.factor, fct_drop)

write.csv(all, "~/Box/BoleraceaLeafScans/data/processed/leaf_phenos.csv", 
          row.names = FALSE)

################################################################################
## Check trait distributions and correlations 
################################################################################

# histograms
all %>% 
  gather(key = measurement, 
         value = value,
         -img, -id, -rep, -leafNum, -morphotype, -pheno) %>% 
  ggplot(., aes(value)) +
  geom_histogram() + 
  facet_wrap(~measurement, scales = "free")

ggsave("~/Box/BoleraceaLeafScans/reports/shape_descriptor_histograms.png")

# correlations
all %>% 
  select(-c(1:6)) %>% 
  ggpairs() +
  theme_bw() +
  theme(strip.text = element_text(size = 14))

ggsave("~/Box/BoleraceaLeafScans/reports/shape_descriptor_correlations.png")

# boxplots
all %>% 
  mutate(morphotype = fct_reorder(morphotype, -area)) %>% 
  gather(key = measurement, 
         value = value,
         -img, -id, -rep, -leafNum, -morphotype, -pheno) %>% 
  ggplot(., aes(morphotype, value, fill = morphotype)) +
  geom_boxplot() + 
  facet_wrap(~measurement, scales = "free") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        strip.text = element_text(size = 14))

ggsave("~/Box/BoleraceaLeafScans/reports/shape_descriptor_boxplots.png")

################################################################################
## Export combined leaf outlines and measurements
################################################################################
save(leaves, leaves_aln, all,
     file = "~/Box/BoleraceaLeafScans/data/processed/aligned_scans.RData")
