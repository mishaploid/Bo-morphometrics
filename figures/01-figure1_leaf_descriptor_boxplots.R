################################################################################
## Bo-morphometrics: Figure 1 Boxplot of traditional shape descriptors
## Authors: Sarah Turner-Hissong & Makenzie Mabry
## Date: 25 November 2019
## Updated: 31 March 2023
################################################################################

library(tidyverse)
library(Momocs)
library(RColorBrewer)
library(cowplot)
library(ggridges)

load('data/processed/02-aligned_scans.RData')


# recode scientific to common names ---------------------------------------

new_names <- c('Lacinato kale' ='palmifolia', 
               'Kohlrabi' = 'gongylodes', 
               'Cabbage' = 'capitata', 
               'Brussels sprouts' = 'gemmifera', 
               'Savoy cabbage' = 'sabauda', 
               'Perpetual kale' = 'ramosa', 
               'Cauliflower' = 'botrytis', 
               'Chinese white kale' = 'alboglabra', 
               'Tronchuda kale' = 'costata', 
               'Broccoli' = 'italica', 
               'B. oleracea' = 'oleracea', 
               'Collard greens' = 'viridis', 
               'B. cretica' = 'cretica', 
               'B. hilarionis' = 'hilarionis', 
               'B. incana' = 'incana', 
               'B. insularis' = 'insularis', 
               'B. macrocarpa' = 'macrocarpa', 
               'B. montana' = 'montana', 
               'B. rupestris' = 'rupestris', 
               'B. villosa' = 'villosa', 
               'Marrow cabbage' = 'medullosa', 
               'Curly kale' = 'sabellica')

all <- all %>% 
  mutate(comName = fct_recode(morphotype, !!!new_names))

# assign colors to morphotypes --------------------------------------------

# morphotype names 
morph_names <- unique(all$comName)
# color palette 
cols <- pal_qual_Dark2(length(morph_names))
# assign names to colors 
names(cols) <- morph_names

# set factor levels for measurements 
all <- all %>% 
  mutate(comName = fct_reorder(comName, -area)) 

all_long <- all %>% 
  # order levels by type of measurement
  gather(key = measurement, 
         value = value,
         -img, -id, -rep, -leafNum, -comName, -morphotype, -pheno) %>% 
  mutate(measurement = fct_relevel(measurement, 
                                   'area', 'width', 'length', 'aspect_ratio', 
                                   'convexity', 'solidity', 'circularity'))


# Fig 1 - boxplot of shape descriptors ------------------------------------

ggplot(all_long, aes(comName, value, fill = comName),
       alpha = 0.8) +
  geom_boxplot(position = position_dodge(0.9)) + 
  theme_bw(base_size = 12) + 
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 60, 
                                   hjust = 1, 
                                   vjust = 1,
                                   size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "none") +
  facet_wrap(measurement ~ ., scales = 'free_y',
             ncol = 2) + 
  xlab('morphotype common name')

ggsave("figures/Fig1-shape_descriptor_boxplots.png",
       height = 7, width = 8)


# Fig S1 - boxplots split by leaf number ----------------------------------

all %>% 
  mutate(comName = fct_reorder(comName, -area)) %>% 
  gather(key = measurement, 
         value = value,
         -img, -id, -rep, -leafNum, -comName, -morphotype, -pheno) %>% 
  ggplot(., aes(comName, value)) +
  geom_boxplot(aes(fill = leafNum),
               alpha = 0.8) +
  facet_grid(measurement ~ ., scales = "free_y") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   vjust = 1)) + 
  scale_fill_viridis_d('Leaf Number') + 
  xlab('morphotype common name')

ggsave("figures/FigS1-shape_descriptor_boxplots_by_leafnum.png", 
       width = 10, height = 10)


# density plots for distinguishing features -------------------------------

# convexity
con_plot <- ggplot(all, aes(x = convexity, 
                y = fct_reorder(comName, -convexity),
                fill = comName)) + 
geom_density_ridges2(jittered_points = TRUE, 
                     scale = .95, 
                     rel_min_height = .01,
                     point_shape = "|", 
                     point_size = 1, 
                     size = 0.25,
                     position = position_points_jitter(height = 0)) +
  theme_ridges(center_axis_labels = TRUE) + 
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        text = element_text(size = 12),
        axis.text = element_text(size = 10)) + 
  scale_fill_manual(values = cols) +
  ylab('Common Name') 

# circularity 
cir_plot <- ggplot(all, aes(x = circularity, 
                y = fct_reorder(comName, -circularity),
                fill = comName)) + 
  geom_density_ridges2(jittered_points = TRUE, 
                       scale = .95, 
                       rel_min_height = .01,
                       point_shape = "|", 
                       point_size = 1, 
                       size = 0.25,
                       position = position_points_jitter(height = 0)) +
  theme_ridges(center_axis_labels = TRUE) + 
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        text = element_text(size = 12),
        axis.text = element_text(size = 10)) + 
  scale_fill_manual(values = cols) +
  ylab('Common Name')

# aspect ratio
ar_plot <- ggplot(all, aes(x = aspect_ratio, 
                y = fct_reorder(comName, -aspect_ratio),
                fill = comName)) + 
  geom_density_ridges2(jittered_points = TRUE, 
                       scale = .95, 
                       rel_min_height = .01,
                       point_shape = "|", 
                       point_size = 1, 
                       size = 0.25,
                       position = position_points_jitter(height = 0)) +
  theme_ridges(center_axis_labels = TRUE) + 
  theme(legend.position = 'none',
        axis.title.y = element_blank(),
        text = element_text(size = 12),
        axis.text = element_text(size = 10)) + 
  scale_fill_manual(values = cols) +
  ylab('Common Name') 

plot_grid(con_plot, cir_plot, ar_plot, nrow = 1)

ggsave('figures/misc-ggridge_plots_for_key_features.png',
       bg = 'white',
       height = 6,
       width = 10)
