################################################################################
## Linear Discriminant Analysis for leaf outlines
## Authors: Sarah Turner-Hissong & Makenzie Mabry
## 10 November 2020
################################################################################

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(Momocs)
library(ggdendro)
library(ggplot2)
library(cowplot)
library(scales)


# setup -------------------------------------------------------------------
load("data/processed/EFDs.RData") # elliptical Fourier descriptors for leaves

# set color palette for morphotype ids (cols)
morph_names <- unique(leaves_aln$fac$morphotype) # extract morphotype ids
cols <- pal_qual_Dark2(length(morph_names)) # vector of unique colors
names(cols) <- morph_names


# PCA  --------------------------------------------------------------------
leaves.p <- leaves_f_sym %>% # symmetric variation only
  Momocs::filter(!pheno %in% "Snogerup") %>% # remove Snogerup images
  Momocs::mutate(id = fct_drop(id)) %>% 
  PCA(.) # principal component analysis

plot.new()

levels(leaves.p$fac$morphotype) <- morph_names


leaves.p %>% 
  plot_PCA(~morphotype,
           points = FALSE,
           chull = FALSE, 
           legend = TRUE,
           axes = c(1,2),
           palette = pal_qual_Dark2,
           title = "symmetric variation only",
           points_transp = 0.5,
           zoom = 1.25) %>%
  layer_points(cex = 1.2,
               transp = 0.3) 

##code below to choose images and put circles on figure
grid(12, 12, lwd = 2, col = "black")

pca_df <- leaves.p %>% 
  as_df() %>%
  select(img, id, rep, leafNum, morphotype, pheno, PC1, PC2, PC3, PC4) 

example_img_ids <- identify(pca_df$PC1, pca_df$PC2)


example_imgs <- pca_df %>% 
  rowid_to_column() %>%
  filter(rowid %in% example_img_ids)

symbols(x=example_imgs$PC1, 
        y=example_imgs$PC2, 
        circles = rep(0.008, nrow(example_imgs)), 
        add = TRUE, 
        inches = FALSE,
        fg = "black",
        lwd = 2)

##plot pc 2v3 and 3v4 for suppl.
pca2_3 <- leaves.p %>% 
  plot_PCA(~morphotype,
           points = FALSE,
           chull = FALSE, 
           legend = TRUE,
           axes = c(2,3),
           palette = pal_qual_Dark2,
           points_transp = 0.5,
           zoom = 1) %>%
  layer_points(cex = 1.2,
               transp = 0.3) 
  
  
pca3_4 <- pca2_3 <- leaves.p %>% 
  plot_PCA(~morphotype,
           points = FALSE,
           chull = FALSE, 
           legend = TRUE,
           axes = c(3,4),
           palette = pal_qual_Dark2,
           points_transp = 0.5,
           zoom = 1.25) %>%
  layer_points(cex = 1.2,
               transp = 0.3) 

plot_grid(pca2_3, pca3_4)

ggsave("~/Box Sync/BoleraceaLeafScans/reports/",
       height = 11,
       width = 8.5)

# pca of traditonal values
all <- read.csv("data/processed/leaf_phenos.csv")
pca_outliers <- c("B_oleracea012_A_1_p1_mask", 
                  "B_oleracea023_B_2_p2_mask",
                  "B_oleracea034_B_leaf2_3_1_p1_mask", 
                  "B_oleracea117_B_2_p2_mask",
                  "B_oleracea157_B_0_p0_mask", 
                  "B_oleracea207_C_0_p0_mask",
                  "B_oleracea218_A_0_p0_mask", 
                  "B_oleracea232_D_0_p0_mask")

all <-  all %>%  filter(!img %in% pca_outliers)

all$id <- as.factor(all$id)
all$leafNum <- as.factor(all$leafNum)
leaves.p$fac$img <- as.factor(leaves.p$fac$img)

##for figure 1
leaves.p %>% 
  as_df() %>% 
  left_join(all, ., by = c("img" = "img")) %>% 
  select(id.x:PC2) %>% 
  gather(key = measurement,
         value = value, 
         -c(area, id.x, id.y, leafNum.y, length, morphotpe_Num, morphotype.x, morphotype.y, pheno.x, pheno.y, rep.x, rep.y, PC1:PC2)) %>% 
         #-c(leafNum.x, convex_area, ar, circularity, solidity, width, PC1:PC2)) %>% 
  group_by(measurement) %>% 
  mutate(percentile = percent_rank(value)) %>% 
  ggplot(., aes(PC1, PC2, color = percentile)) +
  scale_color_viridis_c() +
  #scale_colour_gradientn(colours = c("#1b9e77","#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")) +
  geom_point() +
  facet_wrap(~measurement,
             ncol = 3) +
  theme_bw() +
  theme(strip.text = element_text(size = 14))


##for supplemental figure
leaves.p %>% 
  as_df() %>% 
  left_join(all, ., by = c("img" = "img")) %>% 
  select(id.x:PC2) %>% 
  gather(key = measurement,
         value = value, 
         -c(id.x, id.y, leafNum.y, morphotpe_Num, morphotype.x, morphotype.y, pheno.x, pheno.y, rep.x, rep.y, leafNum.x, convex_area, ar, PC1:PC2)) %>% 
  #-c(leafNum.x, convex_area, ar, PC1:PC2)) %>% 
  group_by(measurement) %>% 
  mutate(percentile = percent_rank(value)) %>% 
  ggplot(., aes(PC1, PC2, color = percentile)) +
  scale_color_viridis_c() +
  #scale_colour_gradientn(colours = c("#1b9e77","#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")) +
  geom_point() +
  facet_wrap(~measurement,
             ncol = 3) +
  theme_bw() +
  theme(strip.text = element_text(size = 14))


# eigenleaf representations for PCs representing variance in EFDs 

gg <- PCcontrib(leaves.p,
                nax = 1:4,
                sd.r = c(-3, -1.5, 0, 1.5, 3))

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

ggsave("reports/figures/Fig3_eigenleaf_PC_variation.png")

### plot Matrix of pairwise combinations of eigenvectors and iso plots

leaves_f_sym %>% 
  MSHAPES(~morphotype) %>% 
  plot_MSHAPES(palette=pal_manual(c("#1B9E77", "#7570B3")))

ggsave("~/Box Sync/BoleraceaLeafScans/reports/figures/SuppFig3_pairwise_eigenleaf_variation.png", 
       height = 6,
       width = 8.5)

mean_species <- MSHAPES(leaves_f_sym,
                        fac = "morphotype")$shp


tps_iso(mean_species$alboglabra, 
        mean_species$sabauda, 
        cont = FALSE, 
        amp = 1,
        grid = TRUE,
        shp.border=pal_qual_Dark2(22)[c(8,5)],
        palette = viridis_pal(option = "D", alpha = 0),
        shp.lwd = c(5,5))


p2 <- leaves.p %>% as_df() %>% 
  filter(morphotype %in% c("macrocarpa", "villosa")) %>%
  ggplot() +
  aes(x = PC1, y = PC2, col = morphotype) + 
  scale_color_manual(values = cols) +
  ggtitle("B. macropcarpa v B. villosa") +
  coord_equal() + 
  geom_point() + 
  geom_density2d() + 
  theme_light() +
  theme(legend.position = "none")

p3 <- leaves.p %>% as_df() %>% 
  filter(morphotype %in% c("palmifolia", "villosa")) %>%
  ggplot() +
  aes(x = PC1, y = PC2, col = morphotype) + 
  scale_color_manual(values = cols) +
  ggtitle("Lacinatio kale v B. villosa") +
  coord_equal() + 
  geom_point() + 
  geom_density2d() + 
  theme_light() +
  theme(legend.position = "none")

p4 <- leaves.p %>% as_df() %>% 
  filter(morphotype %in% c("alboglabra", "sabauda")) %>%
  ggplot() +
  aes(x = PC1, y = PC2, col = morphotype) + 
  scale_color_manual(values = cols) +
  ggtitle("Chinese white kale v savoy cabbage") +
  coord_equal() + 
  geom_point() + 
  geom_density2d() + 
  theme_light() +
  theme(legend.position = "none")


plot_grid(p2, p3, p4, labels = c('A', 'B', 'C'), ncol = 1) ## added the iso leaves to these in inkscape

ggsave("~/Box Sync/BoleraceaLeafScans/reports/Density2D_tpsiso.pdf",
       height = 11,
       width = 8.5)

# LDA  --------------------------------------------------------------------
leaves_f_sym %>%
  LDA(~morphotype) %>%
  classification_metrics() 

leaves.l_id <- LDA(leaves.p,
                'id', 
                retain = 0.99)

leaves.l_morphotype <- LDA(leaves.p,
                'morphotype', 
                retain = 0.99)


LDA_Morph <- plot_LDA(leaves.l_morphotype, legend=TRUE, palette = cols)
LDA_ID <- plot_LDA(leaves.l_id, chull=FALSE)


leaves.cv_morph <- leaves.l_morphotype$CV.tab # table with cross-validation results for morphotype
leaves.cv_id <- leaves.l_id$CV.tab # table with cross-validation results for morphotype

#plot_CV(leaves.l_morphotype, freq = TRUE, rm0 = FALSE) #momocs internal one


# plot confusion matrix for morphotype
counts_morph <- leaves_Snogerup$fac %>%
  group_by(morphotype) %>% 
  count() 

prop_correct_morph <- leaves.cv_morph %>% 
  as_tibble() %>% 
  left_join(., counts_morph, by = c("actual" = "morphotype")) %>% 
  filter(actual == classified) %>%
  mutate(prop_correct_morph = n.x/n.y)

leaves.cv_morph2 <-leaves.cv_morph %>% 
  as_tibble() %>% 
  left_join(., counts_morph, by = c("actual" = "morphotype")) %>% 
  left_join(., prop_correct_morph[,c("actual", "prop_correct_morph")], by = "actual") %>% 
  mutate(prop = n.x/n.y,
         actual = fct_reorder(actual, prop_correct_morph),
         classified = fct_relevel(classified, levels(fct_rev(actual))))


plot_cv_morph <-ggplot(leaves.cv_morph2, aes(classified, actual, fill = prop)) + 
  geom_tile() +
  scale_fill_viridis(option = "mako", direction = -1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 8))

ggsave("/Users/mmabry/OneDrive/Computer/Projects/BoleraceaMorph/reports/lda_confusion_matrix_morphotype_V2.png",
       height = 6,
       width = 8.5)

# plot confusion matrix for id
counts_id <- leaves_Snogerup$fac %>%
  group_by(id) %>% 
  count() 

prop_correct_id <- leaves.cv_id %>% 
  as_tibble() %>% 
  left_join(., counts_id, by = c("actual" = "id")) %>% 
  filter(actual == classified) %>%
  mutate(prop_correct_id = n.x/n.y)

leaves.cv_id2 <- leaves.cv_id %>% 
  as_tibble() %>% 
  left_join(., counts_id, by = c("actual" = "id")) %>% 
  left_join(., prop_correct_id[,c("actual", "prop_correct_id")], by = "actual") %>% 
  mutate(prop = n.x/n.y,
         actual = fct_reorder(actual, -prop_correct_id),
         classified = fct_relevel(classified, levels(actual)))

plot_cv_id <- ggplot(leaves.cv_id2, aes(classified, actual, fill = prop)) + 
  geom_tile() +
  scale_fill_viridis(option = "mako", direction = -1) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 2), axis.text.y = element_text(size = 2),
        text = element_text(size = 8))

ggsave("~/Box Sync/BoleraceaLeafScans/reports/lda_confusion_matrix_id.png",
       height = 6,
       width = 8.5)

#plot both together
plot_grid(plot_cv_morph, plot_cv_id, labels = c('A', 'B'), rel_widths = c(2, 3), ncol = 1)

ggsave("~/Box Sync/BoleraceaLeafScans/reports/lda_confusion_matrix_id_morphotype.png",
       height = 11,
       width = 8.5)


# MANOVA  --------------------------------------------------------------------

res.man <- MANOVA(leaves.p, ~morphotype)
summary_MANOVA <- MANOVA_PW(leaves.p, ~morphotype)

write.csv(summary_MANOVA$summary, "~/Box Sync/BoleraceaLeafScans/reports/MANOVA_summary.csv")
write.csv(summary_MANOVA$stars.tab, "~/Box Sync/BoleraceaLeafScans/reports/MANOVA_summaryStars.csv")

# Hierarchical clustering  ----------------------------------------------------

leaves.hc <- hclust(d  = dist(x = leaves.p$x,
                              method = "euclidean"))

leaves.dendr <- dendro_data(leaves.hc, 
                            type = "rectangle")

morph_labels <- data.frame(label = leaves.p$fac$img,
                           morphotype = leaves.p$fac$morphotype)

leaves.dendr[["labels"]] <- merge(leaves.dendr[["labels"]],
                                  morph_labels,
                                  by = "label")

## to find the best # clusters  (https://statsandr.com/blog/clustering-analysis-k-means-and-hierarchical-clustering-by-hand-and-in-r/#optimal-number-of-clusters-1)
barplot(leaves.hc$height,
        names.arg = (nrow(leaves.p$x) - 1):1 # show the number of cluster below each bars
        )

diff_one_two_clusters <- leaves.hc$height[2493] - leaves.hc$height[2492]  #0.17377

diff_two_three_clusters <- leaves.hc$height[2492] - leaves.hc$height[2491] ##0.1598


#ggplot() + 
#  geom_segment(data=segment(leaves.dendr), 
#               aes(x=x, y=y, xend=xend, yend=yend)) + 
#  geom_text(data=label(leaves.dendr), 
#            aes(x, y, label=label, hjust=0, color=morphotype), 
#            size=2) +
#  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
#  theme(axis.line.y=element_blank(),
#        axis.ticks.y=element_blank(),
#        axis.text.y=element_blank(),
#        axis.title.y=element_blank(),
#        panel.background=element_rect(fill="white"),
#        panel.grid=element_blank())


#ggdendrogram(leaves.dendr)

#http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning

library(ape)
library(ggtree)
library(ggnewscale)

clust3 <- cutree(leaves.hc, 3)
clust2 <- cutree(leaves.hc, 2)

##first get tree in correct format
tree <- as.phylo(leaves.hc)
circ <- ggtree(tree, layout = "circular")

##dataframe for morphotype
df <- data.frame(morphotype = morph_labels$morphotype) 
rownames(df) <- tree$tip.label

## dataframe for clust2 groups
df2 <- as.data.frame(as.factor(clust2)) 

##dataframe for wild/cultivars
df3 <- df %>% mutate(type = case_when(
  morphotype == "palmifolia" ~ "cultivar",
  morphotype == "gongylodes" ~ "cultivar",
  morphotype == "capitata" ~ "cultivar",
  morphotype == "gemmifera"  ~ "cultivar",
  morphotype == "sabauda" ~ "cultivar", 
  morphotype == "ramosa" ~ "cultivar",     
  morphotype == "botrytis" ~ "cultivar",  
  morphotype == "alboglabra" ~ "cultivar",
  morphotype == "costata"  ~ "cultivar", 
  morphotype == "italica"  ~ "cultivar",  
  morphotype == "oleracea" ~ "wild",  
  morphotype == "viridis" ~ "cultivar",   
  morphotype == "cretica"  ~ "wild",  
  morphotype == "hilarionis" ~ "wild",
  morphotype == "incana"  ~ "wild",   
  morphotype == "insularis"  ~ "wild",
  morphotype == "macrocarpa" ~ "wild",
  morphotype == "montana"   ~ "wild",
  morphotype == "rupestris"  ~ "wild",
  morphotype == "villosa" ~ "wild",   
  morphotype == "medullosa"  ~ "cultivar",
  morphotype == "sabellica" ~ "cultivar"
))

df3 <- as.data.frame(df3$type)
rownames(df3) <- tree$tip.label


##dataframe for leaf num
df4 <- data.frame(leafNum = leaves.p$fac$leafNum)
rownames(df4) <- tree$tip.label

##make first plot with tree and morphotype
p1 <- gheatmap(circ, df, width=.1, colnames = FALSE, color = NA) + 
  scale_fill_manual(values = cols, name="morphotype")

##make new scale
p2 <- p1 + new_scale_fill()

##add second plot with wild/clutivars
p3 <- gheatmap(p2, df4, width=.1,
               offset = 0.05, colnames =FALSE, color = NA) +
              scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3"), name = "leafNum")

##add new scale
p4 <- p3 + new_scale_fill()

##add leaf num
p5 <- gheatmap(p4, df3, width=.1,
               offset = 0.1, colnames =FALSE, color = NA) +
              scale_fill_manual(values=c("#1B9E77", "#7570B3"), name = "typw")

##add new scale
p6 <- p5 + new_scale_fill()

##add cluster groups (2)
p7 <- gheatmap(p6, df2, width=.1,
         offset = 0.15, colnames =FALSE, color = NA) +
          scale_fill_manual(values=c("#1B9E77", "#7570B3"))

p7

df5 <- as.data.frame(as.factor(clust3)) 

p8 <- gheatmap(p6, df5, width=.1,
              offset = 0.15, colnames =FALSE, color = NA) +
              scale_fill_manual(values=c("#1B9E77", "#D95F02", "#7570B3"))

p8

#plot both together
plot_grid(p7, p8, labels = c('A', 'B'), ncol = 1)

ggsave("~/Box Sync/BoleraceaLeafScans/reports/denograms_3_2.pdf",
       height = 11,
       width = 8.5)

##getting mean shapes for including in inkscape

library(dplyr)

leaves_f_sym %>%
  Momocs::filter(img %in% rownames(df2)[df2$`as.factor(clust2)`==1]) %>% 
  MSHAPES() %>%
  coo_plot(border = "#1B9E77",
           lwd = 5, 
           centroid = FALSE, 
           xy.axis = FALSE)


leaves_f_sym %>%
  Momocs::filter(img %in% rownames(df2)[df2$`as.factor(clust2)`==2]) %>% 
  MSHAPES() %>% 
  coo_plot(border = "#7570B3",
           lwd = 5, 
           centroid = FALSE, 
           xy.axis = FALSE)

leaves_f_sym %>%
  Momocs::filter(img %in% rownames(df5)[df5$`as.factor(clust3)`==1]) %>% 
  MSHAPES() %>% 
  coo_plot(border = "#1B9E77",
           lwd = 5, 
           centroid = FALSE, 
           xy.axis = FALSE)

leaves_f_sym %>%
  Momocs::filter(img %in% rownames(df5)[df5$`as.factor(clust3)`==2]) %>% 
  MSHAPES() %>% 
  coo_plot(border = "#D95F02",
           lwd = 5, 
           centroid = FALSE, 
           xy.axis = FALSE)

leaves_f_sym %>%
  Momocs::filter(img %in% rownames(df5)[df5$`as.factor(clust3)`==3]) %>% 
  MSHAPES() %>% 
  coo_plot(border = "#7570B3",
           lwd = 5, 
           centroid = FALSE, 
           xy.axis = FALSE)
