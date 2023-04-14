################################################################################
## Import outlines from Snogerup et al. 
## Authors: Sarah Turner-Hissong & Makenzie Mabry
## Date: 25 November 2019
## Updated: 27 May 2020
################################################################################

Snogerup <- list.files("data/Snogerup1990_leaf_scans",
                       pattern = ".jpg",
                       full.names = TRUE,
                       recursive = TRUE) %>%
  import_jpg() 

classifiers <- data.frame(img = names(Snogerup)) %>% 
  mutate(id = img,
         rep = word(id, start = 4, end = 4, sep = "_"),
         leafNum = rep,
         morphotype = word(id, start = 2, end = 2, sep = "_"),
         pheno = "Snogerup") %>%
  mutate_if(is.character, as.factor)

Snogerup_outlines <- Out(Snogerup, fac = classifiers) %>%
  Momocs::arrange(morphotype)

png("figures/misc-panel_Snogerup1990_leaves.png",
    width = 6,
    height = 6,
    units = "in",
    res = 300)
panel(Snogerup_outlines,
      fac = "morphotype",
      names = "morphotype",
      palette = pal_qual_Dark2)
dev.off()

Snogerup_outlines <- def_ldk(Snogerup_outlines, nb.ldk = 4) 

save(Snogerup_outlines, 
     file = "data/intermediate/Snogerup_landmarked.RData")