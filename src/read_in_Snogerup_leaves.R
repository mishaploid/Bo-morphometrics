################################################################################
## test read in of Snogerup et al. outlines
################################################################################

library(Momocs)
library(tidyverse)

spam <- list.files("~/Box/BoleraceaLeafScans/Snogerup_etal_Leaves/Page_10",
                   pattern = ".jpg",
                   full.names = TRUE)[2] %>%
  import_jpg()

eggs <- Out(spam)

pile(eggs)
