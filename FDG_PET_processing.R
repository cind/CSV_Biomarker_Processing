#this calculates pons-corrected meta-roi as described by the UC Berkeley FDG MetaROI methods on LONI
library(dplyr)

fdg_pet_mean <- read.csv("~/UCBERKELEYFDG_8mm_02_17_23_01Feb2024.csv") %>%
  dplyr::filter(!(ROINAME == "Top50PonsVermis")) %>%
  dplyr::mutate(EXAMDATE = as.Date(EXAMDATE),
                RID = as.character(RID)) %>%
  dplyr::rename(Meta_ROI_fdg_mean = MEAN) %>%
  dplyr::select(-ROINAME)

fdg_pons <- read.csv("~/UCBERKELEYFDG_8mm_02_17_23_01Feb2024.csv") %>%
  dplyr::filter(ROINAME == "Top50PonsVermis") %>%
  dplyr::mutate(mean_pons = MEAN) %>%
  dplyr::select(RID, VISCODE, VISCODE2, EXAMDATE, mean_pons, MAX, STDEV, TOTVOX, update_stamp)

fdg_pet <- merge(fdg_pet_mean, fdg_pons, by = c("RID", "VISCODE", "VISCODE2", "EXAMDATE"), all = TRUE) %>%
  dplyr::select(RID, EXAMDATE, Meta_ROI_fdg_mean, mean_pons) %>%
  dplyr::mutate(adjusted_Meta_ROI = Meta_ROI_fdg_mean/mean_pons)
