## table name on LONI: UC Berkeley - Amyloid PET 6mm Res analysis [ADNI1,GO,2,3,4]

library(tidyverse)

amyloid_pet <-readr::read_delim("~/Downloads/UCBERKELEY_AMY_6MM_04Dec2024.csv") ## change to your path

amyloid_pet <- amyloid_pet %>%
  dplyr::mutate(
    amyloid_pet_cs_pos = AMYLOID_STATUS, ## use amyloid_pet_cs_pos for cross-sectional analyses
    amyloid_pet_long_pos = AMYLOID_STATUS_COMPOSITE_REF, ## use amyloid_pet_long_pos for longitudinal analyses
    amyloid_pet_centiloid_pos = as.numeric(CENTILOIDS > 22))

amyloid_pet_all_regions <- amyloid_pet

amyloid_pet <- amyloid_pet %>%
  dplyr::select(RID,VISCODE,SCANDATE,
                amyloid_pet_cs_pos,amyloid_pet_long_pos,amyloid_pet_centiloid_pos,
                CENTILOIDS,SUMMARY_SUVR,SUMMARY_VOLUME,
                COMPOSITE_REF_SUVR,COMPOSITE_REF_VOLUME)

