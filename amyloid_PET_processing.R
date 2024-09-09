## table name on LONI: 	UC Berkeley - Amyloid PET 6mm Res analysis [ADNI1,GO,2,3,4]

#loading in libraries
library(dplyr)

# Reading in data - pay attention to quality control -> qc_flag 
amyloid_pet <- readr::read_delim("~/data/UCBERKELEY_AMY_6MM_31Jul2024.csv")
amyloid_pet <- amyloid_pet %>% dplyr::rename(suvr_summary=SUMMARY_SUVR,
                                             Centiloid=CENTILOIDS,
                                             AmyloidPosPET=AMYLOID_STATUS, # amyloid status is based on SUVR and image resolution - more details explained in methods on LONI
                                             EXAMDATE_pet = SCANDATE)
