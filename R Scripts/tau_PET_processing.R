## table name on LONI: 	UC Berkeley - Amyloid PET 6mm Res analysis [ADNI1,GO,2,3,4]

#loading in libraries
library(dplyr)

# Reading in data - pay attention to quality control -> qc_flag 
tau_pet <- read.csv("~/data/UCBERKELEY_TAU_6MM_17Jun2024.csv")
tau_pet <- tau_pet %>% dplyr::rename(EXAMDATE_pet = SCANDATE)
