## table name: Modified Hachinski Ischemia Scale [ADNI1,GO,2,3]
## also requires vitals_processing.R, egfr_processing.R, demographics_processing.R, and plasma_labdata_processing.R

library(tidyverse)
source("~/vitals_processing.R")
source("~/egfr_processing.R")

## variables for history of hypertension or stroke at baseline from the Hachinski Ischemic Index
hachinski <- readr::read_delim("~/Downloads/MODHACH_03Dec2024.csv") %>%
  dplyr::mutate(hypertension_bl_ind = case_when(
    HMHYPERT == 1 ~ 1,
    TRUE ~ 0
  ),
  stroke_bl_ind = case_when(
    HMSTROKE == 2 ~ 1,
    TRUE ~ 0
  )
  )

medical_covariates <- dplyr::full_join(vitals,hachinski %>%
                                         dplyr::select(RID,hypertension_bl_ind,stroke_bl_ind),
                                       by = c("RID"))

medical_covariates <- dplyr::full_join(medical_covariates,
                                       lab_data,
                                       by = "RID") 
